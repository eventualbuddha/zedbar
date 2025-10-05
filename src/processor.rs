//! Processor lifecycle functions
//!
//! This module provides the processor creation and destruction functions.

use libc::{c_int, c_void, calloc, free};
use std::mem::size_of;
use std::ptr::null_mut;

use crate::{
    error::{ErrInfo, _zbar_err_init, ZBAR_MOD_PROCESSOR},
    img_scanner::{
        zbar_image_scanner_get_results, zbar_image_scanner_recycle_image, zbar_image_scanner_t,
        zbar_symbol_set_t,
    },
    zbar_image_scanner_create, zbar_image_scanner_destroy, zbar_image_t, zbar_scan_image,
};

/// Processor structure - must match processor.h layout exactly
#[repr(C)]
#[allow(non_camel_case_types)]
pub struct zbar_processor_t {
    err: ErrInfo,
    scanner: *mut zbar_image_scanner_t,
    syms: *mut zbar_symbol_set_t,
}

/// Create a new processor instance
///
/// # Parameters
/// - `threaded`: Threading parameter (currently unused)
///
/// # Returns
/// Pointer to new processor or null on allocation failure
#[no_mangle]
pub unsafe extern "C" fn zbar_processor_create(threaded: c_int) -> *mut zbar_processor_t {
    let _ = threaded; // Currently unused

    let proc = calloc(1, size_of::<zbar_processor_t>()) as *mut zbar_processor_t;
    if proc.is_null() {
        return null_mut();
    }

    _zbar_err_init(&mut (*proc).err, ZBAR_MOD_PROCESSOR);

    (*proc).scanner = zbar_image_scanner_create();
    if (*proc).scanner.is_null() {
        free(proc as *mut c_void);
        return null_mut();
    }

    proc
}

/// Destroy a processor instance
///
/// Frees all resources associated with the processor.
///
/// # Parameters
/// - `proc`: Processor to destroy (null-safe)
#[no_mangle]
pub unsafe extern "C" fn zbar_processor_destroy(proc: *mut zbar_processor_t) {
    if proc.is_null() {
        return;
    }

    if !(*proc).syms.is_null() {
        use crate::symbol::zbar_symbol_set_ref;
        zbar_symbol_set_ref((*proc).syms, -1);
        (*proc).syms = null_mut();
    }

    if !(*proc).scanner.is_null() {
        zbar_image_scanner_destroy((*proc).scanner);
        (*proc).scanner = null_mut();
    }

    free(proc as *mut c_void);
}

/// Process an image
///
/// Scans an image for barcodes and stores the results in the processor.
///
/// # Parameters
/// - `proc`: Processor instance
/// - `img`: Image to process
///
/// # Returns
/// Number of symbols found, or -1 on error
#[no_mangle]
pub unsafe extern "C" fn zbar_process_image(
    proc: *mut zbar_processor_t,
    img: *mut zbar_image_t,
) -> c_int {
    if proc.is_null() || img.is_null() {
        return -1;
    }

    // Clean up previous results
    if !(*proc).syms.is_null() {
        use crate::symbol::zbar_symbol_set_ref;
        zbar_symbol_set_ref((*proc).syms, -1);
        (*proc).syms = null_mut();
    }

    // Process the image
    zbar_image_scanner_recycle_image((*proc).scanner, img);
    let nsyms = zbar_scan_image((*proc).scanner, img);

    if nsyms < 0 {
        return nsyms;
    }

    // Store results
    (*proc).syms = zbar_image_scanner_get_results((*proc).scanner);
    if !(*proc).syms.is_null() {
        use crate::symbol::zbar_symbol_set_ref;
        zbar_symbol_set_ref((*proc).syms, 1);
    }

    nsyms
}
