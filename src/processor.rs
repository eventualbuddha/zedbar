//! Processor lifecycle functions
//!
//! This module provides the processor creation and destruction functions.

use libc::{c_int, c_void, calloc, free};
use std::mem::size_of;
use std::ptr::null_mut;

use crate::error::{ErrInfo, _zbar_err_init, ZBAR_MOD_PROCESSOR};

// External C functions
extern "C" {
    fn zbar_image_scanner_create() -> *mut c_void;
    fn zbar_image_scanner_destroy(scanner: *mut c_void);
}

// Forward declaration for symbol set (opaque pointer from C)
#[repr(C)]
#[allow(non_camel_case_types)]
pub struct zbar_symbol_set_t {
    _private: [u8; 0],
}

/// Processor structure - must match processor.h layout exactly
#[repr(C)]
#[allow(non_camel_case_types)]
pub struct zbar_processor_t {
    err: ErrInfo,
    scanner: *mut c_void,
    syms: *const zbar_symbol_set_t,
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
        zbar_symbol_set_ref((*proc).syms as *const c_void, -1);
        (*proc).syms = null_mut();
    }

    if !(*proc).scanner.is_null() {
        zbar_image_scanner_destroy((*proc).scanner);
        (*proc).scanner = null_mut();
    }

    free(proc as *mut c_void);
}
