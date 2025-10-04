use std::{ffi::c_void, ptr::null_mut};

use libc::{c_int, c_uint, c_ulong};

use crate::{line_scanner::zbar_scanner_t, zbar_image_t};

const RECYCLE_BUCKETS: usize = 5;
const NUM_SCN_CFGS: usize = 2; // ZBAR_CFG_Y_DENSITY - ZBAR_CFG_X_DENSITY + 1
const NUM_SYMS: usize = 25; // Number of symbol types

// Forward declarations for opaque C types
#[repr(C)]
#[allow(non_camel_case_types)]
pub struct zbar_symbol_t {
    _private: [u8; 0],
}

#[repr(C)]
#[allow(non_camel_case_types)]
pub struct zbar_decoder_t {
    _private: [u8; 0],
}

#[repr(C)]
#[allow(non_camel_case_types)]
pub struct qr_reader {
    _private: [u8; 0],
}

#[repr(C)]
#[allow(non_camel_case_types)]
pub struct sq_reader {
    _private: [u8; 0],
}

#[repr(C)]
#[allow(non_camel_case_types)]
pub struct zbar_symbol_set_t {
    _private: [u8; 0],
}

#[repr(C)]
#[allow(non_camel_case_types)]
pub struct zbar_image_data_handler_t {
    _private: [u8; 0],
}

#[allow(non_camel_case_types)]
pub struct recycle_bucket_t {
    nsyms: c_int,
    head: *mut zbar_symbol_t,
}

/// image scanner state
#[repr(C)]
#[allow(non_camel_case_types)]
pub struct zbar_image_scanner_t {
    /// associated linear intensity scanner
    scn: *mut zbar_scanner_t,

    /// associated symbol decoder
    dcode: *mut zbar_decoder_t,
    /// QR Code 2D reader
    qr: *mut qr_reader,
    /// SQ Code 2D reader
    sq: *mut sq_reader,

    /// application data
    userdata: *const c_void,
    /// user result callback
    handler: *mut zbar_image_data_handler_t,

    /// scan start time
    time: c_ulong,
    /// currently scanning image *root*
    img: *mut zbar_image_t,
    /// current scan direction
    dx: c_int,
    dy: c_int,
    du: c_int,
    umin: c_int,
    v: c_int,

    /// previous decode results
    syms: *mut zbar_symbol_set_t,
    /// recycled symbols in 4^n size buckets
    recycle: [recycle_bucket_t; RECYCLE_BUCKETS],

    /// current result cache state
    enable_cache: c_int,
    /// inter-image result cache entries
    cache: *mut zbar_symbol_t,

    // configuration settings
    /// config flags
    config: c_uint,
    ean_config: c_uint,
    /// int valued configurations
    configs: [c_int; NUM_SCN_CFGS],
    /// per-symbology configurations
    sym_configs: [[c_int; 1]; NUM_SYMS],
}

/// Get the current set of decoded symbols from the image scanner
///
/// Returns the symbol set containing all symbols detected during the last scan.
/// The returned symbol set is still owned by the scanner and should not be freed.
#[no_mangle]
pub unsafe extern "C" fn zbar_image_scanner_get_results(
    iscn: *const zbar_image_scanner_t,
) -> *const zbar_symbol_set_t {
    if iscn.is_null() {
        return null_mut();
    }
    (*iscn).syms
}
