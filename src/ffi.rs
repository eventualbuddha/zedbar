//! FFI bindings to the C ZBar library
//!
//! This module provides direct bindings to the original C implementation.
//! As modules are converted to Rust, these bindings will be gradually removed.

use libc::{c_char, c_int, c_uint, c_ulong, c_void};

use crate::{
    decoder_types::zbar_decoder_t,
    img_scanner::{zbar_image_scanner_t, zbar_symbol_set_t},
    line_scanner::zbar_scanner_t,
};

#[repr(C)]
#[allow(non_camel_case_types)]
pub struct zbar_image_t {
    pub format: u32,
    pub width: c_uint,
    pub height: c_uint,
    pub data: *mut c_void,
    pub datalen: c_ulong,
    pub userdata: *mut c_void,
    pub cleanup: *mut c_void,
    pub refcnt: c_int,
    pub srcidx: c_int,
    pub next: *mut zbar_image_t,
    pub seq: c_uint,
    pub syms: *mut zbar_symbol_set_t,
}

#[repr(C)]
#[allow(non_camel_case_types)]
pub struct zbar_symbol_t {
    pub symbol_type: c_int,
    pub configs: c_uint,
    pub modifiers: c_uint,
    pub data_alloc: c_uint,
    pub datalen: c_uint,
    pub data: *mut c_char,
    pub pts_alloc: c_uint,
    pub npts: c_uint,
    pub pts: *mut c_void,
    pub orient: c_int,
    pub refcnt: c_int,
    pub next: *mut zbar_symbol_t,
    pub syms: *mut zbar_symbol_set_t,
    pub time: c_ulong,
    pub cache_count: c_int,
    pub quality: c_int,
}

// Reference counting helper
pub(crate) unsafe fn refcnt(cnt: *mut c_int, delta: c_int) -> c_int {
    let rc = *cnt + delta;
    *cnt = rc;
    debug_assert!(rc >= 0);
    rc
}

#[no_mangle]
pub unsafe extern "C" fn _zbar_refcnt(cnt: *mut c_int, delta: c_int) -> c_int {
    refcnt(cnt, delta)
}

#[link(name = "zbar_c", kind = "static")]
extern "C" {
    // Version info
    pub fn zbar_version(major: *mut c_uint, minor: *mut c_uint, patch: *mut c_uint) -> c_int;

    // Error handling
    pub fn zbar_set_verbosity(verbosity: c_int);
    pub fn zbar_increase_verbosity();

    // Symbol name functions
    pub fn zbar_get_symbol_name(sym: c_int) -> *const c_char;
    pub fn zbar_get_config_name(config: c_int) -> *const c_char;

    // Image functions
    pub fn zbar_image_create() -> *mut zbar_image_t;
    pub fn zbar_image_destroy(image: *mut zbar_image_t);
    pub fn zbar_image_ref(image: *mut zbar_image_t, refs: c_int);

    pub fn zbar_image_get_format(image: *const zbar_image_t) -> c_ulong;
    pub fn zbar_image_get_width(image: *const zbar_image_t) -> c_uint;
    pub fn zbar_image_get_height(image: *const zbar_image_t) -> c_uint;
    pub fn zbar_image_get_data(image: *const zbar_image_t) -> *const c_void;

    pub fn zbar_image_set_format(image: *mut zbar_image_t, format: c_ulong);
    pub fn zbar_image_set_size(image: *mut zbar_image_t, width: c_uint, height: c_uint);
    pub fn zbar_image_set_data(
        image: *mut zbar_image_t,
        data: *const c_void,
        data_len: c_ulong,
        cleanup: *const c_void,
    );

    // Scanner functions
    pub fn zbar_scanner_create(dcode: *mut zbar_decoder_t) -> *mut zbar_scanner_t;
    pub fn zbar_image_scanner_create() -> *mut zbar_image_scanner_t;
    pub fn zbar_image_scanner_destroy(scanner: *mut zbar_image_scanner_t);
    pub fn zbar_image_scanner_set_config(
        scanner: *mut zbar_image_scanner_t,
        symbology: c_int,
        config: c_int,
        value: c_int,
    ) -> c_int;
    pub fn zbar_scan_image(scanner: *mut zbar_image_scanner_t, image: *mut zbar_image_t) -> c_int;
    pub fn zbar_image_first_symbol(image: *const zbar_image_t) -> *const c_void;

    // Symbol functions
    pub fn zbar_symbol_get_type(symbol: *const c_void) -> c_int;
    pub fn zbar_symbol_get_data(symbol: *const c_void) -> *const c_char;
    pub fn zbar_symbol_get_data_length(symbol: *const c_void) -> c_uint;
    pub fn zbar_symbol_next(symbol: *const c_void) -> *const c_void;

    // Internal scanner functions (from img_scanner.h)
    pub fn zbar_image_scanner_get_results(scanner: *const zbar_image_scanner_t) -> *mut c_void;

    // From zbar.h
    pub fn zbar_image_scanner_get_config(
        scanner: *const zbar_image_scanner_t,
        symbology: c_int,
        config: c_int,
        value: *mut c_int,
    ) -> c_int;

    // From symbol.h
    pub fn _zbar_symbol_set_create() -> *mut zbar_symbol_set_t;
}
