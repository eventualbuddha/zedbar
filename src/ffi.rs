//! FFI bindings to the C ZBar library
//!
//! This module provides direct bindings to the original C implementation.
//! As modules are converted to Rust, these bindings will be gradually removed.

use libc::{c_char, c_int, c_uint, c_ulong, c_void};

use crate::img_scanner::zbar_symbol_set_t;

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

pub use crate::img_scanner::{
    zbar_image_scanner_create, zbar_image_scanner_destroy, zbar_image_scanner_set_config,
    zbar_scan_image,
};
pub use crate::line_scanner::zbar_scanner_create;
pub use crate::symbol::_zbar_symbol_set_create;
