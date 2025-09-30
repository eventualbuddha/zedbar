//! FFI bindings to the C ZBar library
//! 
//! This module provides direct bindings to the original C implementation.
//! As modules are converted to Rust, these bindings will be gradually removed.

use libc::{c_char, c_int, c_uint, c_void, c_ulong};

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
    pub fn zbar_image_create() -> *mut c_void;
    pub fn zbar_image_destroy(image: *mut c_void);
    pub fn zbar_image_ref(image: *mut c_void, refs: c_int);
    
    pub fn zbar_image_get_format(image: *const c_void) -> c_ulong;
    pub fn zbar_image_get_width(image: *const c_void) -> c_uint;
    pub fn zbar_image_get_height(image: *const c_void) -> c_uint;
    pub fn zbar_image_get_data(image: *const c_void) -> *const c_void;
    
    pub fn zbar_image_set_format(image: *mut c_void, format: c_ulong);
    pub fn zbar_image_set_size(image: *mut c_void, width: c_uint, height: c_uint);
    pub fn zbar_image_set_data(image: *mut c_void, data: *const c_void, data_len: c_ulong, cleanup: *const c_void);
    
    // Scanner functions
    pub fn zbar_image_scanner_create() -> *mut c_void;
    pub fn zbar_image_scanner_destroy(scanner: *mut c_void);
    pub fn zbar_image_scanner_set_config(scanner: *mut c_void, symbology: c_int, config: c_int, value: c_int) -> c_int;
    pub fn zbar_scan_image(scanner: *mut c_void, image: *mut c_void) -> c_int;
    pub fn zbar_image_first_symbol(image: *const c_void) -> *const c_void;
    
    // Symbol functions
    pub fn zbar_symbol_get_type(symbol: *const c_void) -> c_int;
    pub fn zbar_symbol_get_data(symbol: *const c_void) -> *const c_char;
    pub fn zbar_symbol_get_data_length(symbol: *const c_void) -> c_uint;
    pub fn zbar_symbol_next(symbol: *const c_void) -> *const c_void;
}