use std::{ffi::c_void, mem::transmute, ptr::null_mut};

use libc::{c_int, c_uint, c_ulong};

use crate::{line_scanner::zbar_scanner_t, zbar_image_t};

const RECYCLE_BUCKETS: usize = 5;
const NUM_SCN_CFGS: usize = 2; // ZBAR_CFG_Y_DENSITY - ZBAR_CFG_X_DENSITY + 1
const NUM_SYMS: usize = 25; // Number of symbol types

// Import external C function from ffi module
use crate::_zbar_image_scanner_recycle_syms;

// Import functions and constants from symbol module
use crate::symbol::_zbar_get_symbol_hash;

// External C functions
extern "C" {
    fn zbar_decoder_get_config(
        dcode: *mut zbar_decoder_t,
        sym: c_int,
        cfg: c_int,
        val: *mut c_int,
    ) -> c_int;
}

// Config constants (from zbar.h)
const ZBAR_CFG_UNCERTAINTY: c_int = 64;
const ZBAR_CFG_POSITION: c_int = 128;
const ZBAR_CFG_X_DENSITY: c_int = 256;
const ZBAR_CFG_Y_DENSITY: c_int = 257;

// Symbol type constants (from zbar.h)
const ZBAR_PARTIAL: c_int = 1;
const ZBAR_CODE128: c_int = 128;
const ZBAR_COMPOSITE: c_int = 15;

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

#[repr(C)]
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

/// Set the data handler callback for the image scanner
///
/// This function sets a callback that will be invoked when symbols are decoded.
/// Returns the previous handler (or NULL if none was set).
///
/// # Arguments
/// * `iscn` - The image scanner instance
/// * `handler` - The new callback handler (or NULL to disable)
/// * `userdata` - User data pointer to pass to the handler
#[no_mangle]
pub unsafe extern "C" fn zbar_image_scanner_set_data_handler(
    iscn: *mut zbar_image_scanner_t,
    handler: *mut zbar_image_data_handler_t,
    userdata: *const c_void,
) -> *mut zbar_image_data_handler_t {
    if iscn.is_null() {
        return null_mut();
    }

    let result = (*iscn).handler;
    (*iscn).handler = handler;
    (*iscn).userdata = userdata;
    result
}

/// Enable or disable result caching for the image scanner
///
/// When enabled, the scanner caches decoded results to suppress duplicates
/// across consecutive frames. When disabled or when this function is called,
/// all cached symbols are recycled.
///
/// # Arguments
/// * `iscn` - The image scanner instance
/// * `enable` - Non-zero to enable caching, 0 to disable
#[no_mangle]
pub unsafe extern "C" fn zbar_image_scanner_enable_cache(
    iscn: *mut zbar_image_scanner_t,
    enable: c_int,
) {
    if iscn.is_null() {
        return;
    }

    if !(*iscn).cache.is_null() {
        // Recycle all cached symbols
        // Cast to the FFI types expected by the C function
        _zbar_image_scanner_recycle_syms(
            transmute(iscn),
            transmute((*iscn).cache),
        );
        (*iscn).cache = null_mut();
    }
    (*iscn).enable_cache = if enable != 0 { 1 } else { 0 };
}

/// Get configuration value for a specific symbology
///
/// Retrieves the current configuration value for a particular setting
/// of a barcode symbology.
///
/// # Arguments
/// * `iscn` - The image scanner instance
/// * `sym` - The symbology type to query
/// * `cfg` - The configuration parameter to query
/// * `val` - Pointer to store the retrieved value
///
/// # Returns
/// 0 on success, 1 on error
#[no_mangle]
pub unsafe extern "C" fn zbar_image_scanner_get_config(
    iscn: *mut zbar_image_scanner_t,
    sym: c_int,
    cfg: c_int,
    val: *mut c_int,
) -> c_int {
    if iscn.is_null() || val.is_null() {
        return 1;
    }

    // Return error if symbol doesn't have config
    if sym < ZBAR_PARTIAL || sym > ZBAR_CODE128 || sym == ZBAR_COMPOSITE {
        return 1;
    }

    if cfg < ZBAR_CFG_UNCERTAINTY {
        return zbar_decoder_get_config((*iscn).dcode, sym, cfg, val);
    }

    if cfg < ZBAR_CFG_POSITION {
        if sym == ZBAR_PARTIAL {
            return 1;
        }

        let i = _zbar_get_symbol_hash(sym);
        *val = (*iscn).sym_configs[(cfg - ZBAR_CFG_UNCERTAINTY) as usize][i as usize];
        return 0;
    }

    // Image scanner parameters apply only to ZBAR_PARTIAL
    if sym > ZBAR_PARTIAL {
        return 1;
    }

    if cfg < ZBAR_CFG_X_DENSITY {
        *val = if ((*iscn).config & (1 << (cfg - ZBAR_CFG_POSITION))) != 0 {
            1
        } else {
            0
        };
        return 0;
    }

    if cfg <= ZBAR_CFG_Y_DENSITY {
        // CFG macro: ((iscn)->configs[(cfg) - ZBAR_CFG_X_DENSITY])
        *val = (*iscn).configs[(cfg - ZBAR_CFG_X_DENSITY) as usize];
        return 0;
    }

    1
}
