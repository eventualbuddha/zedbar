use std::{ffi::c_void, ptr::null_mut};

use libc::{c_int, c_uint, c_ulong, free};

use crate::{line_scanner::zbar_scanner_t, zbar_image_t};

const RECYCLE_BUCKETS: usize = 5;
const NUM_SCN_CFGS: usize = 2; // ZBAR_CFG_Y_DENSITY - ZBAR_CFG_X_DENSITY + 1
const NUM_SYMS: usize = 25; // Number of symbol types

// Import types and functions from ffi module
use crate::ffi::{refcnt, zbar_symbol_t};

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
    fn _zbar_symbol_set_free(syms: *mut zbar_symbol_set_t);
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
    pub refcnt: c_int,
    pub nsyms: c_int,
    pub head: *mut zbar_symbol_t,
    pub tail: *mut zbar_symbol_t,
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
        _zbar_image_scanner_recycle_syms(iscn as *mut crate::ffi::zbar_image_scanner_t, (*iscn).cache);
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

/// Recycle symbols from the image scanner
///
/// Recursively processes a linked list of symbols, either unlinking referenced
/// symbols or recycling unreferenced ones into the appropriate size bucket.
///
/// # Arguments
/// * `iscn` - The image scanner instance
/// * `sym` - Head of the symbol list to recycle
#[no_mangle]
pub unsafe extern "C" fn _zbar_image_scanner_recycle_syms(
    iscn: *mut crate::ffi::zbar_image_scanner_t,
    mut sym: *mut zbar_symbol_t,
) {
    // Cast to the local type for field access
    let iscn = iscn as *mut zbar_image_scanner_t;
    while !sym.is_null() {
        let next = (*sym).next;

        if (*sym).refcnt != 0 && refcnt(&mut (*sym).refcnt, -1) != 0 {
            // Unlink referenced symbol
            // FIXME handle outstanding component refs (currently unsupported)
            debug_assert!((*sym).data_alloc != 0);
            (*sym).next = null_mut();
        } else {
            // Recycle unreferenced symbol
            if (*sym).data_alloc == 0 {
                (*sym).data = null_mut();
                (*sym).datalen = 0;
            }

            if !(*sym).syms.is_null() {
                let syms = (*sym).syms as *mut zbar_symbol_set_t;
                if refcnt(&mut (*syms).refcnt, -1) != 0 {
                    debug_assert!(false);
                }
                _zbar_image_scanner_recycle_syms(iscn as *mut crate::ffi::zbar_image_scanner_t, (*syms).head);
                (*syms).head = null_mut();
                _zbar_symbol_set_free(syms);
                (*sym).syms = null_mut();
            }

            // Find appropriate bucket based on data allocation size
            let mut i = 0;
            while i < RECYCLE_BUCKETS {
                if ((*sym).data_alloc as c_int) < (1 << (i * 2)) {
                    break;
                }
                i += 1;
            }

            if i == RECYCLE_BUCKETS {
                debug_assert!(!(*sym).data.is_null());
                free((*sym).data as *mut c_void);
                (*sym).data = null_mut();
                (*sym).data_alloc = 0;
                i = 0;
            }

            let bucket = &mut (*iscn).recycle[i];
            // FIXME cap bucket fill
            bucket.nsyms += 1;
            (*sym).next = bucket.head;
            bucket.head = sym;
        }

        sym = next;
    }
}
