use std::{ffi::c_void, mem::size_of, ptr::{copy_nonoverlapping, null_mut}};

use libc::{c_int, c_uint, c_ulong, calloc, free, malloc, size_t};

use crate::line_scanner::zbar_scanner_t;

const RECYCLE_BUCKETS: usize = 5;
const NUM_SCN_CFGS: usize = 2; // ZBAR_CFG_Y_DENSITY - ZBAR_CFG_X_DENSITY + 1
const NUM_SYMS: usize = 25; // Number of symbol types

// Cache timing constants (in milliseconds)
const CACHE_PROXIMITY: c_ulong = 1000;
const CACHE_HYSTERESIS: c_ulong = 2000;
const CACHE_TIMEOUT: c_ulong = CACHE_HYSTERESIS * 2;

// Orientation constant
const ZBAR_ORIENT_UNKNOWN: c_int = -1;

// QR Code finder precision constant
const QR_FINDER_SUBPREC: c_int = 2;

// QR_FIXED macro: ((((v) << 1) + (rnd)) << (QR_FINDER_SUBPREC - 1))
#[inline]
fn qr_fixed(v: c_int, rnd: c_int) -> c_uint {
    (((v as c_uint) << 1) + (rnd as c_uint)) << (QR_FINDER_SUBPREC - 1)
}

// Import types and functions from ffi module
use crate::ffi::{refcnt, zbar_image_t, zbar_symbol_t};

// Import functions and constants from symbol module
use crate::symbol::_zbar_get_symbol_hash;

// qr_point is typedef int qr_point[2] in C
type QrPoint = [c_int; 2];

// qr_finder_line structure from qrcode.h
#[repr(C)]
#[allow(non_camel_case_types)]
pub struct qr_finder_line {
    pub pos: QrPoint,
    pub len: c_int,
    pub boffs: c_int,
    pub eoffs: c_int,
}

// C assert function
extern "C" {
    fn __assert_fail(
        assertion: *const u8,
        file: *const u8,
        line: c_uint,
        function: *const u8,
    ) -> !;
}

// Helper macro to call C assert for compatibility
macro_rules! c_assert {
    ($cond:expr) => {
        if !$cond {
            unsafe {
                __assert_fail(
                    concat!(stringify!($cond), "\0").as_ptr(),
                    concat!(file!(), "\0").as_ptr(),
                    line!(),
                    concat!("", "\0").as_ptr(),
                )
            }
        }
    };
}

// External C functions
extern "C" {
    fn zbar_decoder_get_config(
        dcode: *mut zbar_decoder_t,
        sym: c_int,
        cfg: c_int,
        val: *mut c_int,
    ) -> c_int;
    fn _zbar_symbol_set_free(syms: *mut zbar_symbol_set_t);
    fn zbar_scanner_flush(scn: *mut zbar_scanner_t);
    fn zbar_scanner_new_scan(scn: *mut zbar_scanner_t);
    fn _zbar_decoder_get_qr_finder_line(dcode: *mut zbar_decoder_t) -> *mut qr_finder_line;
    fn _zbar_qr_found_line(qr: *mut qr_reader, direction: c_int, line: *const qr_finder_line) -> c_int;
}

// Import from line_scanner and symbol modules
use crate::line_scanner::zbar_scanner_get_edge;
use crate::symbol::_zbar_symbol_free;

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
    if !(ZBAR_PARTIAL..=ZBAR_CODE128).contains(&sym) || sym == ZBAR_COMPOSITE {
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
            c_assert!((*sym).data_alloc != 0);
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
                    c_assert!(false);
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
                c_assert!(!(*sym).data.is_null());
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

/// Flush scanner pipeline and start new scan
///
/// This function flushes the scanner pipeline twice and then starts a new scan.
/// It's typically called at quiet borders to reset the scanner state.
///
/// # Arguments
/// * `iscn` - The image scanner instance
#[no_mangle]
pub unsafe extern "C" fn _zbar_image_scanner_quiet_border(
    iscn: *mut crate::ffi::zbar_image_scanner_t,
) {
    // Cast to the local type for field access
    let iscn = iscn as *mut zbar_image_scanner_t;
    let scn = (*iscn).scn;

    // Flush scanner pipeline twice
    zbar_scanner_flush(scn);
    zbar_scanner_flush(scn);

    // Start new scan
    zbar_scanner_new_scan(scn);
}

/// Recycle a symbol set
///
/// Decrements the reference count and recycles the symbols if the count reaches zero.
///
/// # Arguments
/// * `iscn` - The image scanner instance
/// * `syms` - The symbol set to recycle
///
/// # Returns
/// 1 if the symbol set is still referenced, 0 if it was recycled
#[no_mangle]
pub unsafe extern "C" fn _zbar_image_scanner_recycle_symbol_set(
    iscn: *mut crate::ffi::zbar_image_scanner_t,
    syms: *mut zbar_symbol_set_t,
) -> c_int {
    // Cast to the local type for field access
    let iscn = iscn as *mut zbar_image_scanner_t;

    if refcnt(&mut (*syms).refcnt, -1) != 0 {
        return 1;
    }

    _zbar_image_scanner_recycle_syms(iscn as *mut crate::ffi::zbar_image_scanner_t, (*syms).head);
    (*syms).head = null_mut();
    (*syms).tail = null_mut();
    (*syms).nsyms = 0;
    0
}

/// Recycle image symbols
///
/// This function recycles symbols from both the scanner's current symbol set
/// and the image's symbol set, managing their lifecycle appropriately.
///
/// # Arguments
/// * `iscn` - The image scanner instance
/// * `img` - The image whose symbols should be recycled
#[no_mangle]
pub unsafe extern "C" fn zbar_image_scanner_recycle_image(
    iscn: *mut crate::ffi::zbar_image_scanner_t,
    img: *mut crate::ffi::zbar_image_t,
) {
    // Cast to the local type for field access
    let iscn = iscn as *mut zbar_image_scanner_t;

    let mut syms = (*iscn).syms;
    if !syms.is_null() && (*syms).refcnt != 0 && _zbar_image_scanner_recycle_symbol_set(iscn as *mut crate::ffi::zbar_image_scanner_t, syms) != 0 {
        (*iscn).syms = null_mut();
    }

    syms = (*img).syms as *mut zbar_symbol_set_t;
    (*img).syms = null_mut();

    if !syms.is_null() && _zbar_image_scanner_recycle_symbol_set(iscn as *mut crate::ffi::zbar_image_scanner_t, syms) != 0 {
        // Symbol set is still referenced
    } else if !syms.is_null() {
        // Select one set to resurrect, destroy the other
        if !(*iscn).syms.is_null() {
            _zbar_symbol_set_free(syms);
        } else {
            (*iscn).syms = syms;
        }
    }
}

/// Allocate a symbol from the recycling pool or allocate a new one
///
/// Attempts to recycle a symbol from the appropriate size bucket,
/// or allocates a new one if no suitable recycled symbol is available.
///
/// # Arguments
/// * `iscn` - The image scanner instance
/// * `sym_type` - The type of symbol to allocate
/// * `datalen` - The length of data to allocate (including null terminator)
///
/// # Returns
/// Pointer to the allocated symbol
#[no_mangle]
pub unsafe extern "C" fn _zbar_image_scanner_alloc_sym(
    iscn: *mut crate::ffi::zbar_image_scanner_t,
    sym_type: c_int,
    datalen: c_int,
) -> *mut zbar_symbol_t {
    // Cast to the local type for field access
    let iscn = iscn as *mut zbar_image_scanner_t;

    // Recycle old or alloc new symbol
    let mut sym: *mut zbar_symbol_t = null_mut();

    // Find appropriate bucket based on datalen
    let mut i: i32 = 0;
    while i < (RECYCLE_BUCKETS - 1) as i32 {
        if datalen <= (1 << (i * 2)) {
            break;
        }
        i += 1;
    }

    // Try to get a symbol from this bucket or larger buckets
    let mut bucket_idx: Option<usize> = None;
    while i >= 0 {
        sym = (*iscn).recycle[i as usize].head;
        if !sym.is_null() {
            bucket_idx = Some(i as usize);
            break;
        }
        i -= 1;
    }

    if let Some(idx) = bucket_idx {
        // Found a recycled symbol
        (*iscn).recycle[idx].head = (*sym).next;
        (*sym).next = null_mut();
        c_assert!((*iscn).recycle[idx].nsyms > 0);
        (*iscn).recycle[idx].nsyms -= 1;
    } else {
        // Allocate a new symbol
        sym = calloc(1, size_of::<zbar_symbol_t>()) as *mut zbar_symbol_t;
    }

    // Initialize the symbol
    (*sym).symbol_type = sym_type;
    (*sym).quality = 1;
    (*sym).npts = 0;
    (*sym).orient = ZBAR_ORIENT_UNKNOWN;
    (*sym).cache_count = 0;
    (*sym).time = (*iscn).time;
    c_assert!((*sym).syms.is_null());

    if datalen > 0 {
        (*sym).datalen = (datalen - 1) as c_uint;
        if (*sym).data_alloc < datalen as c_uint {
            if !(*sym).data.is_null() {
                free((*sym).data as *mut c_void);
            }
            (*sym).data_alloc = datalen as c_uint;
            (*sym).data = malloc(datalen as size_t) as *mut i8;
        }
    } else {
        if !(*sym).data.is_null() {
            free((*sym).data as *mut c_void);
        }
        (*sym).data = null_mut();
        (*sym).datalen = 0;
        (*sym).data_alloc = 0;
    }

    sym
}

/// Handle QR code finder line detection
///
/// Processes a QR code finder line from the decoder and forwards it to the QR reader.
/// Adjusts edge positions based on scanner state and transforms coordinates.
///
/// # Arguments
/// * `iscn` - The image scanner instance
#[no_mangle]
pub unsafe extern "C" fn _zbar_image_scanner_qr_handler(
    iscn: *mut crate::ffi::zbar_image_scanner_t,
) {
    // Cast to the local type for field access
    let iscn = iscn as *mut zbar_image_scanner_t;

    let line = _zbar_decoder_get_qr_finder_line((*iscn).dcode);
    c_assert!(!line.is_null());

    let mut u = zbar_scanner_get_edge((*iscn).scn, (*line).pos[0] as c_uint, QR_FINDER_SUBPREC);
    (*line).boffs = (u as c_int) - zbar_scanner_get_edge((*iscn).scn, (*line).boffs as c_uint, QR_FINDER_SUBPREC) as c_int;
    (*line).len = zbar_scanner_get_edge((*iscn).scn, (*line).len as c_uint, QR_FINDER_SUBPREC) as c_int;
    (*line).eoffs = zbar_scanner_get_edge((*iscn).scn, (*line).eoffs as c_uint, QR_FINDER_SUBPREC) as c_int - (*line).len;
    (*line).len -= u as c_int;

    u = (qr_fixed((*iscn).umin, 0) as i64 + ((*iscn).du as i64) * (u as i64)) as c_uint;
    if (*iscn).du < 0 {
        std::mem::swap(&mut (*line).boffs, &mut (*line).eoffs);
        u = u.wrapping_sub((*line).len as c_uint);
    }

    let vert: c_int = if (*iscn).dx != 0 { 0 } else { 1 };
    (*line).pos[vert as usize] = u as c_int;
    (*line).pos[(1 - vert) as usize] = qr_fixed((*iscn).v, 1) as c_int;

    _zbar_qr_found_line((*iscn).qr, vert, line);
}

/// Look up a symbol in the cache
///
/// Searches for a matching symbol in the cache and recycles stale entries.
///
/// # Arguments
/// * `iscn` - The image scanner instance
/// * `sym` - The symbol to look up
///
/// # Returns
/// Pointer to the matching cache entry, or NULL if not found
#[no_mangle]
pub unsafe extern "C" fn _zbar_image_scanner_cache_lookup(
    iscn: *mut crate::ffi::zbar_image_scanner_t,
    sym: *mut zbar_symbol_t,
) -> *mut zbar_symbol_t {
    // Cast to the local type for field access
    let iscn = iscn as *mut zbar_image_scanner_t;

    let mut entry = &mut (*iscn).cache as *mut *mut zbar_symbol_t;

    while !(*entry).is_null() {
        // Check if this entry matches
        if (*(*entry)).symbol_type == (*sym).symbol_type
            && (*(*entry)).datalen == (*sym).datalen
        {
            // Compare data
            let entry_data = std::slice::from_raw_parts(
                (*(*entry)).data as *const u8,
                (*sym).datalen as usize,
            );
            let sym_data =
                std::slice::from_raw_parts((*sym).data as *const u8, (*sym).datalen as usize);

            if entry_data == sym_data {
                break;
            }
        }

        // Check if entry is stale
        if (*sym).time - (*(*entry)).time > CACHE_TIMEOUT {
            // Recycle stale cache entry
            let next = (*(*entry)).next;
            (*(*entry)).next = null_mut();
            _zbar_image_scanner_recycle_syms(iscn as *mut crate::ffi::zbar_image_scanner_t, *entry);
            *entry = next;
        } else {
            entry = &mut (*(*entry)).next as *mut *mut zbar_symbol_t;
        }
    }

    *entry
}

/// Cache a symbol
///
/// Updates or creates a cache entry for the symbol and sets its cache count.
///
/// # Arguments
/// * `iscn` - The image scanner instance
/// * `sym` - The symbol to cache
#[no_mangle]
pub unsafe extern "C" fn _zbar_image_scanner_cache_sym(
    iscn: *mut crate::ffi::zbar_image_scanner_t,
    sym: *mut zbar_symbol_t,
) {
    // Cast to the local type for field access
    let iscn = iscn as *mut zbar_image_scanner_t;

    if (*iscn).enable_cache != 0 {
        let mut entry = _zbar_image_scanner_cache_lookup(iscn as *mut crate::ffi::zbar_image_scanner_t, sym);

        if entry.is_null() {
            // FIXME reuse sym
            entry = _zbar_image_scanner_alloc_sym(
                iscn as *mut crate::ffi::zbar_image_scanner_t,
                (*sym).symbol_type,
                ((*sym).datalen + 1) as c_int,
            );
            (*entry).configs = (*sym).configs;
            (*entry).modifiers = (*sym).modifiers;
            copy_nonoverlapping(
                (*sym).data as *const u8,
                (*entry).data as *mut u8,
                (*sym).datalen as usize,
            );
            (*entry).time = (*sym).time - CACHE_HYSTERESIS;
            (*entry).cache_count = 0;

            // Add to cache
            (*entry).next = (*iscn).cache;
            (*iscn).cache = entry;
        }

        // Consistency check and hysteresis
        let age = (*sym).time - (*entry).time;
        (*entry).time = (*sym).time;
        let near_thresh = age < CACHE_PROXIMITY;
        let far_thresh = age >= CACHE_HYSTERESIS;
        let dup = (*entry).cache_count >= 0;

        if (!dup && !near_thresh) || far_thresh {
            let sym_type = (*sym).symbol_type;
            let h = _zbar_get_symbol_hash(sym_type);
            (*entry).cache_count = -(*iscn).sym_configs[0][h as usize];
        } else if dup || near_thresh {
            (*entry).cache_count += 1;
        }

        (*sym).cache_count = (*entry).cache_count;
    } else {
        (*sym).cache_count = 0;
    }
}

/// Add a symbol to the scanner's symbol set
///
/// Caches the symbol and adds it to the current symbol set.
///
/// # Arguments
/// * `iscn` - The image scanner instance
/// * `sym` - The symbol to add
#[no_mangle]
pub unsafe extern "C" fn _zbar_image_scanner_add_sym(
    iscn: *mut crate::ffi::zbar_image_scanner_t,
    sym: *mut zbar_symbol_t,
) {
    // Cast to the local type for field access
    let iscn = iscn as *mut zbar_image_scanner_t;

    _zbar_image_scanner_cache_sym(iscn as *mut crate::ffi::zbar_image_scanner_t, sym);

    let syms = (*iscn).syms;

    if (*sym).cache_count != 0 || (*syms).tail.is_null() {
        (*sym).next = (*syms).head;
        (*syms).head = sym;
    } else {
        (*sym).next = (*(*syms).tail).next;
        (*(*syms).tail).next = sym;
    }

    if (*sym).cache_count == 0 {
        (*syms).nsyms += 1;
    } else if (*syms).tail.is_null() {
        (*syms).tail = sym;
    }

    // Increment reference count
    // Note: The condition `&& 1 <= 0` in the original C code is always false,
    // so _zbar_symbol_free is never called
    if refcnt(&mut (*sym).refcnt, 1) == 0 && 1 <= 0 {
        _zbar_symbol_free(sym);
    }
}
