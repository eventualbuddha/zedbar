use std::{
    ffi::c_void,
    ptr::{copy_nonoverlapping, null_mut, NonNull},
};

use libc::{c_int, c_uint, memcmp, size_t};

use crate::{
    decoder::zbar_decoder_get_config,
    decoder_types::{
        zbar_decoder_t, ZBAR_CFG_BINARY, ZBAR_CFG_ENABLE, ZBAR_CFG_POSITION,
        ZBAR_CFG_TEST_INVERTED, ZBAR_CFG_UNCERTAINTY, ZBAR_CFG_X_DENSITY, ZBAR_CFG_Y_DENSITY,
        ZBAR_CODABAR, ZBAR_CODE128, ZBAR_CODE39, ZBAR_CODE93, ZBAR_COMPOSITE, ZBAR_DATABAR,
        ZBAR_DATABAR_EXP, ZBAR_EAN5, ZBAR_ISBN10, ZBAR_ORIENT_UNKNOWN, ZBAR_PARTIAL, ZBAR_QRCODE,
    },
    finder::{decoder_get_qr_finder_line, decoder_get_sq_finder_config},
    image_ffi::zbar_image_t,
    line_scanner::{
        scan_y, scanner_flush, scanner_get_edge, scanner_get_width, scanner_new_scan,
        zbar_scanner_new, zbar_scanner_t,
    },
    qrcode::{
        qr_point,
        qrdec::{
            _zbar_qr_create, _zbar_qr_decode, _zbar_qr_found_line, _zbar_qr_reset,
            qr_finder_lines,
        },
        rs::rs_gf256,
        IsaacCtx,
    },
    sqcode::{SqReader, _zbar_sq_decode},
    symbol::{
        _zbar_get_symbol_hash, _zbar_symbol_add_point, symbol_alloc_zeroed, symbol_clear_data,
        symbol_free, symbol_refcnt, symbol_reserve_data, symbol_set_create, symbol_set_free,
        zbar_symbol_set_ref, zbar_symbol_t,
    },
};

const RECYCLE_BUCKETS: usize = 5;
const NUM_SCN_CFGS: usize = 2; // ZBAR_CFG_Y_DENSITY - ZBAR_CFG_X_DENSITY + 1
const NUM_SYMS: usize = 25; // Number of symbol types

// QR Code finder precision constant
const QR_FINDER_SUBPREC: c_int = 2;

// QR_FIXED macro: ((((v) << 1) + (rnd)) << (QR_FINDER_SUBPREC - 1))
#[inline]
fn qr_fixed(v: c_int, rnd: c_int) -> c_uint {
    (((v as c_uint) << 1) + (rnd as c_uint)) << (QR_FINDER_SUBPREC - 1)
}

// Import types and functions from ffi module
use crate::ffi::refcnt;

// Import functions and constants from symbol module
pub struct qr_finder_line {
    pub pos: qr_point,
    pub len: c_int,
    pub boffs: c_int,
    pub eoffs: c_int,
}

// C assert function
extern "C" {
    fn __assert_fail(assertion: *const u8, file: *const u8, line: c_uint, function: *const u8)
        -> !;
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

// Import from line_scanner, decoder, and symbol modules
use crate::decoder::{
    zbar_decoder_get_configs, zbar_decoder_get_modifiers, zbar_decoder_get_type,
    zbar_decoder_get_userdata, zbar_decoder_set_config, zbar_decoder_set_handler,
    zbar_decoder_set_userdata,
};
use crate::image_ffi::{_zbar_image_copy, _zbar_image_swap_symbols, zbar_image_destroy};

// Helper macros for configuration access
macro_rules! CFG {
    ($iscn:expr, $cfg:expr) => {
        (*$iscn).configs[($cfg - ZBAR_CFG_X_DENSITY) as usize]
    };
}

macro_rules! TEST_CFG {
    ($iscn:expr, $cfg:expr) => {
        (((*$iscn).config >> (($cfg - ZBAR_CFG_POSITION) as u32)) & 1) != 0
    };
}

// FourCC code for image formats
#[inline]
const fn fourcc(a: u8, b: u8, c: u8, d: u8) -> u32 {
    (a as u32) | ((b as u32) << 8) | ((c as u32) << 16) | ((d as u32) << 24)
}

#[allow(non_camel_case_types)]
pub struct qr_reader {
    /// The GF(256) representation used in Reed-Solomon decoding.
    pub gf: rs_gf256,
    /// The random number generator used by RANSAC.
    pub isaac: IsaacCtx,
    ///  current finder state, horizontal and vertical lines
    pub finder_lines: [qr_finder_lines; 2],
}

#[derive(Default)]
#[allow(non_camel_case_types)]
pub struct zbar_symbol_set_t {
    pub refcnt: c_int,
    pub nsyms: c_int,
    pub head: *mut zbar_symbol_t,
    pub tail: *mut zbar_symbol_t,
}

impl Drop for zbar_symbol_set_t {
    fn drop(&mut self) {
        let mut sym = self.head;
        while !sym.is_null() {
            let next = unsafe { (*sym).next };
            unsafe {
                (*sym).next = std::ptr::null_mut();
            }
            unsafe {
                symbol_refcnt(sym, -1);
            }
            sym = next;
        }
        self.head = std::ptr::null_mut();
    }
}

// Function pointer type for image data handler callbacks
pub type zbar_image_data_handler_t = unsafe fn(img: *mut zbar_image_t, userdata: *const c_void);

#[derive(Default)]
#[allow(non_camel_case_types)]
pub struct recycle_bucket_t {
    nsyms: c_int,
    head: *mut zbar_symbol_t,
}

/// image scanner state
#[derive(Default)]
#[allow(non_camel_case_types)]
pub struct zbar_image_scanner_t {
    /// associated linear intensity scanner
    scn: Option<zbar_scanner_t>,

    /// associated symbol decoder
    dcode: Option<zbar_decoder_t>,
    /// QR Code 2D reader
    qr: Option<qr_reader>,
    /// SQ Code 2D reader
    sq: Option<SqReader>,

    /// application data
    userdata: *const c_void,
    /// user result callback
    handler: Option<zbar_image_data_handler_t>,

    /// current scan direction
    dx: c_int,
    dy: c_int,
    du: c_int,
    umin: c_int,
    v: c_int,

    /// previous decode results
    syms: Option<NonNull<zbar_symbol_set_t>>,
    /// recycled symbols in 4^n size buckets
    recycle: [recycle_bucket_t; RECYCLE_BUCKETS],

    // configuration settings
    /// config flags
    config: c_uint,
    ean_config: c_uint,
    /// int valued configurations
    configs: [c_int; NUM_SCN_CFGS],
    /// per-symbology configurations
    sym_configs: [[c_int; NUM_SYMS]; 1],
}

impl zbar_image_scanner_t {
    #[inline]
    pub fn syms_ptr(&self) -> *mut zbar_symbol_set_t {
        self.syms.map_or(null_mut(), NonNull::as_ptr)
    }

    #[inline]
    pub fn set_syms_ptr(&mut self, ptr: *mut zbar_symbol_set_t) {
        self.syms = NonNull::new(ptr);
    }

    #[inline]
    pub fn clear_syms(&mut self) {
        self.syms = None;
    }

    #[inline]
    pub fn take_syms(&mut self) -> Option<NonNull<zbar_symbol_set_t>> {
        self.syms.take()
    }

    // Accessor methods for pointer fields

    #[inline]
    pub(crate) fn scanner(&self) -> Option<&zbar_scanner_t> {
        self.scn.as_ref()
    }

    #[inline]
    pub(crate) fn decoder(&self) -> Option<&zbar_decoder_t> {
        self.dcode.as_ref()
    }

    #[inline]
    pub(crate) fn decoder_mut(&mut self) -> Option<&mut zbar_decoder_t> {
        self.dcode.as_mut()
    }

    #[inline]
    pub(crate) fn set_scanner(&mut self, scn: Option<zbar_scanner_t>) {
        self.scn = scn;
    }

    #[inline]
    pub(crate) fn set_decoder(&mut self, dcode: Option<zbar_decoder_t>) {
        self.dcode = dcode;
    }

    #[inline]
    pub(crate) fn set_qr_reader(&mut self, qr: Option<qr_reader>) {
        self.qr = qr;
    }

    #[inline]
    pub(crate) fn set_sq_reader(&mut self, sq: Option<SqReader>) {
        self.sq = sq;
    }
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
pub unsafe fn zbar_image_scanner_set_data_handler(
    iscn: *mut zbar_image_scanner_t,
    handler: Option<zbar_image_data_handler_t>,
    userdata: *const c_void,
) -> Option<zbar_image_data_handler_t> {
    if iscn.is_null() {
        return None;
    }
    let iscn = &mut *iscn;

    let result = iscn.handler;
    iscn.handler = handler;
    iscn.userdata = userdata;
    result
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
pub unsafe fn zbar_image_scanner_get_config(
    iscn: *mut zbar_image_scanner_t,
    sym: c_int,
    cfg: c_int,
    val: *mut c_int,
) -> c_int {
    if iscn.is_null() || val.is_null() {
        return 1;
    }
    let iscn = &mut *iscn;

    // Return error if symbol doesn't have config
    if !(ZBAR_PARTIAL..=ZBAR_CODE128).contains(&sym) || sym == ZBAR_COMPOSITE {
        return 1;
    }

    if cfg < ZBAR_CFG_UNCERTAINTY {
        let dcode_ptr = iscn.dcode.as_mut().map_or(null_mut(), |d| d as *mut _);
        return zbar_decoder_get_config(dcode_ptr, sym, cfg, val);
    }

    if cfg < ZBAR_CFG_POSITION {
        if sym == ZBAR_PARTIAL {
            return 1;
        }

        let i = _zbar_get_symbol_hash(sym);
        *val = iscn.sym_configs[(cfg - ZBAR_CFG_UNCERTAINTY) as usize][i as usize];
        return 0;
    }

    // Image scanner parameters apply only to ZBAR_PARTIAL
    if sym > ZBAR_PARTIAL {
        return 1;
    }

    if cfg < ZBAR_CFG_X_DENSITY {
        *val = if (iscn.config & (1 << (cfg - ZBAR_CFG_POSITION))) != 0 {
            1
        } else {
            0
        };
        return 0;
    }

    if cfg <= ZBAR_CFG_Y_DENSITY {
        // CFG macro: ((iscn)->configs[(cfg) - ZBAR_CFG_X_DENSITY])
        *val = iscn.configs[(cfg - ZBAR_CFG_X_DENSITY) as usize];
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
pub unsafe fn _zbar_image_scanner_recycle_syms(
    iscn: *mut zbar_image_scanner_t,
    mut sym: *mut zbar_symbol_t,
) {
    if iscn.is_null() {
        return;
    }
    let iscn = &mut *iscn;

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
                symbol_clear_data(sym);
            }

            if !(*sym).syms.is_null() {
                let syms = (*sym).syms;
                if refcnt(&mut (*syms).refcnt, -1) != 0 {
                    c_assert!(false);
                }
                _zbar_image_scanner_recycle_syms(iscn as *mut _, (*syms).head);
                (*syms).head = null_mut();
                {
                    symbol_set_free(syms);
                };
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
                symbol_clear_data(sym);
                i = 0;
            }

            let bucket = &mut iscn.recycle[i];
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
pub unsafe fn _zbar_image_scanner_quiet_border(iscn: *mut zbar_image_scanner_t) {
    if iscn.is_null() {
        return;
    }
    let iscn = &mut *iscn;
    let Some(scn) = iscn.scn.as_mut() else {
        return;
    };

    // Flush scanner pipeline twice
    scanner_flush(scn);
    scanner_flush(scn);

    // Start new scan
    scanner_new_scan(scn);
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
pub(crate) unsafe fn _zbar_image_scanner_alloc_sym(
    iscn: *mut zbar_image_scanner_t,
    sym_type: c_int,
    datalen: c_int,
) -> *mut zbar_symbol_t {
    if iscn.is_null() {
        return null_mut();
    }
    let iscn = &mut *iscn;

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
        sym = iscn.recycle[i as usize].head;
        if !sym.is_null() {
            bucket_idx = Some(i as usize);
            break;
        }
        i -= 1;
    }

    if let Some(idx) = bucket_idx {
        // Found a recycled symbol
        let sym_ref = &mut *sym;
        iscn.recycle[idx].head = sym_ref.next;
        sym_ref.next = null_mut();
        c_assert!(iscn.recycle[idx].nsyms > 0);
        iscn.recycle[idx].nsyms -= 1;
    } else {
        // Allocate a new symbol
        sym = symbol_alloc_zeroed();
    }

    // Initialize the symbol
    let sym_ref = &mut *sym;
    sym_ref.symbol_type = sym_type;
    sym_ref.quality = 1;
    sym_ref.npts = 0;
    sym_ref.orient = ZBAR_ORIENT_UNKNOWN;
    c_assert!(sym_ref.syms.is_null());

    if datalen > 0 {
        if symbol_reserve_data(sym, datalen as usize) {
            sym_ref.datalen = (datalen - 1) as c_uint;
        } else {
            symbol_clear_data(sym);
        }
    } else {
        symbol_clear_data(sym);
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
pub(crate) unsafe fn _zbar_image_scanner_qr_handler(iscn: *mut zbar_image_scanner_t) {
    let iscn = &mut *iscn;
    let Some(dcode) = &mut iscn.dcode else {
        return;
    };
    let line = decoder_get_qr_finder_line(dcode);

    let Some(scn) = iscn.scn.as_ref() else {
        return;
    };
    let mut u = scanner_get_edge(scn, line.pos[0] as c_uint, QR_FINDER_SUBPREC);
    line.boffs =
        (u as c_int) - scanner_get_edge(scn, line.boffs as c_uint, QR_FINDER_SUBPREC) as c_int;
    line.len = scanner_get_edge(scn, line.len as c_uint, QR_FINDER_SUBPREC) as c_int;
    line.eoffs = scanner_get_edge(scn, line.eoffs as c_uint, QR_FINDER_SUBPREC) as c_int - line.len;
    line.len -= u as c_int;

    u = (qr_fixed(iscn.umin, 0) as i64 + (iscn.du as i64) * (u as i64)) as c_uint;
    if iscn.du < 0 {
        std::mem::swap(&mut line.boffs, &mut line.eoffs);
        u = u.wrapping_sub(line.len as c_uint);
    }

    let vert: c_int = if iscn.dx != 0 { 0 } else { 1 };
    line.pos[vert as usize] = u as c_int;
    line.pos[(1 - vert) as usize] = qr_fixed(iscn.v, 1) as c_int;

    if let Some(qr) = &mut iscn.qr {
        _zbar_qr_found_line(qr, vert, line);
    }
}

/// Add a symbol to the scanner's symbol set
///
/// # Arguments
/// * `iscn` - The image scanner instance
/// * `sym` - The symbol to add
pub(crate) unsafe fn _zbar_image_scanner_add_sym(
    iscn: *mut zbar_image_scanner_t,
    sym: *mut zbar_symbol_t,
) {
    let iscn = &mut *iscn;
    let syms = iscn.syms_ptr();
    debug_assert!(!syms.is_null());
    let syms = &mut *syms;
    let sym_ref = &mut *sym;

    // The symbol set takes ownership of the symbol reference.
    // Ensure the reference count reflects this so that recycling
    // and Drop logic can safely release it later.
    let new_refcnt = refcnt(&mut sym_ref.refcnt, 1);
    debug_assert!(new_refcnt > 0);
    let _ = new_refcnt;

    if syms.tail.is_null() {
        sym_ref.next = syms.head;
        syms.head = sym;
    } else {
        let tail = &mut *syms.tail;
        sym_ref.next = tail.next;
        tail.next = sym;
    }

    syms.nsyms += 1;
}

/// SQ code handler - updates SQ reader configuration
///
/// Gets the current SQ finder configuration from the decoder and passes it
/// to the SQ reader for processing.
///
/// # Safety
/// `iscn` must be a valid pointer to a zbar_image_scanner_t
pub(crate) unsafe fn _zbar_image_scanner_sq_handler(iscn: *mut zbar_image_scanner_t) {
    let iscn = &mut *iscn;

    let Some(dcode) = &iscn.dcode else {
        return;
    };
    let config = decoder_get_sq_finder_config(dcode);
    if let Some(sq) = &mut iscn.sq {
        sq.set_enabled(config != 0);
    }
}

#[inline]
unsafe fn image_scanner_alloc_zeroed() -> *mut zbar_image_scanner_t {
    let scanner = Box::new(zbar_image_scanner_t::default());
    Box::into_raw(scanner)
}

#[inline]
unsafe fn image_scanner_free(iscn: *mut zbar_image_scanner_t) {
    drop(Box::from_raw(iscn))
}

/// Create a new image scanner
///
/// Allocates and initializes a new image scanner instance with default configuration.
///
/// # Returns
/// Pointer to new scanner or null on allocation failure
pub(crate) unsafe fn zbar_image_scanner_create() -> *mut zbar_image_scanner_t {
    let iscn = image_scanner_alloc_zeroed();
    if iscn.is_null() {
        return null_mut();
    }

    let Some(decoder) = zbar_decoder_t::new() else {
        return null_mut();
    };

    let iscn_ref = &mut *iscn;
    iscn_ref.set_decoder(Some(decoder));

    // Get a pointer to the decoder for the scanner
    let dcode_ptr = iscn_ref.dcode.as_mut().unwrap() as *mut zbar_decoder_t;
    let scanner = zbar_scanner_new(dcode_ptr as *mut _);
    iscn_ref.set_scanner(Some(scanner));

    if iscn_ref.decoder().is_none() || iscn_ref.scanner().is_none() {
        zbar_image_scanner_destroy(iscn);
        return null_mut();
    }

    if let Some(dcode) = iscn_ref.decoder_mut() {
        zbar_decoder_set_userdata(dcode as *mut _, iscn as *mut c_void);
        zbar_decoder_set_handler(dcode, Some(symbol_handler));
    }

    iscn_ref.set_qr_reader(Some(_zbar_qr_create()));
    iscn_ref.set_sq_reader(Some(SqReader::new()));

    // Apply default configuration
    iscn_ref.configs[0] = 1; // ZBAR_CFG_X_DENSITY
    iscn_ref.configs[1] = 1; // ZBAR_CFG_Y_DENSITY

    zbar_image_scanner_set_config(iscn, 0, ZBAR_CFG_POSITION, 1);
    zbar_image_scanner_set_config(iscn, 0, ZBAR_CFG_UNCERTAINTY, 2);
    zbar_image_scanner_set_config(iscn, 0, 65, 0); // ZBAR_CFG_TEST_INVERTED
    zbar_image_scanner_set_config(iscn, ZBAR_QRCODE, ZBAR_CFG_UNCERTAINTY, 0);
    zbar_image_scanner_set_config(iscn, ZBAR_QRCODE, ZBAR_CFG_BINARY, 0);
    zbar_image_scanner_set_config(iscn, ZBAR_CODE128, ZBAR_CFG_UNCERTAINTY, 0);
    zbar_image_scanner_set_config(iscn, ZBAR_CODE93, ZBAR_CFG_UNCERTAINTY, 0);
    zbar_image_scanner_set_config(iscn, ZBAR_CODE39, ZBAR_CFG_UNCERTAINTY, 0);
    zbar_image_scanner_set_config(iscn, ZBAR_CODABAR, ZBAR_CFG_UNCERTAINTY, 1);
    zbar_image_scanner_set_config(iscn, ZBAR_COMPOSITE, ZBAR_CFG_UNCERTAINTY, 0);

    iscn
}

/// Destroy an image scanner
///
/// Frees all resources associated with the scanner.
///
/// # Parameters
/// - `iscn`: Scanner to destroy (null-safe)
pub unsafe fn zbar_image_scanner_destroy(iscn: *mut zbar_image_scanner_t) {
    if iscn.is_null() {
        return;
    }

    let scanner = &mut *iscn;
    if let Some(syms_handle) = scanner.take_syms() {
        let syms_ptr = syms_handle.as_ptr();
        let syms = &*syms_ptr;
        if syms.refcnt != 0 {
            zbar_symbol_set_ref(syms_ptr, -1);
        } else {
            symbol_set_free(syms_ptr);
        }
    }

    // Scanner is now owned and will be dropped automatically

    // Free recycled symbols
    for i in 0..RECYCLE_BUCKETS {
        let mut sym = scanner.recycle[i].head;
        while !sym.is_null() {
            let sym_ref = &*sym;
            let next = sym_ref.next;
            symbol_free(sym);
            sym = next;
        }
    }

    // qr is owned, so it will be dropped automatically when scanner is dropped
    // No need to call _zbar_qr_destroy

    image_scanner_free(iscn);
}

/// Symbol handler callback for 1D barcode decoding
///
/// This function is called by the decoder when a barcode is successfully decoded.
/// It manages symbol deduplication, position tracking, and adds symbols to the result set.
///
/// # Arguments
/// * `dcode` - The decoder that found the symbol
pub unsafe fn symbol_handler(dcode: *mut zbar_decoder_t) {
    let iscn = zbar_decoder_get_userdata(dcode) as *mut zbar_image_scanner_t;
    let type_ = zbar_decoder_get_type(dcode);
    let mut x = 0;
    let mut y = 0;

    // QR codes are handled separately
    if type_ as c_int == ZBAR_QRCODE {
        _zbar_image_scanner_qr_handler(iscn);
        return;
    }

    // Calculate position if position tracking is enabled
    let iscn_ref = &*iscn;
    if TEST_CFG!(iscn, ZBAR_CFG_POSITION) {
        if let Some(scn) = iscn_ref.scanner() {
            let w = scanner_get_width(scn);
            let u = iscn_ref.umin + iscn_ref.du * scanner_get_edge(scn, w, 0) as c_int;
            if iscn_ref.dx != 0 {
                x = u;
                y = iscn_ref.v;
            } else {
                x = iscn_ref.v;
                y = u;
            }
        }
    }

    // Ignore partial results
    if (type_ as c_int) <= ZBAR_PARTIAL {
        return;
    }

    let data = (*dcode).buffer_slice();
    let datalen = data.len();

    // Check for duplicate symbols
    let syms = iscn_ref.syms_ptr();
    let mut sym = if !syms.is_null() {
        (*syms).head
    } else {
        null_mut()
    };
    while !sym.is_null() {
        let sym_ref = &mut *sym;
        if sym_ref.symbol_type == type_
            && sym_ref.datalen == datalen as c_uint
            && memcmp(
                sym_ref.data as *const c_void,
                data.as_ptr() as *const c_void,
                datalen as size_t,
            ) == 0
        {
            sym_ref.quality += 1;
            if TEST_CFG!(iscn, ZBAR_CFG_POSITION) {
                _zbar_symbol_add_point(sym, x, y);
            }
            return;
        }
        sym = sym_ref.next;
    }

    // Allocate new symbol
    sym = _zbar_image_scanner_alloc_sym(iscn, type_ as c_int, (datalen + 1) as c_int);
    let sym_ref = &mut *sym;
    sym_ref.configs = zbar_decoder_get_configs(dcode, type_);
    sym_ref.modifiers = zbar_decoder_get_modifiers(dcode);

    // Copy data
    copy_nonoverlapping(data.as_ptr(), sym_ref.data as *mut u8, datalen + 1);

    // Initialize position
    if TEST_CFG!(iscn, ZBAR_CFG_POSITION) {
        _zbar_symbol_add_point(sym, x, y);
    }

    // Set orientation
    let dir = (*dcode).direction;
    if dir != 0 {
        sym_ref.orient = (if iscn_ref.dy != 0 { 1 } else { 0 }) + ((iscn_ref.du ^ dir) & 2);
    }

    _zbar_image_scanner_add_sym(iscn, sym);
}

/// Set configuration for image scanner
///
/// Configures various settings for the scanner and its decoders.
///
/// # Arguments
/// * `iscn` - Image scanner instance
/// * `sym` - Symbol type (0 for all symbols, or specific type)
/// * `cfg` - Configuration parameter
/// * `val` - Configuration value
///
/// # Returns
/// 0 on success, 1 on error
pub(crate) unsafe fn zbar_image_scanner_set_config(
    iscn: *mut zbar_image_scanner_t,
    sym: c_int,
    cfg: c_int,
    val: c_int,
) -> c_int {
    let iscn_ref = &mut *iscn;

    // Handle EAN composite configuration
    if (sym == 0 || sym == ZBAR_COMPOSITE) && cfg == ZBAR_CFG_ENABLE {
        iscn_ref.ean_config = if val != 0 { 1 } else { 0 };
        if sym != 0 {
            return 0;
        }
    }

    // Delegate decoder configuration
    if cfg < ZBAR_CFG_UNCERTAINTY {
        if let Some(dcode) = iscn_ref.decoder_mut() {
            return zbar_decoder_set_config(dcode, sym, cfg, val);
        }
        return 1;
    }

    // Handle uncertainty and related configs
    if cfg < ZBAR_CFG_POSITION {
        if cfg > ZBAR_CFG_UNCERTAINTY {
            return 1;
        }
        let c = (cfg - ZBAR_CFG_UNCERTAINTY) as usize;
        if sym > ZBAR_PARTIAL {
            let i = _zbar_get_symbol_hash(sym) as usize;
            iscn_ref.sym_configs[c][i] = val;
        } else {
            for i in 0..NUM_SYMS {
                iscn_ref.sym_configs[c][i] = val;
            }
        }
        return 0;
    }

    // Image scanner parameters apply only to ZBAR_PARTIAL
    if sym > ZBAR_PARTIAL {
        return 1;
    }

    // Handle density configuration
    if (ZBAR_CFG_X_DENSITY..=ZBAR_CFG_Y_DENSITY).contains(&cfg) {
        iscn_ref.configs[(cfg - ZBAR_CFG_X_DENSITY) as usize] = val;
        return 0;
    }

    // Handle position and related configuration flags
    let cfg_bit = cfg - ZBAR_CFG_POSITION;

    if val == 0 {
        iscn_ref.config &= !(1 << cfg_bit);
    } else if val == 1 {
        iscn_ref.config |= 1 << cfg_bit;
    } else {
        return 1;
    }

    0
}

/// Internal image scanning implementation
///
/// Performs the actual barcode scanning on an image with horizontal and vertical passes.
///
/// # Arguments
/// * `iscn` - Image scanner instance
/// * `img` - Image to scan
///
/// # Returns
/// Pointer to symbol set on success, null on error
pub unsafe fn _zbar_scan_image(
    iscn: *mut zbar_image_scanner_t,
    img: *mut zbar_image_t,
) -> *mut zbar_symbol_set_t {
    // Reset QR and SQ decoders
    if let Some(qr) = &mut (*iscn).qr {
        _zbar_qr_reset(qr);
    }
    if let Some(sq) = &mut (*iscn).sq {
        sq.reset();
    }

    // Image must be in grayscale format
    if (*img).format != fourcc(b'Y', b'8', b'0', b'0')
        && (*img).format != fourcc(b'G', b'R', b'E', b'Y')
    {
        return null_mut();
    }

    // Recycle previous scanner and image results
    let scanner = &mut *iscn;
    let mut syms = scanner.syms_ptr();
    if syms.is_null() {
        syms = symbol_set_create();
        scanner.set_syms_ptr(syms);
        zbar_symbol_set_ref(syms, 1);
    } else {
        zbar_symbol_set_ref(syms, 2);
    }
    (*img).set_syms_ptr(syms);

    let w = (*img).width;
    let h = (*img).height;
    let data = (*img).data.as_ptr();

    let Some(scn) = (*iscn).scn.as_mut() else {
        return null_mut();
    };
    scanner_new_scan(scn);

    // Horizontal scanning pass
    let density = CFG!(iscn, ZBAR_CFG_Y_DENSITY);
    if density > 0 {
        let mut p = data;
        let mut x = 0i32;
        let mut y = 0i32;

        let mut border = ((h - 1) % (density as c_uint)).div_ceil(2);
        if border > h / 2 {
            border = h / 2;
        }
        c_assert!(border <= h);
        (*iscn).dy = 0;

        // movedelta(0, border)
        y += border as i32;
        p = p.offset((border * w) as isize);
        (*iscn).v = y;

        while (y as c_uint) < h {
            (*iscn).dx = 1;
            (*iscn).du = 1;
            (*iscn).umin = 0;
            while (x as c_uint) < w {
                let d = *p;
                x += 1;
                p = p.offset(1);
                scan_y(scn, d as c_int);
            }
            _zbar_image_scanner_quiet_border(iscn);

            // movedelta(-1, density)
            x -= 1;
            y += density;
            p = p.offset(-1 + (density * w as i32) as isize);
            (*iscn).v = y;
            if (y as c_uint) >= h {
                break;
            }

            (*iscn).dx = -1;
            (*iscn).du = -1;
            (*iscn).umin = w as c_int;
            while x >= 0 {
                let d = *p;
                x -= 1;
                p = p.offset(-1);
                scan_y(scn, d as c_int);
            }
            _zbar_image_scanner_quiet_border(iscn);

            // movedelta(1, density)
            x += 1;
            y += density;
            p = p.offset(1 + (density * w as i32) as isize);
            (*iscn).v = y;
        }
    }
    (*iscn).dx = 0;

    // Vertical scanning pass
    let density = CFG!(iscn, ZBAR_CFG_X_DENSITY);
    if density > 0 {
        let mut p = data;
        let mut x = 0i32;
        let mut y = 0i32;

        let mut border = ((w - 1) % (density as c_uint)).div_ceil(2);
        if border > w / 2 {
            border = w / 2;
        }
        c_assert!(border <= w);
        // movedelta(border, 0)
        x += border as i32;
        p = p.offset(border as isize);
        (*iscn).v = x;

        while (x as c_uint) < w {
            (*iscn).dy = 1;
            (*iscn).du = 1;
            (*iscn).umin = 0;
            while (y as c_uint) < h {
                let d = *p;
                y += 1;
                p = p.offset(w as isize);
                scan_y(scn, d as c_int);
            }
            _zbar_image_scanner_quiet_border(iscn);

            // movedelta(density, -1)
            x += density;
            y -= 1;
            p = p.offset((density as isize) - (w as isize));
            (*iscn).v = x;
            if (x as c_uint) >= w {
                break;
            }

            (*iscn).dy = -1;
            (*iscn).du = -1;
            (*iscn).umin = h as c_int;
            while y >= 0 {
                let d = *p;
                y -= 1;
                p = p.offset(-(w as isize));
                scan_y(scn, d as c_int);
            }
            _zbar_image_scanner_quiet_border(iscn);

            // movedelta(density, 1)
            x += density;
            y += 1;
            p = p.offset((density as isize) + (w as isize));
            (*iscn).v = x;
        }
    }
    (*iscn).dy = 0;

    // Decode QR and SQ codes
    if let Some(qr) = &mut (*iscn).qr {
        _zbar_qr_decode(qr, iscn, img);
    }

    _zbar_image_scanner_sq_handler(iscn);
    if let Some(sq) = &mut (*iscn).sq {
        _zbar_sq_decode(sq, iscn, img);
    }

    // Filter and merge EAN composite results
    let filter = density == 1 || CFG!(iscn, ZBAR_CFG_Y_DENSITY) == 1;
    let mut nean = 0;
    let mut naddon = 0;

    if (*syms).nsyms != 0 {
        let mut symp = &mut (*syms).head as *mut *mut zbar_symbol_t;
        while !(*symp).is_null() {
            let sym = *symp;
            let sym_type = (*sym).symbol_type as c_int;

            if (sym_type < ZBAR_COMPOSITE && sym_type > ZBAR_PARTIAL)
                || sym_type == ZBAR_DATABAR
                || sym_type == ZBAR_DATABAR_EXP
                || sym_type == ZBAR_CODABAR
            {
                if (sym_type == ZBAR_CODABAR || filter) && (*sym).quality < 4 {
                    // Recycle symbol
                    *symp = (*sym).next;
                    (*syms).nsyms -= 1;
                    (*sym).next = null_mut();
                    _zbar_image_scanner_recycle_syms(iscn, sym);
                    continue;
                } else if sym_type < ZBAR_COMPOSITE && sym_type != ZBAR_ISBN10 {
                    if sym_type > ZBAR_EAN5 {
                        nean += 1;
                    } else {
                        naddon += 1;
                    }
                }
            }
            symp = &mut (*sym).next as *mut *mut zbar_symbol_t;
        }

        // Merge EAN composite if we have exactly one EAN and one add-on
        if nean == 1 && naddon == 1 && (*iscn).ean_config != 0 {
            let mut ean: *mut zbar_symbol_t = null_mut();
            let mut addon: *mut zbar_symbol_t = null_mut();

            // Extract EAN and add-on symbols
            let mut symp = &mut (*syms).head as *mut *mut zbar_symbol_t;
            while !(*symp).is_null() {
                let sym = *symp;
                let sym_type = (*sym).symbol_type as c_int;
                if sym_type < ZBAR_COMPOSITE && sym_type > ZBAR_PARTIAL {
                    // Move to composite
                    *symp = (*sym).next;
                    (*syms).nsyms -= 1;
                    (*sym).next = null_mut();
                    if sym_type <= ZBAR_EAN5 {
                        addon = sym;
                    } else {
                        ean = sym;
                    }
                } else {
                    symp = &mut (*sym).next as *mut *mut zbar_symbol_t;
                }
            }
            c_assert!(!ean.is_null());
            c_assert!(!addon.is_null());

            // Create composite symbol
            let datalen = (*ean).datalen + (*addon).datalen + 1;
            let ean_sym = _zbar_image_scanner_alloc_sym(iscn, ZBAR_COMPOSITE, datalen as c_int);
            (*ean_sym).orient = (*ean).orient;
            (*ean_sym).syms = symbol_set_create();

            // Copy data
            copy_nonoverlapping(
                (*ean).data as *const u8,
                (*ean_sym).data as *mut u8,
                (*ean).datalen as usize,
            );
            copy_nonoverlapping(
                (*addon).data as *const u8,
                ((*ean_sym).data as *mut u8).offset((*ean).datalen as isize),
                ((*addon).datalen + 1) as usize,
            );

            // Link symbols
            (*(*ean_sym).syms).head = ean;
            (*ean).next = addon;
            (*(*ean_sym).syms).nsyms = 2;

            _zbar_image_scanner_add_sym(iscn, ean_sym);
        }
    }

    syms
}

/// Public wrapper for zbar_scan_image
///
/// Scans an image for barcodes, with optional inverted image retry.
///
/// # Arguments
/// * `iscn` - Image scanner instance
/// * `img` - Image to scan
///
/// # Returns
/// Number of symbols found, -1 on error
pub unsafe fn zbar_scan_image(iscn: *mut zbar_image_scanner_t, img: *mut zbar_image_t) -> c_int {
    let mut syms = _zbar_scan_image(iscn, img);
    if syms.is_null() {
        return -1;
    }

    let mut inv: *mut zbar_image_t = null_mut();

    // Try inverted image if no symbols found and TEST_INVERTED is enabled
    if (*syms).nsyms == 0 && TEST_CFG!(iscn, ZBAR_CFG_TEST_INVERTED) {
        inv = _zbar_image_copy(img, 1);
        if !inv.is_null() {
            syms = _zbar_scan_image(iscn, inv);
            _zbar_image_swap_symbols(img, inv);
        }
    }

    // Call user handler if symbols found
    if (*syms).nsyms != 0 {
        if let Some(handler) = (*iscn).handler {
            handler(img, (*iscn).userdata);
        }
    }

    if !inv.is_null() {
        zbar_image_destroy(inv);
    }

    (*syms).nsyms as c_int
}
