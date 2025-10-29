use std::{
    ffi::c_void,
    ptr::{null_mut, NonNull},
};

use libc::{c_int, c_uint};

use crate::{
    decoder::{
        zbar_decoder_t, ZBAR_CFG_BINARY, ZBAR_CFG_ENABLE, ZBAR_CFG_POSITION,
        ZBAR_CFG_TEST_INVERTED, ZBAR_CFG_UNCERTAINTY, ZBAR_CFG_X_DENSITY, ZBAR_CFG_Y_DENSITY,
        ZBAR_ORIENT_UNKNOWN,
    },
    finder::{decoder_get_qr_finder_line, decoder_get_sq_finder_config},
    image_ffi::zbar_image_t,
    line_scanner::{
        scan_y, scanner_flush, scanner_get_edge, scanner_get_width, scanner_new_scan,
        zbar_scanner_new, zbar_scanner_t,
    },
    qrcode::qrdec::{
        _zbar_qr_create, _zbar_qr_found_line, _zbar_qr_reset, qr_decode, qr_finder_lines,
    },
    refcnt,
    sqcode::{sq_decode, SqReader},
    symbol::{
        symbol_alloc_zeroed, symbol_free, symbol_refcnt, symbol_set_create, symbol_set_free,
        zbar_symbol_set_ref, zbar_symbol_t,
    },
    Error, Result, SymbolType,
};

const NUM_SCN_CFGS: usize = 2; // ZBAR_CFG_Y_DENSITY - ZBAR_CFG_X_DENSITY + 1
const NUM_SYMS: usize = 25; // Number of symbol types

// QR Code finder precision constant
const QR_FINDER_SUBPREC: c_int = 2;

// QR_FIXED macro: ((((v) << 1) + (rnd)) << (QR_FINDER_SUBPREC - 1))
#[inline]
fn qr_fixed(v: c_int, rnd: c_int) -> c_uint {
    (((v as c_uint) << 1) + (rnd as c_uint)) << (QR_FINDER_SUBPREC - 1)
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
use crate::image_ffi::_zbar_image_copy;

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

pub(crate) struct qr_reader {
    /// The random number generator used by RANSAC.
    pub(crate) rng: rand_chacha::ChaCha8Rng,
    ///  current finder state, horizontal and vertical lines
    pub(crate) finder_lines: [qr_finder_lines; 2],
}

#[derive(Default)]
pub(crate) struct zbar_symbol_set_t {
    pub(crate) refcnt: c_int,
    pub(crate) nsyms: c_int,
    pub(crate) head: *mut zbar_symbol_t,
    pub(crate) tail: *mut zbar_symbol_t,
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
                symbol_refcnt(&mut *sym, -1);
            }
            sym = next;
        }
        self.head = std::ptr::null_mut();
    }
}

// Function pointer type for image data handler callbacks
pub(crate) type zbar_image_data_handler_t =
    unsafe fn(img: *mut zbar_image_t, userdata: *const c_void);

/// image scanner state
#[derive(Default)]
pub(crate) struct zbar_image_scanner_t {
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
    /// Add a symbol to the scanner's symbol set
    ///
    /// # Arguments
    /// * `sym` - The symbol to add
    pub(crate) unsafe fn add_symbol(&mut self, sym: &mut zbar_symbol_t) {
        let syms = self.syms.expect("syms is set").as_mut();

        // The symbol set takes ownership of the symbol reference.
        // Ensure the reference count reflects this so that recycling
        // and Drop logic can safely release it later.
        let new_refcnt = refcnt!(sym.refcnt, 1);
        debug_assert!(new_refcnt > 0);
        let _ = new_refcnt;

        if syms.tail.is_null() {
            sym.next = syms.head;
            syms.head = sym;
        } else {
            let tail = &mut *syms.tail;
            sym.next = tail.next;
            tail.next = sym;
        }

        syms.nsyms += 1;
    }

    pub(crate) unsafe fn find_duplicate_symbol(
        &mut self,
        symbol_type: SymbolType,
        data: &[u8],
    ) -> *mut zbar_symbol_t {
        let syms = self.syms_ptr();
        let mut sym = if !syms.is_null() {
            (*syms).head
        } else {
            null_mut()
        };
        while !sym.is_null() {
            let sym_ref = &mut *sym;
            if sym_ref.symbol_type == symbol_type && sym_ref.data == data {
                return sym;
            }
            sym = sym_ref.next;
        }
        null_mut()
    }

    #[inline]
    pub(crate) fn syms_ptr(&self) -> *mut zbar_symbol_set_t {
        self.syms.map_or(null_mut(), NonNull::as_ptr)
    }

    #[inline]
    pub(crate) fn set_syms_ptr(&mut self, ptr: *mut zbar_symbol_set_t) {
        self.syms = NonNull::new(ptr);
    }

    #[inline]
    pub(crate) fn take_syms(&mut self) -> Option<NonNull<zbar_symbol_set_t>> {
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

/// Get configuration value for a specific symbology
///
/// Retrieves the current configuration value for a particular setting
/// of a barcode symbology.
///
/// # Arguments
/// * `iscn` - The image scanner instance
/// * `sym` - The symbology type to query
/// * `cfg` - The configuration parameter to query
///
/// # Returns
/// `Ok(value)` on success, `Err(1)` on error
pub(crate) fn zbar_image_scanner_get_config(
    iscn: &mut zbar_image_scanner_t,
    sym: SymbolType,
    cfg: c_int,
) -> Result<c_int, c_int> {
    // Return error if symbol doesn't have config
    if !(SymbolType::Partial..=SymbolType::Code128).contains(&sym) || sym == SymbolType::Composite {
        return Err(1);
    }

    if cfg < ZBAR_CFG_UNCERTAINTY {
        if let Some(dcode) = &mut iscn.dcode.as_mut() {
            return dcode.get_config(sym, cfg);
        } else {
            return Err(1);
        }
    }

    if cfg < ZBAR_CFG_POSITION {
        if sym == SymbolType::Partial {
            return Err(1);
        }

        let i = sym.hash();
        return Ok(iscn.sym_configs[(cfg - ZBAR_CFG_UNCERTAINTY) as usize][i as usize]);
    }

    // Image scanner parameters apply only to ZBAR_PARTIAL
    if sym > SymbolType::Partial {
        return Err(1);
    }

    if cfg < ZBAR_CFG_X_DENSITY {
        return Ok(if (iscn.config & (1 << (cfg - ZBAR_CFG_POSITION))) != 0 {
            1
        } else {
            0
        });
    }

    if cfg <= ZBAR_CFG_Y_DENSITY {
        // CFG macro: ((iscn)->configs[(cfg) - ZBAR_CFG_X_DENSITY])
        return Ok(iscn.configs[(cfg - ZBAR_CFG_X_DENSITY) as usize]);
    }

    Err(1)
}

/// Recycle symbols from the image scanner
///
/// Recursively processes a linked list of symbols, either unlinking referenced
/// symbols or freeing unreferenced ones.
///
/// # Arguments
/// * `iscn` - The image scanner instance
/// * `sym` - Head of the symbol list to process
pub(crate) unsafe fn _zbar_image_scanner_recycle_syms(
    _iscn: &mut zbar_image_scanner_t,
    mut sym: *mut zbar_symbol_t,
) {
    while !sym.is_null() {
        let next = (*sym).next;

        if (*sym).refcnt != 0 && refcnt!((*sym).refcnt, -1) != 0 {
            // Unlink referenced symbol
            // FIXME handle outstanding component refs (currently unsupported)
            c_assert!((*sym).data.capacity() != 0);
            (*sym).next = null_mut();
        } else {
            // Free unreferenced symbol
            if !(*sym).syms.is_null() {
                let syms = (*sym).syms;
                if refcnt!((*syms).refcnt, -1) != 0 {
                    c_assert!(false);
                }
                _zbar_image_scanner_recycle_syms(_iscn, (*syms).head);
                (*syms).head = null_mut();
                symbol_set_free(syms);
                (*sym).syms = null_mut();
            }

            symbol_free(sym);
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
pub(crate) unsafe fn _zbar_image_scanner_quiet_border(iscn: &mut zbar_image_scanner_t) {
    let Some(scn) = iscn.scn.as_mut() else {
        return;
    };

    // Flush scanner pipeline twice
    scanner_flush(scn);
    scanner_flush(scn);

    // Start new scan
    scanner_new_scan(scn);
}

/// Allocate a new symbol
///
/// # Arguments
/// * `iscn` - The image scanner instance
/// * `sym_type` - The type of symbol to allocate
///
/// # Returns
/// Pointer to the allocated symbol
pub(crate) unsafe fn _zbar_image_scanner_alloc_sym(
    _iscn: &mut zbar_image_scanner_t,
    symbol_type: SymbolType,
) -> *mut zbar_symbol_t {
    let sym = symbol_alloc_zeroed();

    // Initialize the symbol
    let sym_ref = &mut *sym;
    sym_ref.symbol_type = symbol_type;
    sym_ref.quality = 1;
    sym_ref.orient = ZBAR_ORIENT_UNKNOWN;
    c_assert!(sym_ref.syms.is_null());

    sym
}

/// Handle QR code finder line detection
///
/// Processes a QR code finder line from the decoder and forwards it to the QR reader.
/// Adjusts edge positions based on scanner state and transforms coordinates.
///
/// # Arguments
/// * `iscn` - The image scanner instance
pub(crate) unsafe fn _zbar_image_scanner_qr_handler(iscn: &mut zbar_image_scanner_t) {
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

/// SQ code handler - updates SQ reader configuration
///
/// Gets the current SQ finder configuration from the decoder and passes it
/// to the SQ reader for processing.
///
/// # Safety
/// `iscn` must be a valid pointer to a zbar_image_scanner_t
pub(crate) unsafe fn _zbar_image_scanner_sq_handler(iscn: &mut zbar_image_scanner_t) {
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
        dcode.set_userdata(iscn as *mut c_void);
        dcode.set_handler(Some(symbol_handler));
    }

    iscn_ref.set_qr_reader(Some(_zbar_qr_create()));
    iscn_ref.set_sq_reader(Some(SqReader::new()));

    // Apply default configuration
    iscn_ref.configs[0] = 1; // ZBAR_CFG_X_DENSITY
    iscn_ref.configs[1] = 1; // ZBAR_CFG_Y_DENSITY

    zbar_image_scanner_set_config(iscn_ref, SymbolType::None, ZBAR_CFG_POSITION, 1);
    zbar_image_scanner_set_config(iscn_ref, SymbolType::None, ZBAR_CFG_UNCERTAINTY, 2);
    zbar_image_scanner_set_config(iscn_ref, SymbolType::None, 65, 0); // ZBAR_CFG_TEST_INVERTED
    zbar_image_scanner_set_config(iscn_ref, SymbolType::QrCode, ZBAR_CFG_UNCERTAINTY, 0);
    zbar_image_scanner_set_config(iscn_ref, SymbolType::QrCode, ZBAR_CFG_BINARY, 0);
    zbar_image_scanner_set_config(iscn_ref, SymbolType::Code128, ZBAR_CFG_UNCERTAINTY, 0);
    zbar_image_scanner_set_config(iscn_ref, SymbolType::Code93, ZBAR_CFG_UNCERTAINTY, 0);
    zbar_image_scanner_set_config(iscn_ref, SymbolType::Code39, ZBAR_CFG_UNCERTAINTY, 0);
    zbar_image_scanner_set_config(iscn_ref, SymbolType::Codabar, ZBAR_CFG_UNCERTAINTY, 1);
    zbar_image_scanner_set_config(iscn_ref, SymbolType::Composite, ZBAR_CFG_UNCERTAINTY, 0);

    iscn
}

/// Destroy an image scanner
///
/// Frees all resources associated with the scanner.
///
/// # Parameters
/// - `iscn`: Scanner to destroy (null-safe)
pub(crate) unsafe fn zbar_image_scanner_destroy(iscn: *mut zbar_image_scanner_t) {
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

    image_scanner_free(iscn);
}

/// Symbol handler callback for 1D barcode decoding
///
/// This function is called by the decoder when a barcode is successfully decoded.
/// It manages symbol deduplication, position tracking, and adds symbols to the result set.
///
/// # Arguments
/// * `dcode` - The decoder that found the symbol
pub(crate) unsafe fn symbol_handler(dcode: *mut zbar_decoder_t) {
    let iscn = (*dcode).get_userdata() as *mut zbar_image_scanner_t;
    let symbol_type = (*dcode).get_type();
    let mut x = 0;
    let mut y = 0;

    // QR codes are handled separately
    if symbol_type == SymbolType::QrCode {
        _zbar_image_scanner_qr_handler(&mut *iscn);
        return;
    }

    // Calculate position if position tracking is enabled
    let iscn_ref = &mut *iscn;
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
    if symbol_type <= SymbolType::Partial {
        return;
    }

    let data = (*dcode).buffer_slice();
    let sym = iscn_ref.find_duplicate_symbol(symbol_type, data);
    if !sym.is_null() {
        (*sym).quality += 1;
        if TEST_CFG!(iscn, ZBAR_CFG_POSITION) {
            (*sym).add_point(x, y);
        }
        return;
    }

    // Allocate new symbol
    let sym = _zbar_image_scanner_alloc_sym(&mut *iscn, symbol_type);
    let sym_ref = &mut *sym;
    sym_ref.configs = (*dcode).get_configs(symbol_type);
    sym_ref.modifiers = (*dcode).get_modifiers();

    // Copy data
    sym_ref.data.extend_from_slice(data);

    // Initialize position
    if TEST_CFG!(iscn, ZBAR_CFG_POSITION) {
        sym_ref.add_point(x, y);
    }

    // Set orientation
    let dir = (*dcode).direction;
    if dir != 0 {
        sym_ref.orient = (if iscn_ref.dy != 0 { 1 } else { 0 }) + ((iscn_ref.du ^ dir) & 2);
    }

    iscn_ref.add_symbol(&mut *sym);
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
    iscn: &mut zbar_image_scanner_t,
    sym: SymbolType,
    cfg: c_int,
    val: c_int,
) -> c_int {
    // Handle EAN composite configuration
    if (sym == SymbolType::None || sym == SymbolType::Composite) && cfg == ZBAR_CFG_ENABLE {
        iscn.ean_config = if val != 0 { 1 } else { 0 };
        if sym != SymbolType::None {
            return 0;
        }
    }

    // Delegate decoder configuration
    if cfg < ZBAR_CFG_UNCERTAINTY {
        if let Some(dcode) = iscn.decoder_mut() {
            return dcode.set_config(sym, cfg, val);
        }
        return 1;
    }

    // Handle uncertainty and related configs
    if cfg < ZBAR_CFG_POSITION {
        if cfg > ZBAR_CFG_UNCERTAINTY {
            return 1;
        }
        let c = (cfg - ZBAR_CFG_UNCERTAINTY) as usize;
        if sym > SymbolType::Partial {
            let i = sym.hash() as usize;
            iscn.sym_configs[c][i] = val;
        } else {
            for i in 0..NUM_SYMS {
                iscn.sym_configs[c][i] = val;
            }
        }
        return 0;
    }

    // Image scanner parameters apply only to ZBAR_PARTIAL
    if sym > SymbolType::Partial {
        return 1;
    }

    // Handle density configuration
    if (ZBAR_CFG_X_DENSITY..=ZBAR_CFG_Y_DENSITY).contains(&cfg) {
        iscn.configs[(cfg - ZBAR_CFG_X_DENSITY) as usize] = val;
        return 0;
    }

    // Handle position and related configuration flags
    let cfg_bit = cfg - ZBAR_CFG_POSITION;

    if val == 0 {
        iscn.config &= !(1 << cfg_bit);
    } else if val == 1 {
        iscn.config |= 1 << cfg_bit;
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
pub(crate) unsafe fn _zbar_scan_image(
    iscn: *mut zbar_image_scanner_t,
    img: &mut zbar_image_t,
) -> *mut zbar_symbol_set_t {
    // Reset QR and SQ decoders
    if let Some(qr) = &mut (*iscn).qr {
        _zbar_qr_reset(qr);
    }
    if let Some(sq) = &mut (*iscn).sq {
        sq.reset();
    }

    // Image must be in grayscale format
    if img.format != fourcc(b'Y', b'8', b'0', b'0') && img.format != fourcc(b'G', b'R', b'E', b'Y')
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
    img.set_syms_ptr(syms);

    let w = img.width;
    let h = img.height;
    let data = img.data.as_ptr();

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
            _zbar_image_scanner_quiet_border(&mut *iscn);

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
            _zbar_image_scanner_quiet_border(&mut *iscn);

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
            _zbar_image_scanner_quiet_border(&mut *iscn);

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
            _zbar_image_scanner_quiet_border(&mut *iscn);

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
        qr_decode(qr, &mut *iscn, img);
    }

    _zbar_image_scanner_sq_handler(&mut *iscn);
    if let Some(sq) = &mut (*iscn).sq {
        sq_decode(sq, &mut *iscn, img);
    }

    // Filter and merge EAN composite results
    let filter = density == 1 || CFG!(iscn, ZBAR_CFG_Y_DENSITY) == 1;
    let mut nean = 0;
    let mut naddon = 0;

    if (*syms).nsyms != 0 {
        let mut symp = &mut (*syms).head as *mut *mut zbar_symbol_t;
        while !(*symp).is_null() {
            let sym = *symp;
            let sym_type = (*sym).symbol_type;

            if (sym_type < SymbolType::Composite && sym_type > SymbolType::Partial)
                || sym_type == SymbolType::Databar
                || sym_type == SymbolType::DatabarExp
                || sym_type == SymbolType::Codabar
            {
                if (sym_type == SymbolType::Codabar || filter) && (*sym).quality < 4 {
                    // Recycle symbol
                    *symp = (*sym).next;
                    (*syms).nsyms -= 1;
                    (*sym).next = null_mut();
                    _zbar_image_scanner_recycle_syms(&mut *iscn, sym);
                    continue;
                } else if sym_type < SymbolType::Composite && sym_type != SymbolType::Isbn10 {
                    if sym_type > SymbolType::Ean5 {
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
                let sym_type = (*sym).symbol_type;
                if sym_type < SymbolType::Composite && sym_type > SymbolType::Partial {
                    // Move to composite
                    *symp = (*sym).next;
                    (*syms).nsyms -= 1;
                    (*sym).next = null_mut();
                    if sym_type <= SymbolType::Ean5 {
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
            let ean_sym = _zbar_image_scanner_alloc_sym(&mut *iscn, SymbolType::Composite);
            (*ean_sym).orient = (*ean).orient;
            (*ean_sym).syms = symbol_set_create();

            // Copy data
            (*ean_sym).data.extend_from_slice(&(*ean).data);
            (*ean_sym).data.extend_from_slice(&(*addon).data);

            // Link symbols
            (*(*ean_sym).syms).head = ean;
            (*ean).next = addon;
            (*(*ean_sym).syms).nsyms = 2;

            (&mut *iscn).add_symbol(&mut *ean_sym);
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
pub(crate) unsafe fn zbar_scan_image(
    iscn: *mut zbar_image_scanner_t,
    img: &mut zbar_image_t,
) -> Result<c_int> {
    let mut syms = _zbar_scan_image(iscn, img);
    if syms.is_null() {
        return Err(Error::Unknown(-1));
    }

    // Try inverted image if no symbols found and TEST_INVERTED is enabled
    if (*syms).nsyms == 0 && TEST_CFG!(iscn, ZBAR_CFG_TEST_INVERTED) {
        if let Some(mut inv) = _zbar_image_copy(&*img, 1) {
            syms = _zbar_scan_image(iscn, &mut inv);
            img.swap_symbols_with(&mut inv);
        }
    }

    // Call user handler if symbols found
    if (*syms).nsyms != 0 {
        if let Some(handler) = (*iscn).handler {
            handler(img, (*iscn).userdata);
        }
    }

    Ok((*syms).nsyms as c_int)
}
