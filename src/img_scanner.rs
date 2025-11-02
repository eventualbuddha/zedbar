use std::ffi::c_void;

use libc::{c_int, c_uint};

use crate::{
    decoder::{
        zbar_decoder_t, ZBAR_CFG_BINARY, ZBAR_CFG_ENABLE, ZBAR_CFG_POSITION,
        ZBAR_CFG_TEST_INVERTED, ZBAR_CFG_UNCERTAINTY, ZBAR_CFG_X_DENSITY, ZBAR_CFG_Y_DENSITY,
    },
    finder::{decoder_get_qr_finder_line, decoder_get_sq_finder_config},
    image_ffi::zbar_image_t,
    line_scanner::zbar_scanner_t,
    qrcode::qrdec::qr_reader,
    sqcode::{sq_decode, SqReader},
    symbol::zbar_symbol_t,
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

#[derive(Default, Clone)]
pub(crate) struct zbar_symbol_set_t {
    pub(crate) symbols: Vec<zbar_symbol_t>,
}

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

    /// current scan direction
    dx: c_int,
    dy: c_int,
    du: c_int,
    umin: c_int,
    v: c_int,

    /// previous decode results
    syms: zbar_symbol_set_t,

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
    /// Create a new image scanner
    ///
    /// Allocates and initializes a new image scanner instance with default configuration.
    ///
    /// # Returns
    /// Pointer to new scanner or null on allocation failure
    pub(crate) unsafe fn new() -> *mut Self {
        let mut iscn = Box::new(Self::default());
        let decoder = zbar_decoder_t::new();

        iscn.set_decoder(Some(decoder));

        // Get a pointer to the decoder for the scanner
        let dcode_ptr = iscn.dcode.as_mut().unwrap() as *mut zbar_decoder_t;
        let scanner = zbar_scanner_t::new(dcode_ptr as *mut _);
        iscn.set_scanner(Some(scanner));

        if iscn.decoder().is_none() || iscn.scanner().is_none() {
            return std::ptr::null_mut();
        }

        iscn.set_qr_reader(Some(qr_reader::new()));
        iscn.set_sq_reader(Some(SqReader::new()));

        // Apply default configuration
        iscn.configs[0] = 1; // ZBAR_CFG_X_DENSITY
        iscn.configs[1] = 1; // ZBAR_CFG_Y_DENSITY

        iscn.set_config(SymbolType::None, ZBAR_CFG_POSITION, 1);
        iscn.set_config(SymbolType::None, ZBAR_CFG_UNCERTAINTY, 2);
        iscn.set_config(SymbolType::None, 65, 0); // ZBAR_CFG_TEST_INVERTED
        iscn.set_config(SymbolType::QrCode, ZBAR_CFG_UNCERTAINTY, 0);
        iscn.set_config(SymbolType::QrCode, ZBAR_CFG_BINARY, 0);
        iscn.set_config(SymbolType::Code128, ZBAR_CFG_UNCERTAINTY, 0);
        iscn.set_config(SymbolType::Code93, ZBAR_CFG_UNCERTAINTY, 0);
        iscn.set_config(SymbolType::Code39, ZBAR_CFG_UNCERTAINTY, 0);
        iscn.set_config(SymbolType::Codabar, ZBAR_CFG_UNCERTAINTY, 1);
        iscn.set_config(SymbolType::Composite, ZBAR_CFG_UNCERTAINTY, 0);

        let iscn_ptr = Box::into_raw(iscn);

        if let Some(dcode) = (*iscn_ptr).decoder_mut() {
            dcode.set_userdata(iscn_ptr as *mut c_void);
            dcode.set_handler(Some(symbol_handler));
        }

        iscn_ptr
    }

    /// Add a symbol to the scanner's symbol set
    ///
    /// # Arguments
    /// * `sym` - The symbol to add
    pub(crate) fn add_symbol(&mut self, sym: zbar_symbol_t) {
        let syms = &mut self.syms;
        syms.symbols.push(sym);
    }

    pub(crate) fn find_duplicate_symbol(
        &mut self,
        symbol_type: SymbolType,
        data: &[u8],
    ) -> Option<&mut zbar_symbol_t> {
        self.syms
            .symbols
            .iter_mut()
            .find(|sym| sym.symbol_type == symbol_type && sym.data == data)
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

    /// Public wrapper for zbar_scan_image
    ///
    /// Scans an image for barcodes, with optional inverted image retry.
    ///
    /// # Arguments
    /// * `img` - Image to scan
    ///
    /// # Returns
    /// Number of symbols found, -1 on error
    pub(crate) unsafe fn scan_image(&mut self, img: &mut zbar_image_t) -> Result<c_int> {
        let Some(syms) = self._zbar_scan_image(img) else {
            return Err(Error::Unknown(-1));
        };

        let nsyms = syms.symbols.len();

        // Try inverted image if no symbols found and TEST_INVERTED is enabled
        if nsyms == 0 && TEST_CFG!(self, ZBAR_CFG_TEST_INVERTED) {
            if let Some(mut inv) = img.copy(true) {
                let _ = self._zbar_scan_image(&mut inv);
                img.swap_symbols_with(&mut inv);
            }
        }

        // Call user handler if symbols found
        let final_nsyms = img.syms().map_or(0, |s| s.symbols.len());
        Ok(final_nsyms as c_int)
    }

    /// Handle QR code finder line detection
    ///
    /// Processes a QR code finder line from the decoder and forwards it to the QR reader.
    /// Adjusts edge positions based on scanner state and transforms coordinates.
    pub(crate) unsafe fn qr_handler(&mut self) {
        let Some(dcode) = &mut self.dcode else {
            return;
        };
        let line = decoder_get_qr_finder_line(dcode);

        let Some(scn) = self.scn.as_ref() else {
            return;
        };
        let mut u = scn.get_edge(line.pos[0] as c_uint, QR_FINDER_SUBPREC);
        line.boffs = (u as c_int) - scn.get_edge(line.boffs as c_uint, QR_FINDER_SUBPREC) as c_int;
        line.len = scn.get_edge(line.len as c_uint, QR_FINDER_SUBPREC) as c_int;
        line.eoffs = scn.get_edge(line.eoffs as c_uint, QR_FINDER_SUBPREC) as c_int - line.len;
        line.len -= u as c_int;

        u = (qr_fixed(self.umin, 0) as i64 + (self.du as i64) * (u as i64)) as c_uint;
        if self.du < 0 {
            std::mem::swap(&mut line.boffs, &mut line.eoffs);
            u = u.wrapping_sub(line.len as c_uint);
        }

        let vert: c_int = if self.dx != 0 { 0 } else { 1 };
        line.pos[vert as usize] = u as c_int;
        line.pos[(1 - vert) as usize] = qr_fixed(self.v, 1) as c_int;

        if let Some(qr) = &mut self.qr {
            qr.found_line(vert, line);
        }
    }

    /// SQ code handler - updates SQ reader configuration
    ///
    /// Gets the current SQ finder configuration from the decoder and passes it
    /// to the SQ reader for processing.
    pub(crate) unsafe fn sq_handler(&mut self) {
        let Some(dcode) = &self.dcode else {
            return;
        };
        let config = decoder_get_sq_finder_config(dcode);
        if let Some(sq) = &mut self.sq {
            sq.set_enabled(config != 0);
        }
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
    pub(crate) unsafe fn set_config(&mut self, sym: SymbolType, cfg: c_int, val: c_int) -> c_int {
        // Handle EAN composite configuration
        if (sym == SymbolType::None || sym == SymbolType::Composite) && cfg == ZBAR_CFG_ENABLE {
            self.ean_config = if val != 0 { 1 } else { 0 };
            if sym != SymbolType::None {
                return 0;
            }
        }

        // Delegate decoder configuration
        if cfg < ZBAR_CFG_UNCERTAINTY {
            if let Some(dcode) = self.decoder_mut() {
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
                self.sym_configs[c][i] = val;
            } else {
                for i in 0..NUM_SYMS {
                    self.sym_configs[c][i] = val;
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
            self.configs[(cfg - ZBAR_CFG_X_DENSITY) as usize] = val;
            return 0;
        }

        // Handle position and related configuration flags
        let cfg_bit = cfg - ZBAR_CFG_POSITION;

        if val == 0 {
            self.config &= !(1 << cfg_bit);
        } else if val == 1 {
            self.config |= 1 << cfg_bit;
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
    unsafe fn _zbar_scan_image(&mut self, img: &mut zbar_image_t) -> Option<zbar_symbol_set_t> {
        // Reset QR and SQ decoders
        if let Some(qr) = &mut self.qr {
            qr.reset();
        }
        if let Some(sq) = &mut self.sq {
            sq.reset();
        }

        // Image must be in grayscale format
        if img.format != fourcc(b'Y', b'8', b'0', b'0')
            && img.format != fourcc(b'G', b'R', b'E', b'Y')
        {
            return None;
        }

        // Clear previous symbols for new scan
        let scanner = &mut *self;
        scanner.syms.symbols.clear();

        // Share the symbol set with the image (we'll clone it at the end)
        // For now, we'll handle this differently - the image will get a clone of the results

        let w = img.width;
        let h = img.height;
        let data = img.data.as_ptr();
        let scn = self.scn.as_mut()?;
        scn.new_scan();

        // Horizontal scanning pass
        let density = CFG!(self, ZBAR_CFG_Y_DENSITY);
        if density > 0 {
            let mut p = data;
            let mut x = 0i32;
            let mut y = 0i32;

            let mut border = ((h - 1) % (density as c_uint)).div_ceil(2);
            if border > h / 2 {
                border = h / 2;
            }
            c_assert!(border <= h);
            self.dy = 0;

            // movedelta(0, border)
            y += border as i32;
            p = p.offset((border * w) as isize);
            self.v = y;

            while (y as c_uint) < h {
                self.dx = 1;
                self.du = 1;
                self.umin = 0;
                while (x as c_uint) < w {
                    let d = *p;
                    x += 1;
                    p = p.offset(1);
                    scn.scan_y(d as c_int);
                }
                scn.quiet_border();

                // movedelta(-1, density)
                x -= 1;
                y += density;
                p = p.offset(-1 + (density * w as i32) as isize);
                self.v = y;
                if (y as c_uint) >= h {
                    break;
                }

                self.dx = -1;
                self.du = -1;
                self.umin = w as c_int;
                while x >= 0 {
                    let d = *p;
                    x -= 1;
                    p = p.offset(-1);
                    scn.scan_y(d as c_int);
                }
                scn.quiet_border();

                // movedelta(1, density)
                x += 1;
                y += density;
                p = p.offset(1 + (density * w as i32) as isize);
                self.v = y;
            }
        }
        self.dx = 0;

        // Vertical scanning pass
        let density = CFG!(self, ZBAR_CFG_X_DENSITY);
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
            self.v = x;

            while (x as c_uint) < w {
                self.dy = 1;
                self.du = 1;
                self.umin = 0;
                while (y as c_uint) < h {
                    let d = *p;
                    y += 1;
                    p = p.offset(w as isize);
                    scn.scan_y(d as c_int);
                }
                scn.quiet_border();

                // movedelta(density, -1)
                x += density;
                y -= 1;
                p = p.offset((density as isize) - (w as isize));
                self.v = x;
                if (x as c_uint) >= w {
                    break;
                }

                self.dy = -1;
                self.du = -1;
                self.umin = h as c_int;
                while y >= 0 {
                    let d = *p;
                    y -= 1;
                    p = p.offset(-(w as isize));
                    scn.scan_y(d as c_int);
                }
                scn.quiet_border();

                // movedelta(density, 1)
                x += density;
                y += 1;
                p = p.offset((density as isize) + (w as isize));
                self.v = x;
            }
        }
        self.dy = 0;

        // Decode QR and SQ codes
        let raw_binary = self
            .get_config(SymbolType::QrCode, ZBAR_CFG_BINARY)
            .unwrap_or(0)
            != 0;
        if let Some(qr) = &mut self.qr {
            let qr_symbols = qr.decode(img, raw_binary);
            for sym in qr_symbols {
                self.add_symbol(sym);
            }
        }

        self.sq_handler();
        if let Some(sq) = &mut self.sq {
            if let Ok(Some(symbol)) = sq_decode(sq, img) {
                self.add_symbol(symbol);
            }
        }

        // Filter and merge EAN composite results
        let filter = density == 1 || CFG!(self, ZBAR_CFG_Y_DENSITY) == 1;
        let scanner = &mut *self;

        // Filter low-quality symbols
        scanner.syms.symbols.retain(|sym| {
            let sym_type = sym.symbol_type;
            if (sym_type < SymbolType::Composite && sym_type > SymbolType::Partial)
                || sym_type == SymbolType::Databar
                || sym_type == SymbolType::DatabarExp
                || sym_type == SymbolType::Codabar
            {
                // Keep if quality >= 4 OR if not Codabar and not filtered
                !((sym_type == SymbolType::Codabar || filter) && sym.quality < 4)
            } else {
                true
            }
        });

        // Count EAN and add-on symbols for potential merging
        let mut nean = 0;
        let mut naddon = 0;
        for sym in &scanner.syms.symbols {
            let sym_type = sym.symbol_type;
            if sym_type < SymbolType::Composite && sym_type != SymbolType::Isbn10 {
                if sym_type > SymbolType::Ean5 {
                    nean += 1;
                } else if sym_type > SymbolType::Partial {
                    naddon += 1;
                }
            }
        }

        // Merge EAN composite if we have exactly one EAN and one add-on
        let ean_config = scanner.ean_config;
        if nean == 1 && naddon == 1 && ean_config != 0 {
            // Extract EAN and add-on symbols
            let mut ean_sym: Option<zbar_symbol_t> = None;
            let mut addon_sym: Option<zbar_symbol_t> = None;
            let mut other_syms = Vec::new();

            for sym in scanner.syms.symbols.drain(..) {
                let sym_type = sym.symbol_type;
                if sym_type < SymbolType::Composite && sym_type > SymbolType::Partial {
                    if sym_type <= SymbolType::Ean5 {
                        addon_sym = Some(sym);
                    } else {
                        ean_sym = Some(sym);
                    }
                } else {
                    other_syms.push(sym);
                }
            }

            if let (Some(ean), Some(addon)) = (ean_sym, addon_sym) {
                // Create composite symbol
                let mut composite = zbar_symbol_t {
                    symbol_type: SymbolType::Composite,
                    orient: ean.orient,
                    quality: 1,
                    ..Default::default()
                };

                // Copy data from both symbols
                composite.data.extend_from_slice(&ean.data);
                composite.data.extend_from_slice(&addon.data);

                // Create component symbol set
                let mut component_set = zbar_symbol_set_t::default();
                component_set.symbols.push(ean);
                component_set.symbols.push(addon);
                composite.components = Some(component_set);

                // Add composite to results
                other_syms.push(composite);
            }

            scanner.syms.symbols = other_syms;
        }

        // Clone the symbol set to the image
        let scanner = &mut *self;
        // Create a clone of the symbol set for the image
        let img_syms = Box::new(zbar_symbol_set_t {
            symbols: scanner.syms.symbols.clone(),
        });
        img.set_syms(img_syms);

        Some(scanner.syms.clone())
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
    fn get_config(&mut self, sym: SymbolType, cfg: c_int) -> Result<c_int, c_int> {
        // Return error if symbol doesn't have config
        if !(SymbolType::Partial..=SymbolType::Code128).contains(&sym)
            || sym == SymbolType::Composite
        {
            return Err(1);
        }

        if cfg < ZBAR_CFG_UNCERTAINTY {
            if let Some(dcode) = &mut self.dcode.as_mut() {
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
            return Ok(self.sym_configs[(cfg - ZBAR_CFG_UNCERTAINTY) as usize][i as usize]);
        }

        // Image scanner parameters apply only to ZBAR_PARTIAL
        if sym > SymbolType::Partial {
            return Err(1);
        }

        if cfg < ZBAR_CFG_X_DENSITY {
            return Ok(if (self.config & (1 << (cfg - ZBAR_CFG_POSITION))) != 0 {
                1
            } else {
                0
            });
        }

        if cfg <= ZBAR_CFG_Y_DENSITY {
            // CFG macro: ((self)->configs[(cfg) - ZBAR_CFG_X_DENSITY])
            return Ok(self.configs[(cfg - ZBAR_CFG_X_DENSITY) as usize]);
        }

        Err(1)
    }
}

/// Symbol handler callback for 1D barcode decoding
///
/// This function is called by the decoder when a barcode is successfully decoded.
/// It manages symbol deduplication, position tracking, and adds symbols to the result set.
///
/// # Arguments
/// * `dcode` - The decoder that found the symbol
pub(crate) unsafe fn symbol_handler(dcode: &mut zbar_decoder_t) {
    let iscn = dcode.get_userdata() as *mut zbar_image_scanner_t;
    let symbol_type = dcode.get_type();
    let mut x = 0;
    let mut y = 0;

    // QR codes are handled separately
    if symbol_type == SymbolType::QrCode {
        (&mut *iscn).qr_handler();
        return;
    }

    // Calculate position if position tracking is enabled
    if TEST_CFG!(iscn, ZBAR_CFG_POSITION) {
        if let Some(scn) = (*iscn).scanner() {
            let w = scn.width();
            let u = (*iscn).umin + (*iscn).du * scn.get_edge(w, 0) as c_int;
            if (*iscn).dx != 0 {
                x = u;
                y = (*iscn).v;
            } else {
                x = (*iscn).v;
                y = u;
            }
        }
    }

    // Ignore partial results
    if symbol_type <= SymbolType::Partial {
        return;
    }

    let data = dcode.buffer_slice();
    if let Some(sym) = (*iscn).find_duplicate_symbol(symbol_type, data) {
        sym.quality += 1;
        if TEST_CFG!(iscn, ZBAR_CFG_POSITION) {
            sym.add_point(x, y);
        }
        return;
    }

    // Allocate new symbol
    let mut sym = zbar_symbol_t::new(symbol_type);
    sym.configs = dcode.get_configs(symbol_type);
    sym.modifiers = dcode.get_modifiers();

    // Copy data
    sym.data.extend_from_slice(data);

    // Initialize position
    if TEST_CFG!(iscn, ZBAR_CFG_POSITION) {
        sym.add_point(x, y);
    }

    // Set orientation
    let dir = dcode.direction;
    if dir != 0 {
        sym.orient = (if (*iscn).dy != 0 { 1 } else { 0 }) + (((*iscn).du ^ dir) & 2);
    }

    (*iscn).add_symbol(sym);
}
