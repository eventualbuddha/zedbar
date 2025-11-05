use std::ffi::c_void;

use libc::{c_int, c_uint};

use crate::{
    config::{internal::DecoderState, DecoderConfig},
    decoder::zbar_decoder_t,
    finder::decoder_get_qr_finder_line,
    image_ffi::zbar_image_t,
    img_scanner_config::ImageScannerConfig,
    qrcode::qrdec::qr_reader,
    sqcode::SqReader,
    symbol::zbar_symbol_t,
    Error, Result, SymbolType,
};

// QR Code finder precision constant
const QR_FINDER_SUBPREC: c_int = 2;

// QR_FIXED macro: ((((v) << 1) + (rnd)) << (QR_FINDER_SUBPREC - 1))
#[inline]
fn qr_fixed(v: c_int, rnd: c_int) -> c_uint {
    (((v as c_uint) << 1) + (rnd as c_uint)) << (QR_FINDER_SUBPREC - 1)
}

// Scanner constants from line_scanner.rs
const ZBAR_FIXED: i32 = 5;
const ROUND: c_uint = 1 << (ZBAR_FIXED - 1); // 16
const ZBAR_SCANNER_THRESH_FADE: u32 = 8;
const ZBAR_SCANNER_THRESH_MIN: c_uint = 4;
const EWMA_WEIGHT: c_uint = 25;
const THRESH_INIT: c_uint = 14;

#[derive(Default, Clone)]
pub(crate) struct zbar_symbol_set_t {
    pub(crate) symbols: Vec<zbar_symbol_t>,
}

/// image scanner state
pub(crate) struct zbar_image_scanner_t {
    /// Scanner state fields (formerly zbar_scanner_t)
    y1_min_thresh: c_uint,
    x: c_uint,
    y0: [i32; 4],
    y1_sign: i32,
    y1_thresh: c_uint,
    cur_edge: c_uint,
    last_edge: c_uint,
    width: c_uint,

    /// associated symbol decoder
    dcode: zbar_decoder_t,
    /// QR Code 2D reader
    qr: qr_reader,
    /// SQ Code 2D reader
    sq: SqReader,

    /// current scan direction
    dx: c_int,
    dy: c_int,
    du: c_int,
    umin: c_int,
    v: c_int,

    /// previous decode results
    syms: zbar_symbol_set_t,

    /// Type-safe scanner configuration
    config: ImageScannerConfig,
}

impl Default for zbar_image_scanner_t {
    /// Create a new image scanner
    ///
    /// Allocates and initializes a new image scanner instance with default configuration.
    ///
    /// # Returns
    /// Pointer to new scanner or null on allocation failure
    fn default() -> Self {
        Self {
            y1_min_thresh: ZBAR_SCANNER_THRESH_MIN,
            x: 0,
            y0: [0; 4],
            y1_sign: 0,
            y1_thresh: 0,
            cur_edge: 0,
            last_edge: 0,
            width: 0,
            dcode: zbar_decoder_t::default(),
            qr: qr_reader::default(),
            sq: SqReader::default(),
            dx: 0,
            dy: 0,
            du: 0,
            umin: 0,
            v: 0,
            syms: zbar_symbol_set_t::default(),
            config: ImageScannerConfig::default(),
        }
    }
}

impl zbar_image_scanner_t {
    /// Create a new image scanner with custom configuration
    ///
    /// This constructor accepts a type-safe `DecoderConfig` and creates a scanner
    /// configured according to those settings.
    ///
    /// # Arguments
    /// * `config` - The decoder configuration to use
    ///
    /// # Returns
    /// A new image scanner instance configured as specified
    pub(crate) fn with_config(config: DecoderConfig) -> Self {
        let decoder_state: DecoderState = (&config).into();

        let mut scanner_config = ImageScannerConfig {
            position_tracking: decoder_state.scanner.position_tracking,
            test_inverted: decoder_state.scanner.test_inverted,
            x_density: decoder_state.scanner.x_density,
            y_density: decoder_state.scanner.y_density,
            ean_composite: decoder_state.is_enabled(SymbolType::Composite),
            uncertainty: Default::default(),
        };

        // Sync uncertainty values from decoder state
        for sym in SymbolType::ALL.iter() {
            if let Some(sym_config) = decoder_state.get(*sym) {
                scanner_config
                    .uncertainty
                    .insert(*sym, sym_config.uncertainty);
            }
        }

        Self {
            y1_min_thresh: ZBAR_SCANNER_THRESH_MIN,
            x: 0,
            y0: [0; 4],
            y1_sign: 0,
            y1_thresh: 0,
            cur_edge: 0,
            last_edge: 0,
            width: 0,
            dcode: zbar_decoder_t::with_config(decoder_state),
            qr: qr_reader::default(),
            sq: SqReader::default(),
            dx: 0,
            dy: 0,
            du: 0,
            umin: 0,
            v: 0,
            syms: zbar_symbol_set_t::default(),
            config: scanner_config,
        }
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

    /// Get the width of the most recent bar or space
    pub(crate) fn width(&self) -> c_uint {
        self.width
    }

    /// Calculate the current threshold for edge detection
    ///
    /// Implements adaptive threshold calculation that slowly fades back to minimum.
    /// This helps with noise rejection while maintaining sensitivity.
    pub(crate) fn calc_thresh(&mut self) -> c_uint {
        // threshold 1st to improve noise rejection
        let thresh = self.y1_thresh;

        if thresh <= self.y1_min_thresh || self.width == 0 {
            return self.y1_min_thresh;
        }

        // slowly return threshold to min
        let dx = (self.x << ZBAR_FIXED) - self.last_edge;
        let mut t = thresh as u64 * dx as u64;
        t /= self.width as u64;
        t /= ZBAR_SCANNER_THRESH_FADE as u64;

        if thresh > t as c_uint {
            let new_thresh = thresh - t as c_uint;
            if new_thresh > self.y1_min_thresh {
                return new_thresh;
            }
        }

        self.y1_thresh = self.y1_min_thresh;
        self.y1_min_thresh
    }

    /// Flush the scanner state
    #[inline]
    pub(crate) fn scanner_flush(&mut self) -> SymbolType {
        if self.y1_sign == 0 {
            return SymbolType::None;
        }

        let x = (self.x << ZBAR_FIXED) + ROUND;

        if self.cur_edge != x || self.y1_sign > 0 {
            let edge = self.process_edge();
            self.cur_edge = x;
            self.y1_sign = -self.y1_sign;
            return edge;
        }

        self.y1_sign = 0;
        self.width = 0;
        unsafe { self.dcode.decode_width(0) }
    }

    /// Start a new scan
    pub(crate) fn new_scan(&mut self) -> SymbolType {
        let mut edge = SymbolType::None;

        while self.y1_sign != 0 {
            let tmp = self.scanner_flush();
            if tmp > edge {
                edge = tmp;
            }
        }

        // reset scanner and associated decoder
        self.x = 0;
        self.y0 = Default::default();
        self.y1_sign = 0;
        self.y1_thresh = self.y1_min_thresh;
        self.cur_edge = 0;
        self.last_edge = 0;
        self.width = 0;

        self.dcode.new_scan();
        edge
    }

    /// Process a single pixel intensity value
    pub(crate) fn scan_y(&mut self, y: c_int) -> SymbolType {
        // retrieve short value history
        let x = self.x;
        let mut y0_1 = self.y0[((x.wrapping_sub(1)) & 3) as usize];
        let mut y0_0 = y0_1;

        if x != 0 {
            // update weighted moving average
            y0_0 += ((y - y0_1) * EWMA_WEIGHT as i32) >> ZBAR_FIXED;
            self.y0[(x & 3) as usize] = y0_0;
        } else {
            y0_0 = y;
            y0_1 = y;
            self.y0[0] = y;
            self.y0[1] = y;
            self.y0[2] = y;
            self.y0[3] = y;
        }

        let y0_2 = self.y0[((x.wrapping_sub(2)) & 3) as usize];
        let y0_3 = self.y0[((x.wrapping_sub(3)) & 3) as usize];

        // 1st differential @ x-1
        let mut y1_1 = y0_1 - y0_2;
        {
            let y1_2 = y0_2 - y0_3;
            if y1_1.abs() < y1_2.abs() && ((y1_1 >= 0) == (y1_2 >= 0)) {
                y1_1 = y1_2;
            }
        }

        // 2nd differentials @ x-1 & x-2
        let y2_1 = y0_0 - (y0_1 * 2) + y0_2;
        let y2_2 = y0_1 - (y0_2 * 2) + y0_3;

        let mut edge = SymbolType::None;

        // 2nd zero-crossing is 1st local min/max - could be edge
        if (y2_1 == 0 || ((y2_1 > 0) == (y2_2 < 0))) && (self.calc_thresh() <= y1_1.unsigned_abs())
        {
            // check for 1st sign change
            let y1_rev = if self.y1_sign > 0 { y1_1 < 0 } else { y1_1 > 0 };

            if y1_rev {
                // intensity change reversal - finalize previous edge
                edge = self.process_edge();
            }

            if y1_rev || (self.y1_sign.abs() < y1_1.abs()) {
                self.y1_sign = y1_1;

                // adaptive thresholding
                // start at multiple of new min/max
                self.y1_thresh =
                    ((y1_1.unsigned_abs() * THRESH_INIT + ROUND) >> ZBAR_FIXED) as c_uint;
                if self.y1_thresh < self.y1_min_thresh {
                    self.y1_thresh = self.y1_min_thresh;
                }

                // update current edge
                let d = y2_1 - y2_2;
                self.cur_edge = 1 << ZBAR_FIXED;
                if d == 0 {
                    self.cur_edge >>= 1;
                } else if y2_1 != 0 {
                    // interpolate zero crossing
                    self.cur_edge -= (((y2_1 << ZBAR_FIXED) + 1) / d) as c_uint;
                }
                self.cur_edge += x << ZBAR_FIXED;
            }
        }

        self.x = x + 1;
        edge
    }

    /// Flush scanner pipeline and start new scan
    ///
    /// This function flushes the scanner pipeline twice and then starts a new scan.
    /// It's typically called at quiet borders to reset the scanner state.
    pub(crate) fn quiet_border(&mut self) {
        // Flush scanner pipeline twice
        self.scanner_flush();
        self.scanner_flush();

        // Start new scan
        self.new_scan();
    }

    /// Get the interpolated position of the last edge
    pub(crate) fn get_edge(&self, offset: c_uint, prec: c_int) -> c_uint {
        let edge = self
            .last_edge
            .wrapping_sub(offset)
            .wrapping_sub(1 << ZBAR_FIXED)
            .wrapping_sub(ROUND);
        let prec = ZBAR_FIXED - prec;

        match prec {
            1.. => edge >> prec,
            0 => edge,
            _ => edge << (-prec),
        }
    }

    /// Process an edge and pass the width to the decoder
    ///
    /// This function is called when an edge (transition) is detected.
    /// It calculates the width of the element and passes it to the decoder.
    fn process_edge(&mut self) -> SymbolType {
        if self.y1_sign == 0 {
            self.last_edge = (1 << ZBAR_FIXED) + ROUND;
            self.cur_edge = (1 << ZBAR_FIXED) + ROUND;
        } else if self.last_edge == 0 {
            self.last_edge = self.cur_edge;
        }

        self.width = self.cur_edge - self.last_edge;
        self.last_edge = self.cur_edge;

        // pass to decoder
        let width = self.width;
        unsafe { self.dcode.decode_width(width) }
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
    pub(crate) fn scan_image(&mut self, img: &mut zbar_image_t) -> Result<c_int> {
        let Some(syms) = self.scan_image_internal(img) else {
            return Err(Error::Unknown(-1));
        };

        let nsyms = syms.symbols.len();

        // Try inverted image if no symbols found and TEST_INVERTED is enabled
        if nsyms == 0 && self.config.test_inverted {
            if let Some(mut inv) = img.copy(true) {
                let _ = self.scan_image_internal(&mut inv);
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
    pub(crate) fn qr_handler(&mut self) {
        let line = decoder_get_qr_finder_line(&mut self.dcode);
        
        // Extract values we need before calling get_edge
        let pos0 = line.pos[0] as c_uint;
        let boffs_in = line.boffs as c_uint;
        let len_in = line.len as c_uint;
        let eoffs_in = line.eoffs as c_uint;
        
        // Now we can safely call get_edge
        let mut u = self.get_edge(pos0, QR_FINDER_SUBPREC);
        let boffs_edge = self.get_edge(boffs_in, QR_FINDER_SUBPREC);
        let len_edge = self.get_edge(len_in, QR_FINDER_SUBPREC);
        let eoffs_edge = self.get_edge(eoffs_in, QR_FINDER_SUBPREC);
        
        // Get line again to modify it
        let line = decoder_get_qr_finder_line(&mut self.dcode);
        line.boffs = (u as c_int) - boffs_edge as c_int;
        line.len = len_edge as c_int;
        line.eoffs = eoffs_edge as c_int - line.len;
        line.len -= u as c_int;

        u = (qr_fixed(self.umin, 0) as i64 + (self.du as i64) * (u as i64)) as c_uint;
        if self.du < 0 {
            std::mem::swap(&mut line.boffs, &mut line.eoffs);
            u = u.wrapping_sub(line.len as c_uint);
        }

        let vert: c_int = if self.dx != 0 { 0 } else { 1 };
        line.pos[vert as usize] = u as c_int;
        line.pos[(1 - vert) as usize] = qr_fixed(self.v, 1) as c_int;

        self.qr.found_line(vert, line);
    }

    /// SQ code handler - updates SQ reader configuration
    ///
    /// Gets the current SQ finder configuration from the decoder and passes it
    /// to the SQ reader for processing.
    pub(crate) fn sq_handler(&mut self) {
        self.sq
            .set_enabled(self.dcode.is_enabled(SymbolType::SqCode));
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
    fn scan_image_internal(&mut self, img: &mut zbar_image_t) -> Option<zbar_symbol_set_t> {
        // Set up decoder's back-pointer to this scanner for symbol callbacks
        let scanner_ptr = self as *mut _ as *mut c_void;
        self.dcode.set_userdata(scanner_ptr);

        self.qr.reset();
        self.sq.reset();

        // Clear previous symbols for new scan
        let scanner = &mut *self;
        scanner.syms.symbols.clear();

        // Share the symbol set with the image (we'll clone it at the end)
        // For now, we'll handle this differently - the image will get a clone of the results

        let w = img.width;
        let h = img.height;
        let data = img.data.as_slice();
        self.new_scan();

        // Horizontal scanning pass
        let density = self.config.y_density;
        if density > 0 {
            let mut p = 0;
            let mut x = 0i32;
            let mut y = 0i32;

            let mut border = ((h - 1) % (density as c_uint)).div_ceil(2);
            if border > h / 2 {
                border = h / 2;
            }
            assert!(border <= h);
            self.dy = 0;

            // movedelta(0, border)
            y += border as i32;
            p += (border * w) as isize;
            self.v = y;

            while (y as c_uint) < h {
                self.dx = 1;
                self.du = 1;
                self.umin = 0;
                while (x as c_uint) < w {
                    let d = data[p as usize];
                    x += 1;
                    p += 1;
                    self.scan_y(d as c_int);
                }
                self.quiet_border();

                // movedelta(-1, density)
                x -= 1;
                y += density as i32;
                p += -1 + (density as i32 * w as i32) as isize;
                self.v = y;
                if (y as c_uint) >= h {
                    break;
                }

                self.dx = -1;
                self.du = -1;
                self.umin = w as c_int;
                while x >= 0 {
                    let d = data[p as usize];
                    x -= 1;
                    p -= 1;
                    self.scan_y(d as c_int);
                }
                self.quiet_border();

                // movedelta(1, density)
                x += 1;
                y += density as i32;
                p += 1 + (density as i32 * w as i32) as isize;
                self.v = y;
            }
        }
        self.dx = 0;

        // Vertical scanning pass
        let density = self.config.x_density;
        if density > 0 {
            let mut p = 0;
            let mut x = 0i32;
            let mut y = 0i32;

            let mut border = ((w - 1) % (density as c_uint)).div_ceil(2);
            if border > w / 2 {
                border = w / 2;
            }
            assert!(border <= w);
            // movedelta(border, 0)
            x += border as i32;
            p += border as isize;
            self.v = x;

            while (x as c_uint) < w {
                self.dy = 1;
                self.du = 1;
                self.umin = 0;
                while (y as c_uint) < h {
                    let d = data[p as usize];
                    y += 1;
                    p += w as isize;
                    self.scan_y(d as c_int);
                }
                self.quiet_border();

                // movedelta(density, -1)
                x += density as i32;
                y -= 1;
                p += (density as isize) - (w as isize);
                self.v = x;
                if (x as c_uint) >= w {
                    break;
                }

                self.dy = -1;
                self.du = -1;
                self.umin = h as c_int;
                while y >= 0 {
                    let d = data[p as usize];
                    y -= 1;
                    p -= w as isize;
                    self.scan_y(d as c_int);
                }
                self.quiet_border();

                // movedelta(density, 1)
                x += density as i32;
                y += 1;
                p += (density as isize) + (w as isize);
                self.v = x;
            }
        }
        self.dy = 0;

        // Decode QR and SQ codes
        let raw_binary = self.dcode.is_binary_mode(SymbolType::QrCode);
        let qr_symbols = self.qr.decode(img, raw_binary);
        for sym in qr_symbols {
            self.add_symbol(sym);
        }

        self.sq_handler();
        if let Ok(Some(symbol)) = self.sq.decode(img) {
            self.add_symbol(symbol);
        }

        // Filter and merge EAN composite results
        let filter = density == 1 || self.config.y_density == 1;
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
        if nean == 1 && naddon == 1 && scanner.config.ean_composite {
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
    if (*iscn).config.position_tracking {
        let w = (*iscn).width();
        let u = (*iscn).umin + (*iscn).du * (*iscn).get_edge(w, 0) as c_int;
        if (*iscn).dx != 0 {
            x = u;
            y = (*iscn).v;
        } else {
            x = (*iscn).v;
            y = u;
        }
    }

    // Ignore partial results
    if symbol_type <= SymbolType::Partial {
        return;
    }

    let data = dcode.buffer_slice();
    if let Some(sym) = (*iscn).find_duplicate_symbol(symbol_type, data) {
        sym.quality += 1;
        if (*iscn).config.position_tracking {
            sym.add_point(x, y);
        }
        return;
    }

    // Allocate new symbol
    let mut sym = zbar_symbol_t::new(symbol_type);
    sym.modifiers = dcode.get_modifiers();

    // Copy data
    sym.data.extend_from_slice(data);

    // Initialize position
    if (*iscn).config.position_tracking {
        sym.add_point(x, y);
    }

    // Set orientation
    let dir = dcode.direction;
    if dir != 0 {
        sym.orient = (if (*iscn).dy != 0 { 1 } else { 0 }) + (((*iscn).du ^ dir) & 2);
    }

    (*iscn).add_symbol(sym);
}
