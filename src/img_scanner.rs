use std::mem::swap;


use crate::{
    color::Color,
    config::{internal::DecoderState, DecoderConfig},
    decoder::DECODE_WINDOW,
    image_ffi::zbar_image_t,
    img_scanner_config::ImageScannerConfig,
    symbol::{Orientation, Symbol},
    Result, SymbolType,
};

#[cfg(feature = "codabar")]
use crate::decoder::codabar_decoder_t;
#[cfg(feature = "code128")]
use crate::decoder::code128_decoder_t;
#[cfg(feature = "code39")]
use crate::decoder::code39_decoder_t;
#[cfg(feature = "code93")]
use crate::decoder::code93_decoder_t;
#[cfg(feature = "databar")]
use crate::decoder::databar_decoder_t;
#[cfg(feature = "i25")]
use crate::decoder::i25_decoder_t;
#[cfg(feature = "qrcode")]
use crate::decoder::qr_finder_t;

#[cfg(feature = "qrcode")]
use crate::finder::find_qr;
#[cfg(feature = "qrcode")]
use crate::qrcode::qrdec::qr_reader;
#[cfg(feature = "sqcode")]
use crate::sqcode::SqReader;

#[cfg(feature = "codabar")]
use crate::decoders::codabar::_zbar_decode_codabar;
#[cfg(feature = "code128")]
use crate::decoders::code128::_zbar_decode_code128;
#[cfg(feature = "code39")]
use crate::decoders::code39::_zbar_decode_code39;
#[cfg(feature = "code93")]
use crate::decoders::code93::_zbar_decode_code93;
#[cfg(feature = "databar")]
use crate::decoders::databar::_zbar_decode_databar;
#[cfg(feature = "ean")]
use crate::decoders::ean::{ean_decoder_t, zbar_decode_ean};
#[cfg(feature = "i25")]
use crate::decoders::i25::_zbar_decode_i25;

// QR Code finder precision constant
const QR_FINDER_SUBPREC: i32 = 2;

// QR_FIXED macro: ((((v) << 1) + (rnd)) << (QR_FINDER_SUBPREC - 1))
fn qr_fixed(v: i32, rnd: i32) -> u32 {
    (((v as u32) << 1) + (rnd as u32)) << (QR_FINDER_SUBPREC - 1)
}

// Scanner constants from line_scanner.rs
const ZBAR_FIXED: i32 = 5;
const ROUND: u32 = 1 << (ZBAR_FIXED - 1); // 16
const ZBAR_SCANNER_THRESH_FADE: u32 = 8;
const ZBAR_SCANNER_THRESH_MIN: u32 = 4;
const EWMA_WEIGHT: u32 = 25;
const THRESH_INIT: u32 = 14;

// Decoder constants from decoder.rs
const BUFFER_MIN: usize = 0x20;
const BUFFER_MAX: usize = 0x100;

/// image scanner state
pub(crate) struct zbar_image_scanner_t {
    /// Scanner state fields (formerly zbar_scanner_t)
    y1_min_thresh: u32,
    x: u32,
    y0: [i32; 4],
    y1_sign: i32,
    y1_thresh: u32,
    cur_edge: u32,
    last_edge: u32,
    width: u32,

    /// Decoder state fields (formerly zbar_decoder_t)
    pub(crate) idx: u8,
    pub(crate) w: [u32; DECODE_WINDOW],
    pub(crate) type_: SymbolType,
    pub(crate) lock: SymbolType,
    pub(crate) modifiers: u32,
    pub(crate) direction: i32,
    pub(crate) s6: u32,

    // Buffer management
    buffer: Vec<u8>,

    // Decoder configuration
    pub(crate) config: DecoderState,

    // Symbology-specific decoders
    #[cfg(feature = "ean")]
    pub(crate) ean: ean_decoder_t,
    #[cfg(feature = "i25")]
    pub(crate) i25: i25_decoder_t,
    #[cfg(feature = "databar")]
    pub(crate) databar: databar_decoder_t,
    #[cfg(feature = "codabar")]
    pub(crate) codabar: codabar_decoder_t,
    #[cfg(feature = "code39")]
    pub(crate) code39: code39_decoder_t,
    #[cfg(feature = "code93")]
    pub(crate) code93: code93_decoder_t,
    #[cfg(feature = "code128")]
    pub(crate) code128: code128_decoder_t,
    #[cfg(feature = "qrcode")]
    pub(crate) qrf: qr_finder_t,

    /// QR Code 2D reader
    #[cfg(feature = "qrcode")]
    qr: qr_reader,
    /// SQ Code 2D reader
    #[cfg(feature = "sqcode")]
    sq: SqReader,

    /// current scan direction
    dx: i32,
    dy: i32,
    du: i32,
    umin: i32,
    v: i32,

    /// previous decode results
    syms: Vec<Symbol>,

    /// Type-safe scanner configuration
    scanner_config: ImageScannerConfig,
}

impl Default for zbar_image_scanner_t {
    /// Create a new image scanner
    ///
    /// Allocates and initializes a new image scanner instance with default configuration.
    ///
    /// # Returns
    /// Pointer to new scanner or null on allocation failure
    fn default() -> Self {
        let mut scanner = Self {
            // Scanner fields
            y1_min_thresh: ZBAR_SCANNER_THRESH_MIN,
            x: 0,
            y0: [0; 4],
            y1_sign: 0,
            y1_thresh: 0,
            cur_edge: 0,
            last_edge: 0,
            width: 0,

            // Decoder fields
            idx: 0,
            w: Default::default(),
            type_: SymbolType::default(),
            lock: SymbolType::default(),
            modifiers: 0,
            direction: 0,
            s6: 0,
            buffer: Vec::with_capacity(BUFFER_MIN),
            config: DecoderState::default(),
            #[cfg(feature = "ean")]
            ean: ean_decoder_t::default(),
            #[cfg(feature = "i25")]
            i25: i25_decoder_t::default(),
            #[cfg(feature = "databar")]
            databar: databar_decoder_t::default(),
            #[cfg(feature = "codabar")]
            codabar: codabar_decoder_t::default(),
            #[cfg(feature = "code39")]
            code39: code39_decoder_t::default(),
            #[cfg(feature = "code93")]
            code93: code93_decoder_t::default(),
            #[cfg(feature = "code128")]
            code128: code128_decoder_t::default(),
            #[cfg(feature = "qrcode")]
            qrf: qr_finder_t::default(),

            // Scanner-specific fields
            #[cfg(feature = "qrcode")]
            qr: qr_reader::default(),
            #[cfg(feature = "sqcode")]
            sq: SqReader::default(),
            dx: 0,
            dy: 0,
            du: 0,
            umin: 0,
            v: 0,
            syms: vec![],
            scanner_config: ImageScannerConfig::default(),
        };

        scanner.sync_config_to_decoders();
        scanner.decoder_reset();

        scanner
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
            upscale_small_images: decoder_state.scanner.upscale_small_images,
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

        let mut scanner = Self {
            config: decoder_state,
            scanner_config,
            ..Default::default()
        };

        scanner.sync_config_to_decoders();
        scanner.decoder_reset();

        scanner
    }

    /// Add a symbol to the scanner's symbol set
    ///
    /// # Arguments
    /// * `sym` - The symbol to add
    pub(crate) fn add_symbol(&mut self, sym: Symbol) {
        self.syms.push(sym);
    }

    pub(crate) fn find_duplicate_symbol(
        &mut self,
        symbol_type: SymbolType,
        data: &[u8],
    ) -> Option<&mut Symbol> {
        self.syms
            .iter_mut()
            .find(|sym| sym.symbol_type() == symbol_type && sym.data == data)
    }

    /// Get the width of the most recent bar or space
    pub(crate) fn width(&self) -> u32 {
        self.width
    }

    /// Calculate the current threshold for edge detection
    ///
    /// Implements adaptive threshold calculation that slowly fades back to minimum.
    /// This helps with noise rejection while maintaining sensitivity.
    pub(crate) fn calc_thresh(&mut self) -> u32 {
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

        if thresh > t as u32 {
            let new_thresh = thresh - t as u32;
            if new_thresh > self.y1_min_thresh {
                return new_thresh;
            }
        }

        self.y1_thresh = self.y1_min_thresh;
        self.y1_min_thresh
    }

    /// Flush the scanner state
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
        self.decode_width(0)
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

        self.decoder_new_scan();
        edge
    }

    /// Process a single pixel intensity value
    pub(crate) fn scan_y(&mut self, y: i32) -> SymbolType {
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
                    (y1_1.unsigned_abs() * THRESH_INIT + ROUND) >> ZBAR_FIXED;
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
                    self.cur_edge -= (((y2_1 << ZBAR_FIXED) + 1) / d) as u32;
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
    pub(crate) fn get_edge(&self, offset: u32, prec: i32) -> u32 {
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
        self.decode_width(width)
    }

    // ========================================================================
    // Decoder methods (formerly zbar_decoder_t)
    // ========================================================================

    /// Sync EAN enable flag from configuration
    fn sync_config_to_decoders(&mut self) {
        #[cfg(feature = "ean")]
        {
            self.ean.enable = self.config.ean_enabled();
        }
    }

    pub(crate) fn color(&self) -> Color {
        self.idx.into()
    }

    /// Get width of a specific element from the decoder's history window
    pub(crate) fn get_width(&self, offset: u8) -> u32 {
        self.w[((self.idx as usize).wrapping_sub(offset as usize)) & (DECODE_WINDOW - 1)]
    }

    /// Get the combined width of two consecutive elements
    pub(crate) fn pair_width(&self, offset: u8) -> u32 {
        self.get_width(offset) + self.get_width(offset + 1)
    }

    /// Calculate sum of n consecutive element widths
    pub(crate) fn calc_s(&self, mut offset: u8, mut n: u8) -> u32 {
        let mut s = 0;
        while n > 0 {
            s += self.get_width(offset);
            offset += 1;
            n -= 1;
        }
        s
    }

    /// Resize the decoder's data buffer if needed
    pub(crate) fn set_buffer_capacity(&mut self, len: usize) -> Result<(), ()> {
        if len <= BUFFER_MIN {
            return Ok(());
        }
        let current_alloc = self.buffer_capacity();
        if len < current_alloc {
            return Ok(());
        }
        if len > BUFFER_MAX {
            return Err(());
        }

        self.buffer.reserve_exact(len - current_alloc);
        Ok(())
    }

    pub(crate) fn buffer_capacity(&self) -> usize {
        self.buffer.capacity()
    }

    pub(crate) fn set_buffer_len(&mut self, len: usize) {
        debug_assert!(len <= self.buffer_capacity());
        self.buffer.resize(len, 0);
    }

    pub(crate) fn buffer_mut_slice(&mut self, len: usize) -> Result<&mut [u8], ()> {
        if len > BUFFER_MAX {
            return Err(());
        }

        self.buffer.resize(len, 0);
        Ok(&mut self.buffer[..len])
    }

    pub(crate) fn buffer_slice(&self) -> &[u8] {
        &self.buffer
    }

    pub(crate) fn write_buffer_byte(&mut self, pos: usize, value: u8) -> Result<(), ()> {
        self.buffer_mut_slice(pos + 1)?[pos] = value;
        Ok(())
    }

    pub(crate) fn truncate_buffer(&mut self, len: usize) {
        self.buffer.truncate(len);
    }

    pub(crate) fn is_enabled(&self, sym: SymbolType) -> bool {
        self.config.is_enabled(sym)
    }

    pub(crate) fn should_emit_checksum(&self, sym: SymbolType) -> bool {
        self.config
            .get(sym)
            .map(|c| c.checksum.emit_check)
            .unwrap_or(false)
    }

    pub(crate) fn should_validate_checksum(&self, sym: SymbolType) -> bool {
        self.config
            .get(sym)
            .map(|c| c.checksum.add_check)
            .unwrap_or(false)
    }

    pub(crate) fn get_length_limits(&self, sym: SymbolType) -> Option<(u32, u32)> {
        self.config
            .get(sym)
            .and_then(|c| c.length_limits)
            .map(|l| (l.min, l.max))
    }

    pub(crate) fn is_binary_mode(&self, sym: SymbolType) -> bool {
        self.config.get(sym).map(|c| c.binary_mode).unwrap_or(false)
    }

    /// Reset decoder to initial state
    pub(crate) fn decoder_reset(&mut self) {
        self.idx = 0;
        self.w.fill(0);
        self.type_ = SymbolType::None;
        self.lock = SymbolType::None;
        self.modifiers = 0;
        self.direction = 0;
        self.s6 = 0;
        self.set_buffer_len(0);

        #[cfg(feature = "ean")]
        self.ean.reset();
        #[cfg(feature = "i25")]
        self.i25.reset();
        #[cfg(feature = "databar")]
        self.databar.reset();
        #[cfg(feature = "codabar")]
        self.codabar.reset();
        #[cfg(feature = "code39")]
        self.code39.reset();
        #[cfg(feature = "code93")]
        self.code93.reset();
        #[cfg(feature = "code128")]
        self.code128.reset();
        #[cfg(feature = "qrcode")]
        self.qrf.reset();
    }

    /// Mark start of a new scan pass (decoder-specific)
    pub(crate) fn decoder_new_scan(&mut self) {
        self.w.fill(0);
        self.lock = SymbolType::None;
        self.idx = 0;
        self.s6 = 0;

        #[cfg(feature = "ean")]
        self.ean.new_scan();
        #[cfg(feature = "i25")]
        self.i25.reset();
        #[cfg(feature = "databar")]
        self.databar.reset();
        #[cfg(feature = "codabar")]
        self.codabar.reset();
        #[cfg(feature = "code39")]
        self.code39.reset();
        #[cfg(feature = "code93")]
        self.code93.reset();
        #[cfg(feature = "code128")]
        self.code128.reset();
        #[cfg(feature = "qrcode")]
        self.qrf.reset();
    }

    pub(crate) fn acquire_lock(&mut self, req: SymbolType) -> bool {
        if self.lock != SymbolType::None {
            return false;
        }
        self.lock = req;
        true
    }

    pub(crate) fn release_lock(&mut self, req: SymbolType) -> bool {
        if self.lock != req {
            return false;
        }
        self.lock = SymbolType::None;
        true
    }

    /// Process next bar/space width from input stream
    pub(crate) fn decode_width(&mut self, w: u32) -> SymbolType {
        let mut sym = SymbolType::None;

        // Store width in circular buffer
        self.w[(self.idx & (DECODE_WINDOW - 1) as u8) as usize] = w;

        // Update shared character width
        self.s6 = self.s6.wrapping_sub(self.get_width(7));
        self.s6 = self.s6.wrapping_add(self.get_width(1));

        // Each decoder processes width stream in parallel
        #[cfg(feature = "qrcode")]
        if self.is_enabled(SymbolType::QrCode) {
            let tmp = find_qr(self);
            if tmp > SymbolType::Partial {
                sym = tmp;
            }
        }

        #[cfg(feature = "ean")]
        if self.ean.enable {
            let tmp = zbar_decode_ean(&mut *self);
            if tmp != SymbolType::None {
                sym = tmp;
            }
        }

        #[cfg(feature = "code39")]
        if self.is_enabled(SymbolType::Code39) {
            let tmp = _zbar_decode_code39(&mut *self);
            if tmp > SymbolType::Partial {
                sym = tmp;
            }
        }

        #[cfg(feature = "code93")]
        if self.is_enabled(SymbolType::Code93) {
            let tmp = _zbar_decode_code93(&mut *self);
            if tmp > SymbolType::Partial {
                sym = tmp;
            }
        }

        #[cfg(feature = "code128")]
        if self.is_enabled(SymbolType::Code128) {
            let tmp = _zbar_decode_code128(&mut *self);
            if tmp > SymbolType::Partial {
                sym = tmp;
            }
        }

        #[cfg(feature = "databar")]
        if self.is_enabled(SymbolType::Databar) || self.is_enabled(SymbolType::DatabarExp) {
            let tmp = _zbar_decode_databar(&mut *self);
            if tmp > SymbolType::Partial {
                sym = tmp;
            }
        }

        #[cfg(feature = "codabar")]
        if self.is_enabled(SymbolType::Codabar) {
            let tmp = _zbar_decode_codabar(&mut *self);
            if tmp > SymbolType::Partial {
                sym = tmp;
            }
        }

        #[cfg(feature = "i25")]
        if self.is_enabled(SymbolType::I25) {
            let tmp = _zbar_decode_i25(self);
            if tmp > SymbolType::Partial {
                sym = tmp;
            }
        }

        self.idx = self.idx.wrapping_add(1);
        self.type_ = sym;

        if sym != SymbolType::None {
            if self.lock != SymbolType::None
                && sym > SymbolType::Partial
                && sym != SymbolType::QrCode
            {
                self.lock = SymbolType::None;
            }

            if sym > SymbolType::Partial {
                self.symbol_handler();
            }
        }

        sym
    }

    pub(crate) fn get_type(&self) -> SymbolType {
        self.type_
    }

    pub(crate) fn get_modifiers(&self) -> u32 {
        self.modifiers
    }

    /// Scans an image for barcodes, with optional inverted image retry.
    pub(crate) fn scan_image(&mut self, img: &mut zbar_image_t) -> Vec<Symbol> {
        let symbols = self.scan_image_internal(img);

        // Try inverted image if no symbols found and TEST_INVERTED is enabled
        if symbols.is_empty()
            && self.scanner_config.test_inverted
            && let Some(mut inv) = img.copy(true)
        {
            let inverted_symbols = self.scan_image_internal(&mut inv);
            if !inverted_symbols.is_empty() {
                return inverted_symbols;
            }
        }

        // Try upscaling small images for QR code detection.
        // Small QR codes (< 200px in either dimension) often have modules that are
        // only 2-3 pixels wide, which is too small for reliable finder pattern detection.
        // Upscaling to 4x improves detection by giving modules more pixels.
        #[cfg(feature = "qrcode")]
        if symbols.is_empty()
            && self.scanner_config.upscale_small_images
            && self.is_enabled(SymbolType::QrCode)
            && (img.width < 200 || img.height < 200)
            && let Some(mut upscaled) = img.upscale(4)
        {
            let upscaled_symbols = self.scan_image_internal(&mut upscaled);
            if !upscaled_symbols.is_empty() {
                return upscaled_symbols;
            }

            // Also try inverted + upscaled
            if self.scanner_config.test_inverted
                && let Some(mut inv_upscaled) = img.copy(true).and_then(|i| i.upscale(4))
            {
                return self.scan_image_internal(&mut inv_upscaled);
            }
        }

        symbols
    }

    /// Handle QR code finder line detection
    ///
    /// Processes a QR code finder line from the decoder and forwards it to the QR reader.
    /// Adjusts edge positions based on scanner state and transforms coordinates.
    #[cfg(feature = "qrcode")]
    pub(crate) fn qr_handler(&mut self) {
        let line = &mut self.qrf.line;

        // Extract values we need before calling get_edge
        let pos0 = line.pos[0] as u32;
        let boffs_in = line.boffs as u32;
        let len_in = line.len as u32;
        let eoffs_in = line.eoffs as u32;

        // Now we can safely call get_edge
        let mut u = self.get_edge(pos0, QR_FINDER_SUBPREC);
        let boffs_edge = self.get_edge(boffs_in, QR_FINDER_SUBPREC);
        let len_edge = self.get_edge(len_in, QR_FINDER_SUBPREC);
        let eoffs_edge = self.get_edge(eoffs_in, QR_FINDER_SUBPREC);

        // Get line again to modify it
        let line = &mut self.qrf.line;
        line.boffs = (u as i32) - boffs_edge as i32;
        line.len = len_edge as i32;
        line.eoffs = eoffs_edge as i32 - line.len;
        line.len -= u as i32;

        u = (qr_fixed(self.umin, 0) as i64 + (self.du as i64) * (u as i64)) as u32;
        if self.du < 0 {
            std::mem::swap(&mut line.boffs, &mut line.eoffs);
            u = u.wrapping_sub(line.len as u32);
        }

        let vert: i32 = if self.dx != 0 { 0 } else { 1 };
        line.pos[vert as usize] = u as i32;
        line.pos[(1 - vert) as usize] = qr_fixed(self.v, 1) as i32;

        self.qr.found_line(vert, line);
    }

    /// SQ code handler - updates SQ reader configuration
    ///
    /// Gets the current SQ finder configuration from the decoder and passes it
    /// to the SQ reader for processing.
    #[cfg(feature = "sqcode")]
    pub(crate) fn sq_handler(&mut self) {
        self.sq.set_enabled(self.is_enabled(SymbolType::SqCode));
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
    fn scan_image_internal(&mut self, img: &mut zbar_image_t) -> Vec<Symbol> {
        #[cfg(feature = "qrcode")]
        self.qr.reset();
        #[cfg(feature = "sqcode")]
        self.sq.reset();

        // Clear previous symbols for new scan
        self.syms.clear();

        // Share the symbol set with the image (we'll clone it at the end)
        // For now, we'll handle this differently - the image will get a clone of the results

        let w = img.width;
        let h = img.height;
        let data = img.data.as_slice();
        self.new_scan();

        // Horizontal scanning pass
        let density = self.scanner_config.y_density;
        if density > 0 {
            let mut p = 0;
            let mut x = 0i32;
            let mut y = 0i32;

            let mut border = ((h - 1) % density).div_ceil(2);
            if border > h / 2 {
                border = h / 2;
            }
            assert!(border <= h);
            self.dy = 0;

            // movedelta(0, border)
            y += border as i32;
            p += (border * w) as isize;
            self.v = y;

            while (y as u32) < h {
                self.dx = 1;
                self.du = 1;
                self.umin = 0;
                while (x as u32) < w {
                    let d = data[p as usize];
                    x += 1;
                    p += 1;
                    self.scan_y(d as i32);
                }
                self.quiet_border();

                // movedelta(-1, density)
                x -= 1;
                y += density as i32;
                p += -1 + (density as i32 * w as i32) as isize;
                self.v = y;
                if (y as u32) >= h {
                    break;
                }

                self.dx = -1;
                self.du = -1;
                self.umin = w as i32;
                while x >= 0 {
                    let d = data[p as usize];
                    x -= 1;
                    p -= 1;
                    self.scan_y(d as i32);
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
        let density = self.scanner_config.x_density;
        if density > 0 {
            let mut p = 0;
            let mut x = 0i32;
            let mut y = 0i32;

            let mut border = ((w - 1) % density).div_ceil(2);
            if border > w / 2 {
                border = w / 2;
            }
            assert!(border <= w);
            // movedelta(border, 0)
            x += border as i32;
            p += border as isize;
            self.v = x;

            while (x as u32) < w {
                self.dy = 1;
                self.du = 1;
                self.umin = 0;
                while (y as u32) < h {
                    let d = data[p as usize];
                    y += 1;
                    p += w as isize;
                    self.scan_y(d as i32);
                }
                self.quiet_border();

                // movedelta(density, -1)
                x += density as i32;
                y -= 1;
                p += (density as isize) - (w as isize);
                self.v = x;
                if (x as u32) >= w {
                    break;
                }

                self.dy = -1;
                self.du = -1;
                self.umin = h as i32;
                while y >= 0 {
                    let d = data[p as usize];
                    y -= 1;
                    p -= w as isize;
                    self.scan_y(d as i32);
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
        #[cfg(feature = "qrcode")]
        {
            let raw_binary = self.is_binary_mode(SymbolType::QrCode);
            let qr_symbols = self.qr.decode(img, raw_binary);
            for sym in qr_symbols {
                self.add_symbol(sym);
            }
        }

        #[cfg(feature = "sqcode")]
        {
            self.sq_handler();
            if let Ok(Some(symbol)) = self.sq.decode(img) {
                self.add_symbol(symbol);
            }
        }

        // Filter and merge EAN composite results
        let filter = density == 1 || self.scanner_config.y_density == 1;
        let scanner = &mut *self;

        // Filter low-quality symbols
        scanner.syms.retain(|sym| {
            let sym_type = sym.symbol_type();
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
        for sym in &scanner.syms {
            let sym_type = sym.symbol_type();
            if sym_type < SymbolType::Composite && sym_type != SymbolType::Isbn10 {
                if sym_type > SymbolType::Ean5 {
                    nean += 1;
                } else if sym_type > SymbolType::Partial {
                    naddon += 1;
                }
            }
        }

        // Merge EAN composite if we have exactly one EAN and one add-on
        if nean == 1 && naddon == 1 && scanner.scanner_config.ean_composite {
            // Extract EAN and add-on symbols
            let mut ean_sym: Option<Symbol> = None;
            let mut addon_sym: Option<Symbol> = None;
            let mut other_syms = Vec::new();

            for sym in scanner.syms.drain(..) {
                let sym_type = sym.symbol_type();
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
                other_syms.push(Symbol::composite(ean, addon));
            }

            scanner.syms = other_syms;
        }

        let mut symbols = vec![];
        swap(&mut symbols, &mut scanner.syms);
        symbols
    }

    /// Symbol handler callback for 1D barcode decoding
    ///
    /// This function is called by the decoder when a barcode is successfully decoded.
    /// It manages symbol deduplication, position tracking, and adds symbols to the result set.
    ///
    /// # Arguments
    /// * `dcode` - The decoder that found the symbol
    pub(crate) fn symbol_handler(&mut self) {
        let symbol_type = self.get_type();
        let mut x = 0;
        let mut y = 0;

        // QR codes are handled separately
        #[cfg(feature = "qrcode")]
        if symbol_type == SymbolType::QrCode {
            self.qr_handler();
            return;
        }

        #[cfg(not(feature = "qrcode"))]
        if symbol_type == SymbolType::QrCode {
            return;
        }

        // Calculate position if position tracking is enabled
        let position_tracking = self.scanner_config.position_tracking;
        if position_tracking {
            let w = self.width();
            let u = self.umin + self.du * self.get_edge(w, 0) as i32;
            if self.dx != 0 {
                x = u;
                y = self.v;
            } else {
                x = self.v;
                y = u;
            }
        }

        // Ignore partial results
        if symbol_type <= SymbolType::Partial {
            return;
        }

        // Check for duplicate - we need to copy the data first to avoid borrowing issues
        let data_vec = self.buffer_slice().to_vec();

        if let Some(sym) = self.find_duplicate_symbol(symbol_type, &data_vec) {
            sym.quality += 1;
            if position_tracking {
                sym.add_point(x, y);
            }
            return;
        }

        // Allocate new symbol
        let mut sym = Symbol::new(symbol_type);
        sym.modifiers = self.get_modifiers();

        // Copy data
        sym.data.extend_from_slice(&data_vec);

        // Initialize position
        if position_tracking {
            sym.add_point(x, y);
        }

        // Set orientation
        let dir = self.direction;
        if dir != 0 {
            let base = if self.dy != 0 { 1 } else { 0 };
            let offset = (self.du ^ dir) & 2;
            sym.orientation = match base + offset {
                0 => Orientation::Up,
                1 => Orientation::Down,
                2 => Orientation::Right,
                3 => Orientation::Left,
                _ => Orientation::Unknown,
            };
        }

        self.add_symbol(sym);
    }
}
