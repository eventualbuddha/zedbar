//! Decoder type definitions
//!
//! This module contains Rust definitions for all the barcode decoder types
//! that mirror the C struct layouts exactly for FFI compatibility.

use libc::{c_char, c_int, c_short, c_uint, c_void};
use std::ptr;

use crate::{
    decoder::{decoder_alloc_databar_segments, decoder_free_databar_segments},
    line_scanner::zbar_color_t,
};

/// Number of integer configs (ZBAR_CFG_MAX_LEN - ZBAR_CFG_MIN_LEN + 1)
pub const NUM_CFGS: usize = 2;

/// Window size for bar width history (must be power of 2)
pub const DECODE_WINDOW: usize = 16;

// Forward declare zbar_symbol_type_t as c_int
#[allow(non_camel_case_types)]
pub type zbar_symbol_type_t = c_int;

// Macro equivalents
#[inline]
unsafe fn cfg_set(configs: &mut [c_int; 2], cfg: c_int, val: c_int) {
    configs[(cfg - ZBAR_CFG_MIN_LEN) as usize] = val;
}

// ============================================================================
// Symbol type constants
// ============================================================================

pub const ZBAR_NONE: zbar_symbol_type_t = 0;
pub const ZBAR_PARTIAL: zbar_symbol_type_t = 1;
pub const ZBAR_EAN2: zbar_symbol_type_t = 2;
pub const ZBAR_EAN5: zbar_symbol_type_t = 5;
pub const ZBAR_EAN8: zbar_symbol_type_t = 8;
pub const ZBAR_UPCE: zbar_symbol_type_t = 9;
pub const ZBAR_ISBN10: zbar_symbol_type_t = 10;
pub const ZBAR_UPCA: zbar_symbol_type_t = 12;
pub const ZBAR_EAN13: zbar_symbol_type_t = 13;
pub const ZBAR_ISBN13: zbar_symbol_type_t = 14;
pub const ZBAR_COMPOSITE: zbar_symbol_type_t = 15;
pub const ZBAR_I25: zbar_symbol_type_t = 25;
pub const ZBAR_DATABAR: zbar_symbol_type_t = 34;
pub const ZBAR_DATABAR_EXP: zbar_symbol_type_t = 35;
pub const ZBAR_CODABAR: zbar_symbol_type_t = 38;
pub const ZBAR_CODE39: zbar_symbol_type_t = 39;
pub const ZBAR_QRCODE: zbar_symbol_type_t = 64;
pub const ZBAR_SQCODE: zbar_symbol_type_t = 80;
pub const ZBAR_CODE93: zbar_symbol_type_t = 93;
pub const ZBAR_CODE128: zbar_symbol_type_t = 128;

pub const ZBAR_SYMBOL: zbar_symbol_type_t = 0x00ff;

// ============================================================================
// Configuration constants
// ============================================================================

pub const ZBAR_CFG_ENABLE: c_int = 0;
pub const ZBAR_CFG_ADD_CHECK: c_int = 1;
pub const ZBAR_CFG_EMIT_CHECK: c_int = 2;
pub const ZBAR_CFG_ASCII: c_int = 3;
pub const ZBAR_CFG_BINARY: c_int = 4;
pub const ZBAR_CFG_NUM: c_int = 5;
pub const ZBAR_CFG_MIN_LEN: c_int = 0x20;
pub const ZBAR_CFG_MAX_LEN: c_int = 0x21;
pub const ZBAR_CFG_UNCERTAINTY: c_int = 64;
pub const ZBAR_CFG_POSITION: c_int = 128;
pub const ZBAR_CFG_TEST_INVERTED: c_int = 129;
pub const ZBAR_CFG_X_DENSITY: c_int = 256;
pub const ZBAR_CFG_Y_DENSITY: c_int = 257;

// ============================================================================
// Modifier constants
// ============================================================================

pub const ZBAR_MOD_GS1: c_int = 0;
pub const ZBAR_MOD_AIM: c_int = 1;

// ============================================================================
// Orientation constants
// ============================================================================

pub const ZBAR_ORIENT_UNKNOWN: c_int = -1;
pub const ZBAR_ORIENT_UP: c_int = 0;
pub const ZBAR_ORIENT_RIGHT: c_int = 1;
pub const ZBAR_ORIENT_DOWN: c_int = 2;
pub const ZBAR_ORIENT_LEFT: c_int = 3;

// ============================================================================
// Buffer size constants
// ============================================================================

pub const BUFFER_MIN: c_uint = 0x20;
pub const BUFFER_MAX: c_uint = 0x100;
pub const BUFFER_INCR: c_uint = 0x10;

// ============================================================================
// Simple decoder types
// ============================================================================

/// Interleaved 2 of 5 decoder state
#[derive(Default)]
#[allow(non_camel_case_types)]
pub struct i25_decoder_t {
    // Bitfields packed into first 32 bits:
    // direction: 1 bit, element: 4 bits, character: 12 bits = 17 bits used
    // We'll use a u32 and provide accessor methods
    bitfields: c_uint,
    pub s10: c_uint,
    pub width: c_uint,
    pub buf: [u8; 4],
    pub config: c_uint,
    pub configs: [c_int; NUM_CFGS],
}

impl i25_decoder_t {
    #[inline]
    pub fn direction(&self) -> bool {
        (self.bitfields & 0x1) != 0
    }

    #[inline]
    pub fn set_direction(&mut self, val: bool) {
        self.bitfields = (self.bitfields & !0x1) | (val as c_uint);
    }

    #[inline]
    pub fn element(&self) -> u8 {
        ((self.bitfields >> 1) & 0xF) as u8
    }

    #[inline]
    pub fn set_element(&mut self, val: u8) {
        self.bitfields = (self.bitfields & !(0xF << 1)) | ((val as c_uint & 0xF) << 1);
    }

    #[inline]
    pub fn character(&self) -> i16 {
        // Sign extend the 12-bit value
        let val = ((self.bitfields >> 5) & 0xFFF) as i16;
        // If the sign bit (bit 11) is set, extend it
        if val & 0x800 != 0 {
            val | !0xFFF
        } else {
            val
        }
    }

    #[inline]
    pub fn set_character(&mut self, val: i16) {
        self.bitfields = (self.bitfields & !(0xFFF << 5)) | (((val as c_uint) & 0xFFF) << 5);
    }

    /// Reset i25 decoder state
    pub(crate) unsafe fn reset(&mut self) {
        self.set_direction(false);
        self.set_element(0);
        self.set_character(-1);
        self.s10 = 0;
    }
}

/// Code 39 decoder state
#[derive(Default)]
#[allow(non_camel_case_types)]
pub struct code39_decoder_t {
    // Bitfields: direction: 1, element: 4, character: 12
    bitfields: c_uint,
    pub s9: c_uint,
    pub width: c_uint,
    pub config: c_uint,
    pub configs: [c_int; NUM_CFGS],
}

impl code39_decoder_t {
    #[inline]
    pub fn direction(&self) -> bool {
        (self.bitfields & 0x1) != 0
    }

    #[inline]
    pub fn set_direction(&mut self, val: bool) {
        self.bitfields = (self.bitfields & !0x1) | (val as c_uint);
    }

    #[inline]
    pub fn element(&self) -> u8 {
        ((self.bitfields >> 1) & 0xF) as u8
    }

    #[inline]
    pub fn set_element(&mut self, val: u8) {
        self.bitfields = (self.bitfields & !(0xF << 1)) | ((val as c_uint & 0xF) << 1);
    }

    #[inline]
    pub fn character(&self) -> i16 {
        // Sign extend the 12-bit value
        let val = ((self.bitfields >> 5) & 0xFFF) as i16;
        // If the sign bit (bit 11) is set, extend it
        if val & 0x800 != 0 {
            val | !0xFFF
        } else {
            val
        }
    }

    #[inline]
    pub fn set_character(&mut self, val: i16) {
        self.bitfields = (self.bitfields & !(0xFFF << 5)) | (((val as c_uint) & 0xFFF) << 5);
    }

    /// Reset code39 decoder state
    pub(crate) unsafe fn reset(&mut self) {
        self.set_direction(false);
        self.set_element(0);
        self.set_character(-1);
        self.s9 = 0;
    }
}

/// Code 93 decoder state
#[derive(Default)]
#[allow(non_camel_case_types)]
pub struct code93_decoder_t {
    // Bitfields: direction: 1, element: 3, character: 12
    bitfields: c_uint,
    pub width: c_uint,
    pub buf: u8,
    pub config: c_uint,
    pub configs: [c_int; NUM_CFGS],
}

impl code93_decoder_t {
    #[inline]
    pub fn direction(&self) -> bool {
        (self.bitfields & 0x1) != 0
    }

    #[inline]
    pub fn set_direction(&mut self, val: bool) {
        self.bitfields = (self.bitfields & !0x1) | (val as c_uint);
    }

    #[inline]
    pub fn element(&self) -> u8 {
        ((self.bitfields >> 1) & 0x7) as u8
    }

    #[inline]
    pub fn set_element(&mut self, val: u8) {
        self.bitfields = (self.bitfields & !(0x7 << 1)) | ((val as c_uint & 0x7) << 1);
    }

    #[inline]
    pub fn character(&self) -> i16 {
        // Sign extend the 12-bit value
        let val = ((self.bitfields >> 4) & 0xFFF) as i16;
        // If the sign bit (bit 11) is set, extend it
        if val & 0x800 != 0 {
            val | !0xFFF
        } else {
            val
        }
    }

    #[inline]
    pub fn set_character(&mut self, val: i16) {
        self.bitfields = (self.bitfields & !(0xFFF << 4)) | (((val as c_uint) & 0xFFF) << 4);
    }

    /// Reset code93 decoder state
    pub(crate) unsafe fn reset(&mut self) {
        self.set_direction(false);
        self.set_element(0);
        self.set_character(-1);
    }
}

/// Codabar decoder state
#[derive(Default)]
#[allow(non_camel_case_types)]
pub struct codabar_decoder_t {
    // Bitfields: direction: 1, element: 4, character: 12
    bitfields: c_uint,
    pub s7: c_uint,
    pub width: c_uint,
    pub buf: [u8; 6],
    pub config: c_uint,
    pub configs: [c_int; NUM_CFGS],
}

impl codabar_decoder_t {
    #[inline]
    pub fn direction(&self) -> bool {
        (self.bitfields & 0x1) != 0
    }

    #[inline]
    pub fn set_direction(&mut self, val: bool) {
        self.bitfields = (self.bitfields & !0x1) | (val as c_uint);
    }

    #[inline]
    pub fn element(&self) -> u8 {
        ((self.bitfields >> 1) & 0xF) as u8
    }

    #[inline]
    pub fn set_element(&mut self, val: u8) {
        self.bitfields = (self.bitfields & !(0xF << 1)) | ((val as c_uint & 0xF) << 1);
    }

    #[inline]
    pub fn character(&self) -> i16 {
        // Sign extend the 12-bit value
        let val = ((self.bitfields >> 5) & 0xFFF) as i16;
        // If the sign bit (bit 11) is set, extend it
        if val & 0x800 != 0 {
            val | !0xFFF
        } else {
            val
        }
    }

    #[inline]
    pub fn set_character(&mut self, val: i16) {
        self.bitfields = (self.bitfields & !(0xFFF << 5)) | (((val as c_uint) & 0xFFF) << 5);
    }

    /// Reset codabar decoder state
    pub(crate) unsafe fn reset(&mut self) {
        self.set_direction(false);
        self.set_element(0);
        self.set_character(-1);
        self.s7 = 0;
    }
}

/// Code 128 decoder state
#[derive(Default)]
#[allow(non_camel_case_types)]
pub struct code128_decoder_t {
    // Bitfields: direction: 1, element: 3, character: 12 (16 bits)
    // start: 8 bits - packed into same u32
    // Total: 24 bits used in first u32
    bitfields_and_start: c_uint,
    pub s6: c_uint,
    pub width: c_uint,
    pub config: c_uint,
    pub configs: [c_int; NUM_CFGS],
}

impl code128_decoder_t {
    #[inline]
    pub fn direction(&self) -> u8 {
        (self.bitfields_and_start & 0x1) as u8
    }

    #[inline]
    pub fn set_direction(&mut self, val: u8) {
        self.bitfields_and_start = (self.bitfields_and_start & !0x1) | (val as c_uint & 0x1);
    }

    #[inline]
    pub fn element(&self) -> u8 {
        ((self.bitfields_and_start >> 1) & 0x7) as u8
    }

    #[inline]
    pub fn set_element(&mut self, val: u8) {
        self.bitfields_and_start =
            (self.bitfields_and_start & !(0x7 << 1)) | ((val as c_uint & 0x7) << 1);
    }

    #[inline]
    pub fn character(&self) -> i16 {
        // Sign extend the 12-bit value
        let val = ((self.bitfields_and_start >> 4) & 0xFFF) as i16;
        // If the sign bit (bit 11) is set, extend it
        if val & 0x800 != 0 {
            val | !0xFFF
        } else {
            val
        }
    }

    #[inline]
    pub fn set_character(&mut self, val: i16) {
        self.bitfields_and_start =
            (self.bitfields_and_start & !(0xFFF << 4)) | (((val as c_uint) & 0xFFF) << 4);
    }

    #[inline]
    pub fn start(&self) -> u8 {
        ((self.bitfields_and_start >> 16) & 0xFF) as u8
    }

    #[inline]
    pub fn set_start(&mut self, val: u8) {
        self.bitfields_and_start =
            (self.bitfields_and_start & !(0xFF << 16)) | ((val as c_uint) << 16);
    }

    /// Reset code128 decoder state
    pub(crate) unsafe fn reset(&mut self) {
        self.set_direction(0);
        self.set_element(0);
        self.set_character(-1);
        self.s6 = 0;
    }
}

/// SQ Code finder state (already defined but include here for completeness)
#[derive(Default)]
#[allow(non_camel_case_types)]
pub struct sq_finder_t {
    pub config: c_uint,
}

// ============================================================================
// Complex decoder types
// ============================================================================

/// DataBar segment (partial)
#[allow(non_camel_case_types)]
pub struct databar_segment_t {
    // First 32 bits of bitfields:
    // finder: 5, exp: 1, color: 1, side: 1,
    // partial: 1, count: 7, epoch: 8, check: 8 = 32 bits
    bitfields: c_uint,
    pub data: c_short,
    pub width: c_short,
}

impl databar_segment_t {
    #[inline]
    pub fn finder(&self) -> i8 {
        // finder is a signed 5-bit field (bits 0-4)
        let val = (self.bitfields & 0x1F) as i8;
        // Sign extend from 5 bits to 8 bits
        if val & 0x10 != 0 {
            val | 0xE0_u8 as i8
        } else {
            val
        }
    }

    #[inline]
    pub fn set_finder(&mut self, val: i8) {
        self.bitfields = (self.bitfields & !0x1F) | ((val as c_uint) & 0x1F);
    }

    #[inline]
    pub fn partial(&self) -> bool {
        // partial is bit 8
        (self.bitfields & (1 << 8)) != 0
    }

    #[inline]
    pub fn set_partial(&mut self, val: bool) {
        if val {
            self.bitfields |= 1 << 8;
        } else {
            self.bitfields &= !(1 << 8);
        }
    }

    #[inline]
    pub fn exp(&self) -> bool {
        // exp is bit 5
        (self.bitfields & (1 << 5)) != 0
    }

    #[inline]
    pub fn set_exp(&mut self, val: bool) {
        if val {
            self.bitfields |= 1 << 5;
        } else {
            self.bitfields &= !(1 << 5);
        }
    }

    #[inline]
    pub fn color(&self) -> zbar_color_t {
        // color is bit 6
        (((self.bitfields >> 6) & 1) as u8).into()
    }

    #[inline]
    pub fn set_color(&mut self, val: zbar_color_t) {
        self.bitfields = (self.bitfields & !(1 << 6)) | ((val as c_uint) << 6);
    }

    #[inline]
    pub fn side(&self) -> u8 {
        // side is bit 7
        ((self.bitfields >> 7) & 1) as u8
    }

    #[inline]
    pub fn set_side(&mut self, val: u8) {
        self.bitfields = (self.bitfields & !(1 << 7)) | (((val & 1) as c_uint) << 7);
    }

    #[inline]
    pub fn count(&self) -> u8 {
        // count is bits 9-15 (7 bits)
        ((self.bitfields >> 9) & 0x7F) as u8
    }

    #[inline]
    pub fn set_count(&mut self, val: u8) {
        self.bitfields = (self.bitfields & !(0x7F << 9)) | (((val & 0x7F) as c_uint) << 9);
    }

    #[inline]
    pub fn epoch(&self) -> u8 {
        // epoch is bits 16-23 (8 bits)
        ((self.bitfields >> 16) & 0xFF) as u8
    }

    #[inline]
    pub fn set_epoch(&mut self, val: u8) {
        self.bitfields = (self.bitfields & !(0xFF << 16)) | ((val as c_uint) << 16);
    }

    #[inline]
    pub fn check(&self) -> u8 {
        // check is bits 24-31 (8 bits)
        ((self.bitfields >> 24) & 0xFF) as u8
    }

    #[inline]
    pub fn set_check(&mut self, val: u8) {
        self.bitfields = (self.bitfields & !(0xFF << 24)) | ((val as c_uint) << 24);
    }
}

/// DataBar decoder state
#[derive(Default)]
#[allow(non_camel_case_types)]
pub struct databar_decoder_t {
    pub config: c_uint,
    pub config_exp: c_uint,
    // Bitfields: csegs: 8, epoch: 8 (16 bits in a u32)
    pub bitfields: c_uint,
    pub segs: *mut databar_segment_t,
    pub chars: [c_char; 16],
}

impl databar_decoder_t {
    #[inline]
    pub fn csegs(&self) -> u8 {
        (self.bitfields & 0xFF) as u8
    }

    #[inline]
    pub fn set_csegs(&mut self, val: u8) {
        self.bitfields = (self.bitfields & !0xFF) | (val as c_uint);
    }

    #[inline]
    pub fn epoch(&self) -> u8 {
        ((self.bitfields >> 8) & 0xFF) as u8
    }

    #[inline]
    pub fn set_epoch(&mut self, val: u8) {
        self.bitfields = (self.bitfields & !(0xFF << 8)) | ((val as c_uint) << 8);
    }

    /// Reset DataBar decoder state
    pub(crate) unsafe fn reset(&mut self) {
        let n = self.csegs() as isize;
        self.new_scan();
        for i in 0..n {
            let seg = (self.segs).offset(i);
            (*seg).set_finder(-1);
        }
    }

    /// Prepare DataBar decoder for new scan
    pub(crate) unsafe fn new_scan(&mut self) {
        for i in 0..16 {
            if self.chars[i] >= 0 {
                let seg = (self.segs).offset(self.chars[i] as isize);
                if (*seg).partial() {
                    (*seg).set_finder(-1);
                }
                self.chars[i] = -1;
            }
        }
    }
}

/// EAN pass state
#[derive(Default)]
#[allow(non_camel_case_types)]
pub struct ean_pass_t {
    pub state: c_char,
    pub width: c_uint,
    pub raw: [u8; 7],
}

/// EAN/UPC decoder state
#[derive(Default)]
#[allow(non_camel_case_types)]
pub struct ean_decoder_t {
    pub pass: [ean_pass_t; 4],
    pub left: zbar_symbol_type_t,
    pub right: zbar_symbol_type_t,
    pub direction: c_int,
    pub s4: c_uint,
    pub width: c_uint,
    pub buf: [c_char; 18],
    pub enable: c_char,
    pub ean13_config: c_uint,
    pub ean8_config: c_uint,
    pub upca_config: c_uint,
    pub upce_config: c_uint,
    pub isbn10_config: c_uint,
    pub isbn13_config: c_uint,
    pub ean5_config: c_uint,
    pub ean2_config: c_uint,
}

impl ean_decoder_t {
    /// Prepare EAN decoder for new scan
    pub(crate) unsafe fn new_scan(&mut self) {
        self.pass[0].state = -1;
        self.pass[1].state = -1;
        self.pass[2].state = -1;
        self.pass[3].state = -1;
        self.s4 = 0;
    }

    /// Reset EAN decoder state
    pub(crate) unsafe fn reset(&mut self) {
        self.new_scan();
        self.left = 0; // ZBAR_NONE
        self.right = 0; // ZBAR_NONE
    }
}

/// QR finder line (from qrcode.h)
#[derive(Default, Copy, Clone)]
#[allow(non_camel_case_types)]
pub struct qr_finder_line {
    pub pos: [c_int; 2], // qr_point
    pub len: c_int,
    pub boffs: c_int,
    pub eoffs: c_int,
}

/// QR Code finder state
#[derive(Default)]
#[allow(non_camel_case_types)]
pub struct qr_finder_t {
    pub s5: c_uint,
    pub line: qr_finder_line,
    pub config: c_uint,
}

impl qr_finder_t {
    /// Reset QR finder state
    pub(crate) fn reset(&mut self) {
        self.s5 = 0;
    }
}

// ============================================================================
// Main decoder struct
// ============================================================================

/// Decoder handler callback
#[allow(non_camel_case_types)]
pub type zbar_decoder_handler_t = unsafe fn(*mut zbar_decoder_t);

/// Main barcode decoder state
#[derive(Default)]
#[allow(non_camel_case_types)]
pub struct zbar_decoder_t {
    // Basic decoder state
    pub idx: u8,
    pub w: [c_uint; DECODE_WINDOW],
    pub type_: zbar_symbol_type_t,
    pub lock: zbar_symbol_type_t,
    pub modifiers: c_uint,
    pub direction: c_int,
    pub s6: c_uint,

    // Buffer management (everything above here is reset)
    buffer: Vec<u8>,
    pub userdata: *mut c_void,
    pub handler: Option<zbar_decoder_handler_t>,

    // Symbology-specific decoders
    pub ean: ean_decoder_t,
    pub i25: i25_decoder_t,
    pub databar: databar_decoder_t,
    pub codabar: codabar_decoder_t,
    pub code39: code39_decoder_t,
    pub code93: code93_decoder_t,
    pub code128: code128_decoder_t,
    pub qrf: qr_finder_t,
    pub sqf: sq_finder_t,
}

impl zbar_decoder_t {
    /// Create a new decoder instance
    pub unsafe fn new() -> Option<Self> {
        let mut decoder = Self::default();

        decoder.buffer = Vec::with_capacity(BUFFER_MIN as usize);

        // Initialize default configs
        decoder.ean.enable = 1;
        decoder.ean.ean13_config = (1 << ZBAR_CFG_ENABLE) | (1 << ZBAR_CFG_EMIT_CHECK);
        decoder.ean.ean8_config = (1 << ZBAR_CFG_ENABLE) | (1 << ZBAR_CFG_EMIT_CHECK);
        decoder.ean.upca_config = 1 << ZBAR_CFG_EMIT_CHECK;
        decoder.ean.upce_config = 1 << ZBAR_CFG_EMIT_CHECK;
        decoder.ean.isbn10_config = 1 << ZBAR_CFG_EMIT_CHECK;
        decoder.ean.isbn13_config = 1 << ZBAR_CFG_EMIT_CHECK;
        // FIXME_ADDON_SYNC not defined, skip ean2/ean5 config

        decoder.i25.config = 1 << ZBAR_CFG_ENABLE;
        cfg_set(&mut decoder.i25.configs, ZBAR_CFG_MIN_LEN, 6);

        decoder.databar.config = (1 << ZBAR_CFG_ENABLE) | (1 << ZBAR_CFG_EMIT_CHECK);
        decoder.databar.config_exp = (1 << ZBAR_CFG_ENABLE) | (1 << ZBAR_CFG_EMIT_CHECK);
        decoder.databar.set_csegs(4);
        decoder.databar.segs = decoder_alloc_databar_segments(4);
        if decoder.databar.segs.is_null() {
            return None;
        }

        decoder.codabar.config = 1 << ZBAR_CFG_ENABLE;
        cfg_set(&mut decoder.codabar.configs, ZBAR_CFG_MIN_LEN, 4);

        decoder.code39.config = 1 << ZBAR_CFG_ENABLE;
        cfg_set(&mut decoder.code39.configs, ZBAR_CFG_MIN_LEN, 1);

        decoder.code93.config = 1 << ZBAR_CFG_ENABLE;
        decoder.code128.config = 1 << ZBAR_CFG_ENABLE;
        decoder.qrf.config = 1 << ZBAR_CFG_ENABLE;
        decoder.sqf.config = 1 << ZBAR_CFG_ENABLE;

        decoder.reset();
        Some(decoder)
    }

    pub(crate) fn color(&self) -> zbar_color_t {
        self.idx.into()
    }

    /// Get width of a specific element from the decoder's history window
    pub(crate) unsafe fn get_width(&self, offset: u8) -> c_uint {
        self.w[((self.idx as usize).wrapping_sub(offset as usize)) & (DECODE_WINDOW - 1)]
    }

    /// Get the combined width of two consecutive elements
    pub(crate) unsafe fn pair_width(&self, offset: u8) -> c_uint {
        self.get_width(offset) + self.get_width(offset + 1)
    }

    /// Calculate sum of n consecutive element widths
    pub(crate) unsafe fn calc_s(&self, mut offset: u8, mut n: u8) -> c_uint {
        let mut s = 0;
        while n > 0 {
            s += self.get_width(offset);
            offset += 1;
            n -= 1;
        }
        s
    }

    /// Resize the decoder's data buffer if needed.
    ///
    /// # Errors
    /// Returns `Err` on allocation failure or if max size exceeded, `Ok` on success.
    pub(crate) unsafe fn set_buffer_capacity(&mut self, len: c_uint) -> Result<(), ()> {
        if len <= BUFFER_MIN {
            return Ok(());
        }
        let current_alloc = self.buffer_capacity();
        if len < current_alloc {
            // FIXME: size reduction heuristic?
            return Ok(());
        }
        if len > BUFFER_MAX {
            return Err(());
        }

        self.buffer.reserve_exact((len - current_alloc) as usize);
        Ok(())
    }

    #[inline]
    pub(crate) fn buffer_capacity(&self) -> c_uint {
        self.buffer.capacity() as c_uint
    }

    #[inline]
    pub(crate) fn buffer_len(&self) -> c_uint {
        self.buffer.len() as c_uint
    }

    #[inline]
    pub(crate) fn set_buffer_len(&mut self, len: c_uint) {
        debug_assert!(len <= self.buffer_capacity());
        let len = len as usize;
        let current = self.buffer.len();
        if len <= current {
            self.buffer.truncate(len);
        } else {
            unsafe {
                self.buffer.set_len(len);
            }
        }
    }

    #[inline]
    pub(crate) fn buffer_ptr(&self) -> *const c_char {
        self.buffer.as_ptr() as *const c_char
    }

    #[inline]
    pub(crate) fn buffer_mut_ptr(&mut self) -> *mut c_char {
        self.buffer.as_mut_ptr() as *mut c_char
    }

    /// Acquire a decoder lock for a specific symbology type
    /// Returns 1 if already locked, 0 if lock acquired
    pub(crate) unsafe fn _zbar_decoder_acquire_lock(&mut self, req: zbar_symbol_type_t) -> c_char {
        if self.lock != 0 {
            return 1;
        }
        self.lock = req;
        0
    }

    /// Release a decoder lock
    /// Returns 0 on success
    pub(crate) unsafe fn _zbar_decoder_release_lock(&mut self, req: zbar_symbol_type_t) -> c_char {
        debug_assert_eq!(self.lock, req, "lock={} req={}", self.lock, req);
        self.lock = 0;
        0
    }

    /// Reset decoder to initial state
    pub unsafe fn reset(&mut self) {
        self.idx = 0;
        self.w.fill(0);
        self.type_ = ZBAR_NONE;
        self.lock = ZBAR_NONE;
        self.modifiers = 0;
        self.direction = 0;
        self.s6 = 0;
        self.set_buffer_len(0);

        self.ean.reset();
        self.i25.reset();
        self.databar.reset();
        self.codabar.reset();
        self.code39.reset();
        self.code93.reset();
        self.code128.reset();
        self.qrf.reset();
    }

    /// Mark start of a new scan pass
    ///
    /// Clears any intra-symbol state and resets color to ZBAR_SPACE.
    /// Any partially decoded symbol state is retained.
    pub(crate) unsafe fn new_scan(&mut self) {
        // Soft reset decoder
        self.w.fill(0);
        self.lock = ZBAR_NONE;
        self.idx = 0;
        self.s6 = 0;

        self.ean.new_scan();
        self.i25.reset();
        self.databar.reset();
        self.codabar.reset();
        self.code39.reset();
        self.code93.reset();
        self.code128.reset();
        self.qrf.reset();
    }
}

impl Drop for zbar_decoder_t {
    fn drop(&mut self) {
        unsafe {
            if !self.databar.segs.is_null() {
                decoder_free_databar_segments(self.databar.segs);
                self.databar.segs = ptr::null_mut();
            }
        }
    }
}
