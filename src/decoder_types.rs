//! Decoder type definitions
//!
//! This module contains Rust definitions for all the barcode decoder types
//! that mirror the C struct layouts exactly for FFI compatibility.

use libc::{c_char, c_int, c_short, c_uint, c_void};

/// Number of integer configs (ZBAR_CFG_MAX_LEN - ZBAR_CFG_MIN_LEN + 1)
pub const NUM_CFGS: usize = 2;

/// Window size for bar width history (must be power of 2)
pub const DECODE_WINDOW: usize = 16;

// Forward declare zbar_symbol_type_t as c_int
#[allow(non_camel_case_types)]
pub type zbar_symbol_type_t = c_int;

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
pub const ZBAR_CFG_X_DENSITY: c_int = 256;
pub const ZBAR_CFG_Y_DENSITY: c_int = 257;

// ============================================================================
// Color constants
// ============================================================================

pub const ZBAR_SPACE: u8 = 0;
pub const ZBAR_BAR: u8 = 1;

// ============================================================================
// Modifier constants
// ============================================================================

pub const ZBAR_MOD_GS1: c_int = 0;
pub const ZBAR_MOD_AIM: c_int = 1;

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
#[repr(C)]
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
}

/// Code 39 decoder state
#[repr(C)]
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
}

/// Code 93 decoder state
#[repr(C)]
#[allow(non_camel_case_types)]
pub struct code93_decoder_t {
    // Bitfields: direction: 1, element: 3, character: 12
    bitfields: c_uint,
    pub width: c_uint,
    pub buf: u8,
    _padding: [u8; 3],
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
}

/// Codabar decoder state
#[repr(C)]
#[allow(non_camel_case_types)]
pub struct codabar_decoder_t {
    // Bitfields: direction: 1, element: 4, character: 12
    bitfields: c_uint,
    pub s7: c_uint,
    pub width: c_uint,
    pub buf: [u8; 6],
    _padding: [u8; 2],
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
}

/// Code 128 decoder state
#[repr(C)]
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
}

/// SQ Code finder state (already defined but include here for completeness)
#[repr(C)]
#[allow(non_camel_case_types)]
pub struct sq_finder_t {
    pub config: c_uint,
}

// ============================================================================
// Complex decoder types
// ============================================================================

/// DataBar segment (partial)
#[repr(C)]
#[allow(non_camel_case_types)]
pub struct databar_segment_t {
    // First 32 bits of bitfields:
    // finder: 5, exp: 1, color: 1, side: 1,
    // partial: 1, count: 7, epoch: 8, check: 8 = 32 bits
    bitfields: c_uint,
    pub data: c_short,
    pub width: c_short,
}

/// DataBar decoder state
#[repr(C)]
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
}

/// EAN pass state
#[repr(C)]
#[allow(non_camel_case_types)]
pub struct ean_pass_t {
    pub state: c_char,
    _padding: [u8; 3],
    pub width: c_uint,
    pub raw: [u8; 7],
    _padding2: u8,
}

/// EAN/UPC decoder state
#[repr(C)]
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
    _padding: u8, // Align to 4 bytes
    pub ean13_config: c_uint,
    pub ean8_config: c_uint,
    pub upca_config: c_uint,
    pub upce_config: c_uint,
    pub isbn10_config: c_uint,
    pub isbn13_config: c_uint,
    pub ean5_config: c_uint,
    pub ean2_config: c_uint,
}

/// QR finder line (from qrcode.h)
#[repr(C)]
#[allow(non_camel_case_types)]
pub struct qr_finder_line {
    pub pos: [c_int; 2], // qr_point
    pub len: c_int,
    pub boffs: c_int,
    pub eoffs: c_int,
}

/// QR Code finder state
#[repr(C)]
#[allow(non_camel_case_types)]
pub struct qr_finder_t {
    pub s5: c_uint,
    pub line: qr_finder_line,
    pub config: c_uint,
}

// ============================================================================
// Main decoder struct
// ============================================================================

/// Decoder handler callback
#[allow(non_camel_case_types)]
pub type zbar_decoder_handler_t = unsafe extern "C" fn(*mut zbar_decoder_t);

/// Main barcode decoder state
#[repr(C)]
#[allow(non_camel_case_types)]
pub struct zbar_decoder_t {
    // Basic decoder state
    pub idx: u8,
    _padding1: [u8; 3], // Align to 4 bytes
    pub w: [c_uint; DECODE_WINDOW],
    pub type_: zbar_symbol_type_t,
    pub lock: zbar_symbol_type_t,
    pub modifiers: c_uint,
    pub direction: c_int,
    pub s6: c_uint,

    // Buffer management (everything above here is reset)
    pub buf_alloc: c_uint,
    pub buflen: c_uint,
    pub buf: *mut c_char,
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
