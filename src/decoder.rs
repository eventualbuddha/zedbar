//! Decoder type definitions
//!
//! This module contains Rust definitions for all the barcode decoder types
//! that mirror the C struct layouts exactly for FFI compatibility.

use std::{ffi::c_void, ptr::null_mut};

use libc::{c_char, c_int, c_short, c_uint};

use crate::{
    config::internal::DecoderState,
    decoders::{
        codabar::_zbar_decode_codabar,
        code128::_zbar_decode_code128,
        code39::_zbar_decode_code39,
        code93::_zbar_decode_code93,
        databar::_zbar_decode_databar,
        ean::{ean_decoder_t, zbar_decode_ean},
        i25::_zbar_decode_i25,
    },
    finder::find_qr,
    img_scanner::symbol_handler,
    line_scanner::zbar_color_t,
    Result, SymbolType,
};

/// Window size for bar width history (must be power of 2)
pub(crate) const DECODE_WINDOW: usize = 16;

// Assertion macro
macro_rules! zassert {
    ($condition:expr, $retval:expr, $($arg:tt)*) => {
        if !$condition {
            return $retval;
        }
    };
}

/// Decode element width into a discrete value
/// Returns -1 if the element width is invalid
pub(crate) fn _zbar_decoder_decode_e(e: c_uint, s: c_uint, n: c_uint) -> c_int {
    let big_e = ((e * n * 2 + 1) / s).wrapping_sub(3) / 2;
    if big_e >= n - 3 {
        -1
    } else {
        big_e as c_int
    }
}

// ============================================================================
// Configuration parameters
// ============================================================================

// ============================================================================
// Modifier constants
// ============================================================================

pub(crate) const ZBAR_MOD_GS1: c_int = 0;
pub(crate) const ZBAR_MOD_AIM: c_int = 1;

// ============================================================================
// Orientation constants
// ============================================================================

pub(crate) const ZBAR_ORIENT_UNKNOWN: c_int = -1;

// ============================================================================
// Buffer size constants
// ============================================================================

pub(crate) const BUFFER_MIN: c_uint = 0x20;
pub(crate) const BUFFER_MAX: c_uint = 0x100;

// ============================================================================
// Simple decoder types
// ============================================================================

/// Interleaved 2 of 5 decoder state
#[derive(Default)]
pub(crate) struct i25_decoder_t {
    // Bitfields packed into first 32 bits:
    // direction: 1 bit, element: 4 bits, character: 12 bits = 17 bits used
    // We'll use a u32 and provide accessor methods
    bitfields: c_uint,
    pub(crate) s10: c_uint,
    pub(crate) width: c_uint,
    pub(crate) buffer: Vec<u8>,
}

impl i25_decoder_t {
    #[inline]
    pub(crate) fn direction(&self) -> bool {
        (self.bitfields & 0x1) != 0
    }

    #[inline]
    pub(crate) fn set_direction(&mut self, val: bool) {
        self.bitfields = (self.bitfields & !0x1) | (val as c_uint);
    }

    #[inline]
    pub(crate) fn element(&self) -> u8 {
        ((self.bitfields >> 1) & 0xF) as u8
    }

    #[inline]
    pub(crate) fn set_element(&mut self, val: u8) {
        self.bitfields = (self.bitfields & !(0xF << 1)) | ((val as c_uint & 0xF) << 1);
    }

    #[inline]
    pub(crate) fn character(&self) -> i16 {
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
    pub(crate) fn set_character(&mut self, val: i16) {
        self.bitfields = (self.bitfields & !(0xFFF << 5)) | (((val as c_uint) & 0xFFF) << 5);
    }

    /// Reset i25 decoder state
    pub(crate) fn reset(&mut self) {
        self.set_direction(false);
        self.set_element(0);
        self.set_character(-1);
        self.s10 = 0;
    }

    pub(crate) fn set_byte(&mut self, index: usize, value: u8) {
        if self.buffer.len() <= index {
            self.buffer.resize(index + 1, 0);
        }
        self.buffer[index] = value;
    }
}

/// Code 39 decoder state
#[derive(Default)]
pub(crate) struct code39_decoder_t {
    // Bitfields: direction: 1, element: 4, character: 12
    bitfields: c_uint,
    pub s9: c_uint,
    pub width: c_uint,
    pub buffer: Vec<u8>,
}

impl code39_decoder_t {
    #[inline]
    pub(crate) fn direction(&self) -> bool {
        (self.bitfields & 0x1) != 0
    }

    #[inline]
    pub(crate) fn set_direction(&mut self, val: bool) {
        self.bitfields = (self.bitfields & !0x1) | (val as c_uint);
    }

    #[inline]
    pub(crate) fn element(&self) -> u8 {
        ((self.bitfields >> 1) & 0xF) as u8
    }

    #[inline]
    pub(crate) fn set_element(&mut self, val: u8) {
        self.bitfields = (self.bitfields & !(0xF << 1)) | ((val as c_uint & 0xF) << 1);
    }

    #[inline]
    pub(crate) fn character(&self) -> i16 {
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
    pub(crate) fn set_character(&mut self, val: i16) {
        self.bitfields = (self.bitfields & !(0xFFF << 5)) | (((val as c_uint) & 0xFFF) << 5);
    }

    /// Reset code39 decoder state
    pub(crate) fn reset(&mut self) {
        self.set_direction(false);
        self.set_element(0);
        self.set_character(-1);
        self.s9 = 0;
    }

    pub(crate) fn set_byte(&mut self, index: usize, value: u8) {
        if self.buffer.len() <= index {
            self.buffer.resize(index + 1, 0);
        }
        self.buffer[index] = value;
    }
}

/// Code 93 decoder state
#[derive(Default)]
pub(crate) struct code93_decoder_t {
    // Bitfields: direction: 1, element: 3, character: 12
    bitfields: c_uint,
    pub(crate) width: c_uint,
    pub(crate) buf: u8,
}

impl code93_decoder_t {
    #[inline]
    pub(crate) fn direction(&self) -> bool {
        (self.bitfields & 0x1) != 0
    }

    #[inline]
    pub(crate) fn set_direction(&mut self, val: bool) {
        self.bitfields = (self.bitfields & !0x1) | (val as c_uint);
    }

    #[inline]
    pub(crate) fn element(&self) -> u8 {
        ((self.bitfields >> 1) & 0x7) as u8
    }

    #[inline]
    pub(crate) fn set_element(&mut self, val: u8) {
        self.bitfields = (self.bitfields & !(0x7 << 1)) | ((val as c_uint & 0x7) << 1);
    }

    #[inline]
    pub(crate) fn character(&self) -> i16 {
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
    pub(crate) fn set_character(&mut self, val: i16) {
        self.bitfields = (self.bitfields & !(0xFFF << 4)) | (((val as c_uint) & 0xFFF) << 4);
    }

    /// Reset code93 decoder state
    pub(crate) fn reset(&mut self) {
        self.set_direction(false);
        self.set_element(0);
        self.set_character(-1);
    }
}

/// Codabar decoder state
#[derive(Default)]
pub(crate) struct codabar_decoder_t {
    // Bitfields: direction: 1, element: 4, character: 12
    bitfields: c_uint,
    pub(crate) s7: c_uint,
    pub(crate) width: c_uint,
    pub(crate) buf: [u8; 6],
}

impl codabar_decoder_t {
    #[inline]
    pub(crate) fn direction(&self) -> bool {
        (self.bitfields & 0x1) != 0
    }

    #[inline]
    pub(crate) fn set_direction(&mut self, val: bool) {
        self.bitfields = (self.bitfields & !0x1) | (val as c_uint);
    }

    #[inline]
    pub(crate) fn element(&self) -> u8 {
        ((self.bitfields >> 1) & 0xF) as u8
    }

    #[inline]
    pub(crate) fn set_element(&mut self, val: u8) {
        self.bitfields = (self.bitfields & !(0xF << 1)) | ((val as c_uint & 0xF) << 1);
    }

    #[inline]
    pub(crate) fn character(&self) -> i16 {
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
    pub(crate) fn set_character(&mut self, val: i16) {
        self.bitfields = (self.bitfields & !(0xFFF << 5)) | (((val as c_uint) & 0xFFF) << 5);
    }

    /// Reset codabar decoder state
    pub(crate) fn reset(&mut self) {
        self.set_direction(false);
        self.set_element(0);
        self.set_character(-1);
        self.s7 = 0;
    }
}

/// Code 128 decoder state
#[derive(Default)]
pub(crate) struct code128_decoder_t {
    // Bitfields: direction: 1, element: 3, character: 12 (16 bits)
    // start: 8 bits - packed into same u32
    // Total: 24 bits used in first u32
    bitfields_and_start: c_uint,
    pub(crate) s6: c_uint,
    pub(crate) width: c_uint,
}

impl code128_decoder_t {
    #[inline]
    pub(crate) fn direction(&self) -> u8 {
        (self.bitfields_and_start & 0x1) as u8
    }

    #[inline]
    pub(crate) fn set_direction(&mut self, val: u8) {
        self.bitfields_and_start = (self.bitfields_and_start & !0x1) | (val as c_uint & 0x1);
    }

    #[inline]
    pub(crate) fn element(&self) -> u8 {
        ((self.bitfields_and_start >> 1) & 0x7) as u8
    }

    #[inline]
    pub(crate) fn set_element(&mut self, val: u8) {
        self.bitfields_and_start =
            (self.bitfields_and_start & !(0x7 << 1)) | ((val as c_uint & 0x7) << 1);
    }

    #[inline]
    pub(crate) fn character(&self) -> i16 {
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
    pub(crate) fn set_character(&mut self, val: i16) {
        self.bitfields_and_start =
            (self.bitfields_and_start & !(0xFFF << 4)) | (((val as c_uint) & 0xFFF) << 4);
    }

    #[inline]
    pub(crate) fn start(&self) -> u8 {
        ((self.bitfields_and_start >> 16) & 0xFF) as u8
    }

    #[inline]
    pub(crate) fn set_start(&mut self, val: u8) {
        self.bitfields_and_start =
            (self.bitfields_and_start & !(0xFF << 16)) | ((val as c_uint) << 16);
    }

    /// Reset code128 decoder state
    pub(crate) fn reset(&mut self) {
        self.set_direction(0);
        self.set_element(0);
        self.set_character(-1);
        self.s6 = 0;
    }
}

/// SQ Code finder state (already defined but include here for completeness)
#[derive(Default)]
pub(crate) struct sq_finder_t {}

// ============================================================================
// Complex decoder types
// ============================================================================

/// DataBar segment (partial)
#[derive(Clone)]
pub(crate) struct databar_segment_t {
    // First 32 bits of bitfields:
    // finder: 5, exp: 1, color: 1, side: 1,
    // partial: 1, count: 7, epoch: 8, check: 8 = 32 bits
    bitfields: c_uint,
    pub(crate) data: c_short,
    pub(crate) width: c_short,
}

impl databar_segment_t {
    #[inline]
    pub(crate) fn finder(&self) -> i8 {
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
    pub(crate) fn set_finder(&mut self, val: i8) {
        self.bitfields = (self.bitfields & !0x1F) | ((val as c_uint) & 0x1F);
    }

    #[inline]
    pub(crate) fn partial(&self) -> bool {
        // partial is bit 8
        (self.bitfields & (1 << 8)) != 0
    }

    #[inline]
    pub(crate) fn set_partial(&mut self, val: bool) {
        if val {
            self.bitfields |= 1 << 8;
        } else {
            self.bitfields &= !(1 << 8);
        }
    }

    #[inline]
    pub(crate) fn exp(&self) -> bool {
        // exp is bit 5
        (self.bitfields & (1 << 5)) != 0
    }

    #[inline]
    pub(crate) fn set_exp(&mut self, val: bool) {
        if val {
            self.bitfields |= 1 << 5;
        } else {
            self.bitfields &= !(1 << 5);
        }
    }

    #[inline]
    pub(crate) fn color(&self) -> zbar_color_t {
        // color is bit 6
        (((self.bitfields >> 6) & 1) as u8).into()
    }

    #[inline]
    pub(crate) fn set_color(&mut self, val: zbar_color_t) {
        self.bitfields = (self.bitfields & !(1 << 6)) | ((val as c_uint) << 6);
    }

    #[inline]
    pub(crate) fn side(&self) -> u8 {
        // side is bit 7
        ((self.bitfields >> 7) & 1) as u8
    }

    #[inline]
    pub(crate) fn set_side(&mut self, val: u8) {
        self.bitfields = (self.bitfields & !(1 << 7)) | (((val & 1) as c_uint) << 7);
    }

    #[inline]
    pub(crate) fn count(&self) -> u8 {
        // count is bits 9-15 (7 bits)
        ((self.bitfields >> 9) & 0x7F) as u8
    }

    #[inline]
    pub(crate) fn set_count(&mut self, val: u8) {
        self.bitfields = (self.bitfields & !(0x7F << 9)) | (((val & 0x7F) as c_uint) << 9);
    }

    #[inline]
    pub(crate) fn epoch(&self) -> u8 {
        // epoch is bits 16-23 (8 bits)
        ((self.bitfields >> 16) & 0xFF) as u8
    }

    #[inline]
    pub(crate) fn set_epoch(&mut self, val: u8) {
        self.bitfields = (self.bitfields & !(0xFF << 16)) | ((val as c_uint) << 16);
    }

    #[inline]
    pub(crate) fn check(&self) -> u8 {
        // check is bits 24-31 (8 bits)
        ((self.bitfields >> 24) & 0xFF) as u8
    }

    #[inline]
    pub(crate) fn set_check(&mut self, val: u8) {
        self.bitfields = (self.bitfields & !(0xFF << 24)) | ((val as c_uint) << 24);
    }

    #[inline]
    pub(crate) fn segment_index(&self) -> i32 {
        ((self.finder() as i32) << 2)
            | ((self.color() as i32) << 1)
            | (((self.color() as u8 ^ self.side()) as i32) & 1)
    }
}

impl Default for databar_segment_t {
    fn default() -> Self {
        let mut seg = Self {
            bitfields: 0,
            data: 0,
            width: 0,
        };
        seg.set_finder(-1);
        seg
    }
}

/// DataBar decoder state
pub(crate) struct databar_decoder_t {
    epoch: u8,
    pub(crate) segs: Vec<databar_segment_t>,
    chars: [c_char; 16],
}

impl Default for databar_decoder_t {
    fn default() -> Self {
        Self {
            epoch: 0,
            segs: Vec::new(),
            chars: [-1; 16],
        }
    }
}

impl databar_decoder_t {
    pub(crate) fn csegs(&self) -> usize {
        self.segs.len()
    }

    pub(crate) fn seg(&self, index: usize) -> &databar_segment_t {
        &self.segs[index]
    }

    pub(crate) fn seg_mut(&mut self, index: usize) -> &mut databar_segment_t {
        &mut self.segs[index]
    }

    pub(crate) fn char(&self, index: usize) -> c_char {
        self.chars[index]
    }

    pub(crate) fn set_char(&mut self, index: usize, value: c_char) {
        self.chars[index] = value;
    }

    pub(crate) fn resize_segs(&mut self, size: usize) {
        self.segs.resize_with(size, databar_segment_t::default);
    }

    pub(crate) fn epoch(&self) -> u8 {
        self.epoch
    }

    pub(crate) fn set_epoch(&mut self, val: u8) {
        self.epoch = val;
    }

    /// Reset DataBar decoder state
    pub(crate) fn reset(&mut self) {
        let n = self.segs.len();
        self.new_scan();
        for i in 0..n {
            let seg = self.seg_mut(i);
            seg.set_finder(-1);
        }
    }

    /// Prepare DataBar decoder for new scan
    pub(crate) fn new_scan(&mut self) {
        for i in 0..16 {
            if self.chars[i] >= 0 {
                let seg = &mut self.segs[self.chars[i] as usize];
                if seg.partial() {
                    seg.set_finder(-1);
                }
                self.chars[i] = -1;
            }
        }
    }
}

/// QR finder line (from qrcode.h)
#[derive(Default, Copy, Clone)]
pub(crate) struct qr_finder_line {
    pub(crate) pos: [c_int; 2], // qr_point
    pub(crate) len: c_int,
    pub(crate) boffs: c_int,
    pub(crate) eoffs: c_int,
}

/// QR Code finder state
#[derive(Default)]
pub(crate) struct qr_finder_t {
    pub(crate) s5: c_uint,
    pub(crate) line: qr_finder_line,
}

impl qr_finder_t {
    /// Reset QR finder state
    pub(crate) fn reset(&mut self) {
        self.s5 = 0;
    }
}
