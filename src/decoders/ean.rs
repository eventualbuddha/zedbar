//! EAN/UPC barcode decoder
//!
//! This module implements decoding for EAN-8, EAN-13, UPC-A, UPC-E,
//! ISBN-10, ISBN-13, EAN-2, and EAN-5 barcodes.

use crate::{
    config::internal::DecoderState, finder::decode_e, img_scanner::zbar_image_scanner_t,
    line_scanner::zbar_color_t, SymbolType,
};
use libc::{c_char, c_uint};

// State constants for ean_pass_t
const STATE_REV: i8 = -0x80; // 0x80 as signed
const STATE_ADDON: i8 = 0x40;
const STATE_IDX: i8 = 0x3f;

// Assertion macro
macro_rules! zassert {
    ($condition:expr, $retval:expr, $($arg:tt)*) => {
        if !$condition {
            return $retval;
        }
    };
}

// ============================================================================
// Lookup tables
// ============================================================================

/// Convert compact encoded D2E1E2 to character (bit4 is parity)
static DIGITS: [u8; 20] = [
    // E1   E2
    0x06, 0x10, 0x04, 0x13, //  2  2-5
    0x19, 0x08, 0x11, 0x05, //  3  2-5 (d2 <= thr)
    0x09, 0x12, 0x07, 0x15, //  4  2-5 (d2 <= thr)
    0x16, 0x00, 0x14, 0x03, //  5  2-5
    0x18, 0x01, 0x02, 0x17, // E1E2=43,44,33,34 (d2 > thr)
];

/// Parity decoding for UPC-E check digit and EAN-13 leading digit
static PARITY_DECODE: [u8; 64] = [
    0xf0, // [00] [xx] BBBBBB = RIGHT half EAN-13
    // UPC-E check digit encoding
    0xff, 0xff, 0x0f, // [01-03] [07] BBBAAA = 0
    0xff, 0x1f, // [04-05] [0b] BBABAA = 1
    0x2f, // [06] [0d] BBAABA = 2
    0xf3, // [07] [0e] BBAAAB = 3
    0xff, 0x4f, // [08-09] [13] BABBAA = 4
    0x7f, // [0a] [15] BABABA = 7
    0xf8, // [0b] [16] BABAAB = 8
    0x5f, // [0c] [19] BAABBA = 5
    0xf9, // [0d] [1a] BAABAB = 9
    0xf6, // [0e] [1c] BAAABB = 6
    0xff, // [0f]
    // LEFT half EAN-13 leading digit
    0xff, 0x6f, // [10-11] [23] ABBBAA = 6
    0x9f, // [12] [25] ABBABA = 9
    0xf5, // [13] [26] ABBAAB = 5
    0x8f, // [14] [29] ABABBA = 8
    0xf7, // [15] [2a] ABABAB = 7
    0xf4, // [16] [2c] ABAABB = 4
    0xff, 0x3f, // [17-18] [31] AABBBA = 3
    0xf2, // [19] [32] AABBAB = 2
    0xf1, // [1a] [34] AABABB = 1
    0xff, 0xff, 0xff, 0xff, 0x0f, // [1b-1f] [3f] AAAAAA = 0
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, // [20-27] padding
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, // [28-2f] padding
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, // [30-37] padding
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, // [38-3f] padding
];

// ============================================================================
// Helper functions
// ============================================================================

/// Calculate total character width "s"
#[inline]
fn calc_s(dcode: &zbar_image_scanner_t, mut offset: u8, mut n: u8) -> c_uint {
    let mut s = 0;
    while n > 0 {
        s += dcode.get_width(offset);
        offset += 1;
        n -= 1;
    }
    s
}

/// Check if two widths are within tolerance
#[inline]
fn check_width(w0: u32, w1: u32) -> c_uint {
    let dw0 = w0;
    let w0_scaled = w0 * 8;
    let w1_scaled = w1 * 8;
    if w0_scaled >= dw0 && w0_scaled.wrapping_sub(dw0) <= w1_scaled && w1_scaled <= w0_scaled + dw0
    {
        1
    } else {
        0
    }
}

/// Evaluate previous N (>= 2) widths as auxiliary pattern,
/// using preceding 4 as character width
fn aux_end(dcode: &zbar_image_scanner_t, fwd: u8) -> i8 {
    // reference width from previous character
    let s = calc_s(dcode, 4 + fwd, 4);

    // check quiet zone
    let qz = dcode.get_width(0);
    if fwd == 0 && qz != 0 && qz <= s * 3 / 4 {
        return -1;
    }

    let mut code: i8 = 0;
    for i in (1 - fwd)..(3 + fwd) {
        let e = dcode.get_width(i) + dcode.get_width(i + 1);
        let e_code = decode_e(e, s, 7);
        if e_code < 0 {
            return -1;
        }
        code = (code << 2) | (e_code as i8);
    }
    code
}

/// Determine possible auxiliary pattern using current 4 as possible character
fn aux_start(dcode: &zbar_image_scanner_t) -> i8 {
    // FIXME NB add-on has no guard in reverse
    let e2 = dcode.get_width(5) + dcode.get_width(6);

    if dcode.ean.s4 < 6 {
        return -1;
    }

    let e2_code = decode_e(e2, dcode.ean.s4, 7);
    if e2_code != 0 {
        return -1;
    }

    let e1 = dcode.get_width(4) + dcode.get_width(5);
    let e1_code = decode_e(e1, dcode.ean.s4, 7);

    if dcode.color() == zbar_color_t::ZBAR_BAR {
        // check for quiet-zone
        let qz = dcode.get_width(7);
        if qz == 0 || qz > dcode.ean.s4 * 3 / 4 {
            if e1_code == 0 {
                return 0; // normal symbol start
            } else if e1_code == 1 {
                return STATE_ADDON; // add-on symbol start
            }
        }
        return -1;
    }

    if e1_code == 0 {
        // attempting decode from SPACE => validate center guard
        let e3 = dcode.get_width(6) + dcode.get_width(7);
        let e4 = dcode.get_width(7) + dcode.get_width(8);
        if decode_e(e3, dcode.ean.s4, 7) == 0 && decode_e(e4, dcode.ean.s4, 7) == 0 {
            return 0; // start after center guard
        }
    }
    -1
}

/// Check addon delimiter using current 4 as character
#[inline]
fn aux_mid(dcode: &zbar_image_scanner_t) -> i8 {
    let e = dcode.get_width(4) + dcode.get_width(5);
    decode_e(e, dcode.ean.s4, 7) as i8
}

/// Attempt to decode previous 4 widths (2 bars and 2 spaces) as a character
fn decode4(dcode: &zbar_image_scanner_t) -> i8 {
    // calculate similar edge measurements
    let e1 = if dcode.color() == zbar_color_t::ZBAR_BAR {
        dcode.get_width(0) + dcode.get_width(1)
    } else {
        dcode.get_width(2) + dcode.get_width(3)
    };
    let e2 = dcode.get_width(1) + dcode.get_width(2);

    if dcode.ean.s4 < 6 {
        return -1;
    }

    // create compacted encoding for direct lookup
    let e1_code = decode_e(e1, dcode.ean.s4, 7);
    let e2_code = decode_e(e2, dcode.ean.s4, 7);
    if e1_code < 0 || e2_code < 0 {
        return -1;
    }
    let mut code = ((e1_code << 2) | e2_code) as i8;

    // 4 combinations require additional determinant (D2)
    // E1E2 == 34 (0110)
    // E1E2 == 43 (1001)
    // E1E2 == 33 (0101)
    // E1E2 == 44 (1010)
    if ((1 << code) & 0x0660) != 0 {
        // use sum of bar widths
        let d2 = if dcode.color() == zbar_color_t::ZBAR_BAR {
            dcode.get_width(0) + dcode.get_width(2)
        } else {
            dcode.get_width(1) + dcode.get_width(3)
        } * 7;

        let mid = if ((1 << code) & 0x0420) != 0 {
            3 // E1E2 in 33,44
        } else {
            4 // E1E2 in 34,43
        };
        let alt = d2 > (mid * dcode.ean.s4);
        if alt {
            code = ((code >> 1) & 3) | 0x10; // compress code space
        }
    }

    code
}

// ============================================================================
// Partial decode functions
// ============================================================================

/// EAN pass state
#[derive(Default)]
pub(crate) struct ean_pass_t {
    pub(crate) state: c_char,
    pub(crate) width: u32,
    pub(crate) raw: [u8; 7],
}

/// EAN/UPC decoder state
#[derive(Default)]
pub(crate) struct ean_decoder_t {
    pub(crate) pass: [ean_pass_t; 4],
    pub(crate) left: SymbolType,
    pub(crate) right: SymbolType,
    pub(crate) direction: i32,
    pub(crate) s4: u32,
    pub(crate) width: u32,
    pub(crate) buf: [c_char; 18],
    pub(crate) enable: bool,
}

impl ean_decoder_t {
    /// Prepare EAN decoder for new scan
    pub(crate) fn new_scan(&mut self) {
        self.pass[0].state = -1;
        self.pass[1].state = -1;
        self.pass[2].state = -1;
        self.pass[3].state = -1;
        self.s4 = 0;
    }

    /// Reset EAN decoder state
    pub(crate) fn reset(&mut self) {
        self.new_scan();
        self.left = SymbolType::None;
        self.right = SymbolType::None;
    }

    // Removed ean_get_config - now use dcode.is_enabled() and related helpers

    /// Handle EAN-13/UPC-E partial
    fn ean_part_end7(
        &mut self,
        config: &DecoderState,
        pass_index: usize,
        fwd: u8,
    ) -> PartialSymbolType {
        // calculate parity index
        let par = if fwd != 0 {
            ((self.pass[pass_index].raw[1] & 0x10) << 1)
                | (self.pass[pass_index].raw[2] & 0x10)
                | ((self.pass[pass_index].raw[3] & 0x10) >> 1)
                | ((self.pass[pass_index].raw[4] & 0x10) >> 2)
                | ((self.pass[pass_index].raw[5] & 0x10) >> 3)
                | ((self.pass[pass_index].raw[6] & 0x10) >> 4)
        } else {
            ((self.pass[pass_index].raw[1] & 0x10) >> 4)
                | ((self.pass[pass_index].raw[2] & 0x10) >> 3)
                | ((self.pass[pass_index].raw[3] & 0x10) >> 2)
                | ((self.pass[pass_index].raw[4] & 0x10) >> 1)
                | (self.pass[pass_index].raw[5] & 0x10)
                | ((self.pass[pass_index].raw[6] & 0x10) << 1)
        };

        // lookup parity combination
        self.pass[pass_index].raw[0] = PARITY_DECODE[(par >> 1) as usize];
        if (par & 1) != 0 {
            self.pass[pass_index].raw[0] >>= 4;
        }
        self.pass[pass_index].raw[0] &= 0xf;

        if self.pass[pass_index].raw[0] == 0xf {
            // invalid parity combination
            return PartialSymbolType::None;
        }

        if (par != 0) != (fwd != 0) {
            self.pass[pass_index].state |= STATE_REV;
            // reverse sampled digits
            for i in 1..4 {
                self.pass[pass_index].raw.swap(i, 7 - i);
            }
        }

        if config.is_enabled(SymbolType::Ean13) {
            if par == 0 {
                return PartialSymbolType::Ean13(Side::Right);
            }
            if (par & 0x20) != 0 {
                return PartialSymbolType::Ean13(Side::Left);
            }
        }

        if par != 0 && (par & 0x20) == 0 {
            return PartialSymbolType::Upce;
        }

        PartialSymbolType::None
    }

    /// Handle EAN-5 addon
    fn ean_part_end5(&self, config: &DecoderState, pass_index: usize) -> PartialSymbolType {
        if !config.is_enabled(SymbolType::Ean5) {
            return PartialSymbolType::None;
        }

        // extract parity bits
        let par = (self.pass[pass_index].raw[1] & 0x10)
            | ((self.pass[pass_index].raw[2] & 0x10) >> 1)
            | ((self.pass[pass_index].raw[3] & 0x10) >> 2)
            | ((self.pass[pass_index].raw[4] & 0x10) >> 3)
            | ((self.pass[pass_index].raw[5] & 0x10) >> 4);

        // calculate checksum
        let chk = ((((self.pass[pass_index].raw[1] & 0x0f) as c_uint
            + (self.pass[pass_index].raw[2] & 0x0f) as c_uint * 3
            + (self.pass[pass_index].raw[3] & 0x0f) as c_uint
            + (self.pass[pass_index].raw[4] & 0x0f) as c_uint * 3
            + (self.pass[pass_index].raw[5] & 0x0f) as c_uint)
            * 3)
            % 10) as u8;

        let mut parchk = PARITY_DECODE[(par >> 1) as usize];
        if (par & 1) != 0 {
            parchk >>= 4;
        }
        parchk &= 0xf;

        if parchk != chk {
            return PartialSymbolType::None;
        }

        PartialSymbolType::Ean5
    }

    /// EAN checksum verification
    fn ean_verify_checksum(&self, n: usize) -> i8 {
        let mut chk: u8 = 0;
        for i in 0..n {
            let d = self.buf[i] as u8;
            zassert!(d < 10, -1, "i={:x} d={:x} chk={:x}", i, d, chk);
            chk = chk.wrapping_add(d);
            if (i ^ n) & 1 != 0 {
                chk = chk.wrapping_add(d << 1);
                if chk >= 20 {
                    chk -= 20;
                }
            }
            if chk >= 10 {
                chk -= 10;
            }
        }
        zassert!(chk < 10, -1, "chk={:x} n={:x}", chk, n);
        if chk != 0 {
            chk = 10 - chk;
        }
        let d = self.buf[n] as u8;
        zassert!(d < 10, -1, "n={:x} d={:x} chk={:x}", n, d, chk);
        if chk != d {
            return -1;
        }
        0
    }

    /// Calculate ISBN-10 checksum
    fn isbn10_calc_checksum(&self) -> c_char {
        let mut chk: u32 = 0;
        for w in (2..=10).rev() {
            let d = self.buf[13 - w] as u8;
            zassert!(d < 10, b'?' as c_char, "w={:x} d={:x} chk={:x}", w, d, chk);
            chk += d as c_uint * w as c_uint;
        }
        chk %= 11;
        if chk == 0 {
            return b'0' as c_char;
        }
        chk = 11 - chk;
        if chk < 10 {
            return (chk + b'0' as c_uint) as c_char;
        }
        b'X' as c_char
    }

    /// Expand UPC-E to UPC-A
    fn ean_expand_upce(&mut self, pass_index: usize) {
        let mut i = 0;

        // parity encoded digit is checksum
        self.buf[12] = self.pass[pass_index].raw[i] as c_char;
        i += 1;

        let decode = self.pass[pass_index].raw[6] & 0xf;
        self.buf[0] = 0;
        self.buf[1] = 0;
        self.buf[2] = (self.pass[pass_index].raw[i] & 0xf) as c_char;
        i += 1;
        self.buf[3] = (self.pass[pass_index].raw[i] & 0xf) as c_char;
        i += 1;
        self.buf[4] = if decode < 3 {
            decode as c_char
        } else {
            (self.pass[pass_index].raw[i] & 0xf) as c_char
        };
        if decode >= 3 {
            i += 1;
        }
        self.buf[5] = if decode < 4 {
            0
        } else {
            (self.pass[pass_index].raw[i] & 0xf) as c_char
        };
        if decode >= 4 {
            i += 1;
        }
        self.buf[6] = if decode < 5 {
            0
        } else {
            (self.pass[pass_index].raw[i] & 0xf) as c_char
        };
        if decode >= 5 {
            i += 1;
        }
        self.buf[7] = 0;
        self.buf[8] = 0;
        self.buf[9] = if decode < 3 {
            (self.pass[pass_index].raw[i] & 0xf) as c_char
        } else {
            0
        };
        if decode < 3 {
            i += 1;
        }
        self.buf[10] = if decode < 4 {
            (self.pass[pass_index].raw[i] & 0xf) as c_char
        } else {
            0
        };
        if decode < 4 {
            i += 1;
        }
        self.buf[11] = if decode < 5 {
            (self.pass[pass_index].raw[i] & 0xf) as c_char
        } else {
            decode as c_char
        };
    }

    /// Integrate partial decode results
    fn integrate_partial(
        &mut self,
        config: &DecoderState,
        pass_index: usize,
        mut part: PartialSymbolType,
    ) -> SymbolType {
        // if same partial is not consistent, reset others
        if (self.left != SymbolType::None && (SymbolType::from(part) != self.left))
            || (self.right != SymbolType::None && (SymbolType::from(part) != self.right))
        {
            // partial mismatch - reset collected parts
            self.left = SymbolType::None;
            self.right = SymbolType::None;
        }

        if (self.left != SymbolType::None || self.right != SymbolType::None)
            && check_width(self.width, self.pass[pass_index].width) == 0
        {
            self.left = SymbolType::None;
            self.right = SymbolType::None;
        }

        if part.side() == Some(Side::Right) {
            let mut j = i32::from(part) - 1;
            let mut i = i32::from(part) >> 1;
            while i > 0 {
                let digit = (self.pass[pass_index].raw[i as usize] & 0xf) as c_char;
                if self.right != SymbolType::None && self.buf[j as usize] != digit {
                    // partial mismatch - reset collected parts
                    self.left = SymbolType::None;
                    self.right = SymbolType::None;
                }
                self.buf[j as usize] = digit;
                i -= 1;
                j -= 1;
            }
            self.right = part.into();
            part = part
                .replace_symbol_type(self.left)
                .expect("left and right should be compatible"); // FIXME!?
        } else if matches!(
            part,
            PartialSymbolType::Ean13(_) | PartialSymbolType::Ean8(_)
        ) {
            // EAN_LEFT
            let mut j = (i32::from(part) - 1) >> 1;
            let mut i = i32::from(part) >> 1;
            while j >= 0 {
                let digit = (self.pass[pass_index].raw[i as usize] & 0xf) as c_char;
                if self.left != SymbolType::None && self.buf[j as usize] != digit {
                    // partial mismatch - reset collected parts
                    self.left = SymbolType::None;
                    self.right = SymbolType::None;
                }
                self.buf[j as usize] = digit;
                i -= 1;
                j -= 1;
            }
            self.left = part.into();
            part = part
                .replace_symbol_type(self.right)
                .expect("left and right should be compatible"); // FIXME!?
        } else if part != PartialSymbolType::Upce {
            // add-ons
            for i in (1..=(i32::from(part) as usize)).rev() {
                self.buf[i - 1] = (self.pass[pass_index].raw[i] & 0xf) as c_char;
            }
            self.left = part.into();
        } else {
            self.ean_expand_upce(pass_index);
        }

        self.width = self.pass[pass_index].width;

        // Initialize symbol_type from part, then override for special cases
        let mut symbol_type = SymbolType::from(part);

        if part == PartialSymbolType::None {
            symbol_type = SymbolType::Partial;
        }

        if (matches!(part, PartialSymbolType::Ean13(_) | PartialSymbolType::Upce)
            && self.ean_verify_checksum(12) != 0)
            || (matches!(part, PartialSymbolType::Ean8(_)) && self.ean_verify_checksum(7) != 0)
        {
            // invalid checksum
            if self.right != SymbolType::None {
                self.left = SymbolType::None;
            } else {
                self.right = SymbolType::None;
            }
            symbol_type = SymbolType::None;
        }

        if matches!(part, PartialSymbolType::Ean13(_)) {
            // special case ean-13 subsets
            if self.buf[0] == 0 && config.is_enabled(SymbolType::Upca) {
                symbol_type = SymbolType::Upca;
            } else if self.buf[0] == 9 && self.buf[1] == 7 {
                if (self.buf[2] == 8 || self.buf[2] == 9) && config.is_enabled(SymbolType::Isbn13) {
                    symbol_type = SymbolType::Isbn13;
                } else if self.buf[2] == 8 && config.is_enabled(SymbolType::Isbn10) {
                    symbol_type = SymbolType::Isbn10;
                }
            }
        } else if matches!(part, PartialSymbolType::Upce) {
            if config.is_enabled(SymbolType::Upce) {
                // UPC-E was decompressed for checksum verification,
                // but user requested compressed result
                self.buf[0] = 0;
                self.buf[1] = 0;
                for i in 2..8 {
                    self.buf[i] = (self.pass[pass_index].raw[i - 1] & 0xf) as c_char;
                }
                self.buf[8] = (self.pass[pass_index].raw[0] & 0xf) as c_char;
            } else if config.is_enabled(SymbolType::Upca) {
                // UPC-E reported as UPC-A has priority over ean-13
                symbol_type = SymbolType::Upca;
            } else if config.is_enabled(SymbolType::Ean13) {
                symbol_type = SymbolType::Ean13;
            } else {
                symbol_type = SymbolType::None;
            }
        }

        symbol_type
    }

    /// Handle EAN-2 addon
    fn ean_part_end2(&mut self, config: &DecoderState, pass_index: usize) -> PartialSymbolType {
        if !config.is_enabled(SymbolType::Ean2) {
            return PartialSymbolType::None;
        }

        // extract parity bits
        let par = ((self.pass[pass_index].raw[1] & 0x10) >> 3)
            | ((self.pass[pass_index].raw[2] & 0x10) >> 4);
        // calculate "checksum"
        let chk = (!(((self.pass[pass_index].raw[1] & 0xf) * 10)
            + (self.pass[pass_index].raw[2] & 0xf)))
            & 0x3;
        if par != chk {
            return PartialSymbolType::None;
        }
        PartialSymbolType::Ean2
    }
}

/// Handle EAN-8 partial (left or right half)
fn ean_part_end4(pass: &mut ean_pass_t, fwd: u8) -> PartialSymbolType {
    // extract parity bits
    let par = ((pass.raw[1] & 0x10) >> 1)
        | ((pass.raw[2] & 0x10) >> 2)
        | ((pass.raw[3] & 0x10) >> 3)
        | ((pass.raw[4] & 0x10) >> 4);

    if par != 0 && par != 0xf {
        // invalid parity combination
        return PartialSymbolType::None;
    }

    if (par == 0) != (fwd != 0) {
        // reverse sampled digits
        pass.state |= STATE_REV;
        pass.raw.swap(1, 4);
        pass.raw.swap(2, 3);
    }

    if par == 0 {
        PartialSymbolType::Ean8(Side::Right)
    } else {
        PartialSymbolType::Ean8(Side::Left)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Side {
    Left,
    Right,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PartialSymbolType {
    Ean2,
    Ean5,
    Ean8(Side),
    Ean13(Side),
    Upce,
    None,
}

impl PartialSymbolType {
    fn side(self) -> Option<Side> {
        match self {
            Self::Ean8(side) | Self::Ean13(side) => Some(side),
            Self::Ean2 | Self::Ean5 | Self::Upce | Self::None => None,
        }
    }

    fn replace_symbol_type(self, symbol_type: SymbolType) -> Option<Self> {
        let side = self.side();
        match symbol_type {
            SymbolType::Ean2 => Some(Self::Ean2),
            SymbolType::Ean5 => Some(Self::Ean5),
            SymbolType::Ean8 => side.map(Self::Ean8),
            SymbolType::Ean13 => side.map(Self::Ean13),
            SymbolType::Upce => Some(Self::Upce),
            SymbolType::None => Some(Self::None),
            _ => None,
        }
    }
}

impl From<PartialSymbolType> for SymbolType {
    fn from(value: PartialSymbolType) -> Self {
        match value {
            PartialSymbolType::Ean2 => Self::Ean2,
            PartialSymbolType::Ean5 => Self::Ean5,
            PartialSymbolType::Ean8(_) => Self::Ean8,
            PartialSymbolType::Ean13(_) => Self::Ean13,
            PartialSymbolType::Upce => Self::Upce,
            PartialSymbolType::None => Self::None,
        }
    }
}

impl From<PartialSymbolType> for i32 {
    fn from(value: PartialSymbolType) -> Self {
        SymbolType::from(value).into()
    }
}

/// Copy result to output buffer
fn postprocess(dcode: &mut zbar_image_scanner_t, sym: SymbolType) {
    let mut j: usize = 0;
    let new_direction;
    let base;
    let mut i: usize = 0;
    let mut temp_buf = [0i8; 18]; // Max EAN buffer size
    let mut buf_len = 0;
    let needs_isbn10_check;
    let isbn10_check_digit;

    {
        let ean = &dcode.ean;
        base = sym;

        if base > SymbolType::Partial {
            match base {
                SymbolType::Upca | SymbolType::Upce => i = 1,
                SymbolType::Isbn10 => i = 3,
                _ => {}
            }

            let mut calc_base = base;
            match base {
                SymbolType::Isbn13 => calc_base = SymbolType::Ean13,
                SymbolType::Upce => calc_base = SymbolType::Ean8,
                _ => {}
            }

            if base == SymbolType::Isbn10
                || (calc_base > SymbolType::Ean5 && !dcode.should_emit_checksum(sym))
            {
                calc_base = (calc_base as i32 - 1).into();
            }

            // Copy to temp buffer
            while j < calc_base as usize && j < 18 && ean.buf[i] >= 0 {
                temp_buf[j] = ean.buf[i];
                i += 1;
                j += 1;
            }
            buf_len = j;

            needs_isbn10_check = base == SymbolType::Isbn10
                && j == 9
                && dcode.should_emit_checksum(SymbolType::Isbn10);
            isbn10_check_digit = if needs_isbn10_check {
                ean.isbn10_calc_checksum() as u8
            } else {
                0
            };
        } else {
            needs_isbn10_check = false;
            isbn10_check_digit = 0;
        }

        new_direction = 1 - 2 * ean.direction;
    }

    if base > SymbolType::Partial {
        // Get mutable slice for writing
        let buffer = match dcode.buffer_mut_slice(buf_len + if needs_isbn10_check { 2 } else { 1 })
        {
            Ok(buf) => buf,
            Err(_) => return,
        };

        // Copy from temp buffer
        for k in 0..buf_len {
            buffer[k] = (temp_buf[k] + b'0' as c_char) as u8;
        }

        if needs_isbn10_check {
            // Use pre-calculated ISBN-10 check digit
            buffer[buf_len] = isbn10_check_digit;
            buf_len += 1;
        }

        // Add null terminator
        buffer[buf_len] = 0;
    }

    dcode.truncate_buffer(buf_len);
    dcode.direction = new_direction;
    dcode.modifiers = 0;
}

/// Update state for one of 4 parallel passes
fn decode_pass(dcode: &mut zbar_image_scanner_t, pass_index: usize) -> PartialSymbolType {
    dcode.ean.pass[pass_index].state = dcode.ean.pass[pass_index].state.wrapping_add(1);
    let idx = dcode.ean.pass[pass_index].state & STATE_IDX;
    let fwd = (dcode.ean.pass[pass_index].state & 1) as u8;

    if dcode.color() == zbar_color_t::ZBAR_SPACE {
        if (dcode.ean.pass[pass_index].state & STATE_ADDON) != 0 {
            if idx == 0x09 || idx == 0x21 {
                let qz = dcode.get_width(0);
                let s = calc_s(dcode, 1, 4);
                let part = if qz == 0 || qz >= s * 3 / 4 {
                    if idx == 0x09 {
                        dcode.ean.ean_part_end2(&dcode.config, pass_index)
                    } else {
                        dcode.ean.ean_part_end5(&dcode.config, pass_index)
                    }
                } else {
                    PartialSymbolType::None
                };

                if part != PartialSymbolType::None || idx == 0x21 {
                    let ean = &mut dcode.ean;
                    ean.direction = 0;
                    dcode.ean.pass[pass_index].state = -1;
                    return part;
                }
            }
            if (idx & 7) == 1 {
                dcode.ean.pass[pass_index].state += 2;
            }
        } else if (idx == 0x10 || idx == 0x11)
            && dcode.is_enabled(SymbolType::Ean8)
            && aux_end(dcode, fwd) == 0
        {
            let part = ean_part_end4(&mut dcode.ean.pass[pass_index], fwd);
            if part != PartialSymbolType::None {
                dcode.ean.direction = if (dcode.ean.pass[pass_index].state & STATE_REV) != 0 {
                    1
                } else {
                    0
                };
            }
            dcode.ean.pass[pass_index].state = -1;
            return part;
        } else if idx == 0x18 || idx == 0x19 {
            let mut part = PartialSymbolType::None;
            if aux_end(dcode, fwd) == 0 && dcode.ean.pass[pass_index].raw[5] != 0xff {
                part = dcode.ean.ean_part_end7(&dcode.config, pass_index, fwd);
            }
            if part != PartialSymbolType::None {
                dcode.ean.direction = if (dcode.ean.pass[pass_index].state & STATE_REV) != 0 {
                    1
                } else {
                    0
                };
            }
            dcode.ean.pass[pass_index].state = -1;
            return part;
        }
    }

    let mut idx = idx;
    if (dcode.ean.pass[pass_index].state & STATE_ADDON) != 0 {
        idx >>= 1;
    }

    if (idx & 0x03) == 0 && idx <= 0x14 {
        let mut code: i8 = -1;
        let mut w = dcode.ean.pass[pass_index].width;

        if dcode.ean.s4 == 0 {
            return PartialSymbolType::None;
        }

        // validate guard bars before decoding first char of symbol
        if dcode.ean.pass[pass_index].state == 0 {
            dcode.ean.pass[pass_index].state = aux_start(dcode);
            dcode.ean.pass[pass_index].width = dcode.ean.s4;
            if dcode.ean.pass[pass_index].state < 0 {
                return PartialSymbolType::None;
            }
            idx = dcode.ean.pass[pass_index].state & STATE_IDX;
        } else {
            w = check_width(w, dcode.ean.s4);
            if w != 0 {
                dcode.ean.pass[pass_index].width =
                    (dcode.ean.pass[pass_index].width + dcode.ean.s4 * 3) / 4;
            }
        }

        if w != 0 {
            code = decode4(dcode);
        }

        if (code < 0 && idx != 0x10)
            || (idx > 0
                && (dcode.ean.pass[pass_index].state & STATE_ADDON) != 0
                && aux_mid(dcode) != 0)
        {
            dcode.ean.pass[pass_index].state = -1;
        } else if code < 0 {
            dcode.ean.pass[pass_index].raw[5] = 0xff;
        } else {
            let raw_idx = ((idx >> 2) + 1) as usize;
            dcode.ean.pass[pass_index].raw[raw_idx] = DIGITS[code as usize];
        }
    }

    PartialSymbolType::None
}

// ============================================================================
// Main decoder function
// ============================================================================

/// Main EAN/UPC decoder entry point
pub(crate) fn zbar_decode_ean(dcode: &mut zbar_image_scanner_t) -> SymbolType {
    // process up to 4 separate passes
    let mut sym = SymbolType::None;
    let pass_idx = (dcode.idx & 3) as usize;

    // update latest character width
    dcode.ean.s4 = dcode.ean.s4.wrapping_sub(dcode.get_width(4));
    dcode.ean.s4 = dcode.ean.s4.wrapping_add(dcode.get_width(0));

    for i in 0..4 {
        if dcode.ean.pass[i].state >= 0 || i == pass_idx {
            let part = decode_pass(dcode, i);
            if part != PartialSymbolType::None {
                // update accumulated data from new partial decode
                sym = dcode.ean.integrate_partial(&dcode.config, i, part);
                if sym != SymbolType::None {
                    // this pass valid => _reset_ all passes
                    dcode.ean.pass[0].state = -1;
                    dcode.ean.pass[1].state = -1;
                    dcode.ean.pass[2].state = -1;
                    dcode.ean.pass[3].state = -1;
                    if sym > SymbolType::Partial {
                        if dcode.acquire_lock(sym) {
                            postprocess(dcode, sym);
                        } else {
                            sym = SymbolType::Partial;
                        }
                    }
                }
            }
        }
    }

    sym
}
