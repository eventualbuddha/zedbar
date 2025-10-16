//! EAN/UPC barcode decoder
//!
//! This module implements decoding for EAN-8, EAN-13, UPC-A, UPC-E,
//! ISBN-10, ISBN-13, EAN-2, and EAN-5 barcodes.

use crate::{
    decoder_types::{
        ean_decoder_t, ean_pass_t, zbar_decoder_t, zbar_symbol_type_t, DECODE_WINDOW,
        ZBAR_CFG_EMIT_CHECK, ZBAR_CFG_ENABLE, ZBAR_EAN13, ZBAR_EAN2, ZBAR_EAN5, ZBAR_EAN8,
        ZBAR_ISBN10, ZBAR_ISBN13, ZBAR_NONE, ZBAR_PARTIAL, ZBAR_SYMBOL, ZBAR_UPCA, ZBAR_UPCE,
    },
    line_scanner::zbar_color_t,
};
use libc::{c_char, c_int, c_uint};

// State constants for ean_pass_t
const STATE_REV: i8 = -0x80; // 0x80 as signed
const STATE_ADDON: i8 = 0x40;
const STATE_IDX: i8 = 0x3f;

// Partial decode symbol location
const EAN_LEFT: zbar_symbol_type_t = 0x0000;
const EAN_RIGHT: zbar_symbol_type_t = 0x1000;

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

/// Retrieve i-th previous element width
#[inline]
fn get_width(dcode: &zbar_decoder_t, offset: u8) -> c_uint {
    dcode.w[((dcode.idx.wrapping_sub(offset)) & (DECODE_WINDOW as u8 - 1)) as usize]
}

/// Calculate total character width "s"
#[inline]
fn calc_s(dcode: &zbar_decoder_t, mut offset: u8, mut n: u8) -> c_uint {
    let mut s = 0;
    while n > 0 {
        s += get_width(dcode, offset);
        offset += 1;
        n -= 1;
    }
    s
}

/// Fixed character width decode assist
#[inline]
fn decode_e(e: c_uint, s: c_uint, n: c_uint) -> i8 {
    let e_val = ((e * n * 2 + 1) / s).wrapping_sub(3) / 2;
    if e_val >= n - 3 {
        -1
    } else {
        e_val as i8
    }
}

/// Check if two widths are within tolerance
#[inline]
fn check_width(w0: c_uint, w1: c_uint) -> c_uint {
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

/// Test if a config flag is enabled
#[inline]
fn test_cfg(config: c_uint, cfg: c_int) -> bool {
    ((config >> cfg) & 1) != 0
}

/// Get config for a specific symbol type
#[inline]
fn ean_get_config(ean: &ean_decoder_t, sym: zbar_symbol_type_t) -> c_uint {
    match sym {
        ZBAR_EAN2 => ean.ean2_config,
        ZBAR_EAN5 => ean.ean5_config,
        ZBAR_EAN8 => ean.ean8_config,
        ZBAR_UPCE => ean.upce_config,
        ZBAR_ISBN10 => ean.isbn10_config,
        ZBAR_UPCA => ean.upca_config,
        ZBAR_EAN13 => ean.ean13_config,
        ZBAR_ISBN13 => ean.isbn13_config,
        _ => 0,
    }
}

/// Evaluate previous N (>= 2) widths as auxiliary pattern,
/// using preceding 4 as character width
fn aux_end(dcode: &zbar_decoder_t, fwd: u8) -> i8 {
    // reference width from previous character
    let s = calc_s(dcode, 4 + fwd, 4);

    // check quiet zone
    let qz = get_width(dcode, 0);
    if fwd == 0 && qz != 0 && qz <= s * 3 / 4 {
        return -1;
    }

    let mut code: i8 = 0;
    for i in (1 - fwd)..(3 + fwd) {
        let e = get_width(dcode, i) + get_width(dcode, i + 1);
        let e_code = decode_e(e, s, 7);
        if e_code < 0 {
            return -1;
        }
        code = (code << 2) | e_code;
    }
    code
}

/// Determine possible auxiliary pattern using current 4 as possible character
fn aux_start(dcode: &zbar_decoder_t) -> i8 {
    // FIXME NB add-on has no guard in reverse
    let e2 = get_width(dcode, 5) + get_width(dcode, 6);

    if dcode.ean.s4 < 6 {
        return -1;
    }

    let e2_code = decode_e(e2, dcode.ean.s4, 7);
    if e2_code != 0 {
        return -1;
    }

    let e1 = get_width(dcode, 4) + get_width(dcode, 5);
    let e1_code = decode_e(e1, dcode.ean.s4, 7);

    if dcode.color() == zbar_color_t::ZBAR_BAR {
        // check for quiet-zone
        let qz = get_width(dcode, 7);
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
        let e3 = get_width(dcode, 6) + get_width(dcode, 7);
        let e4 = get_width(dcode, 7) + get_width(dcode, 8);
        if decode_e(e3, dcode.ean.s4, 7) == 0 && decode_e(e4, dcode.ean.s4, 7) == 0 {
            return 0; // start after center guard
        }
    }
    -1
}

/// Check addon delimiter using current 4 as character
#[inline]
fn aux_mid(dcode: &zbar_decoder_t) -> i8 {
    let e = get_width(dcode, 4) + get_width(dcode, 5);
    decode_e(e, dcode.ean.s4, 7)
}

/// Attempt to decode previous 4 widths (2 bars and 2 spaces) as a character
fn decode4(dcode: &zbar_decoder_t) -> i8 {
    // calculate similar edge measurements
    let e1 = if dcode.color() == zbar_color_t::ZBAR_BAR {
        get_width(dcode, 0) + get_width(dcode, 1)
    } else {
        get_width(dcode, 2) + get_width(dcode, 3)
    };
    let e2 = get_width(dcode, 1) + get_width(dcode, 2);

    if dcode.ean.s4 < 6 {
        return -1;
    }

    // create compacted encoding for direct lookup
    let e1_code = decode_e(e1, dcode.ean.s4, 7);
    let e2_code = decode_e(e2, dcode.ean.s4, 7);
    if e1_code < 0 || e2_code < 0 {
        return -1;
    }
    let mut code = (e1_code << 2) | e2_code;

    // 4 combinations require additional determinant (D2)
    // E1E2 == 34 (0110)
    // E1E2 == 43 (1001)
    // E1E2 == 33 (0101)
    // E1E2 == 44 (1010)
    if ((1 << code) & 0x0660) != 0 {
        // use sum of bar widths
        let d2 = if dcode.color() == zbar_color_t::ZBAR_BAR {
            get_width(dcode, 0) + get_width(dcode, 2)
        } else {
            get_width(dcode, 1) + get_width(dcode, 3)
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

/// Handle EAN-2 addon
fn ean_part_end2(ean: &ean_decoder_t, pass: &ean_pass_t) -> zbar_symbol_type_t {
    if !test_cfg(ean.ean2_config, ZBAR_CFG_ENABLE) {
        return ZBAR_NONE;
    }

    // extract parity bits
    let par = ((pass.raw[1] & 0x10) >> 3) | ((pass.raw[2] & 0x10) >> 4);
    // calculate "checksum"
    let chk = (!(((pass.raw[1] & 0xf) * 10) + (pass.raw[2] & 0xf))) & 0x3;
    if par != chk {
        return ZBAR_NONE;
    }
    ZBAR_EAN2
}

/// Handle EAN-8 partial (left or right half)
fn ean_part_end4(pass: &mut ean_pass_t, fwd: u8) -> zbar_symbol_type_t {
    // extract parity bits
    let par = ((pass.raw[1] & 0x10) >> 1)
        | ((pass.raw[2] & 0x10) >> 2)
        | ((pass.raw[3] & 0x10) >> 3)
        | ((pass.raw[4] & 0x10) >> 4);

    if par != 0 && par != 0xf {
        // invalid parity combination
        return ZBAR_NONE;
    }

    if (par == 0) != (fwd != 0) {
        // reverse sampled digits
        pass.state |= STATE_REV;
        pass.raw.swap(1, 4);
        pass.raw.swap(2, 3);
    }

    if par == 0 {
        ZBAR_EAN8 | EAN_RIGHT
    } else {
        ZBAR_EAN8 | EAN_LEFT
    }
}

/// Handle EAN-5 addon
fn ean_part_end5(ean: &ean_decoder_t, pass: &ean_pass_t) -> zbar_symbol_type_t {
    if !test_cfg(ean.ean5_config, ZBAR_CFG_ENABLE) {
        return ZBAR_NONE;
    }

    // extract parity bits
    let par = (pass.raw[1] & 0x10)
        | ((pass.raw[2] & 0x10) >> 1)
        | ((pass.raw[3] & 0x10) >> 2)
        | ((pass.raw[4] & 0x10) >> 3)
        | ((pass.raw[5] & 0x10) >> 4);

    // calculate checksum
    let chk = ((((pass.raw[1] & 0x0f) as c_uint
        + (pass.raw[2] & 0x0f) as c_uint * 3
        + (pass.raw[3] & 0x0f) as c_uint
        + (pass.raw[4] & 0x0f) as c_uint * 3
        + (pass.raw[5] & 0x0f) as c_uint)
        * 3)
        % 10) as u8;

    let mut parchk = PARITY_DECODE[(par >> 1) as usize];
    if (par & 1) != 0 {
        parchk >>= 4;
    }
    parchk &= 0xf;

    if parchk != chk {
        return ZBAR_NONE;
    }

    ZBAR_EAN5
}

/// Handle EAN-13/UPC-E partial
fn ean_part_end7(ean: &ean_decoder_t, pass: &mut ean_pass_t, fwd: u8) -> zbar_symbol_type_t {
    // calculate parity index
    let par = if fwd != 0 {
        ((pass.raw[1] & 0x10) << 1)
            | (pass.raw[2] & 0x10)
            | ((pass.raw[3] & 0x10) >> 1)
            | ((pass.raw[4] & 0x10) >> 2)
            | ((pass.raw[5] & 0x10) >> 3)
            | ((pass.raw[6] & 0x10) >> 4)
    } else {
        ((pass.raw[1] & 0x10) >> 4)
            | ((pass.raw[2] & 0x10) >> 3)
            | ((pass.raw[3] & 0x10) >> 2)
            | ((pass.raw[4] & 0x10) >> 1)
            | (pass.raw[5] & 0x10)
            | ((pass.raw[6] & 0x10) << 1)
    };

    // lookup parity combination
    pass.raw[0] = PARITY_DECODE[(par >> 1) as usize];
    if (par & 1) != 0 {
        pass.raw[0] >>= 4;
    }
    pass.raw[0] &= 0xf;

    if pass.raw[0] == 0xf {
        // invalid parity combination
        return ZBAR_NONE;
    }

    if (par != 0) != (fwd != 0) {
        pass.state |= STATE_REV;
        // reverse sampled digits
        for i in 1..4 {
            pass.raw.swap(i, 7 - i);
        }
    }

    if test_cfg(ean.ean13_config, ZBAR_CFG_ENABLE) {
        if par == 0 {
            return ZBAR_EAN13 | EAN_RIGHT;
        }
        if (par & 0x20) != 0 {
            return ZBAR_EAN13 | EAN_LEFT;
        }
    }

    if par != 0 && (par & 0x20) == 0 {
        return ZBAR_UPCE;
    }

    ZBAR_NONE
}

/// Acquire lock for a symbol
#[inline]
fn acquire_lock(dcode: &mut zbar_decoder_t, sym: zbar_symbol_type_t) -> bool {
    if dcode.lock != 0 {
        return true;
    }
    dcode.lock = sym;
    false
}

/// EAN checksum verification
fn ean_verify_checksum(ean: &ean_decoder_t, n: usize) -> i8 {
    let mut chk: u8 = 0;
    for i in 0..n {
        let d = ean.buf[i] as u8;
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
    let d = ean.buf[n] as u8;
    zassert!(d < 10, -1, "n={:x} d={:x} chk={:x}", n, d, chk);
    if chk != d {
        return -1;
    }
    0
}

/// Calculate ISBN-10 checksum
fn isbn10_calc_checksum(ean: &ean_decoder_t) -> c_char {
    let mut chk: c_uint = 0;
    for w in (2..=10).rev() {
        let d = ean.buf[13 - w] as u8;
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
fn ean_expand_upce(ean: &mut ean_decoder_t, pass: &ean_pass_t) {
    let mut i = 0;

    // parity encoded digit is checksum
    ean.buf[12] = pass.raw[i] as c_char;
    i += 1;

    let decode = pass.raw[6] & 0xf;
    ean.buf[0] = 0;
    ean.buf[1] = 0;
    ean.buf[2] = (pass.raw[i] & 0xf) as c_char;
    i += 1;
    ean.buf[3] = (pass.raw[i] & 0xf) as c_char;
    i += 1;
    ean.buf[4] = if decode < 3 {
        decode as c_char
    } else {
        (pass.raw[i] & 0xf) as c_char
    };
    if decode >= 3 {
        i += 1;
    }
    ean.buf[5] = if decode < 4 {
        0
    } else {
        (pass.raw[i] & 0xf) as c_char
    };
    if decode >= 4 {
        i += 1;
    }
    ean.buf[6] = if decode < 5 {
        0
    } else {
        (pass.raw[i] & 0xf) as c_char
    };
    if decode >= 5 {
        i += 1;
    }
    ean.buf[7] = 0;
    ean.buf[8] = 0;
    ean.buf[9] = if decode < 3 {
        (pass.raw[i] & 0xf) as c_char
    } else {
        0
    };
    if decode < 3 {
        i += 1;
    }
    ean.buf[10] = if decode < 4 {
        (pass.raw[i] & 0xf) as c_char
    } else {
        0
    };
    if decode < 4 {
        i += 1;
    }
    ean.buf[11] = if decode < 5 {
        (pass.raw[i] & 0xf) as c_char
    } else {
        decode as c_char
    };
}

/// Integrate partial decode results
fn integrate_partial(
    ean: &mut ean_decoder_t,
    pass: &mut ean_pass_t,
    mut part: zbar_symbol_type_t,
) -> zbar_symbol_type_t {
    // if same partial is not consistent, reset others
    if (ean.left != 0 && ((part & ZBAR_SYMBOL) != ean.left))
        || (ean.right != 0 && ((part & ZBAR_SYMBOL) != ean.right))
    {
        // partial mismatch - reset collected parts
        ean.left = ZBAR_NONE;
        ean.right = ZBAR_NONE;
    }

    if (ean.left != 0 || ean.right != 0) && check_width(ean.width, pass.width) == 0 {
        ean.left = ZBAR_NONE;
        ean.right = ZBAR_NONE;
    }

    if (part & EAN_RIGHT) != 0 {
        part &= ZBAR_SYMBOL;
        let mut j = part - 1;
        let mut i = part >> 1;
        while i > 0 {
            let digit = (pass.raw[i as usize] & 0xf) as c_char;
            if ean.right != 0 && ean.buf[j as usize] != digit {
                // partial mismatch - reset collected parts
                ean.left = ZBAR_NONE;
                ean.right = ZBAR_NONE;
            }
            ean.buf[j as usize] = digit;
            i -= 1;
            j -= 1;
        }
        ean.right = part;
        part &= ean.left; // FIXME!?
    } else if part == ZBAR_EAN13 || part == ZBAR_EAN8 {
        // EAN_LEFT
        let mut j = (part - 1) >> 1;
        let mut i = part >> 1;
        while j >= 0 {
            let digit = (pass.raw[i as usize] & 0xf) as c_char;
            if ean.left != 0 && ean.buf[j as usize] != digit {
                // partial mismatch - reset collected parts
                ean.left = ZBAR_NONE;
                ean.right = ZBAR_NONE;
            }
            ean.buf[j as usize] = digit;
            i -= 1;
            j -= 1;
        }
        ean.left = part;
        part &= ean.right; // FIXME!?
    } else if part != ZBAR_UPCE {
        // add-ons
        for i in (1..=(part as usize)).rev() {
            ean.buf[i - 1] = (pass.raw[i] & 0xf) as c_char;
        }
        ean.left = part;
    } else {
        ean_expand_upce(ean, pass);
    }

    ean.width = pass.width;

    if part == 0 {
        part = ZBAR_PARTIAL;
    }

    if ((part == ZBAR_EAN13 || part == ZBAR_UPCE) && ean_verify_checksum(ean, 12) != 0)
        || (part == ZBAR_EAN8 && ean_verify_checksum(ean, 7) != 0)
    {
        // invalid checksum
        if ean.right != 0 {
            ean.left = ZBAR_NONE;
        } else {
            ean.right = ZBAR_NONE;
        }
        part = ZBAR_NONE;
    }

    if part == ZBAR_EAN13 {
        // special case EAN-13 subsets
        if ean.buf[0] == 0 && test_cfg(ean.upca_config, ZBAR_CFG_ENABLE) {
            part = ZBAR_UPCA;
        } else if ean.buf[0] == 9 && ean.buf[1] == 7 {
            if (ean.buf[2] == 8 || ean.buf[2] == 9) && test_cfg(ean.isbn13_config, ZBAR_CFG_ENABLE)
            {
                part = ZBAR_ISBN13;
            } else if ean.buf[2] == 8 && test_cfg(ean.isbn10_config, ZBAR_CFG_ENABLE) {
                part = ZBAR_ISBN10;
            }
        }
    } else if part == ZBAR_UPCE {
        if test_cfg(ean.upce_config, ZBAR_CFG_ENABLE) {
            // UPC-E was decompressed for checksum verification,
            // but user requested compressed result
            ean.buf[0] = 0;
            ean.buf[1] = 0;
            for i in 2..8 {
                ean.buf[i] = (pass.raw[i - 1] & 0xf) as c_char;
            }
            ean.buf[8] = (pass.raw[0] & 0xf) as c_char;
        } else if test_cfg(ean.upca_config, ZBAR_CFG_ENABLE) {
            // UPC-E reported as UPC-A has priority over EAN-13
            part = ZBAR_UPCA;
        } else if test_cfg(ean.ean13_config, ZBAR_CFG_ENABLE) {
            part = ZBAR_EAN13;
        } else {
            part = ZBAR_NONE;
        }
    }

    part
}

/// Copy result to output buffer
fn postprocess(dcode: &mut zbar_decoder_t, sym: zbar_symbol_type_t) {
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

        if base > ZBAR_PARTIAL {
            match base {
                ZBAR_UPCA | ZBAR_UPCE => i = 1,
                ZBAR_ISBN10 => i = 3,
                _ => {}
            }

            let mut calc_base = base;
            match base {
                ZBAR_ISBN13 => calc_base = ZBAR_EAN13,
                ZBAR_UPCE => calc_base -= 1,
                _ => {}
            }

            if base == ZBAR_ISBN10
                || (calc_base > ZBAR_EAN5
                    && !test_cfg(ean_get_config(ean, sym), ZBAR_CFG_EMIT_CHECK))
            {
                calc_base -= 1;
            }

            // Copy to temp buffer
            while j < calc_base as usize && j < 18 && ean.buf[i] >= 0 {
                temp_buf[j] = ean.buf[i];
                i += 1;
                j += 1;
            }
            buf_len = j;

            needs_isbn10_check =
                base == ZBAR_ISBN10 && j == 9 && test_cfg(ean.isbn10_config, ZBAR_CFG_EMIT_CHECK);
            isbn10_check_digit = if needs_isbn10_check {
                isbn10_calc_checksum(ean) as u8
            } else {
                0
            };
        } else {
            needs_isbn10_check = false;
            isbn10_check_digit = 0;
        }

        new_direction = 1 - 2 * ean.direction;
    }

    if base > ZBAR_PARTIAL {
        // Get mutable slice for writing
        let buffer = unsafe {
            match dcode.buffer_mut_slice(buf_len + if needs_isbn10_check { 2 } else { 1 }) {
                Ok(buf) => buf,
                Err(_) => return,
            }
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
fn decode_pass(dcode: &mut zbar_decoder_t, pass: &mut ean_pass_t) -> zbar_symbol_type_t {
    pass.state = pass.state.wrapping_add(1);
    let idx = pass.state & STATE_IDX;
    let fwd = (pass.state & 1) as u8;

    if dcode.color() == zbar_color_t::ZBAR_SPACE {
        if (pass.state & STATE_ADDON) != 0 {
            if idx == 0x09 || idx == 0x21 {
                let qz = get_width(dcode, 0);
                let s = calc_s(dcode, 1, 4);
                let part = if qz == 0 || qz >= s * 3 / 4 {
                    if idx == 0x09 {
                        ean_part_end2(&dcode.ean, pass)
                    } else {
                        ean_part_end5(&dcode.ean, pass)
                    }
                } else {
                    ZBAR_NONE
                };

                if part != 0 || idx == 0x21 {
                    let ean = &mut dcode.ean;
                    ean.direction = 0;
                    pass.state = -1;
                    return part;
                }
            }
            if (idx & 7) == 1 {
                pass.state += 2;
            }
        } else if (idx == 0x10 || idx == 0x11)
            && test_cfg(dcode.ean.ean8_config, ZBAR_CFG_ENABLE)
            && aux_end(dcode, fwd) == 0
        {
            let part = ean_part_end4(pass, fwd);
            if part != 0 {
                let ean = &mut dcode.ean;
                ean.direction = if (pass.state & STATE_REV) != 0 { 1 } else { 0 };
            }
            pass.state = -1;
            return part;
        } else if idx == 0x18 || idx == 0x19 {
            let mut part = ZBAR_NONE;
            if aux_end(dcode, fwd) == 0 && pass.raw[5] != 0xff {
                part = ean_part_end7(&dcode.ean, pass, fwd);
            }
            if part != 0 {
                let ean = &mut dcode.ean;
                ean.direction = if (pass.state & STATE_REV) != 0 { 1 } else { 0 };
            }
            pass.state = -1;
            return part;
        }
    }

    let mut idx = idx;
    if (pass.state & STATE_ADDON) != 0 {
        idx >>= 1;
    }

    if (idx & 0x03) == 0 && idx <= 0x14 {
        let mut code: i8 = -1;
        let mut w = pass.width;

        if dcode.ean.s4 == 0 {
            return 0;
        }

        // validate guard bars before decoding first char of symbol
        if pass.state == 0 {
            pass.state = aux_start(dcode);
            pass.width = dcode.ean.s4;
            if pass.state < 0 {
                return 0;
            }
            idx = pass.state & STATE_IDX;
        } else {
            w = check_width(w, dcode.ean.s4);
            if w != 0 {
                pass.width = (pass.width + dcode.ean.s4 * 3) / 4;
            }
        }

        if w != 0 {
            code = decode4(dcode);
        }

        if (code < 0 && idx != 0x10)
            || (idx > 0 && (pass.state & STATE_ADDON) != 0 && aux_mid(dcode) != 0)
        {
            pass.state = -1;
        } else if code < 0 {
            pass.raw[5] = 0xff;
        } else {
            let raw_idx = ((idx >> 2) + 1) as usize;
            pass.raw[raw_idx] = DIGITS[code as usize];
        }
    }

    0
}

// ============================================================================
// Main decoder function
// ============================================================================

/// Main EAN/UPC decoder entry point
pub unsafe fn _zbar_decode_ean(dcode: *mut zbar_decoder_t) -> zbar_symbol_type_t {
    let dcode = &mut *dcode;

    // process up to 4 separate passes
    let mut sym = ZBAR_NONE;
    let pass_idx = (dcode.idx & 3) as usize;

    // update latest character width
    dcode.ean.s4 = dcode.ean.s4.wrapping_sub(get_width(dcode, 4));
    dcode.ean.s4 = dcode.ean.s4.wrapping_add(get_width(dcode, 0));

    for i in 0..4 {
        let pass = unsafe { &mut *(&mut dcode.ean.pass[i] as *mut ean_pass_t) };
        if pass.state >= 0 || i == pass_idx {
            let part = decode_pass(dcode, pass);
            if part != 0 {
                // update accumulated data from new partial decode
                let ean = unsafe { &mut *(&mut dcode.ean as *mut ean_decoder_t) };
                let pass = unsafe { &mut *(&mut dcode.ean.pass[i] as *mut ean_pass_t) };
                sym = integrate_partial(ean, pass, part);
                if sym != 0 {
                    // this pass valid => _reset_ all passes
                    dcode.ean.pass[0].state = -1;
                    dcode.ean.pass[1].state = -1;
                    dcode.ean.pass[2].state = -1;
                    dcode.ean.pass[3].state = -1;
                    if sym > ZBAR_PARTIAL {
                        if !acquire_lock(dcode, sym) {
                            postprocess(dcode, sym);
                        } else {
                            sym = ZBAR_PARTIAL;
                        }
                    }
                }
            }
        }
    }

    sym
}
