//! Code 93 barcode decoder
//!
//! This module implements decoding for Code 93 barcodes.

use crate::{
    decoder_types::{
        code93_decoder_t, zbar_decoder_t, zbar_symbol_type_t, DECODE_WINDOW, ZBAR_CFG_MAX_LEN,
        ZBAR_CFG_MIN_LEN, ZBAR_CODE93, ZBAR_NONE, ZBAR_PARTIAL,
    },
    line_scanner::zbar_color_t,
};
use libc::{c_int, c_uint};

// Checksum constant
const CHKMOD: i32 = 47;

// Assertion macro
macro_rules! zassert {
    ($condition:expr, $retval:expr, $($arg:tt)*) => {
        if !$condition {
            return $retval;
        }
    };
}

// ============================================================================
// Code 93 lookup table
// ============================================================================

static CODE93_HASH: [i8; 0x40] = [
    0x0f, 0x2b, 0x30, 0x38, 0x13, 0x1b, 0x11, 0x2a, 0x0a, -1, 0x2f, 0x0f, 0x38, 0x38, 0x2f, 0x37,
    0x24, 0x3a, 0x1b, 0x36, 0x18, 0x26, 0x02, 0x2c, 0x2b, 0x05, 0x21, 0x3b, 0x04, 0x15, 0x12, 0x0c,
    0x00, 0x26, 0x23, 0x00, -1, 0x2e, 0x3f, 0x13, 0x2e, 0x36, -1, 0x08, 0x09, -1, 0x15, 0x14, -1,
    0x00, 0x21, 0x3b, -1, 0x33, 0x00, -1, 0x2d, 0x0c, 0x1b, 0x0a, 0x3f, 0x3f, 0x29, 0x1c,
];

static CODE93_GRAPH: &[u8; 7] = b"-. $/+%";
static CODE93_S2: &[u8; 26] = b"\x1b\x1c\x1d\x1e\x1f;<=>?[\\]^_{|}~\x7f\x00\x40`\x7f\x7f\x7f";

// ============================================================================
// Helper functions from decoder.h
// ============================================================================

/// Retrieve i-th previous element width
#[inline]
fn get_width(dcode: &zbar_decoder_t, offset: u8) -> c_uint {
    dcode.w[((dcode.idx.wrapping_sub(offset)) & (DECODE_WINDOW as u8 - 1)) as usize]
}

/// Retrieve i-th previous pair width
#[inline]
fn pair_width(dcode: &zbar_decoder_t, offset: u8) -> c_uint {
    get_width(dcode, offset) + get_width(dcode, offset + 1)
}

/// Fixed character width decode assist
#[inline]
fn decode_e(e: c_uint, s: c_uint, n: c_uint) -> c_uint {
    ((e * n * 2 + 1) / s).wrapping_sub(3) / 2
}

/// Acquire shared state lock
#[inline]
fn acquire_lock(dcode: &mut zbar_decoder_t, req: zbar_symbol_type_t) -> bool {
    if dcode.lock != 0 {
        return true;
    }
    dcode.lock = req;
    false
}

/// Check and release shared state lock
#[inline]
fn release_lock(dcode: &mut zbar_decoder_t, req: zbar_symbol_type_t) -> i8 {
    zassert!(dcode.lock == req, 1, "lock={} req={}\n", dcode.lock, req);
    dcode.lock = 0;
    0
}

/// Access config value by index
#[inline]
fn cfg(decoder: &code93_decoder_t, cfg: c_int) -> c_int {
    decoder.configs[(cfg - ZBAR_CFG_MIN_LEN) as usize]
}

// ============================================================================
// Code 93 Decoder functions
// ============================================================================

/// Check width variance
#[inline]
fn check_width(cur: c_uint, prev: c_uint) -> bool {
    let dw = prev.abs_diff(cur);
    let dw = dw * 4;
    dw > prev
}

/// Build edge signature of character
#[inline]
fn encode6(dcode: &zbar_decoder_t) -> i32 {
    let s = dcode.s6;
    if s < 9 {
        return -1;
    }

    let mut sig: i32 = 0;
    for i in (1..6).rev() {
        let c = decode_e(pair_width(dcode, i), s, 9);
        if c > 3 {
            return -1;
        }
        sig = (sig << 2) | (c as i32);
    }

    sig
}

/// Validate signature
#[inline]
fn validate_sig(mut sig: i32) -> i32 {
    let mut sum = 0;
    let mut emin = 0;
    let mut sig0 = 0;
    let mut sig1 = 0;

    for i in (0..3).rev() {
        let e = sig & 3;
        sig >>= 2;
        sum = e - sum;
        sig1 <<= 4;
        sig1 += sum;
        if i == 0 {
            break;
        }

        let e = sig & 3;
        sig >>= 2;
        sum = e - sum;
        sig0 <<= 4;
        if emin > sum {
            emin = sum;
        }
        sig0 += sum;
    }

    emin = emin + (emin << 4) + (emin << 8);
    sig0 -= emin;
    sig1 += emin;
    (sig0 | sig1) & 0x888
}

/// Decode 6 elements
#[inline]
fn decode6(dcode: &zbar_decoder_t) -> i32 {
    let mut sig = encode6(dcode);
    if sig < 0 {
        return -1;
    }

    if (sig & 0x3) + ((sig >> 4) & 0x3) + ((sig >> 8) & 0x3) != 3 || validate_sig(sig) != 0 {
        return -1;
    }

    if dcode.code93.direction() {
        // Reverse signature
        let tmp = sig & 0x030;
        sig = ((sig & 0x3c0) >> 6) | ((sig & 0x00f) << 6);
        sig = ((sig & 0x30c) >> 2) | ((sig & 0x0c3) << 2) | tmp;
    }

    let g0 = CODE93_HASH[((sig - (sig >> 4)) & 0x3f) as usize];
    let g1 = CODE93_HASH[(((sig >> 2) - (sig >> 7)) & 0x3f) as usize];

    zassert!(
        g0 >= 0 && g1 >= 0,
        -1,
        "dir={:x} sig={:03x} g0={:03x} g1={:03x}\n",
        dcode.code93.direction() as i32,
        sig,
        g0,
        g1
    );

    let c = (g0 + g1) & 0x3f;
    c as i32
}

/// Decode start pattern
#[inline]
fn decode_start(dcode: &mut zbar_decoder_t) -> zbar_symbol_type_t {
    let s = dcode.s6;
    let c = encode6(dcode);
    if c < 0 || (c != 0x00f && c != 0x0f0) {
        return ZBAR_NONE;
    }

    let dir = (c >> 7) != 0;

    let qz = if dir {
        if decode_e(pair_width(dcode, 0), s, 9) != 0 {
            return ZBAR_NONE;
        }
        get_width(dcode, 8)
    } else {
        get_width(dcode, 7)
    };

    if qz != 0 && qz < (s * 3) / 4 {
        return ZBAR_NONE;
    }

    // Decoded valid start/stop - initialize state
    dcode.code93.set_direction(dir);
    dcode.code93.set_element(if !dir { 0 } else { 7 });
    dcode.code93.set_character(0);
    dcode.code93.width = s;
    ZBAR_PARTIAL
}

/// Abort decoding
#[inline]
fn decode_abort(dcode: &mut zbar_decoder_t) -> zbar_symbol_type_t {
    if dcode.code93.character() > 1 {
        release_lock(dcode, ZBAR_CODE93);
    }
    dcode.code93.set_character(-1);
    ZBAR_NONE
}

/// Check stop pattern
#[inline]
fn check_stop(dcode: &zbar_decoder_t) -> bool {
    let n = dcode.code93.character() as i32;
    let s = dcode.s6;
    let max_len = cfg(&dcode.code93, ZBAR_CFG_MAX_LEN);

    if n < 2 || n < cfg(&dcode.code93, ZBAR_CFG_MIN_LEN) || (max_len != 0 && n > max_len) {
        return false;
    }

    if dcode.code93.direction() {
        let qz = get_width(dcode, 0);
        if qz != 0 && qz < (s * 3) / 4 {
            return false;
        }
    } else if decode_e(pair_width(dcode, 0), s, 9) != 0 {
        return false;
    }

    true
}

/// Plus modulo 47
#[inline]
fn plusmod47(mut acc: i32, add: i32) -> i32 {
    acc += add;
    if acc >= CHKMOD {
        acc -= CHKMOD;
    }
    acc
}

/// Validate checksums
#[inline]
unsafe fn validate_checksums(dcode: &zbar_decoder_t) -> bool {
    let n = dcode.code93.character() as usize;
    let buf = dcode.buffer_slice();

    // Ensure buffer has enough data
    if buf.len() < n {
        return false;
    }

    let mut sum_c = 0;
    let mut acc_c = 0;
    let mut i_c = ((n - 2) % 20) as i32;
    let mut sum_k = 0;
    let mut acc_k = 0;
    let mut i_k = ((n - 1) % 15) as i32;

    for i in 0..(n - 2) {
        let d = if dcode.code93.direction() {
            buf[n - 1 - i] as i32
        } else {
            buf[i] as i32
        };

        i_c -= 1;
        if i_c < 0 {
            acc_c = 0;
            i_c = 19;
        }
        acc_c = plusmod47(acc_c, d);
        sum_c = plusmod47(sum_c, acc_c);

        i_k -= 1;
        if i_k < 0 {
            acc_k = 0;
            i_k = 14;
        }
        acc_k = plusmod47(acc_k, d);
        sum_k = plusmod47(sum_k, acc_k);
    }

    let d = if dcode.code93.direction() {
        buf[1] as i32
    } else {
        buf[n - 2] as i32
    };
    if d != sum_c {
        return false;
    }

    acc_k = plusmod47(acc_k, sum_c);
    sum_k = plusmod47(sum_k, acc_k);
    let d = if dcode.code93.direction() {
        buf[0] as i32
    } else {
        buf[n - 1] as i32
    };
    if d != sum_k {
        return false;
    }

    true
}

/// Resolve scan direction and convert to ASCII
#[inline]
unsafe fn postprocess(dcode: &mut zbar_decoder_t) -> bool {
    let n = dcode.code93.character() as usize;
    let direction = dcode.code93.direction();

    dcode.direction = 1 - 2 * (direction as c_int);

    // Get mutable slice for processing
    let buffer = match dcode.buffer_mut_slice(n) {
        Ok(buf) => buf,
        Err(_) => return true,
    };

    if direction {
        // Reverse buffer
        buffer[..n].reverse();
    }

    let n = n - 2;
    let mut i = 0;
    let mut j = 0;
    while i < n {
        let mut d = buffer[i];
        i += 1;

        if d < 0xa {
            d += b'0';
        } else if d < 0x24 {
            d = b'A' + d - 0xa;
        } else if d < 0x2b {
            d = CODE93_GRAPH[(d - 0x24) as usize];
        } else {
            let shift = d;
            zassert!(shift < 0x2f, true, "shift={:02x}\n", shift);
            d = buffer[i];
            i += 1;
            if !(0xa..0x24).contains(&d) {
                return true;
            }
            d -= 0xa;
            match shift {
                0x2b => d += 1,
                0x2c => d = CODE93_S2[d as usize],
                0x2d => d += 0x21,
                0x2e => d += 0x61,
                _ => return true,
            }
        }
        buffer[j] = d;
        j += 1;
    }

    // Add null terminator and truncate
    buffer[j] = 0;
    dcode.truncate_buffer(j);

    dcode.modifiers = 0;
    false
}

/// Main Code 93 decode function
pub unsafe fn _zbar_decode_code93(dcode: *mut zbar_decoder_t) -> zbar_symbol_type_t {
    let dcode = &mut *dcode;

    if dcode.code93.character() < 0 {
        if dcode.color() != zbar_color_t::ZBAR_BAR {
            return ZBAR_NONE;
        }
        return decode_start(dcode);
    }

    // Process every 6th element of active symbol
    let element = dcode.code93.element() + 1;
    dcode.code93.set_element(element);

    if element != 6 || dcode.color() as u8 == (dcode.code93.direction() as u8) {
        return ZBAR_NONE;
    }

    dcode.code93.set_element(0);

    if check_width(dcode.s6, dcode.code93.width) {
        return decode_abort(dcode);
    }

    let c = decode6(dcode);
    if c < 0 {
        return decode_abort(dcode);
    }

    if c == 0x2f {
        if !check_stop(dcode) {
            return ZBAR_NONE;
        }
        if !validate_checksums(dcode) {
            return decode_abort(dcode);
        }
        if postprocess(dcode) {
            return decode_abort(dcode);
        }

        dcode.code93.set_character(-1);
        return ZBAR_CODE93;
    }

    let character = dcode.code93.character();
    if dcode
        .set_buffer_capacity((character + 1) as c_uint)
        .is_err()
    {
        return decode_abort(dcode);
    }

    dcode.code93.width = dcode.s6;

    if character == 1 {
        // Lock shared resources
        if acquire_lock(dcode, ZBAR_CODE93) {
            return decode_abort(dcode);
        }
        // Copy from holding buffer
        if dcode.write_buffer_byte(0, dcode.code93.buf).is_err() {
            return decode_abort(dcode);
        }
    }

    if character == 0 {
        dcode.code93.buf = c as u8;
    } else {
        if dcode
            .write_buffer_byte(character as usize, c as u8)
            .is_err()
        {
            return decode_abort(dcode);
        }
    }
    dcode.code93.set_character(character + 1);
    ZBAR_NONE
}
