//! Interleaved 2 of 5 barcode decoder
//!
//! This module implements decoding for Interleaved 2 of 5 (I25) barcodes.

use crate::decoder_types::{
    i25_decoder_t, zbar_decoder_t, zbar_symbol_type_t, BUFFER_INCR, BUFFER_MAX, BUFFER_MIN,
    DECODE_WINDOW, ZBAR_BAR, ZBAR_CFG_MAX_LEN, ZBAR_CFG_MIN_LEN, ZBAR_I25, ZBAR_NONE, ZBAR_PARTIAL,
};
use libc::{c_char, c_int, c_uint};

// Assertion macro
macro_rules! zassert {
    ($condition:expr, $retval:expr, $($arg:tt)*) => {
        if !$condition {
            return $retval;
        }
    };
}

// ============================================================================
// Helper functions from decoder.h
// ============================================================================

/// Return current element color
#[inline]
fn get_color(dcode: &zbar_decoder_t) -> u8 {
    dcode.idx & 1
}

/// Retrieve i-th previous element width
#[inline]
fn get_width(dcode: &zbar_decoder_t, offset: u8) -> c_uint {
    dcode.w[((dcode.idx.wrapping_sub(offset)) & (DECODE_WINDOW as u8 - 1)) as usize]
}

/// Fixed character width decode assist
#[inline]
fn decode_e(e: c_uint, s: c_uint, n: c_uint) -> i32 {
    let e_val = ((e * n * 2 + 1) / s).wrapping_sub(3) / 2;
    if e_val >= n - 3 {
        -1
    } else {
        e_val as i32
    }
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

/// Ensure output buffer has sufficient allocation for request
#[inline]
unsafe fn size_buf(dcode: &mut zbar_decoder_t, mut len: c_uint) -> i8 {
    if len <= BUFFER_MIN {
        return 0;
    }
    if len < dcode.buf_alloc {
        return 0;
    }
    if len > BUFFER_MAX {
        return 1;
    }
    if len < dcode.buf_alloc + BUFFER_INCR {
        len = dcode.buf_alloc + BUFFER_INCR;
        if len > BUFFER_MAX {
            len = BUFFER_MAX;
        }
    }
    let new_buf = libc::realloc(dcode.buf as *mut libc::c_void, len as usize) as *mut c_char;
    if new_buf.is_null() {
        return 1;
    }
    dcode.buf = new_buf;
    dcode.buf_alloc = len;
    0
}

/// Access config value by index
#[inline]
fn cfg(decoder: &i25_decoder_t, cfg: c_int) -> c_int {
    decoder.configs[(cfg - ZBAR_CFG_MIN_LEN) as usize]
}

// ============================================================================
// I25 Decoder functions
// ============================================================================

/// Decode a single element
#[inline]
fn i25_decode1(enc: u8, e: c_uint, s: c_uint) -> u8 {
    let e_val = decode_e(e, s, 45);
    if e_val > 7 {
        return 0xff;
    }
    let mut enc = enc << 1;
    if e_val > 2 {
        enc |= 1;
    }
    enc
}

/// Decode 10 elements into a digit (5 bars + 5 spaces)
#[inline]
fn i25_decode10(dcode: &zbar_decoder_t, offset: u8) -> u8 {
    let dcode25 = &dcode.i25;

    if dcode25.s10 < 10 {
        return 0xff;
    }

    let mut enc: u8 = 0;
    let mut par: u8 = 0;

    // Threshold bar width ratios
    let mut i: i8 = 8;
    while i >= 0 {
        let j = offset
            + if dcode25.direction() {
                i as u8
            } else {
                8 - i as u8
            };
        enc = i25_decode1(enc, get_width(dcode, j), dcode25.s10);
        if enc == 0xff {
            return 0xff;
        }
        if enc & 1 != 0 {
            par += 1;
        }
        i -= 2;
    }

    // Parity check
    if par != 2 {
        return 0xff;
    }

    // Decode binary weights
    enc &= 0xf;
    if enc & 8 != 0 {
        if enc == 12 {
            enc = 0;
        } else {
            enc -= 1;
            if enc > 9 {
                return 0xff;
            }
        }
    }

    enc
}

/// Decode start pattern
#[inline]
fn i25_decode_start(dcode: &mut zbar_decoder_t) -> zbar_symbol_type_t {
    let s10 = dcode.i25.s10;

    if s10 < 10 {
        return ZBAR_NONE;
    }

    let mut enc: u8 = 0;
    let mut i: u8 = 10;

    enc = i25_decode1(enc, get_width(dcode, i), s10);
    i += 1;
    enc = i25_decode1(enc, get_width(dcode, i), s10);
    i += 1;
    enc = i25_decode1(enc, get_width(dcode, i), s10);
    i += 1;

    let valid = if get_color(dcode) == ZBAR_BAR {
        enc == 4
    } else {
        enc = i25_decode1(enc, get_width(dcode, i), s10);
        i += 1;
        enc == 0
    };

    if !valid {
        return ZBAR_NONE;
    }

    // Check leading quiet zone - spec is 10n(?)
    // we require 5.25n for w=2n to 6.75n for w=3n
    let quiet = get_width(dcode, i);
    if quiet != 0 && quiet < s10 * 3 / 8 {
        return ZBAR_NONE;
    }

    dcode.i25.set_direction(get_color(dcode) != 0);
    dcode.i25.set_element(1);
    dcode.i25.set_character(0);
    ZBAR_PARTIAL
}

/// Acquire lock and copy holding buffer
#[inline]
unsafe fn i25_acquire_lock(dcode: &mut zbar_decoder_t) -> bool {
    // Lock shared resources
    if acquire_lock(dcode, ZBAR_I25) {
        dcode.i25.set_character(-1);
        return true;
    }

    // Copy holding buffer
    for i in 0..4 {
        *dcode.buf.offset(i) = dcode.i25.buf[i as usize] as c_char;
    }
    false
}

/// Decode end pattern and validate
#[inline]
unsafe fn i25_decode_end(dcode: &mut zbar_decoder_t) -> zbar_symbol_type_t {
    let width = dcode.i25.width;
    let direction = dcode.i25.direction();
    let character = dcode.i25.character();

    // Check trailing quiet zone
    let quiet = get_width(dcode, 0);
    if (quiet != 0 && quiet < width * 3 / 8)
        || decode_e(get_width(dcode, 1), width, 45) > 2
        || decode_e(get_width(dcode, 2), width, 45) > 2
    {
        return ZBAR_NONE;
    }

    // Check exit condition
    let e = decode_e(get_width(dcode, 3), width, 45);
    let valid = if !direction {
        (e as u32).wrapping_sub(3) <= 4
    } else {
        e <= 2 && decode_e(get_width(dcode, 4), width, 45) <= 2
    };

    if !valid {
        return ZBAR_NONE;
    }

    if character <= 4 && i25_acquire_lock(dcode) {
        return ZBAR_PARTIAL;
    }

    dcode.direction = 1 - 2 * (direction as c_int);
    if direction {
        // Reverse buffer
        let char_count = character as isize;
        for i in 0..(char_count / 2) {
            let j = char_count - 1 - i;
            let c = *dcode.buf.offset(i);
            *dcode.buf.offset(i) = *dcode.buf.offset(j);
            *dcode.buf.offset(j) = c;
        }
    }

    let min_len = cfg(&dcode.i25, ZBAR_CFG_MIN_LEN);
    let max_len = cfg(&dcode.i25, ZBAR_CFG_MAX_LEN);

    if character < min_len as i16 || (max_len > 0 && character > max_len as i16) {
        release_lock(dcode, ZBAR_I25);
        dcode.i25.set_character(-1);
        return ZBAR_NONE;
    }

    let char_count = character as c_uint;
    zassert!(
        char_count < dcode.buf_alloc,
        ZBAR_NONE,
        "i={:02x}\n",
        char_count
    );

    dcode.buflen = char_count;
    *dcode.buf.offset(char_count as isize) = 0;
    dcode.modifiers = 0;
    dcode.i25.set_character(-1);
    ZBAR_I25
}

/// Main I25 decode function
#[no_mangle]
pub unsafe extern "C" fn _zbar_decode_i25(dcode: *mut zbar_decoder_t) -> zbar_symbol_type_t {
    let dcode = &mut *dcode;

    // Update latest character width
    let w10 = get_width(dcode, 10);
    let w0 = get_width(dcode, 0);
    dcode.i25.s10 = dcode.i25.s10.wrapping_sub(w10).wrapping_add(w0);

    if dcode.i25.character() < 0 && i25_decode_start(dcode) == ZBAR_NONE {
        return ZBAR_NONE;
    }

    let element = dcode.i25.element();
    dcode.i25.set_element(element.wrapping_sub(1));

    let check_element = 6u8.wrapping_sub(dcode.i25.direction() as u8);
    if dcode.i25.element() == check_element {
        return i25_decode_end(dcode);
    } else if dcode.i25.element() != 0 {
        return ZBAR_NONE;
    }

    // FIXME check current character width against previous
    dcode.i25.width = dcode.i25.s10;

    let _direction = dcode.i25.direction();
    let _character = dcode.i25.character();
    let _element = dcode.i25.element();

    if dcode.i25.character() == 4 && i25_acquire_lock(dcode) {
        return ZBAR_PARTIAL;
    }

    let c = i25_decode10(dcode, 1);
    if c > 9 {
        // goto reset
        if dcode.i25.character() >= 4 {
            release_lock(dcode, ZBAR_I25);
        }
        dcode.i25.set_character(-1);
        return ZBAR_NONE;
    }

    if size_buf(dcode, dcode.i25.character() as c_uint + 3) != 0 {
        // goto reset
        if dcode.i25.character() >= 4 {
            release_lock(dcode, ZBAR_I25);
        }
        dcode.i25.set_character(-1);
        return ZBAR_NONE;
    }

    let character = dcode.i25.character();
    let buf = if character >= 4 {
        dcode.buf
    } else {
        dcode.i25.buf.as_mut_ptr() as *mut c_char
    };

    *buf.offset(character as isize) = (c + b'0') as c_char;
    dcode.i25.set_character(character + 1);

    let c = i25_decode10(dcode, 0);
    if c > 9 {
        // goto reset
        if dcode.i25.character() >= 4 {
            release_lock(dcode, ZBAR_I25);
        }
        dcode.i25.set_character(-1);
        return ZBAR_NONE;
    }

    let character = dcode.i25.character();
    *buf.offset(character as isize) = (c + b'0') as c_char;
    dcode.i25.set_character(character + 1);
    dcode.i25.set_element(10);

    if dcode.i25.character() == 2 {
        ZBAR_PARTIAL
    } else {
        ZBAR_NONE
    }
}
