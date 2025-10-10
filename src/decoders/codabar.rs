//! Codabar barcode decoder
//!
//! This module implements decoding for Codabar barcodes.

use crate::decoder_types::{
    codabar_decoder_t, zbar_decoder_t, zbar_symbol_type_t, BUFFER_INCR, BUFFER_MAX, BUFFER_MIN,
    DECODE_WINDOW, ZBAR_CFG_ADD_CHECK, ZBAR_CFG_EMIT_CHECK, ZBAR_CFG_MAX_LEN, ZBAR_CFG_MIN_LEN,
    ZBAR_CODABAR, ZBAR_NONE, ZBAR_PARTIAL, ZBAR_SPACE,
};
use libc::{c_char, c_int, c_uint};

// Buffer constants
const NIBUF: usize = 6;

// Assertion macro
macro_rules! zassert {
    ($condition:expr, $retval:expr, $($arg:tt)*) => {
        if !$condition {
            return $retval;
        }
    };
}

// ============================================================================
// Codabar lookup tables
// ============================================================================

static CODABAR_LO: [i8; 12] = [0x0, 0x1, 0x4, 0x5, 0x2, 0xa, 0xb, 0x9, 0x6, 0x7, 0x8, 0x3];

static CODABAR_HI: [u8; 8] = [0x1, 0x4, 0x7, 0x6, 0x2, 0x3, 0x0, 0x5];

static CODABAR_CHARACTERS: &[u8; 20] = b"0123456789-$:/.+ABCD";

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

/// Sort 3 like-colored elements and return ordering
#[inline]
fn decode_sort3(dcode: &zbar_decoder_t, i0: u8) -> c_uint {
    let w0 = get_width(dcode, i0);
    let w2 = get_width(dcode, i0 + 2);
    let w4 = get_width(dcode, i0 + 4);

    if w0 < w2 {
        if w2 < w4 {
            ((i0 as c_uint) << 8) | (((i0 + 2) as c_uint) << 4) | ((i0 + 4) as c_uint)
        } else if w0 < w4 {
            ((i0 as c_uint) << 8) | (((i0 + 4) as c_uint) << 4) | ((i0 + 2) as c_uint)
        } else {
            (((i0 + 4) as c_uint) << 8) | ((i0 as c_uint) << 4) | ((i0 + 2) as c_uint)
        }
    } else if w4 < w2 {
        (((i0 + 4) as c_uint) << 8) | (((i0 + 2) as c_uint) << 4) | (i0 as c_uint)
    } else if w0 < w4 {
        (((i0 + 2) as c_uint) << 8) | ((i0 as c_uint) << 4) | ((i0 + 4) as c_uint)
    } else {
        (((i0 + 2) as c_uint) << 8) | (((i0 + 4) as c_uint) << 4) | (i0 as c_uint)
    }
}

/// Sort N like-colored elements and return ordering
#[inline]
fn decode_sortn(dcode: &zbar_decoder_t, n: i32, i0: u8) -> c_uint {
    let mut mask: c_uint = 0;
    let mut sort: c_uint = 0;

    for _i in (0..n).rev() {
        let mut wmin = c_uint::MAX;
        let mut jmin: i32 = -1;

        for j in (0..n).rev() {
            if (mask >> j) & 1 != 0 {
                continue;
            }
            let w = get_width(dcode, i0 + (j as u8) * 2);
            if wmin >= w {
                wmin = w;
                jmin = j;
            }
        }
        zassert!(jmin >= 0, 0, "sortn({},{}) jmin={}", n, i0, jmin);
        sort <<= 4;
        mask |= 1 << jmin;
        sort |= (i0 as c_uint) + (jmin as c_uint) * 2;
    }
    sort
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
fn cfg(decoder: &codabar_decoder_t, cfg: c_int) -> c_int {
    decoder.configs[(cfg - ZBAR_CFG_MIN_LEN) as usize]
}

/// Test config bit
#[inline]
fn test_cfg(config: c_uint, cfg: c_int) -> bool {
    ((config >> cfg) & 1) != 0
}

// ============================================================================
// Codabar Decoder functions
// ============================================================================

/// Check width against reference
#[inline]
fn check_width(ref_width: c_uint, w: c_uint) -> bool {
    let dref = ref_width;
    let ref_4 = ref_width * 4;
    let w_4 = w * 4;
    ref_4.wrapping_sub(dref) <= w_4 && w_4 <= ref_4.wrapping_add(dref)
}

/// Decode 7 elements into a character
#[inline]
fn codabar_decode7(dcode: &zbar_decoder_t) -> i8 {
    let codabar = &dcode.codabar;
    let s = codabar.s7;
    if s < 7 {
        return -1;
    }

    // Check width
    if !check_width(codabar.width, s) {
        return -1;
    }

    // Extract min/max bar
    let ibar = decode_sortn(dcode, 4, 1);

    let wbmax = get_width(dcode, (ibar & 0xf) as u8);
    let wbmin = get_width(dcode, (ibar >> 12) as u8);
    if 8 * wbmin < wbmax || 3 * wbmin > 2 * wbmax {
        return -1;
    }

    let wb1 = get_width(dcode, ((ibar >> 8) & 0xf) as u8);
    let wb2 = get_width(dcode, ((ibar >> 4) & 0xf) as u8);
    let b0b3 = wbmin * wbmax;
    let b1b2 = wb1 * wb2;

    let ibar = if b1b2 + b1b2 / 8 < b0b3 {
        // Single wide bar combinations
        if 8 * wbmin < 5 * wb1
            || 8 * wb1 < 5 * wb2
            || 4 * wb2 > 3 * wbmax
            || wb2 * wb2 >= wb1 * wbmax
        {
            return -1;
        }
        (ibar >> 1) & 0x3
    } else if b1b2 > b0b3 + b0b3 / 8 {
        // Three wide bars, no wide spaces
        if 4 * wbmin > 3 * wb1
            || 8 * wb1 < 5 * wb2
            || 8 * wb2 < 5 * wbmax
            || wbmin * wb2 >= wb1 * wb1
        {
            return -1;
        }
        (ibar >> 13) + 4
    } else {
        return -1;
    };

    let ispc = decode_sort3(dcode, 2);

    let wsmax = get_width(dcode, (ispc & 0xf) as u8);
    let wsmid = get_width(dcode, ((ispc >> 4) & 0xf) as u8);
    let wsmin = get_width(dcode, ((ispc >> 8) & 0xf) as u8);

    if ibar >> 2 != 0 {
        // Verify no wide spaces
        if 8 * wsmin < wsmax || 8 * wsmin < 5 * wsmid || 8 * wsmid < 5 * wsmax {
            return -1;
        }
        let mut ibar = ibar & 0x3;
        if codabar.direction() {
            ibar = 3 - ibar;
        }
        let c = (0xfcde >> (ibar << 2)) & 0xf;
        return c as i8;
    }
    if 8 * wsmin < wsmax || 3 * wsmin > 2 * wsmax {
        return -1;
    }

    let s0s2 = wsmin * wsmax;
    let s1s1 = wsmid * wsmid;

    if s1s1 + s1s1 / 8 < s0s2 {
        // Single wide space
        if 8 * wsmin < 5 * wsmid || 4 * wsmid > 3 * wsmax {
            return -1;
        }
        let ispc = (((ispc & 0xf) as u8) >> 1) - 1;
        let mut ic = ((ispc as c_uint) << 2) | ibar;
        if codabar.direction() {
            ic = 11 - ic;
        }

        CODABAR_LO[ic as usize]
    } else if s1s1 > s0s2 + s0s2 / 8 {
        // Two wide spaces, check start/stop
        if 4 * wsmin > 3 * wsmid || 8 * wsmid < 5 * wsmax {
            return -1;
        }
        if (ispc >> 8) == 4 {
            return -1;
        }
        let ispc = ispc >> 10;
        let ic = ispc * 4 + ibar;
        zassert!(ic < 8, -1, "ic={} ispc={} ibar={}", ic, ispc, ibar);
        let c = CODABAR_HI[ic as usize];
        if (c >> 2) != (codabar.direction() as u8) {
            return -1;
        }
        let c = (c & 0x3) | 0x10;
        return c as i8;
    } else {
        return -1;
    }
}

/// Decode start pattern
#[inline]
fn codabar_decode_start(dcode: &mut zbar_decoder_t) -> zbar_symbol_type_t {
    let s = dcode.codabar.s7;
    if s < 8 {
        return ZBAR_NONE;
    }

    // Check leading quiet zone - spec is 10x
    let qz = get_width(dcode, 8);
    if (qz != 0 && qz * 2 < s) || 4 * get_width(dcode, 0) > 3 * s {
        return ZBAR_NONE;
    }

    // Check space ratios first
    let ispc = decode_sort3(dcode, 2);
    if (ispc >> 8) == 4 {
        return ZBAR_NONE;
    }

    // Require 2 wide and 1 narrow spaces
    let wsmax = get_width(dcode, (ispc & 0xf) as u8);
    let wsmin = get_width(dcode, (ispc >> 8) as u8);
    let wsmid = get_width(dcode, ((ispc >> 4) & 0xf) as u8);
    if 8 * wsmin < wsmax
        || 3 * wsmin > 2 * wsmax
        || 4 * wsmin > 3 * wsmid
        || 8 * wsmid < 5 * wsmax
        || wsmid * wsmid <= wsmax * wsmin
    {
        return ZBAR_NONE;
    }
    let ispc = ispc >> 10;

    // Check bar ratios
    let ibar = decode_sortn(dcode, 4, 1);

    let wbmax = get_width(dcode, (ibar & 0xf) as u8);
    let wbmin = get_width(dcode, (ibar >> 12) as u8);
    if 8 * wbmin < wbmax || 3 * wbmin > 2 * wbmax {
        return ZBAR_NONE;
    }

    // Require 1 wide & 3 narrow bars
    let wb1 = get_width(dcode, ((ibar >> 8) & 0xf) as u8);
    let wb2 = get_width(dcode, ((ibar >> 4) & 0xf) as u8);
    if 8 * wbmin < 5 * wb1
        || 8 * wb1 < 5 * wb2
        || 4 * wb2 > 3 * wbmax
        || wb1 * wb2 >= wbmin * wbmax
        || wb2 * wb2 >= wb1 * wbmax
    {
        return ZBAR_NONE;
    }
    let ibar = (((ibar & 0xf) as u8) - 1) >> 1;

    // Decode combination
    let ic = ispc * 4 + (ibar as c_uint);
    zassert!(ic < 8, ZBAR_NONE, "ic={} ispc={} ibar={}", ic, ispc, ibar);
    let c = CODABAR_HI[ic as usize];
    dcode.codabar.buf[0] = (c & 0x3) | 0x10;

    // Set character direction
    dcode.codabar.set_direction((c >> 2) != 0);

    dcode.codabar.set_element(4);
    dcode.codabar.set_character(1);
    dcode.codabar.width = dcode.codabar.s7;
    ZBAR_PARTIAL
}

/// Calculate checksum
#[inline]
unsafe fn codabar_checksum(dcode: &zbar_decoder_t, n: usize) -> bool {
    let mut chk: c_uint = 0;
    let buf = dcode.buf;
    for i in 0..n {
        chk += *buf.add(i) as c_uint;
    }
    (chk & 0xf) != 0
}

/// Post-process decoded buffer
#[inline]
unsafe fn codabar_postprocess(dcode: &mut zbar_decoder_t) -> zbar_symbol_type_t {
    let dir = dcode.codabar.direction();
    dcode.direction = 1 - 2 * (dir as c_int);
    let mut n = dcode.codabar.character() as usize;

    // Copy from holding buffer to main buffer
    for i in 0..NIBUF {
        *dcode.buf.add(i) = dcode.codabar.buf[i] as c_char;
    }

    if dir {
        // Reverse buffer
        for i in 0..(n / 2) {
            let j = n - 1 - i;
            let code = *dcode.buf.add(i);
            *dcode.buf.add(i) = *dcode.buf.add(j);
            *dcode.buf.add(j) = code;
        }
    }

    if test_cfg(dcode.codabar.config, ZBAR_CFG_ADD_CHECK) {
        // Validate checksum
        if codabar_checksum(dcode, n) {
            return ZBAR_NONE;
        }
        if !test_cfg(dcode.codabar.config, ZBAR_CFG_EMIT_CHECK) {
            *dcode.buf.add(n - 2) = *dcode.buf.add(n - 1);
            n -= 1;
        }
    }

    for i in 0..n {
        let c = *dcode.buf.add(i) as u8;
        *dcode.buf.add(i) = if (c as usize) < 0x14 {
            CODABAR_CHARACTERS[c as usize] as c_char
        } else {
            b'?' as c_char
        };
    }
    dcode.buflen = n as c_uint;
    *dcode.buf.add(n) = 0;
    dcode.modifiers = 0;

    dcode.codabar.set_character(-1);
    ZBAR_CODABAR
}

/// Main Codabar decode function
#[no_mangle]
pub unsafe extern "C" fn _zbar_decode_codabar(dcode: *mut zbar_decoder_t) -> zbar_symbol_type_t {
    let dcode = &mut *dcode;

    // Update latest character width
    let w8 = get_width(dcode, 8);
    let w1 = get_width(dcode, 1);
    dcode.codabar.s7 = dcode.codabar.s7.wrapping_sub(w8).wrapping_add(w1);

    if get_color(dcode) != ZBAR_SPACE {
        return ZBAR_NONE;
    }

    if dcode.codabar.character() < 0 {
        return codabar_decode_start(dcode);
    }

    if dcode.codabar.character() < 2 && codabar_decode_start(dcode) != ZBAR_NONE {
        return ZBAR_PARTIAL;
    }

    let element = dcode.codabar.element() - 1;
    dcode.codabar.set_element(element);
    if element != 0 {
        return ZBAR_NONE;
    }
    dcode.codabar.set_element(4);

    let c = codabar_decode7(dcode);
    if c < 0 {
        // goto reset
        let character = dcode.codabar.character();
        if character >= NIBUF as i16 {
            release_lock(dcode, ZBAR_CODABAR);
        }
        dcode.codabar.set_character(-1);
        return ZBAR_NONE;
    }

    let character = dcode.codabar.character();
    if character < NIBUF as i16 {
        dcode.codabar.buf[character as usize] = c as u8;
    } else {
        if character >= BUFFER_MIN as i16 && size_buf(dcode, (character + 1) as c_uint) != 0 {
            // goto reset
            release_lock(dcode, ZBAR_CODABAR);
            dcode.codabar.set_character(-1);
            return ZBAR_NONE;
        }
        *dcode.buf.offset(character as isize) = c as c_char;
    }
    dcode.codabar.set_character(character + 1);

    // Lock shared resources
    if character + 1 == NIBUF as i16 && acquire_lock(dcode, ZBAR_CODABAR) {
        dcode.codabar.set_character(-1);
        return ZBAR_PARTIAL;
    }

    let s = dcode.codabar.s7;
    if (c & 0x10) != 0 {
        let qz = get_width(dcode, 0);
        if qz != 0 && qz * 2 < s {
            // goto reset
            let character = dcode.codabar.character();
            if character >= NIBUF as i16 {
                release_lock(dcode, ZBAR_CODABAR);
            }
            dcode.codabar.set_character(-1);
            return ZBAR_NONE;
        }

        let n = dcode.codabar.character();
        let min_len = cfg(&dcode.codabar, ZBAR_CFG_MIN_LEN);
        let max_len = cfg(&dcode.codabar, ZBAR_CFG_MAX_LEN);
        if n < min_len as i16 || (max_len > 0 && n > max_len as i16) {
            // goto reset
            let character = dcode.codabar.character();
            if character >= NIBUF as i16 {
                release_lock(dcode, ZBAR_CODABAR);
            }
            dcode.codabar.set_character(-1);
            return ZBAR_NONE;
        }

        if dcode.codabar.character() < NIBUF as i16 && acquire_lock(dcode, ZBAR_CODABAR) {
            dcode.codabar.set_character(-1);
            return ZBAR_PARTIAL;
        }

        let sym = codabar_postprocess(dcode);
        if sym <= ZBAR_PARTIAL {
            release_lock(dcode, ZBAR_CODABAR);
            dcode.codabar.set_character(-1);
        }
        return sym;
    } else if 4 * get_width(dcode, 0) > 3 * s {
        // goto reset
        let character = dcode.codabar.character();
        if character >= NIBUF as i16 {
            release_lock(dcode, ZBAR_CODABAR);
        }
        dcode.codabar.set_character(-1);
        return ZBAR_NONE;
    }

    ZBAR_NONE
}
