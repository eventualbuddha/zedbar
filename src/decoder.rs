//! Low-level barcode decoder

use libc::{c_char, c_int, c_uint};
use crate::decoder_types::{zbar_decoder_t, zbar_symbol_type_t, DECODE_WINDOW};

// Buffer allocation constants
const BUFFER_MIN: c_uint = 0x20;
const BUFFER_MAX: c_uint = 0x100;
const BUFFER_INCR: c_uint = 0x10;

/// Get the color (bar/space) of the current element
#[no_mangle]
pub unsafe extern "C" fn _zbar_decoder_get_color(dcode: *const zbar_decoder_t) -> c_char {
    ((*dcode).idx & 1) as c_char
}

/// Get width of a specific element from the decoder's history window
#[no_mangle]
pub unsafe extern "C" fn _zbar_decoder_get_width(
    dcode: *const zbar_decoder_t,
    offset: u8,
) -> c_uint {
    (*dcode).w[(((*dcode).idx as usize).wrapping_sub(offset as usize)) & (DECODE_WINDOW - 1)]
}

/// Get the combined width of two consecutive elements
#[no_mangle]
pub unsafe extern "C" fn _zbar_decoder_pair_width(
    dcode: *const zbar_decoder_t,
    offset: u8,
) -> c_uint {
    _zbar_decoder_get_width(dcode, offset) + _zbar_decoder_get_width(dcode, offset + 1)
}

/// Calculate sum of n consecutive element widths
#[no_mangle]
pub unsafe extern "C" fn _zbar_decoder_calc_s(
    dcode: *const zbar_decoder_t,
    mut offset: u8,
    mut n: u8,
) -> c_uint {
    let mut s = 0;
    while n > 0 {
        s += _zbar_decoder_get_width(dcode, offset);
        offset += 1;
        n -= 1;
    }
    s
}

/// Decode element width into a discrete value
/// Returns -1 if the element width is invalid
#[no_mangle]
pub unsafe extern "C" fn _zbar_decoder_decode_e(e: c_uint, s: c_uint, n: c_uint) -> c_int {
    let big_e = ((e * n * 2 + 1) / s).wrapping_sub(3) / 2;
    if big_e >= n - 3 {
        -1
    } else {
        big_e as c_int
    }
}

/// Sort 3 elements by width and return their indices
#[no_mangle]
pub unsafe extern "C" fn _zbar_decoder_decode_sort3(
    dcode: *mut zbar_decoder_t,
    i0: c_int,
) -> c_uint {
    let w0 = _zbar_decoder_get_width(dcode, i0 as u8);
    let w2 = _zbar_decoder_get_width(dcode, (i0 + 2) as u8);
    let w4 = _zbar_decoder_get_width(dcode, (i0 + 4) as u8);

    if w0 < w2 {
        if w2 < w4 {
            ((i0 << 8) | ((i0 + 2) << 4) | (i0 + 4)) as c_uint
        } else if w0 < w4 {
            ((i0 << 8) | ((i0 + 4) << 4) | (i0 + 2)) as c_uint
        } else {
            (((i0 + 4) << 8) | (i0 << 4) | (i0 + 2)) as c_uint
        }
    } else {
        if w4 < w2 {
            (((i0 + 4) << 8) | ((i0 + 2) << 4) | i0) as c_uint
        } else if w0 < w4 {
            (((i0 + 2) << 8) | (i0 << 4) | (i0 + 4)) as c_uint
        } else {
            (((i0 + 2) << 8) | ((i0 + 4) << 4) | i0) as c_uint
        }
    }
}

/// Sort n elements by width and return their indices packed into a value
#[no_mangle]
pub unsafe extern "C" fn _zbar_decoder_decode_sortn(
    dcode: *mut zbar_decoder_t,
    n: c_int,
    i0: c_int,
) -> c_uint {
    let mut mask: c_uint = 0;
    let mut sort: c_uint = 0;

    for _ in (0..n).rev() {
        let mut wmin = c_uint::MAX;
        let mut jmin: c_int = -1;

        for j in (0..n).rev() {
            if (mask >> j) & 1 != 0 {
                continue;
            }
            let w = _zbar_decoder_get_width(dcode, (i0 + j * 2) as u8);
            if wmin >= w {
                wmin = w;
                jmin = j;
            }
        }

        debug_assert!(jmin >= 0, "sortn({}, {}) jmin={}", n, i0, jmin);
        sort <<= 4;
        mask |= 1 << jmin;
        sort |= (i0 + jmin * 2) as c_uint;
    }

    sort
}

/// Acquire a decoder lock for a specific symbology type
/// Returns 1 if already locked, 0 if lock acquired
#[no_mangle]
pub unsafe extern "C" fn _zbar_decoder_acquire_lock(
    dcode: *mut zbar_decoder_t,
    req: zbar_symbol_type_t,
) -> c_char {
    if (*dcode).lock != 0 {
        return 1;
    }
    (*dcode).lock = req;
    0
}

/// Release a decoder lock
/// Returns 0 on success
#[no_mangle]
pub unsafe extern "C" fn _zbar_decoder_release_lock(
    dcode: *mut zbar_decoder_t,
    req: zbar_symbol_type_t,
) -> c_char {
    debug_assert_eq!(
        (*dcode).lock, req,
        "lock={} req={}",
        (*dcode).lock, req
    );
    (*dcode).lock = 0;
    0
}

/// Resize the decoder's data buffer if needed
/// Returns 1 on allocation failure or if max size exceeded, 0 on success
#[no_mangle]
pub unsafe extern "C" fn _zbar_decoder_size_buf(
    dcode: *mut zbar_decoder_t,
    len: c_uint,
) -> c_char {
    if len <= BUFFER_MIN {
        return 0;
    }
    if len < (*dcode).buf_alloc {
        // FIXME: size reduction heuristic?
        return 0;
    }
    if len > BUFFER_MAX {
        return 1;
    }

    let mut new_len = len;
    if len < (*dcode).buf_alloc + BUFFER_INCR {
        new_len = (*dcode).buf_alloc + BUFFER_INCR;
        if new_len > BUFFER_MAX {
            new_len = BUFFER_MAX;
        }
    }

    let buf = libc::realloc((*dcode).buf as *mut libc::c_void, new_len as usize) as *mut c_char;
    if buf.is_null() {
        return 1;
    }

    (*dcode).buf = buf;
    (*dcode).buf_alloc = new_len;
    0
}

/// Low-level decoder for processing bar/space width streams
pub struct Decoder {
    _ptr: *mut std::ffi::c_void,
}

impl Decoder {
    /// Create a new decoder
    pub fn new() -> Self {
        // Note: This is a placeholder - the C library doesn't expose the decoder directly
        // In a full conversion, we'd implement the decoder logic in Rust
        todo!("Decoder needs to be implemented as part of the Rust conversion")
    }
}

impl Default for Decoder {
    fn default() -> Self {
        Self::new()
    }
}