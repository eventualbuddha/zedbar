//! Low-level barcode line scanner
//!
//! This module provides the line scanner functionality that processes 1D scan lines
//! to detect bar/space transitions and widths.

use libc::{c_int, c_uint, c_void, free, malloc, memset};
use std::ptr;

// Constants from scanner.c
const ZBAR_FIXED: i32 = 5;
const ROUND: c_uint = 1 << (ZBAR_FIXED - 1); // 16
const ZBAR_SCANNER_THRESH_MIN: c_uint = 4;
const ZBAR_SCANNER_THRESH_FADE: c_uint = 8;

// Calculated constants
// THRESH_INIT = (0.44 * (1 << (ZBAR_FIXED + 1)) + 1) / 2 = (0.44 * 64 + 1) / 2 = 14
const THRESH_INIT: c_uint = 14;
// EWMA_WEIGHT = (0.78 * (1 << (ZBAR_FIXED + 1)) + 1) / 2 = (0.78 * 64 + 1) / 2 = 25
const EWMA_WEIGHT: c_uint = 25;

/// Color of element: bar or space
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(non_camel_case_types)]
pub enum zbar_color_t {
    ZBAR_SPACE = 0, // light area or space between bars
    ZBAR_BAR = 1,   // dark area or colored bar segment
}

// Forward declaration for decoder type (opaque pointer from C)
#[repr(C)]
#[allow(non_camel_case_types)]
pub struct zbar_decoder_t {
    _private: [u8; 0],
}

/// Scanner state structure - must match the C layout exactly
#[repr(C)]
#[allow(non_camel_case_types)]
pub struct zbar_scanner_t {
    decoder: *mut zbar_decoder_t,
    y1_min_thresh: c_uint,
    x: c_uint,
    y0: [i32; 4],
    y1_sign: i32,
    y1_thresh: c_uint,
    cur_edge: c_uint,
    last_edge: c_uint,
    width: c_uint,
}

/// Get the width of the most recent bar or space
///
/// This function returns the width of the last decoded element (bar or space)
/// in the barcode scan line.
#[no_mangle]
pub unsafe extern "C" fn zbar_scanner_get_width(scn: *const zbar_scanner_t) -> c_uint {
    if scn.is_null() {
        return 0;
    }

    (*scn).width
}

/// Get the interpolated position of the last edge
///
/// Returns the interpolated position of the last processed edge, adjusted by
/// the specified offset and precision.
#[no_mangle]
pub unsafe extern "C" fn zbar_scanner_get_edge(
    scn: *const zbar_scanner_t,
    offset: c_uint,
    prec: c_int,
) -> c_uint {
    if scn.is_null() {
        return 0;
    }

    let edge = (*scn)
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

/// Get the current color of the scanner
///
/// Returns ZBAR_SPACE or ZBAR_BAR depending on whether the scanner is
/// currently processing a space (light area) or bar (dark area).
#[no_mangle]
pub unsafe extern "C" fn zbar_scanner_get_color(scn: *const zbar_scanner_t) -> zbar_color_t {
    if scn.is_null() {
        return zbar_color_t::ZBAR_SPACE;
    }

    if (*scn).y1_sign <= 0 {
        zbar_color_t::ZBAR_SPACE
    } else {
        zbar_color_t::ZBAR_BAR
    }
}

/// Destroy a scanner instance
///
/// Frees all resources associated with the scanner.
#[no_mangle]
pub extern "C" fn zbar_scanner_destroy(scn: *mut zbar_scanner_t) {
    if !scn.is_null() {
        unsafe {
            free(scn as *mut c_void);
        }
    }
}

// Symbol type constants (from zbar.h)
#[allow(non_camel_case_types)]
pub type zbar_symbol_type_t = c_int;
pub const ZBAR_NONE: zbar_symbol_type_t = 0;
pub const ZBAR_PARTIAL: zbar_symbol_type_t = 1;

// External C functions we need to call
extern "C" {
    fn zbar_decoder_reset(decoder: *mut zbar_decoder_t) -> zbar_symbol_type_t;
    fn zbar_decode_width(decoder: *mut zbar_decoder_t, width: c_uint) -> zbar_symbol_type_t;
    fn zbar_decoder_new_scan(decoder: *mut zbar_decoder_t);
}

/// Create a new scanner instance
///
/// Allocates and initializes a new scanner with the given decoder.
#[no_mangle]
pub unsafe extern "C" fn zbar_scanner_create(
    dcode: *mut zbar_decoder_t,
) -> *mut zbar_scanner_t {
    let scn = malloc(std::mem::size_of::<zbar_scanner_t>()) as *mut zbar_scanner_t;
    if scn.is_null() {
        return ptr::null_mut();
    }

    (*scn).decoder = dcode;
    (*scn).y1_min_thresh = ZBAR_SCANNER_THRESH_MIN;
    zbar_scanner_reset(scn);
    scn
}

/// Reset scanner state
///
/// Clears all scanner state and resets the threshold. Also resets the
/// associated decoder if present.
#[no_mangle]
pub unsafe extern "C" fn zbar_scanner_reset(scn: *mut zbar_scanner_t) -> zbar_symbol_type_t {
    if scn.is_null() {
        return ZBAR_NONE;
    }

    // memset everything from 'x' field to end of struct
    let offset = std::mem::offset_of!(zbar_scanner_t, x);
    let size = std::mem::size_of::<zbar_scanner_t>() - offset;
    let ptr = (scn as *mut u8).add(offset);
    memset(ptr as *mut c_void, 0, size);

    (*scn).y1_thresh = (*scn).y1_min_thresh;

    if !(*scn).decoder.is_null() {
        zbar_decoder_reset((*scn).decoder);
    }

    ZBAR_NONE
}

/// Calculate adaptive threshold
///
/// Returns the threshold for edge detection, which fades back to minimum
/// over distance from the last edge.
#[inline]
unsafe fn calc_thresh(scn: *mut zbar_scanner_t) -> c_uint {
    let thresh = (*scn).y1_thresh;

    if thresh <= (*scn).y1_min_thresh || (*scn).width == 0 {
        return (*scn).y1_min_thresh;
    }

    // Slowly return threshold to min
    let dx = ((*scn).x << ZBAR_FIXED) - (*scn).last_edge;
    let t = (thresh as u64 * dx as u64) / (*scn).width as u64;
    let t = t / ZBAR_SCANNER_THRESH_FADE as u64;

    if thresh as u64 > t {
        let new_thresh = (thresh as u64 - t) as c_uint;
        if new_thresh > (*scn).y1_min_thresh {
            return new_thresh;
        }
    }

    (*scn).y1_thresh = (*scn).y1_min_thresh;
    (*scn).y1_min_thresh
}

/// Process detected edge
///
/// Updates edge positions and passes the width to the decoder.
#[inline]
unsafe fn process_edge(scn: *mut zbar_scanner_t, _y1: i32) -> zbar_symbol_type_t {
    if (*scn).y1_sign == 0 {
        (*scn).last_edge = (1 << ZBAR_FIXED) + ROUND;
        (*scn).cur_edge = (*scn).last_edge;
    } else if (*scn).last_edge == 0 {
        (*scn).last_edge = (*scn).cur_edge;
    }

    (*scn).width = (*scn).cur_edge - (*scn).last_edge;
    (*scn).last_edge = (*scn).cur_edge;

    // Pass to decoder
    if !(*scn).decoder.is_null() {
        return zbar_decode_width((*scn).decoder, (*scn).width);
    }
    ZBAR_PARTIAL
}

/// Flush pending edge
///
/// Forces processing of any pending edge data.
#[no_mangle]
pub unsafe extern "C" fn zbar_scanner_flush(scn: *mut zbar_scanner_t) -> zbar_symbol_type_t {
    if scn.is_null() {
        return ZBAR_NONE;
    }

    if (*scn).y1_sign == 0 {
        return ZBAR_NONE;
    }

    let x = ((*scn).x << ZBAR_FIXED) + ROUND;

    if (*scn).cur_edge != x || (*scn).y1_sign > 0 {
        let edge = process_edge(scn, -(*scn).y1_sign);
        (*scn).cur_edge = x;
        (*scn).y1_sign = -(*scn).y1_sign;
        return edge;
    }

    (*scn).y1_sign = 0;
    (*scn).width = 0;

    if !(*scn).decoder.is_null() {
        return zbar_decode_width((*scn).decoder, 0);
    }
    ZBAR_PARTIAL
}

/// Start a new scan line
///
/// Flushes any pending data and resets the scanner for a new scan line.
#[no_mangle]
pub unsafe extern "C" fn zbar_scanner_new_scan(scn: *mut zbar_scanner_t) -> zbar_symbol_type_t {
    if scn.is_null() {
        return ZBAR_NONE;
    }

    let mut edge = ZBAR_NONE;
    while (*scn).y1_sign != 0 {
        let tmp = zbar_scanner_flush(scn);
        if tmp < 0 || tmp > edge {
            edge = tmp;
        }
    }

    // Reset scanner state
    let offset = std::mem::offset_of!(zbar_scanner_t, x);
    let size = std::mem::size_of::<zbar_scanner_t>() - offset;
    let ptr = (scn as *mut u8).add(offset);
    memset(ptr as *mut c_void, 0, size);

    (*scn).y1_thresh = (*scn).y1_min_thresh;

    if !(*scn).decoder.is_null() {
        zbar_decoder_new_scan((*scn).decoder);
    }

    edge
}

/// Process a single intensity sample
///
/// This is the main scanning function that processes each pixel's intensity
/// value to detect edges and measure bar/space widths.
#[no_mangle]
pub unsafe extern "C" fn zbar_scan_y(scn: *mut zbar_scanner_t, y: c_int) -> zbar_symbol_type_t {
    if scn.is_null() {
        return ZBAR_NONE;
    }

    let x = (*scn).x;
    let y0_1 = (*scn).y0[((x.wrapping_sub(1)) & 3) as usize];
    let mut y0_0 = y0_1;

    if x != 0 {
        // Update weighted moving average
        y0_0 += ((y - y0_1) * EWMA_WEIGHT as i32) >> ZBAR_FIXED;
        (*scn).y0[(x & 3) as usize] = y0_0;
    } else {
        y0_0 = y;
        (*scn).y0[0] = y;
        (*scn).y0[1] = y;
        (*scn).y0[2] = y;
        (*scn).y0[3] = y;
    }

    let y0_2 = (*scn).y0[((x.wrapping_sub(2)) & 3) as usize];
    let y0_3 = (*scn).y0[((x.wrapping_sub(3)) & 3) as usize];

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

    let mut edge = ZBAR_NONE;

    // 2nd zero-crossing is 1st local min/max - could be edge
    if (y2_1 == 0 || ((y2_1 > 0) != (y2_2 >= 0)))
        && (calc_thresh(scn) <= y1_1.abs() as c_uint)
    {
        // Check for 1st sign change
        let y1_rev = if (*scn).y1_sign > 0 {
            y1_1 < 0
        } else {
            y1_1 > 0
        };

        if y1_rev {
            // Intensity change reversal - finalize previous edge
            edge = process_edge(scn, y1_1);
        }

        if y1_rev || (*scn).y1_sign.abs() < y1_1.abs() {
            (*scn).y1_sign = y1_1;

            // Adaptive thresholding - start at multiple of new min/max
            (*scn).y1_thresh = (y1_1.abs() as c_uint * THRESH_INIT + ROUND) >> ZBAR_FIXED;
            if (*scn).y1_thresh < (*scn).y1_min_thresh {
                (*scn).y1_thresh = (*scn).y1_min_thresh;
            }

            // Update current edge
            let d = y2_1 - y2_2;
            (*scn).cur_edge = 1 << ZBAR_FIXED;
            if d == 0 {
                (*scn).cur_edge >>= 1;
            } else if y2_1 != 0 {
                // Interpolate zero crossing
                let adjustment = ((y2_1 << ZBAR_FIXED) + 1) / d;
                (*scn).cur_edge = ((*scn).cur_edge as i32 - adjustment) as c_uint;
            }
            (*scn).cur_edge += x << ZBAR_FIXED;
        }
    }

    (*scn).x = x + 1;
    edge
}
