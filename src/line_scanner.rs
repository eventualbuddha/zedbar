//! Low-level barcode line scanner
//!
//! This module provides the line scanner functionality that processes 1D scan lines
//! to detect bar/space transitions and widths.

use libc::{c_int, c_uint, c_void, free};

// Constants from scanner.c
const ZBAR_FIXED: i32 = 5;
const ROUND: c_uint = 1 << (ZBAR_FIXED - 1); // 16

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
