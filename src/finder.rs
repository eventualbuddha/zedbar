//! Barcode finder pattern detection
//!
//! This module implements finder pattern detection for QR codes and SQ codes.

use crate::decoder_types::{
    qr_finder_line, qr_finder_t, zbar_decoder_t, zbar_symbol_type_t, DECODE_WINDOW, ZBAR_QRCODE,
    ZBAR_SPACE,
};
use libc::c_uint;

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

/// Retrieve bar+space pair width starting at offset i
#[inline]
fn pair_width(dcode: &zbar_decoder_t, offset: u8) -> c_uint {
    get_width(dcode, offset) + get_width(dcode, offset + 1)
}

/// Fixed character width decode assist
///
/// Bar+space width are compared as a fraction of the reference dimension "x"
/// - +/- 1/2 x tolerance
/// - measured total character width (s) compared to symbology baseline (n)
/// - bar+space *pair width* "e" is used to factor out bad "exposures"
///
/// Returns encoded number of units - 2 (for use as zero based index)
/// or -1 if invalid
#[inline]
fn decode_e(e: c_uint, s: c_uint, n: c_uint) -> i32 {
    let e_val = ((e * n * 2 + 1) / s).wrapping_sub(3) / 2;
    if e_val >= n - 3 {
        -1
    } else {
        e_val as i32
    }
}

// ============================================================================
// SQ Finder functions
// ============================================================================

#[no_mangle]
pub unsafe extern "C" fn _zbar_decoder_get_sq_finder_config(dcode: *mut zbar_decoder_t) -> c_uint {
    (*dcode).sqf.config
}

// ============================================================================
// QR Finder functions
// ============================================================================

/// Get pointer to QR finder line state
///
/// At this point lengths are all decode unit offsets from the decode edge.
/// Note: owned by finder
#[no_mangle]
pub unsafe extern "C" fn _zbar_decoder_get_qr_finder_line(
    dcode: *mut zbar_decoder_t,
) -> *mut qr_finder_line {
    &mut (*dcode).qrf.line
}

/// Find QR code finder pattern
///
/// Searches for the 1:1:3:1:1 ratio pattern characteristic of QR code finders.
#[no_mangle]
pub unsafe extern "C" fn _zbar_find_qr(dcode: *mut zbar_decoder_t) -> zbar_symbol_type_t {
    let dcode_ref = &mut *dcode;
    let qrf: *mut qr_finder_t = &mut dcode_ref.qrf;

    // Update latest finder pattern width
    (*qrf).s5 -= get_width(dcode_ref, 6);
    (*qrf).s5 += get_width(dcode_ref, 1);
    let s = (*qrf).s5;

    // TODO: The 2005 standard allows reflectance-reversed codes (light on dark
    // instead of dark on light).
    // If we find finder patterns with the opposite polarity, we should invert
    // the final binarized image and use them to search for QR codes in that.
    if get_color(dcode_ref) != ZBAR_SPACE || s < 7 {
        return 0;
    }

    // Check for 1:1:3:1:1 ratio pattern
    let mut ei = decode_e(pair_width(dcode_ref, 1), s, 7);
    if ei != 0 {
        return 0;
    }

    ei = decode_e(pair_width(dcode_ref, 2), s, 7);
    if ei != 2 {
        return 0;
    }

    ei = decode_e(pair_width(dcode_ref, 3), s, 7);
    if ei != 2 {
        return 0;
    }

    ei = decode_e(pair_width(dcode_ref, 4), s, 7);
    if ei != 0 {
        return 0;
    }

    // Valid QR finder symbol - mark positions needed by decoder
    let qz = get_width(dcode_ref, 0);
    let mut w = get_width(dcode_ref, 1);
    (*qrf).line.eoffs = (qz + w.div_ceil(2)) as i32;
    (*qrf).line.len = (qz + w + get_width(dcode_ref, 2)) as i32;
    (*qrf).line.pos[0] = ((*qrf).line.len + get_width(dcode_ref, 3) as i32) as i32;
    (*qrf).line.pos[1] = (*qrf).line.pos[0];
    w = get_width(dcode_ref, 5);
    (*qrf).line.boffs =
        ((*qrf).line.pos[0] + get_width(dcode_ref, 4) as i32 + w.div_ceil(2) as i32) as i32;

    dcode_ref.direction = 0;
    dcode_ref.buflen = 0;

    ZBAR_QRCODE
}
