//! Barcode finder pattern detection
//!
//! This module implements finder pattern detection for QR codes and SQ codes.

use crate::{
    decoder::{qr_finder_line, zbar_decoder_t},
    line_scanner::zbar_color_t,
    SymbolType,
};
use libc::c_uint;

// ============================================================================
// Helper functions from decoder.h
// ============================================================================

/// Retrieve bar+space pair width starting at offset i
#[inline]
fn pair_width(dcode: &zbar_decoder_t, offset: u8) -> c_uint {
    dcode.get_width(offset) + dcode.get_width(offset + 1)
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

pub(crate) fn decoder_get_sq_finder_config(dcode: &zbar_decoder_t) -> c_uint {
    dcode.sqf.config
}

// ============================================================================
// QR Finder functions
// ============================================================================

/// Get mutable reference to QR finder line state
///
/// At this point lengths are all decode unit offsets from the decode edge.
pub(crate) fn decoder_get_qr_finder_line(dcode: &mut zbar_decoder_t) -> &mut qr_finder_line {
    &mut dcode.qrf.line
}

/// Find QR code finder pattern
///
/// Searches for the 1:1:3:1:1 ratio pattern characteristic of QR code finders.
pub(crate) fn find_qr(dcode: &mut zbar_decoder_t) -> SymbolType {
    // Update latest finder pattern width
    dcode.qrf.s5 -= dcode.get_width(6);
    dcode.qrf.s5 += dcode.get_width(1);
    let s = dcode.qrf.s5;

    // TODO: The 2005 standard allows reflectance-reversed codes (light on dark
    // instead of dark on light).
    // If we find finder patterns with the opposite polarity, we should invert
    // the final binarized image and use them to search for QR codes in that.
    if dcode.color() != zbar_color_t::ZBAR_SPACE || s < 7 {
        return SymbolType::None;
    }

    // Check for 1:1:3:1:1 ratio pattern
    let mut ei = decode_e(pair_width(dcode, 1), s, 7);
    if ei != 0 {
        return SymbolType::None;
    }

    ei = decode_e(pair_width(dcode, 2), s, 7);
    if ei != 2 {
        return SymbolType::None;
    }

    ei = decode_e(pair_width(dcode, 3), s, 7);
    if ei != 2 {
        return SymbolType::None;
    }

    ei = decode_e(pair_width(dcode, 4), s, 7);
    if ei != 0 {
        return SymbolType::None;
    }

    // Valid QR finder symbol - mark positions needed by decoder
    // Calculate all values first before modifying the line structure
    let qz = dcode.get_width(0);
    let w1 = dcode.get_width(1);
    let w2 = dcode.get_width(2);
    let w3 = dcode.get_width(3);
    let w4 = dcode.get_width(4);
    let w5 = dcode.get_width(5);

    let eoffs = (qz + w1.div_ceil(2)) as i32;
    let len = (qz + w1 + w2) as i32;
    let pos = len + w3 as i32;
    let boffs = pos + w4 as i32 + w5.div_ceil(2) as i32;

    // Now update the line structure
    dcode.qrf.line.eoffs = eoffs;
    dcode.qrf.line.len = len;
    dcode.qrf.line.pos[0] = pos;
    dcode.qrf.line.pos[1] = pos;
    dcode.qrf.line.boffs = boffs;

    dcode.direction = 0;
    dcode.set_buffer_len(0);

    SymbolType::QrCode
}
