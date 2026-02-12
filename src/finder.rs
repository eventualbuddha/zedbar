//! Barcode finder pattern detection
//!
//! This module implements finder pattern detection for QR codes and SQ codes.

use crate::{color::Color, decoder::decode_e, img_scanner::ImageScanner, SymbolType};

// ============================================================================
// QR Finder functions
// ============================================================================

/// Find QR code finder pattern
///
/// Searches for the 1:1:3:1:1 ratio pattern characteristic of QR code finders.
pub(crate) fn find_qr(dcode: &mut ImageScanner) -> SymbolType {
    // Update latest finder pattern width
    dcode.qrf.s5 -= dcode.get_width(6);
    dcode.qrf.s5 += dcode.get_width(1);
    let s = dcode.qrf.s5;

    // TODO: The 2005 standard allows reflectance-reversed codes (light on dark
    // instead of dark on light).
    // If we find finder patterns with the opposite polarity, we should invert
    // the final binarized image and use them to search for QR codes in that.
    if dcode.color() != Color::Space || s < 7 {
        return SymbolType::None;
    }

    // Check for 1:1:3:1:1 ratio pattern
    let mut ei = decode_e(dcode.pair_width(1), s, 7);
    if ei != 0 {
        return SymbolType::None;
    }

    ei = decode_e(dcode.pair_width(2), s, 7);
    if ei != 2 {
        return SymbolType::None;
    }

    ei = decode_e(dcode.pair_width(3), s, 7);
    if ei != 2 {
        return SymbolType::None;
    }

    ei = decode_e(dcode.pair_width(4), s, 7);
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
