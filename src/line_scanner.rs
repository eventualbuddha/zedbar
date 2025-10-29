//! Low-level barcode line scanner
//!
//! This module provides the line scanner functionality that processes 1D scan lines
//! to detect bar/space transitions and widths.

use libc::{c_int, c_uint};

use crate::decoder::zbar_decoder_t;
use crate::SymbolType;

// Constants from scanner.c
const ZBAR_FIXED: i32 = 5;
const ROUND: c_uint = 1 << (ZBAR_FIXED - 1); // 16
const ZBAR_SCANNER_THRESH_FADE: u32 = 8;
const ZBAR_SCANNER_THRESH_MIN: c_uint = 4;

// EWMA_WEIGHT = (unsigned)((0.78 * (1 << 6) + 1) / 2) = 25
const EWMA_WEIGHT: c_uint = 25;

// THRESH_INIT = (unsigned)((0.44 * (1 << 6) + 1) / 2) = 14
const THRESH_INIT: c_uint = 14;

/// Calculate the current threshold for edge detection
///
/// Implements adaptive threshold calculation that slowly fades back to minimum.
/// This helps with noise rejection while maintaining sensitivity.
pub fn calc_thresh(scn: &mut zbar_scanner_t) -> c_uint {
    // threshold 1st to improve noise rejection
    let thresh = scn.y1_thresh;

    if thresh <= scn.y1_min_thresh || scn.width == 0 {
        return scn.y1_min_thresh;
    }

    // slowly return threshold to min
    let dx = (scn.x << ZBAR_FIXED) - scn.last_edge;
    let mut t = thresh as u64 * dx as u64;
    t /= scn.width as u64;
    t /= ZBAR_SCANNER_THRESH_FADE as u64;

    if thresh > t as c_uint {
        let new_thresh = thresh - t as c_uint;
        if new_thresh > scn.y1_min_thresh {
            return new_thresh;
        }
    }

    scn.y1_thresh = scn.y1_min_thresh;
    scn.y1_min_thresh
}

/// Color of element: bar or space
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum zbar_color_t {
    ZBAR_SPACE = 0, // light area or space between bars
    ZBAR_BAR = 1,   // dark area or colored bar segment
}

impl From<u8> for zbar_color_t {
    fn from(value: u8) -> Self {
        if value & 1 == 1 {
            Self::ZBAR_BAR
        } else {
            Self::ZBAR_SPACE
        }
    }
}

/// Scanner state structure - must match the C layout exactly
#[derive(Default)]
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

impl zbar_scanner_t {
    /// Get a mutable reference to the decoder (if present)
    #[inline]
    pub(crate) fn decoder_mut(&mut self) -> Option<&mut zbar_decoder_t> {
        unsafe { self.decoder.as_mut() }
    }

    /// Set the decoder pointer
    #[inline]
    pub(crate) fn set_decoder(&mut self, decoder: *mut zbar_decoder_t) {
        self.decoder = decoder;
    }
}

// ============================================================================
// Safe reference-based APIs
// ============================================================================

/// Get the width of the most recent bar or space
#[inline]
pub fn scanner_get_width(scn: &zbar_scanner_t) -> c_uint {
    scn.width
}

/// Get the interpolated position of the last edge
#[inline]
pub fn scanner_get_edge(scn: &zbar_scanner_t, offset: c_uint, prec: c_int) -> c_uint {
    let edge = scn
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

/// Create a new scanner instance (owned version)
///
/// Initializes a new scanner with the specified decoder.
pub(crate) fn zbar_scanner_new(dcode: *mut zbar_decoder_t) -> zbar_scanner_t {
    let mut scn = zbar_scanner_t::default();
    scn.set_decoder(dcode);
    scn.y1_min_thresh = ZBAR_SCANNER_THRESH_MIN;
    scn
}

/// Process an edge and pass the width to the decoder
///
/// This function is called when an edge (transition) is detected.
/// It calculates the width of the element and passes it to the decoder.
fn process_edge(scn: &mut zbar_scanner_t, _y1: i32) -> SymbolType {
    if scn.y1_sign == 0 {
        scn.last_edge = (1 << ZBAR_FIXED) + ROUND;
        scn.cur_edge = (1 << ZBAR_FIXED) + ROUND;
    } else if scn.last_edge == 0 {
        scn.last_edge = scn.cur_edge;
    }

    scn.width = scn.cur_edge - scn.last_edge;
    scn.last_edge = scn.cur_edge;

    // pass to decoder
    let width = scn.width;
    if let Some(decoder) = scn.decoder_mut() {
        unsafe { decoder.zbar_decode_width(width) }
    } else {
        SymbolType::Partial
    }
}

// ============================================================================
// Safe reference-based APIs (main implementations)
// ============================================================================

/// Flush the scanner state
#[inline]
pub fn scanner_flush(scn: &mut zbar_scanner_t) -> SymbolType {
    if scn.y1_sign == 0 {
        return SymbolType::None;
    }

    let x = (scn.x << ZBAR_FIXED) + ROUND;

    if scn.cur_edge != x || scn.y1_sign > 0 {
        let edge = process_edge(scn, -scn.y1_sign);
        scn.cur_edge = x;
        scn.y1_sign = -scn.y1_sign;
        return edge;
    }

    scn.y1_sign = 0;
    scn.width = 0;
    if let Some(decoder) = scn.decoder_mut() {
        unsafe { decoder.zbar_decode_width(0) }
    } else {
        SymbolType::Partial
    }
}

/// Start a new scan
pub fn scanner_new_scan(scn: &mut zbar_scanner_t) -> SymbolType {
    let mut edge = SymbolType::None;

    while scn.y1_sign != 0 {
        let tmp = scanner_flush(scn);
        if tmp > edge {
            edge = tmp;
        }
    }

    // reset scanner and associated decoder
    scn.x = 0;
    scn.y0 = Default::default();
    scn.y1_sign = 0;
    scn.y1_thresh = scn.y1_min_thresh;
    scn.cur_edge = 0;
    scn.last_edge = 0;
    scn.width = 0;

    if let Some(decoder) = scn.decoder_mut() {
        unsafe { decoder.new_scan() };
    }
    edge
}

/// Process a single pixel intensity value
pub fn scan_y(scn: &mut zbar_scanner_t, y: c_int) -> SymbolType {
    // retrieve short value history
    let x = scn.x;
    let mut y0_1 = scn.y0[((x.wrapping_sub(1)) & 3) as usize];
    let mut y0_0 = y0_1;

    if x != 0 {
        // update weighted moving average
        y0_0 += ((y - y0_1) * EWMA_WEIGHT as i32) >> ZBAR_FIXED;
        scn.y0[(x & 3) as usize] = y0_0;
    } else {
        y0_0 = y;
        y0_1 = y;
        scn.y0[0] = y;
        scn.y0[1] = y;
        scn.y0[2] = y;
        scn.y0[3] = y;
    }

    let y0_2 = scn.y0[((x.wrapping_sub(2)) & 3) as usize];
    let y0_3 = scn.y0[((x.wrapping_sub(3)) & 3) as usize];

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

    let mut edge = SymbolType::None;

    // 2nd zero-crossing is 1st local min/max - could be edge
    if (y2_1 == 0 || ((y2_1 > 0) == (y2_2 < 0))) && (calc_thresh(scn) <= y1_1.unsigned_abs()) {
        // check for 1st sign change
        let y1_rev = if scn.y1_sign > 0 { y1_1 < 0 } else { y1_1 > 0 };

        if y1_rev {
            // intensity change reversal - finalize previous edge
            edge = process_edge(scn, y1_1);
        }

        if y1_rev || (scn.y1_sign.abs() < y1_1.abs()) {
            scn.y1_sign = y1_1;

            // adaptive thresholding
            // start at multiple of new min/max
            scn.y1_thresh = ((y1_1.unsigned_abs() * THRESH_INIT + ROUND) >> ZBAR_FIXED) as c_uint;
            if scn.y1_thresh < scn.y1_min_thresh {
                scn.y1_thresh = scn.y1_min_thresh;
            }

            // update current edge
            let d = y2_1 - y2_2;
            scn.cur_edge = 1 << ZBAR_FIXED;
            if d == 0 {
                scn.cur_edge >>= 1;
            } else if y2_1 != 0 {
                // interpolate zero crossing
                scn.cur_edge -= (((y2_1 << ZBAR_FIXED) + 1) / d) as c_uint;
            }
            scn.cur_edge += x << ZBAR_FIXED;
        }
    }

    scn.x = x + 1;
    edge
}
