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

/// Color of element: bar or space
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum zbar_color_t {
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
pub(crate) struct zbar_scanner_t {
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
    /// Create a new scanner instance (owned version)
    ///
    /// Initializes a new scanner with the specified decoder.
    pub(crate) fn new() -> Self {
        Self {
            y1_min_thresh: ZBAR_SCANNER_THRESH_MIN,
            ..Self::default()
        }
    }

    /// Get the width of the most recent bar or space
    pub(crate) fn width(&self) -> c_uint {
        self.width
    }

    /// Calculate the current threshold for edge detection
    ///
    /// Implements adaptive threshold calculation that slowly fades back to minimum.
    /// This helps with noise rejection while maintaining sensitivity.
    pub(crate) fn calc_thresh(&mut self) -> c_uint {
        // threshold 1st to improve noise rejection
        let thresh = self.y1_thresh;

        if thresh <= self.y1_min_thresh || self.width == 0 {
            return self.y1_min_thresh;
        }

        // slowly return threshold to min
        let dx = (self.x << ZBAR_FIXED) - self.last_edge;
        let mut t = thresh as u64 * dx as u64;
        t /= self.width as u64;
        t /= ZBAR_SCANNER_THRESH_FADE as u64;

        if thresh > t as c_uint {
            let new_thresh = thresh - t as c_uint;
            if new_thresh > self.y1_min_thresh {
                return new_thresh;
            }
        }

        self.y1_thresh = self.y1_min_thresh;
        self.y1_min_thresh
    }

    /// Flush the scanner state
    #[inline]
    pub(crate) fn scanner_flush(&mut self, decoder: &mut zbar_decoder_t) -> SymbolType {
        if self.y1_sign == 0 {
            return SymbolType::None;
        }

        let x = (self.x << ZBAR_FIXED) + ROUND;

        if self.cur_edge != x || self.y1_sign > 0 {
            let edge = self.process_edge(decoder);
            self.cur_edge = x;
            self.y1_sign = -self.y1_sign;
            return edge;
        }

        self.y1_sign = 0;
        self.width = 0;
        unsafe { decoder.decode_width(0) }
    }

    /// Start a new scan
    pub(crate) fn new_scan(&mut self, decoder: &mut zbar_decoder_t) -> SymbolType {
        let mut edge = SymbolType::None;

        while self.y1_sign != 0 {
            let tmp = self.scanner_flush(decoder);
            if tmp > edge {
                edge = tmp;
            }
        }

        // reset scanner and associated decoder
        self.x = 0;
        self.y0 = Default::default();
        self.y1_sign = 0;
        self.y1_thresh = self.y1_min_thresh;
        self.cur_edge = 0;
        self.last_edge = 0;
        self.width = 0;

        decoder.new_scan();
        edge
    }

    /// Process a single pixel intensity value
    pub(crate) fn scan_y(&mut self, y: c_int, decoder: &mut zbar_decoder_t) -> SymbolType {
        // retrieve short value history
        let x = self.x;
        let mut y0_1 = self.y0[((x.wrapping_sub(1)) & 3) as usize];
        let mut y0_0 = y0_1;

        if x != 0 {
            // update weighted moving average
            y0_0 += ((y - y0_1) * EWMA_WEIGHT as i32) >> ZBAR_FIXED;
            self.y0[(x & 3) as usize] = y0_0;
        } else {
            y0_0 = y;
            y0_1 = y;
            self.y0[0] = y;
            self.y0[1] = y;
            self.y0[2] = y;
            self.y0[3] = y;
        }

        let y0_2 = self.y0[((x.wrapping_sub(2)) & 3) as usize];
        let y0_3 = self.y0[((x.wrapping_sub(3)) & 3) as usize];

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
        if (y2_1 == 0 || ((y2_1 > 0) == (y2_2 < 0))) && (self.calc_thresh() <= y1_1.unsigned_abs())
        {
            // check for 1st sign change
            let y1_rev = if self.y1_sign > 0 { y1_1 < 0 } else { y1_1 > 0 };

            if y1_rev {
                // intensity change reversal - finalize previous edge
                edge = self.process_edge(decoder);
            }

            if y1_rev || (self.y1_sign.abs() < y1_1.abs()) {
                self.y1_sign = y1_1;

                // adaptive thresholding
                // start at multiple of new min/max
                self.y1_thresh =
                    ((y1_1.unsigned_abs() * THRESH_INIT + ROUND) >> ZBAR_FIXED) as c_uint;
                if self.y1_thresh < self.y1_min_thresh {
                    self.y1_thresh = self.y1_min_thresh;
                }

                // update current edge
                let d = y2_1 - y2_2;
                self.cur_edge = 1 << ZBAR_FIXED;
                if d == 0 {
                    self.cur_edge >>= 1;
                } else if y2_1 != 0 {
                    // interpolate zero crossing
                    self.cur_edge -= (((y2_1 << ZBAR_FIXED) + 1) / d) as c_uint;
                }
                self.cur_edge += x << ZBAR_FIXED;
            }
        }

        self.x = x + 1;
        edge
    }

    /// Flush scanner pipeline and start new scan
    ///
    /// This function flushes the scanner pipeline twice and then starts a new scan.
    /// It's typically called at quiet borders to reset the scanner state.
    pub(crate) fn quiet_border(&mut self, decoder: &mut zbar_decoder_t) {
        // Flush scanner pipeline twice
        self.scanner_flush(decoder);
        self.scanner_flush(decoder);

        // Start new scan
        self.new_scan(decoder);
    }

    /// Get the interpolated position of the last edge
    pub(crate) fn get_edge(&self, offset: c_uint, prec: c_int) -> c_uint {
        let edge = self
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

    /// Process an edge and pass the width to the decoder
    ///
    /// This function is called when an edge (transition) is detected.
    /// It calculates the width of the element and passes it to the decoder.
    fn process_edge(&mut self, decoder: &mut zbar_decoder_t) -> SymbolType {
        if self.y1_sign == 0 {
            self.last_edge = (1 << ZBAR_FIXED) + ROUND;
            self.cur_edge = (1 << ZBAR_FIXED) + ROUND;
        } else if self.last_edge == 0 {
            self.last_edge = self.cur_edge;
        }

        self.width = self.cur_edge - self.last_edge;
        self.last_edge = self.cur_edge;

        // pass to decoder
        let width = self.width;
        unsafe { decoder.decode_width(width) }
    }
}
