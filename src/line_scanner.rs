//! Low-level barcode line scanner
//!
//! This module provides color definitions used during 1D barcode decoding.

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
