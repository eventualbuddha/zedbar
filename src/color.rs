//! Low-level barcode line scanner
//!
//! This module provides color definitions used during 1D barcode decoding.

/// Color of element: bar or space
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Color {
    Space = 0, // light area or space between bars
    Bar = 1,   // dark area or colored bar segment
}

impl From<u8> for Color {
    fn from(value: u8) -> Self {
        if value & 1 == 1 {
            Self::Bar
        } else {
            Self::Space
        }
    }
}
