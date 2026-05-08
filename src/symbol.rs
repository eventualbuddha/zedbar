//! Decoded barcode symbols and types
//!
//! This module provides types for working with decoded barcode data:
//! - [`Symbol`] - A decoded barcode symbol with its data and metadata
//! - [`SymbolType`] - The barcode format (QR Code, EAN-13, etc.)
//! - [`Orientation`] - The orientation of the barcode in the image
//!
//! # Example
//!
//! ```no_run
//! use zedbar::{Image, Scanner};
//!
//! # let data = vec![0u8; 640 * 480];
//! # let mut image = Image::from_gray(&data, 640, 480).unwrap();
//! let mut scanner = Scanner::new();
//! let symbols = scanner.scan(&mut image);
//!
//! for symbol in symbols {
//!     println!("Found {:?}", symbol.symbol_type());
//!     if let Some(text) = symbol.data_string() {
//!         println!("Data: {}", text);
//!     }
//!     println!("Orientation: {:?}", symbol.orientation());
//! }
//! ```

use std::{fmt::Display, str::from_utf8};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Hash)]
pub enum Orientation {
    #[default]
    Unknown,
    Up,
    Right,
    Down,
    Left,
}

/// A 2D pixel coordinate in image space.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Point {
    pub x: i32,
    pub y: i32,
}

impl Point {
    pub const fn new(x: i32, y: i32) -> Self {
        Self { x, y }
    }
}

impl From<[i32; 2]> for Point {
    fn from([x, y]: [i32; 2]) -> Self {
        Self { x, y }
    }
}

impl From<Point> for [i32; 2] {
    fn from(p: Point) -> Self {
        [p.x, p.y]
    }
}

#[derive(Default, Clone)]
pub struct Symbol {
    symbol_type: SymbolType,
    pub(crate) modifiers: u32,
    pub(crate) data: Vec<u8>,
    pub(crate) raw_data: Option<Vec<u8>>,
    pub(crate) pts: Vec<Point>,
    pub(crate) orientation: Orientation,
    pub(crate) components: Vec<Symbol>,
    pub(crate) quality: i32,
}

// Internal API
impl Symbol {
    /// Create a new symbol
    pub(crate) fn new(symbol_type: SymbolType) -> Self {
        Self {
            symbol_type,
            quality: 1,
            ..Default::default()
        }
    }

    /// Create a composite symbol from two symbols.
    pub(crate) fn composite(self, other: Self) -> Self {
        Self {
            symbol_type: SymbolType::Composite,
            orientation: self.orientation,
            data: self.data.iter().chain(&other.data).cloned().collect(),
            components: vec![self, other],
            ..Default::default()
        }
    }

    pub(crate) fn add_point(&mut self, x: i32, y: i32) {
        self.pts.push(Point::new(x, y));
    }
}

// Public API
impl Symbol {
    /// Get the symbol type of this symbol.
    pub fn symbol_type(&self) -> SymbolType {
        self.symbol_type
    }

    /// Get the decoded data as bytes.
    ///
    /// For QR codes, this is the text-decoded output after encoding detection
    /// (UTF-8, Shift-JIS, Windows-1252, etc.). For SQ codes, this is the
    /// base64-encoded payload. Use [`raw_data`](Self::raw_data) to get the
    /// original bytes before conversion.
    ///
    /// For linear barcodes, this is the raw data (which is always ASCII text).
    pub fn data(&self) -> &[u8] {
        &self.data
    }

    /// Get the raw bytes before encoding conversion, if available.
    ///
    /// For QR codes, this returns the original bytes as encoded in the barcode,
    /// before any text encoding detection or conversion. In mixed-mode QR codes,
    /// numeric and alphanumeric segments appear as their ASCII representation
    /// (since those modes encode text, not arbitrary bytes), while byte-mode
    /// segments are the original uninterpreted bytes.
    ///
    /// For SQ codes, this returns the raw bytes before base64 encoding.
    ///
    /// Returns `None` for linear barcodes (use [`data`](Self::data) instead,
    /// which is already raw ASCII).
    pub fn raw_data(&self) -> Option<&[u8]> {
        self.raw_data.as_deref()
    }

    /// Get the decoded data as a string, if decodable.
    pub fn data_string(&self) -> Option<&str> {
        from_utf8(&self.data).ok()
    }

    /// Get the orientation of the barcode
    pub fn orientation(&self) -> Orientation {
        self.orientation
    }

    /// Get the component symbols (for composite symbols like EAN+add-on or QR structured append)
    pub fn components(&self) -> &[Symbol] {
        &self.components
    }

    /// Get the recorded position points of this symbol in image coordinates.
    ///
    /// For QR codes, this is the four corner points of the QR's bounding
    /// rectangle (in an implementation-defined order). For linear barcodes,
    /// this is one or more touchpoints accumulated as the symbol was
    /// scanned, suitable for approximating the symbol's extent.
    ///
    /// Returns an empty slice when no points were recorded — for example,
    /// when position tracking is disabled for linear barcodes, or when the
    /// symbol came from a code path that does not populate points.
    pub fn points(&self) -> &[Point] {
        &self.pts
    }

    /// Get the axis-aligned bounding rectangle of this symbol in image
    /// coordinates, computed from [`points`](Self::points).
    ///
    /// `width` and `height` are reported as `max - min` of the recorded
    /// points (i.e. the horizontal and vertical extent between the
    /// outermost points), not as a pixel count.
    ///
    /// Returns `None` when no points were recorded.
    pub fn bounds(&self) -> Option<Bounds> {
        let mut iter = self.pts.iter().copied();
        let Point { x: x0, y: y0 } = iter.next()?;
        let (mut min_x, mut max_x) = (x0, x0);
        let (mut min_y, mut max_y) = (y0, y0);
        for Point { x, y } in iter {
            if x < min_x {
                min_x = x;
            } else if x > max_x {
                max_x = x;
            }
            if y < min_y {
                min_y = y;
            } else if y > max_y {
                max_y = y;
            }
        }
        let width = (max_x as i64 - min_x as i64).max(0).min(u32::MAX as i64) as u32;
        let height = (max_y as i64 - min_y as i64).max(0).min(u32::MAX as i64) as u32;
        Some(Bounds {
            x: min_x,
            y: min_y,
            width,
            height,
        })
    }
}

/// Axis-aligned bounding rectangle of a [`Symbol`] in image coordinates.
///
/// `x` and `y` are the top-left corner. `width` and `height` are the
/// horizontal and vertical extent of the symbol's recorded points, computed
/// as `max - min` (so a single recorded point yields `width == height == 0`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Bounds {
    pub x: i32,
    pub y: i32,
    pub width: u32,
    pub height: u32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, PartialOrd, Ord, Hash)]
pub enum SymbolType {
    #[default]
    None = 0,
    Partial = 1,
    Ean2 = 2,
    Ean5 = 5,
    Ean8 = 8,
    Upce = 9,
    Isbn10 = 10,
    Upca = 12,
    Ean13 = 13,
    Isbn13 = 14,
    Composite = 15,
    I25 = 25,
    Databar = 34,
    DatabarExp = 35,
    Codabar = 38,
    Code39 = 39,
    QrCode = 64,
    SqCode = 80,
    Code93 = 93,
    Code128 = 128,
}

impl SymbolType {
    pub(crate) const ALL: [Self; 17] = [
        Self::Ean13,
        Self::Ean2,
        Self::Ean5,
        Self::Ean8,
        Self::Upca,
        Self::Upce,
        Self::Isbn10,
        Self::Isbn13,
        Self::I25,
        Self::Databar,
        Self::DatabarExp,
        Self::Codabar,
        Self::Code39,
        Self::Code93,
        Self::Code128,
        Self::QrCode,
        Self::SqCode,
    ];
}

impl Display for SymbolType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::None => "None",
                Self::Partial => "Partial",
                Self::Ean2 => "EAN-2",
                Self::Ean5 => "EAN-5",
                Self::Ean8 => "EAN-8",
                Self::Upce => "UPC-E",
                Self::Isbn10 => "ISBN-10",
                Self::Upca => "UPC-A",
                Self::Ean13 => "EAN-13",
                Self::Isbn13 => "ISBN-13",
                Self::Composite => "COMPOSITE",
                Self::I25 => "I2/5",
                Self::Databar => "DataBar",
                Self::DatabarExp => "DataBar-Exp",
                Self::Codabar => "Codabar",
                Self::Code39 => "CODE-39",
                Self::Code93 => "CODE-93",
                Self::Code128 => "CODE-128",
                Self::QrCode => "QR-Code",
                Self::SqCode => "SQ-Code",
            }
        )
    }
}

impl From<SymbolType> for i32 {
    fn from(value: SymbolType) -> Self {
        value as i32
    }
}

impl From<i32> for SymbolType {
    fn from(value: i32) -> Self {
        match value {
            0 => Self::None,
            1 => Self::Partial,
            2 => Self::Ean2,
            5 => Self::Ean5,
            8 => Self::Ean8,
            9 => Self::Upce,
            10 => Self::Isbn10,
            12 => Self::Upca,
            13 => Self::Ean13,
            14 => Self::Isbn13,
            15 => Self::Composite,
            25 => Self::I25,
            34 => Self::Databar,
            35 => Self::DatabarExp,
            38 => Self::Codabar,
            39 => Self::Code39,
            64 => Self::QrCode,
            80 => Self::SqCode,
            93 => Self::Code93,
            128 => Self::Code128,
            _ => Self::None,
        }
    }
}

/// Returned when a string does not match any [`SymbolType`] display name.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseSymbolTypeError(pub String);

impl Display for ParseSymbolTypeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "unknown symbology: {:?}", self.0)
    }
}

impl std::error::Error for ParseSymbolTypeError {}

impl std::str::FromStr for SymbolType {
    type Err = ParseSymbolTypeError;

    /// Parse a symbology by its [`Display`] name (e.g. `"QR-Code"`,
    /// `"EAN-13"`, `"CODE-128"`). Recognises every variant in
    /// [`SymbolType::ALL`].
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::ALL
            .iter()
            .copied()
            .find(|sym| sym.to_string() == s)
            .ok_or_else(|| ParseSymbolTypeError(s.to_string()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_points_empty_when_no_points_recorded() {
        let sym = Symbol::new(SymbolType::QrCode);
        assert!(sym.points().is_empty());
    }

    #[test]
    fn test_bounds_returns_none_when_no_points_recorded() {
        let sym = Symbol::new(SymbolType::QrCode);
        assert_eq!(sym.bounds(), None);
    }

    #[test]
    fn test_bounds_with_single_point_has_zero_extent() {
        let mut sym = Symbol::new(SymbolType::Code128);
        sym.add_point(7, 11);
        assert_eq!(
            sym.bounds(),
            Some(Bounds {
                x: 7,
                y: 11,
                width: 0,
                height: 0,
            })
        );
    }

    #[test]
    fn test_bounds_aabb_of_qr_corners() {
        let mut sym = Symbol::new(SymbolType::QrCode);
        // Corners pushed in qrdec's [0, 2, 3, 1] order; bounds() must be
        // independent of insertion order.
        sym.add_point(88, 1996);
        sym.add_point(211, 2121);
        sym.add_point(88, 2121);
        sym.add_point(211, 1996);
        assert_eq!(sym.points().len(), 4);
        assert_eq!(
            sym.bounds(),
            Some(Bounds {
                x: 88,
                y: 1996,
                width: 123,
                height: 125,
            })
        );
    }

    #[test]
    fn test_bounds_with_negative_coordinates() {
        let mut sym = Symbol::new(SymbolType::Code128);
        sym.add_point(-3, -7);
        sym.add_point(2, 4);
        assert_eq!(
            sym.bounds(),
            Some(Bounds {
                x: -3,
                y: -7,
                width: 5,
                height: 11,
            })
        );
    }
}
