//! Symbol management module
//!
//! Handles symbol lifecycle, reference counting, and data access.
//!
//! Rust port based on C code from the ZBar library.
//! Original C code copyright (C) 2007-2010 Jeff Brown <spadix@users.sourceforge.net>
//! Licensed under LGPL 3.0 or later

use crate::{img_scanner::zbar_symbol_set_t, qrcode::qrdec::qr_point};
use std::fmt::Display;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Hash)]
pub enum Orientation {
    #[default]
    Unknown,
    Up,
    Right,
    Down,
    Left,
}

#[derive(Default, Clone)]
pub(crate) struct zbar_symbol_t {
    pub(crate) symbol_type: SymbolType,
    pub(crate) modifiers: u32,
    pub(crate) data: Vec<u8>,
    pub(crate) pts: Vec<qr_point>,
    pub(crate) orientation: Orientation,
    pub(crate) components: Option<zbar_symbol_set_t>,
    pub(crate) quality: i32,
}

impl zbar_symbol_t {
    /// Create a new symbol
    pub(crate) fn new(symbol_type: SymbolType) -> zbar_symbol_t {
        zbar_symbol_t {
            symbol_type,
            quality: 1,
            ..Default::default()
        }
    }

    pub(crate) fn add_point(&mut self, x: i32, y: i32) {
        self.pts.push([x, y]);
    }
}

// High-level Rust API types

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
        SymbolType::Ean13,
        SymbolType::Ean2,
        SymbolType::Ean5,
        SymbolType::Ean8,
        SymbolType::Upca,
        SymbolType::Upce,
        SymbolType::Isbn10,
        SymbolType::Isbn13,
        SymbolType::I25,
        SymbolType::Databar,
        SymbolType::DatabarExp,
        SymbolType::Codabar,
        SymbolType::Code39,
        SymbolType::Code93,
        SymbolType::Code128,
        SymbolType::QrCode,
        SymbolType::SqCode,
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

/// A reference to a decoded barcode symbol
pub struct Symbol<'a> {
    inner: &'a zbar_symbol_t,
}

impl<'a> Symbol<'a> {
    pub(crate) fn from_ref(sym: &'a zbar_symbol_t) -> Self {
        Symbol { inner: sym }
    }

    /// Get the symbol type
    pub fn symbol_type(&self) -> SymbolType {
        self.inner.symbol_type
    }

    /// Get the decoded data as bytes
    pub fn data(&self) -> &[u8] {
        &self.inner.data
    }

    /// Get the decoded data as a string (if valid UTF-8)
    pub fn data_string(&self) -> Option<&str> {
        std::str::from_utf8(self.data()).ok()
    }

    /// Get the orientation of the barcode
    pub fn orientation(&self) -> Orientation {
        self.inner.orientation
    }

    /// Get the component symbols (for composite symbols like EAN+add-on or QR structured append)
    pub fn components(&self) -> Option<SymbolSet<'_>> {
        self.inner.components.as_ref().map(|set| SymbolSet {
            symbols: &set.symbols,
        })
    }
}

/// Iterator over symbols
pub struct SymbolIterator<'a> {
    iter: std::slice::Iter<'a, zbar_symbol_t>,
}

impl<'a> Iterator for SymbolIterator<'a> {
    type Item = Symbol<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(Symbol::from_ref)
    }
}

/// Collection of decoded symbols
pub struct SymbolSet<'a> {
    symbols: &'a [zbar_symbol_t],
}

impl<'a> SymbolSet<'a> {
    pub(crate) fn from_slice(symbols: &'a [zbar_symbol_t]) -> Self {
        Self { symbols }
    }

    /// Get an iterator over the symbols
    pub fn iter(&self) -> SymbolIterator<'a> {
        SymbolIterator {
            iter: self.symbols.iter(),
        }
    }

    /// Check if there are any symbols
    pub fn is_empty(&self) -> bool {
        self.symbols.is_empty()
    }

    /// Get the number of symbols
    pub fn len(&self) -> usize {
        self.symbols.len()
    }
}

impl<'a> IntoIterator for SymbolSet<'a> {
    type Item = Symbol<'a>;
    type IntoIter = SymbolIterator<'a>;

    fn into_iter(self) -> Self::IntoIter {
        SymbolIterator {
            iter: self.symbols.iter(),
        }
    }
}
