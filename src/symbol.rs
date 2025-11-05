//! Symbol management module
//!
//! Handles symbol lifecycle, reference counting, and data access.
//!
//! Rust port based on C code from the ZBar library.
//! Original C code copyright (C) 2007-2010 Jeff Brown <spadix@users.sourceforge.net>
//! Licensed under LGPL 3.0 or later

use crate::qrcode::qrdec::qr_point;
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

#[derive(Default, Clone)]
pub struct Symbol {
    symbol_type: SymbolType,
    pub(crate) modifiers: u32,
    pub(crate) data: Vec<u8>,
    pub(crate) pts: Vec<qr_point>,
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
        self.pts.push([x, y]);
    }
}

// Public API
impl Symbol {
    /// Get the symbol type of this symbol.
    pub fn symbol_type(&self) -> SymbolType {
        self.symbol_type
    }

    /// Get the decoded data as bytes
    pub fn data(&self) -> &[u8] {
        &self.data
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
