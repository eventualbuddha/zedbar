//! Symbol management module
//!
//! Copyright 2007-2010 (c) Jeff Brown <spadix@users.sourceforge.net>
//! Rust port based on the C implementation
//!
//! Handles symbol lifecycle, reference counting, and data access.

use crate::{
    decoder_types::{
        ZBAR_CFG_ADD_CHECK, ZBAR_CFG_ASCII, ZBAR_CFG_BINARY, ZBAR_CFG_EMIT_CHECK, ZBAR_CFG_ENABLE,
        ZBAR_CFG_MAX_LEN, ZBAR_CFG_MIN_LEN, ZBAR_CFG_POSITION, ZBAR_CFG_UNCERTAINTY,
        ZBAR_CFG_X_DENSITY, ZBAR_CFG_Y_DENSITY, ZBAR_MOD_AIM, ZBAR_MOD_GS1, ZBAR_ORIENT_DOWN,
        ZBAR_ORIENT_LEFT, ZBAR_ORIENT_RIGHT, ZBAR_ORIENT_UP,
    },
    img_scanner::zbar_symbol_set_t,
    qrcode::qr_point,
    refcnt,
};
use libc::{c_int, c_uint};
use std::{fmt::Display, ptr};

/// Get config name
pub fn get_config_name(cfg: i32) -> &'static str {
    match cfg {
        ZBAR_CFG_ENABLE => "ENABLE",
        ZBAR_CFG_ADD_CHECK => "ADD_CHECK",
        ZBAR_CFG_EMIT_CHECK => "EMIT_CHECK",
        ZBAR_CFG_ASCII => "ASCII",
        ZBAR_CFG_BINARY => "BINARY",
        ZBAR_CFG_MIN_LEN => "MIN_LEN",
        ZBAR_CFG_MAX_LEN => "MAX_LEN",
        ZBAR_CFG_UNCERTAINTY => "UNCERTAINTY",
        ZBAR_CFG_POSITION => "POSITION",
        ZBAR_CFG_X_DENSITY => "X_DENSITY",
        ZBAR_CFG_Y_DENSITY => "Y_DENSITY",
        _ => "",
    }
}

/// Get modifier name
pub fn get_modifier_name(mod_type: i32) -> &'static str {
    match mod_type {
        ZBAR_MOD_GS1 => "GS1",
        ZBAR_MOD_AIM => "AIM",
        _ => "",
    }
}

/// Get orientation name
pub fn get_orientation_name(orient: i32) -> &'static str {
    match orient {
        ZBAR_ORIENT_UP => "UP",
        ZBAR_ORIENT_RIGHT => "RIGHT",
        ZBAR_ORIENT_DOWN => "DOWN",
        ZBAR_ORIENT_LEFT => "LEFT",
        _ => "UNKNOWN",
    }
}

#[derive(Default)]
pub(crate) struct zbar_symbol_t {
    pub(crate) symbol_type: SymbolType,
    pub(crate) configs: c_uint,
    pub(crate) modifiers: c_uint,
    pub(crate) data: Vec<u8>,
    pub(crate) pts: Vec<qr_point>,
    pub(crate) orient: c_int,
    pub(crate) refcnt: c_int,
    pub(crate) next: *mut zbar_symbol_t,
    pub(crate) syms: *mut zbar_symbol_set_t,
    pub(crate) quality: c_int,
}

impl zbar_symbol_t {
    pub(crate) fn add_point(&mut self, x: c_int, y: c_int) {
        self.pts.push([x, y]);
    }
}

impl Drop for zbar_symbol_t {
    fn drop(&mut self) {
        if !self.syms.is_null() {
            unsafe {
                zbar_symbol_set_ref(self.syms, -1);
            }
            self.syms = ptr::null_mut();
        }
    }
}

/// Free a symbol
///
/// # Safety
///
/// The symbol pointer must be valid and not previously freed.
pub(crate) unsafe fn symbol_free(sym: *mut zbar_symbol_t) {
    drop(Box::from_raw(sym));
}

/// Adjust symbol reference count
pub(crate) unsafe fn symbol_refcnt(sym: &mut zbar_symbol_t, delta: c_int) {
    if refcnt!(sym.refcnt, delta) == 0 && delta <= 0 {
        symbol_free(sym as *mut _);
    }
}

/// Create a new symbol set
///
/// # Safety
///
/// Allocates memory that must be freed with `symbol_set_free`.
pub(crate) unsafe fn symbol_set_create() -> *mut zbar_symbol_set_t {
    let mut symbol_set = Box::new(zbar_symbol_set_t::default());
    refcnt!(symbol_set.refcnt, 1);
    Box::into_raw(symbol_set)
}

/// Free a symbol set
///
/// # Safety
///
/// The symbol set pointer must be valid and not previously freed.
pub(crate) unsafe fn symbol_set_free(syms: *mut zbar_symbol_set_t) {
    drop(Box::from_raw(syms));
}

/// Allocate a zeroed symbol instance suitable for initialization.
pub(crate) unsafe fn symbol_alloc_zeroed() -> *mut zbar_symbol_t {
    let symbol = Box::new(zbar_symbol_t::default());
    Box::into_raw(symbol)
}

pub(crate) unsafe fn zbar_symbol_next(sym: *const zbar_symbol_t) -> *const zbar_symbol_t {
    if sym.is_null() {
        ptr::null()
    } else {
        let sym = &*sym;
        sym.next
    }
}

pub(crate) unsafe fn zbar_symbol_set_ref(syms: *mut zbar_symbol_set_t, delta: c_int) {
    if syms.is_null() {
        return;
    }
    let syms_ref = &mut *syms;

    if refcnt!(syms_ref.refcnt, delta) == 0 && delta <= 0 {
        symbol_set_free(syms);
    }
}

// High-level Rust API types

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, PartialOrd, Ord)]
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

    pub(crate) fn hash(self) -> i8 {
        match self {
            Self::SqCode => 1,
            Self::Code128 => 2,
            Self::Ean13 => 3,
            Self::Upca => 4,
            Self::Ean8 => 5,
            Self::Upce => 6,
            Self::Isbn13 => 7,
            Self::Isbn10 => 8,
            Self::Code39 => 9,
            Self::I25 => 10,
            Self::QrCode => 12,
            Self::Databar => 13,
            Self::DatabarExp => 14,
            Self::Code93 => 15,
            Self::Ean2 => 16,
            Self::Ean5 => 17,
            Self::Composite => 18,
            Self::Codabar => 19,
            _ => -1,
        }
    }
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

/// A decoded barcode symbol
pub struct Symbol {
    ptr: *const zbar_symbol_t,
}

impl Symbol {
    pub(crate) unsafe fn from_ptr(ptr: *const zbar_symbol_t) -> Option<Self> {
        if ptr.is_null() {
            None
        } else {
            Some(Symbol { ptr })
        }
    }

    /// Get the symbol type
    pub fn symbol_type(&self) -> SymbolType {
        unsafe { (*self.ptr).symbol_type }
    }

    /// Get the decoded data as bytes
    pub fn data(&self) -> &[u8] {
        unsafe { &(*self.ptr).data }
    }

    /// Get the decoded data as a string (if valid UTF-8)
    pub fn data_string(&self) -> Option<&str> {
        std::str::from_utf8(self.data()).ok()
    }

    /// Get the next symbol in the result set
    pub fn next(&self) -> Option<Symbol> {
        let next_ptr = unsafe { zbar_symbol_next(self.ptr) };
        unsafe { Symbol::from_ptr(next_ptr) }
    }
}

/// Iterator over symbols
pub struct SymbolIterator {
    current: Option<Symbol>,
}

impl SymbolIterator {
    pub(crate) fn new(first: Option<Symbol>) -> Self {
        SymbolIterator { current: first }
    }
}

impl Iterator for SymbolIterator {
    type Item = Symbol;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(current) = self.current.take() {
            self.current = current.next();
            Some(current)
        } else {
            None
        }
    }
}

/// Collection of decoded symbols
pub struct SymbolSet {
    symbols: Option<Symbol>,
}

impl SymbolSet {
    pub(crate) fn new(first_symbol: Option<Symbol>) -> Self {
        Self {
            symbols: first_symbol,
        }
    }

    /// Get an iterator over the symbols
    pub fn iter(&self) -> SymbolIterator {
        SymbolIterator::new(self.symbols.as_ref().map(|s| Symbol { ptr: s.ptr }))
    }

    /// Check if there are any symbols
    pub fn is_empty(&self) -> bool {
        self.symbols.is_none()
    }
}

impl IntoIterator for SymbolSet {
    type Item = Symbol;
    type IntoIter = SymbolIterator;

    fn into_iter(self) -> Self::IntoIter {
        SymbolIterator::new(self.symbols)
    }
}
