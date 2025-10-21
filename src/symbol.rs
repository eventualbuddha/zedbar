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
        ZBAR_CFG_X_DENSITY, ZBAR_CFG_Y_DENSITY, ZBAR_CODABAR, ZBAR_CODE128, ZBAR_CODE39,
        ZBAR_CODE93, ZBAR_COMPOSITE, ZBAR_DATABAR, ZBAR_DATABAR_EXP, ZBAR_EAN13, ZBAR_EAN2,
        ZBAR_EAN5, ZBAR_EAN8, ZBAR_I25, ZBAR_ISBN10, ZBAR_ISBN13, ZBAR_MOD_AIM, ZBAR_MOD_GS1,
        ZBAR_ORIENT_DOWN, ZBAR_ORIENT_LEFT, ZBAR_ORIENT_RIGHT, ZBAR_ORIENT_UP, ZBAR_PARTIAL,
        ZBAR_QRCODE, ZBAR_SQCODE, ZBAR_SYMBOL, ZBAR_UPCA, ZBAR_UPCE,
    },
    ffi::refcnt,
    img_scanner::zbar_symbol_set_t,
};
use libc::{c_char, c_int, c_uint, c_ulong, c_void};
use std::{mem::size_of, ptr};

const NUM_SYMS: usize = 20;

// Symbol hash table for fast lookup
static SYMBOL_HASH: [i8; (ZBAR_CODE128 + 1) as usize] = {
    let mut hash = [-1i8; (ZBAR_CODE128 + 1) as usize];
    hash[ZBAR_SQCODE as usize] = 1;
    hash[ZBAR_CODE128 as usize] = 2;
    hash[ZBAR_EAN13 as usize] = 3;
    hash[ZBAR_UPCA as usize] = 4;
    hash[ZBAR_EAN8 as usize] = 5;
    hash[ZBAR_UPCE as usize] = 6;
    hash[ZBAR_ISBN13 as usize] = 7;
    hash[ZBAR_ISBN10 as usize] = 8;
    hash[ZBAR_CODE39 as usize] = 9;
    hash[ZBAR_I25 as usize] = 10;
    hash[ZBAR_QRCODE as usize] = 12;
    hash[ZBAR_DATABAR as usize] = 13;
    hash[ZBAR_DATABAR_EXP as usize] = 14;
    hash[ZBAR_CODE93 as usize] = 15;
    hash[ZBAR_EAN2 as usize] = 16;
    hash[ZBAR_EAN5 as usize] = 17;
    hash[ZBAR_COMPOSITE as usize] = 18;
    hash[ZBAR_CODABAR as usize] = 19;
    hash
};

/// Get symbol type name
pub fn get_symbol_name(sym: i32) -> &'static str {
    match sym & ZBAR_SYMBOL {
        ZBAR_EAN2 => "EAN-2",
        ZBAR_EAN5 => "EAN-5",
        ZBAR_EAN8 => "EAN-8",
        ZBAR_UPCE => "UPC-E",
        ZBAR_ISBN10 => "ISBN-10",
        ZBAR_UPCA => "UPC-A",
        ZBAR_EAN13 => "EAN-13",
        ZBAR_ISBN13 => "ISBN-13",
        ZBAR_COMPOSITE => "COMPOSITE",
        ZBAR_I25 => "I2/5",
        ZBAR_DATABAR => "DataBar",
        ZBAR_DATABAR_EXP => "DataBar-Exp",
        ZBAR_CODABAR => "Codabar",
        ZBAR_CODE39 => "CODE-39",
        ZBAR_CODE93 => "CODE-93",
        ZBAR_CODE128 => "CODE-128",
        ZBAR_QRCODE => "QR-Code",
        ZBAR_SQCODE => "SQ-Code",
        _ => "UNKNOWN",
    }
}

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

/// Get symbol hash index
pub fn get_symbol_hash(sym: i32) -> i32 {
    debug_assert!((ZBAR_PARTIAL..=ZBAR_CODE128).contains(&sym));
    let h = SYMBOL_HASH[sym as usize];
    debug_assert!(h >= 0 && (h as usize) < NUM_SYMS);
    h as i32
}

#[derive(Default)]
#[allow(non_camel_case_types)]
pub struct zbar_symbol_t {
    pub symbol_type: c_int,
    pub configs: c_uint,
    pub modifiers: c_uint,
    pub data_alloc: c_uint,
    pub datalen: c_uint,
    pub data: *mut c_char,
    pub pts_alloc: c_uint,
    pub npts: c_uint,
    pub pts: *mut c_void,
    pub orient: c_int,
    pub refcnt: c_int,
    pub next: *mut zbar_symbol_t,
    pub syms: *mut zbar_symbol_set_t,
    pub time: c_ulong,
    pub quality: c_int,
}

impl Drop for zbar_symbol_t {
    fn drop(&mut self) {
        if !self.syms.is_null() {
            unsafe {
                zbar_symbol_set_ref(self.syms, -1);
            }
            self.syms = ptr::null_mut();
        }
        if !self.pts.is_null() {
            unsafe {
                libc::free(self.pts);
            }
        }
        if self.data_alloc != 0 && !self.data.is_null() {
            unsafe {
                libc::free(self.data as *mut c_void);
            }
        }
    }
}

/// Free a symbol
///
/// # Safety
///
/// The symbol pointer must be valid and not previously freed.
pub unsafe fn symbol_free(sym: *mut zbar_symbol_t) {
    drop(Box::from_raw(sym));
}

/// Adjust symbol reference count
///
/// # Safety
///
/// The symbol pointer must be valid.
pub unsafe fn symbol_refcnt(sym: *mut zbar_symbol_t, delta: c_int) {
    if sym.is_null() {
        return;
    }
    let sym = &mut *sym;

    if refcnt(&mut sym.refcnt, delta) == 0 && delta <= 0 {
        symbol_free(sym as *mut _);
    }
}

/// Create a new symbol set
///
/// # Safety
///
/// Allocates memory that must be freed with `symbol_set_free`.
pub unsafe fn symbol_set_create() -> *mut zbar_symbol_set_t {
    let mut symbol_set = Box::new(zbar_symbol_set_t::default());
    refcnt(&mut symbol_set.refcnt, 1);
    Box::into_raw(symbol_set)
}

/// Free a symbol set
///
/// # Safety
///
/// The symbol set pointer must be valid and not previously freed.
pub unsafe fn symbol_set_free(syms: *mut zbar_symbol_set_t) {
    drop(Box::from_raw(syms));
}

/// Allocate a zeroed symbol instance suitable for initialization.
pub unsafe fn symbol_alloc_zeroed() -> *mut zbar_symbol_t {
    let symbol = Box::new(zbar_symbol_t::default());
    Box::into_raw(symbol)
}

/// Release any allocated symbol data buffer and reset metadata.
pub unsafe fn symbol_clear_data(sym: *mut zbar_symbol_t) {
    if sym.is_null() {
        return;
    }
    let sym = &mut *sym;

    if !sym.data.is_null() {
        libc::free(sym.data as *mut c_void);
    }
    sym.data = ptr::null_mut();
    sym.datalen = 0;
    sym.data_alloc = 0;
}

/// Ensure the symbol data buffer can hold at least `capacity` bytes (including null terminator).
/// Returns `true` on success and leaves the buffer unchanged on failure.
pub unsafe fn symbol_reserve_data(sym: *mut zbar_symbol_t, capacity: usize) -> bool {
    if sym.is_null() {
        return false;
    }
    let sym = &mut *sym;

    if capacity == 0 {
        symbol_clear_data(sym as *mut _);
        return true;
    }
    if sym.data_alloc as usize >= capacity {
        return true;
    }

    let new_ptr = libc::realloc(sym.data as *mut c_void, capacity) as *mut c_char;
    if new_ptr.is_null() {
        return false;
    }

    sym.data = new_ptr;
    sym.data_alloc = capacity as c_uint;
    true
}

/// Ensure the symbol point array can store at least `capacity` points.
/// Returns `true` on success and leaves the buffer unchanged on failure.
pub unsafe fn symbol_reserve_points(sym: *mut zbar_symbol_t, capacity: u32) -> bool {
    if sym.is_null() {
        return false;
    }
    let sym = &mut *sym;

    if capacity == 0 {
        if !sym.pts.is_null() {
            libc::free(sym.pts);
            sym.pts = ptr::null_mut();
        }
        sym.pts_alloc = 0;
        return true;
    }
    if sym.pts_alloc >= capacity {
        return true;
    }

    let new_size = capacity as usize * size_of::<Point>();
    let new_ptr = libc::realloc(sym.pts, new_size);
    if new_ptr.is_null() {
        return false;
    }

    sym.pts = new_ptr;
    sym.pts_alloc = capacity;
    true
}

// C FFI exports

pub unsafe fn zbar_get_symbol_name(sym: c_int) -> *const c_char {
    get_symbol_name(sym).as_ptr() as *const c_char
}

pub unsafe fn _zbar_get_symbol_hash(sym: c_int) -> c_int {
    get_symbol_hash(sym)
}

pub unsafe fn zbar_symbol_next(sym: *const zbar_symbol_t) -> *const zbar_symbol_t {
    if sym.is_null() {
        ptr::null()
    } else {
        let sym = &*sym;
        sym.next
    }
}

pub unsafe fn zbar_symbol_set_ref(syms: *mut zbar_symbol_set_t, delta: c_int) {
    if syms.is_null() {
        return;
    }
    let syms_ref = &mut *syms;

    if refcnt(&mut syms_ref.refcnt, delta) == 0 && delta <= 0 {
        symbol_set_free(syms);
    }
}

// Point storage for location data (two consecutive c_int values)
type Point = [c_int; 2];

pub unsafe fn _zbar_symbol_add_point(sym: *mut zbar_symbol_t, x: c_int, y: c_int) {
    if sym.is_null() {
        return;
    }
    let sym = &mut *sym;

    let i = sym.npts as usize;
    sym.npts += 1;

    if sym.npts >= sym.pts_alloc {
        let new_capacity = sym.pts_alloc + 1;
        if !symbol_reserve_points(sym as *mut _, new_capacity) {
            sym.npts -= 1;
            return;
        }
    }

    let pts = sym.pts as *mut Point;
    (*pts.add(i))[0] = x;
    (*pts.add(i))[1] = y;
}

// High-level Rust API types

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SymbolType {
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

impl From<i32> for SymbolType {
    fn from(value: i32) -> Self {
        match value {
            0 => SymbolType::None,
            1 => SymbolType::Partial,
            2 => SymbolType::Ean2,
            5 => SymbolType::Ean5,
            8 => SymbolType::Ean8,
            9 => SymbolType::Upce,
            10 => SymbolType::Isbn10,
            12 => SymbolType::Upca,
            13 => SymbolType::Ean13,
            14 => SymbolType::Isbn13,
            15 => SymbolType::Composite,
            25 => SymbolType::I25,
            34 => SymbolType::Databar,
            35 => SymbolType::DatabarExp,
            38 => SymbolType::Codabar,
            39 => SymbolType::Code39,
            64 => SymbolType::QrCode,
            80 => SymbolType::SqCode,
            93 => SymbolType::Code93,
            128 => SymbolType::Code128,
            _ => SymbolType::None,
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
        let type_code = unsafe { (*self.ptr).symbol_type };
        SymbolType::from(type_code)
    }

    /// Get the decoded data as bytes
    pub fn data(&self) -> &[u8] {
        unsafe {
            let data_ptr = (*self.ptr).data;
            let data_len = (*self.ptr).datalen as usize;
            if data_ptr.is_null() || data_len == 0 {
                &[]
            } else {
                std::slice::from_raw_parts(data_ptr as *const u8, data_len)
            }
        }
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

    pub fn add_point(&mut self, x: c_int, y: c_int) {
        unsafe {
            _zbar_symbol_add_point(self.ptr as *mut _, x, y);
        }
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
