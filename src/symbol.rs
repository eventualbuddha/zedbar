//! Symbol management module
//!
//! Copyright 2007-2010 (c) Jeff Brown <spadix@users.sourceforge.net>
//! Rust port based on the C implementation
//!
//! Handles symbol lifecycle, reference counting, and data access.

use crate::{ffi::zbar_symbol_t, refcnt};
use libc::{c_char, c_int, c_uint, c_void};
use std::ptr;

// Symbol type constants (from zbar.h)
const ZBAR_PARTIAL: i32 = 1;
const ZBAR_EAN2: i32 = 2;
const ZBAR_EAN5: i32 = 5;
const ZBAR_EAN8: i32 = 8;
const ZBAR_UPCE: i32 = 9;
const ZBAR_ISBN10: i32 = 10;
const ZBAR_UPCA: i32 = 12;
const ZBAR_EAN13: i32 = 13;
const ZBAR_ISBN13: i32 = 14;
const ZBAR_COMPOSITE: i32 = 15;
const ZBAR_I25: i32 = 25;
const ZBAR_DATABAR: i32 = 34;
const ZBAR_DATABAR_EXP: i32 = 35;
const ZBAR_CODABAR: i32 = 38;
const ZBAR_CODE39: i32 = 39;
const ZBAR_CODE93: i32 = 93;
const ZBAR_CODE128: i32 = 128;
const ZBAR_QRCODE: i32 = 64;
const ZBAR_SQCODE: i32 = 80;
const ZBAR_SYMBOL: i32 = 0x00ff;

// Config type constants
const ZBAR_CFG_ENABLE: i32 = 0;
const ZBAR_CFG_ADD_CHECK: i32 = 1;
const ZBAR_CFG_EMIT_CHECK: i32 = 2;
const ZBAR_CFG_ASCII: i32 = 3;
const ZBAR_CFG_BINARY: i32 = 4;
const ZBAR_CFG_MIN_LEN: i32 = 32;
const ZBAR_CFG_MAX_LEN: i32 = 33;
const ZBAR_CFG_UNCERTAINTY: i32 = 64;
const ZBAR_CFG_POSITION: i32 = 128;
const ZBAR_CFG_X_DENSITY: i32 = 256;
const ZBAR_CFG_Y_DENSITY: i32 = 257;
// Modifier constants
const ZBAR_MOD_GS1: i32 = 0;
const ZBAR_MOD_AIM: i32 = 1;

// Orientation constants
const ZBAR_ORIENT_UP: i32 = 0;
const ZBAR_ORIENT_RIGHT: i32 = 1;
const ZBAR_ORIENT_DOWN: i32 = 2;
const ZBAR_ORIENT_LEFT: i32 = 3;

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

/// Free a symbol
///
/// # Safety
///
/// The symbol pointer must be valid and not previously freed.
pub unsafe fn symbol_free(sym: *mut zbar_symbol_t) {
    if !(*sym).syms.is_null() {
        zbar_symbol_set_ref((*sym).syms as *const _, -1);
        (*sym).syms = ptr::null_mut();
    }
    if !(*sym).pts.is_null() {
        libc::free((*sym).pts);
    }
    if (*sym).data_alloc != 0 && !(*sym).data.is_null() {
        libc::free((*sym).data as *mut c_void);
    }
    libc::free(sym as *mut c_void);
}

/// Adjust symbol reference count
///
/// # Safety
///
/// The symbol pointer must be valid.
pub unsafe fn symbol_refcnt(sym: *mut zbar_symbol_t, delta: c_int) {
    if refcnt(&mut (*sym).refcnt, delta) == 0 && delta <= 0 {
        symbol_free(sym);
    }
}

/// Create a new symbol set
///
/// # Safety
///
/// Allocates memory that must be freed with `symbol_set_free`.
pub unsafe fn symbol_set_create() -> *mut c_void {
    let syms = libc::calloc(1, std::mem::size_of::<CSymbolSet>()) as *mut CSymbolSet;
    refcnt(&mut (*syms).refcnt, 1);
    syms as *mut c_void
}

/// Free a symbol set
///
/// # Safety
///
/// The symbol set pointer must be valid and not previously freed.
pub unsafe fn symbol_set_free(syms: *mut c_void) {
    let syms = syms as *mut CSymbolSet;
    let mut sym = (*syms).head;
    while !sym.is_null() {
        let next = (*sym).next;
        (*sym).next = ptr::null_mut();
        symbol_refcnt(sym, -1);
        sym = next;
    }
    (*syms).head = ptr::null_mut();
    libc::free(syms as *mut c_void);
}

// Internal symbol set structure (C FFI)
#[repr(C)]
pub struct CSymbolSet {
    pub refcnt: c_int,
    pub nsyms: c_int,
    pub head: *mut zbar_symbol_t,
    pub tail: *mut zbar_symbol_t,
}

// C FFI exports

#[no_mangle]
pub unsafe extern "C" fn zbar_get_symbol_name(sym: c_int) -> *const c_char {
    get_symbol_name(sym).as_ptr() as *const c_char
}

#[no_mangle]
pub unsafe extern "C" fn zbar_get_config_name(cfg: c_int) -> *const c_char {
    get_config_name(cfg).as_ptr() as *const c_char
}

#[no_mangle]
pub unsafe extern "C" fn zbar_get_modifier_name(mod_type: c_int) -> *const c_char {
    get_modifier_name(mod_type).as_ptr() as *const c_char
}

#[no_mangle]
pub unsafe extern "C" fn zbar_get_orientation_name(orient: c_int) -> *const c_char {
    get_orientation_name(orient).as_ptr() as *const c_char
}

#[no_mangle]
pub unsafe extern "C" fn _zbar_get_symbol_hash(sym: c_int) -> c_int {
    get_symbol_hash(sym)
}

#[no_mangle]
pub unsafe extern "C" fn _zbar_symbol_free(sym: *mut zbar_symbol_t) {
    symbol_free(sym);
}

#[no_mangle]
pub unsafe extern "C" fn zbar_symbol_ref(sym: *const zbar_symbol_t, refs: c_int) {
    let sym = sym as *mut zbar_symbol_t;
    symbol_refcnt(sym, refs);
}

#[no_mangle]
pub unsafe extern "C" fn zbar_symbol_get_type(sym: *const zbar_symbol_t) -> c_int {
    (*sym).symbol_type
}

#[no_mangle]
pub unsafe extern "C" fn zbar_symbol_get_configs(sym: *const zbar_symbol_t) -> c_uint {
    (*sym).configs
}

#[no_mangle]
pub unsafe extern "C" fn zbar_symbol_get_modifiers(sym: *const zbar_symbol_t) -> c_uint {
    (*sym).modifiers
}

#[no_mangle]
pub unsafe extern "C" fn zbar_symbol_get_data(sym: *const zbar_symbol_t) -> *const c_char {
    (*sym).data
}

#[no_mangle]
pub unsafe extern "C" fn zbar_symbol_get_data_length(sym: *const zbar_symbol_t) -> c_uint {
    (*sym).datalen
}

#[no_mangle]
pub unsafe extern "C" fn zbar_symbol_get_count(sym: *const zbar_symbol_t) -> c_int {
    (*sym).cache_count
}

#[no_mangle]
pub unsafe extern "C" fn zbar_symbol_get_quality(sym: *const zbar_symbol_t) -> c_int {
    (*sym).quality
}

#[no_mangle]
pub unsafe extern "C" fn zbar_symbol_get_loc_size(sym: *const zbar_symbol_t) -> c_uint {
    (*sym).npts
}

#[no_mangle]
pub unsafe extern "C" fn zbar_symbol_get_loc_x(sym: *const zbar_symbol_t, idx: c_uint) -> c_int {
    if idx < (*sym).npts {
        let pts = (*sym).pts as *const Point;
        (*pts.add(idx as usize)).x
    } else {
        -1
    }
}

#[no_mangle]
pub unsafe extern "C" fn zbar_symbol_get_loc_y(sym: *const zbar_symbol_t, idx: c_uint) -> c_int {
    if idx < (*sym).npts {
        let pts = (*sym).pts as *const Point;
        (*pts.add(idx as usize)).y
    } else {
        -1
    }
}

#[no_mangle]
pub unsafe extern "C" fn zbar_symbol_get_orientation(sym: *const zbar_symbol_t) -> c_int {
    (*sym).orient
}

#[no_mangle]
pub unsafe extern "C" fn zbar_symbol_next(sym: *const zbar_symbol_t) -> *const zbar_symbol_t {
    if sym.is_null() {
        ptr::null()
    } else {
        (*sym).next
    }
}

#[no_mangle]
pub unsafe extern "C" fn zbar_symbol_get_components(sym: *const zbar_symbol_t) -> *const c_void {
    (*sym).syms
}

#[no_mangle]
pub unsafe extern "C" fn zbar_symbol_first_component(
    sym: *const zbar_symbol_t,
) -> *const zbar_symbol_t {
    if !sym.is_null() && !(*sym).syms.is_null() {
        let syms = (*sym).syms as *const CSymbolSet;
        (*syms).head
    } else {
        ptr::null()
    }
}

#[no_mangle]
pub unsafe extern "C" fn _zbar_symbol_set_create() -> *mut c_void {
    symbol_set_create()
}

#[no_mangle]
pub unsafe extern "C" fn _zbar_symbol_set_free(syms: *mut c_void) {
    symbol_set_free(syms);
}

#[no_mangle]
pub unsafe extern "C" fn zbar_symbol_set_ref(syms: *const c_void, delta: c_int) {
    let syms = syms as *mut CSymbolSet;
    if refcnt(&mut (*syms).refcnt, delta) == 0 && delta <= 0 {
        symbol_set_free(syms as *mut c_void);
    }
}

#[no_mangle]
pub unsafe extern "C" fn zbar_symbol_set_get_size(syms: *const c_void) -> c_int {
    let syms = syms as *const CSymbolSet;
    (*syms).nsyms
}

#[no_mangle]
pub unsafe extern "C" fn zbar_symbol_set_first_symbol(syms: *const c_void) -> *const zbar_symbol_t {
    let syms = syms as *const CSymbolSet;
    let sym = (*syms).tail;
    if !sym.is_null() {
        (*sym).next
    } else {
        (*syms).head
    }
}

#[no_mangle]
pub unsafe extern "C" fn zbar_symbol_set_first_unfiltered(
    syms: *const c_void,
) -> *const zbar_symbol_t {
    let syms = syms as *const CSymbolSet;
    (*syms).head
}

// Point structure for location data
#[repr(C)]
struct Point {
    x: c_int,
    y: c_int,
}

// Note: XML serialization (zbar_symbol_xml) and base64_encode are complex and
// rarely used in the core scanning functionality. They can be added later if needed.

// High-level Rust API types
use crate::ffi;

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
    ptr: *const std::ffi::c_void,
}

impl Symbol {
    pub(crate) unsafe fn from_ptr(ptr: *const std::ffi::c_void) -> Option<Self> {
        if ptr.is_null() {
            None
        } else {
            Some(Symbol { ptr })
        }
    }

    /// Get the symbol type
    pub fn symbol_type(&self) -> SymbolType {
        let type_code = unsafe { ffi::zbar_symbol_get_type(self.ptr) };
        SymbolType::from(type_code)
    }

    /// Get the decoded data as bytes
    pub fn data(&self) -> &[u8] {
        unsafe {
            let data_ptr = ffi::zbar_symbol_get_data(self.ptr);
            let data_len = ffi::zbar_symbol_get_data_length(self.ptr) as usize;
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
        let next_ptr = unsafe { ffi::zbar_symbol_next(self.ptr) };
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
        SymbolSet {
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
