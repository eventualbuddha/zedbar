//! Symbol types and decoding results

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