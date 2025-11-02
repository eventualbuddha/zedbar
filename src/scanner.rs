//! Image scanner for finding barcodes in 2D images

use crate::image::Image;
use crate::img_scanner::{zbar_image_scanner_create, zbar_image_scanner_t};
use crate::{Error, Result, SymbolType};

/// Configuration options for barcode scanning
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Config {
    Enable = 0,
    AddCheck = 1,
    EmitCheck = 2,
    Ascii = 3,
    Binary = 4,
    MinLen = 0x20,
    MaxLen = 0x21,
    Uncertainty = 0x40,
    Position = 0x80,
    TestInverted = 0x81,
    XDensity = 0x100,
    YDensity = 0x101,
}

/// Image scanner that can find barcodes in 2D images
pub struct Scanner {
    ptr: *mut zbar_image_scanner_t,
}

impl Scanner {
    /// Create a new image scanner
    pub fn new() -> Self {
        let ptr = unsafe { zbar_image_scanner_create() };
        Scanner { ptr }
    }

    /// Configure the scanner for a specific symbology
    pub fn set_config(&mut self, symbology: SymbolType, config: Config, value: i32) -> Result<()> {
        let result = unsafe { (&mut *self.ptr).set_config(symbology, config as i32, value) };

        if result == 0 {
            Ok(())
        } else {
            Err(Error::Invalid)
        }
    }

    /// Scan an image for barcodes
    pub fn scan(&mut self, image: &mut Image) -> Result<i32> {
        unsafe { (*self.ptr).scan_image(image.as_mut_image()) }
    }
}

impl Drop for Scanner {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe { drop(Box::from_raw(self.ptr)) }
        }
    }
}

impl Default for Scanner {
    fn default() -> Self {
        Self::new()
    }
}
