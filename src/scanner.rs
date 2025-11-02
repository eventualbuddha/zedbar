//! Image scanner for finding barcodes in 2D images

use crate::image::Image;
use crate::img_scanner::zbar_image_scanner_t;
use crate::{Result, SymbolType};

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
    scanner: zbar_image_scanner_t,
}

impl Scanner {
    /// Create a new image scanner
    pub fn new() -> Self {
        Scanner {
            scanner: zbar_image_scanner_t::default(),
        }
    }

    /// Configure the scanner for a specific symbology
    pub fn set_config(&mut self, symbology: SymbolType, config: Config, value: i32) -> Result<()> {
        self.scanner.set_config(symbology, config as i32, value)
    }

    /// Scan an image for barcodes
    pub fn scan(&mut self, image: &mut Image) -> Result<i32> {
        self.scanner.scan_image(image.as_mut_image())
    }
}

impl Default for Scanner {
    fn default() -> Self {
        Self::new()
    }
}
