//! Image scanner configuration
//!
//! Type-safe configuration structures for the image scanner, replacing C-style
//! bitfields and integer arrays with idiomatic Rust types.

use crate::SymbolType;
use std::collections::HashMap;

/// Configuration options for the image scanner
#[derive(Debug, Clone)]
pub(crate) struct ImageScannerConfig {
    /// Position tracking enabled
    pub(crate) position_tracking: bool,

    /// Test inverted images
    pub(crate) test_inverted: bool,

    /// Horizontal scan density
    pub(crate) x_density: u32,

    /// Vertical scan density
    pub(crate) y_density: u32,

    /// EAN composite symbology enabled
    pub(crate) ean_composite: bool,

    /// Per-symbology uncertainty thresholds
    pub(crate) uncertainty: HashMap<SymbolType, u32>,
}

impl Default for ImageScannerConfig {
    fn default() -> Self {
        let mut config = Self {
            position_tracking: true,
            test_inverted: false,
            x_density: 1,
            y_density: 1,
            ean_composite: false,
            uncertainty: HashMap::new(),
        };

        // Set default uncertainty values
        config.uncertainty.insert(SymbolType::Ean2, 2);
        config.uncertainty.insert(SymbolType::Ean5, 2);
        config.uncertainty.insert(SymbolType::Ean8, 2);
        config.uncertainty.insert(SymbolType::Ean13, 2);
        config.uncertainty.insert(SymbolType::Upca, 2);
        config.uncertainty.insert(SymbolType::Upce, 2);
        config.uncertainty.insert(SymbolType::Isbn10, 2);
        config.uncertainty.insert(SymbolType::Isbn13, 2);
        config.uncertainty.insert(SymbolType::I25, 2);
        config.uncertainty.insert(SymbolType::Databar, 2);
        config.uncertainty.insert(SymbolType::DatabarExp, 2);
        config.uncertainty.insert(SymbolType::Codabar, 1);
        config.uncertainty.insert(SymbolType::Code39, 0);
        config.uncertainty.insert(SymbolType::Code93, 0);
        config.uncertainty.insert(SymbolType::Code128, 0);
        config.uncertainty.insert(SymbolType::QrCode, 0);
        config.uncertainty.insert(SymbolType::SqCode, 0);
        config.uncertainty.insert(SymbolType::Composite, 0);

        config
    }
}
