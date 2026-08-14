//! Image scanner configuration
//!
//! Type-safe configuration structures for the image scanner, replacing C-style
//! bitfields and integer arrays with idiomatic Rust types.

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
}

impl Default for ImageScannerConfig {
    fn default() -> Self {
        Self {
            position_tracking: true,
            test_inverted: false,
            x_density: 1,
            y_density: 1,
            ean_composite: false,
        }
    }
}
