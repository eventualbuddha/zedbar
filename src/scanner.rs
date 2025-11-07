//! Image scanner for finding barcodes in 2D images
//!
//! The [`Scanner`] is the main entry point for barcode scanning operations.
//! It processes grayscale images to detect and decode barcodes.
//!
//! # Example
//!
//! ```no_run
//! use zedbar::{Image, Scanner};
//!
//! // Create a scanner with default settings
//! let mut scanner = Scanner::new();
//!
//! // Or with custom configuration
//! use zedbar::config::*;
//! let config = DecoderConfig::new()
//!     .enable(QrCode)
//!     .enable(Ean13);
//! let mut scanner = Scanner::with_config(config);
//!
//! // Scan an image
//! # let data = vec![0u8; 640 * 480];
//! let mut image = Image::from_gray(&data, 640, 480).unwrap();
//! let symbols = scanner.scan(&mut image);
//!
//! for symbol in symbols {
//!     println!("{:?}: {:?}", symbol.symbol_type(), symbol.data_string());
//! }
//! ```

use crate::config::DecoderConfig;
use crate::image::Image;
use crate::img_scanner::zbar_image_scanner_t;
use crate::symbol::Symbol;

/// Image scanner that can find barcodes in 2D images
///
/// # Example
/// ```no_run
/// use zedbar::config::*;
/// use zedbar::{Scanner, DecoderConfig, Image};
///
/// // Create scanner with type-safe configuration
/// let config = DecoderConfig::new()
///     .enable(Ean13)
///     .enable(QrCode)
///     .set_binary(QrCode, true)
///     .position_tracking(true)
///     .scan_density(1, 1);
///
/// let mut scanner = Scanner::with_config(config);
///
/// // Scan an image
/// let data = vec![0u8; 640 * 480];
/// let mut image = Image::from_gray(&data, 640, 480).unwrap();
/// let symbols = scanner.scan(&mut image);
/// ```
pub struct Scanner {
    scanner: zbar_image_scanner_t,
}

impl Scanner {
    /// Create a new image scanner with default configuration
    ///
    /// For more control over the configuration, use [`Scanner::with_config()`].
    pub fn new() -> Self {
        Scanner {
            scanner: zbar_image_scanner_t::default(),
        }
    }

    /// Create a new image scanner with custom configuration
    ///
    /// This is the recommended way to create a scanner with specific settings.
    ///
    /// # Example
    /// ```no_run
    /// use zedbar::config::*;
    /// use zedbar::{Scanner, DecoderConfig};
    ///
    /// let config = DecoderConfig::new()
    ///     .enable(Ean13)
    ///     .enable(Code39)
    ///     .set_length_limits(Code39, 4, 20)
    ///     .position_tracking(true);
    ///
    /// let scanner = Scanner::with_config(config);
    /// ```
    pub fn with_config(config: DecoderConfig) -> Self {
        Scanner {
            scanner: zbar_image_scanner_t::with_config(config),
        }
    }

    /// Scan an image for barcodes
    pub fn scan(&mut self, image: &mut Image) -> Vec<Symbol> {
        self.scanner.scan_image(image.as_mut_image())
    }
}

impl Default for Scanner {
    fn default() -> Self {
        Self::new()
    }
}
