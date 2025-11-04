//! Image scanner for finding barcodes in 2D images

use crate::config::DecoderConfig;
use crate::image::Image;
use crate::img_scanner::zbar_image_scanner_t;
use crate::Result;

/// Image scanner that can find barcodes in 2D images
///
/// # Example
/// ```no_run
/// use zbar::config::*;
/// use zbar::{Scanner, DecoderConfig, Image};
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
/// let num_symbols = scanner.scan(&mut image).unwrap();
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
    /// use zbar::config::*;
    /// use zbar::{Scanner, DecoderConfig};
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
    ///
    /// Returns the number of symbols found in the image.
    pub fn scan(&mut self, image: &mut Image) -> Result<i32> {
        self.scanner.scan_image(image.as_mut_image())
    }
}

impl Default for Scanner {
    fn default() -> Self {
        Self::new()
    }
}
