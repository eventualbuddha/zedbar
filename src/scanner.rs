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
//! let result = scanner.scan(&mut image);
//!
//! for symbol in &result {
//!     println!("{:?}: {:?}", symbol.symbol_type(), symbol.data_string());
//! }
//! ```

use crate::config::DecoderConfig;
use crate::image::Image;
use crate::img_scanner::ImageScanner;
use crate::symbol::Symbol;

/// A region where QR finder patterns were detected but decoding failed.
///
/// The bounding box describes the area in the original image where
/// finder pattern lines were found. To attempt decoding, crop the
/// image to this region (with padding for the quiet zone) and
/// upscale before re-scanning.
///
/// # Example
///
/// ```no_run
/// # use zedbar::{Image, Scanner};
/// # let data = vec![0u8; 800 * 600];
/// # let mut image = Image::from_gray(&data, 800, 600).unwrap();
/// # let mut scanner = Scanner::new();
/// let result = scanner.scan(&mut image);
///
/// if let Some(region) = result.finder_region() {
///     let pad = region.width.max(region.height) / 2;
///     let x = region.x.saturating_sub(pad);
///     let y = region.y.saturating_sub(pad);
///     let w = (region.width + 2 * pad).min(image.width() - x);
///     let h = (region.height + 2 * pad).min(image.height() - y);
///
///     if let Some(cropped) = image.crop(x, y, w, h) {
///         if let Some(mut upscaled) = cropped.upscale(4) {
///             let retry = scanner.scan(&mut upscaled);
///             for symbol in retry.symbols() {
///                 println!("Recovered: {}", symbol.data_string().unwrap_or(""));
///             }
///         }
///     }
/// }
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FinderRegion {
    /// X coordinate of the top-left corner (pixels).
    pub x: u32,
    /// Y coordinate of the top-left corner (pixels).
    pub y: u32,
    /// Width of the bounding box (pixels).
    pub width: u32,
    /// Height of the bounding box (pixels).
    pub height: u32,
}

/// Result of scanning an image for barcodes.
///
/// Contains decoded symbols and metadata about regions where
/// QR finder patterns were detected but decoding failed.
///
/// `ScanResult` implements [`Deref<Target = [Symbol]>`](std::ops::Deref) and
/// [`IntoIterator`], so most code that previously worked with `Vec<Symbol>`
/// continues to work unchanged.
pub struct ScanResult {
    symbols: Vec<Symbol>,
    finder_region: Option<FinderRegion>,
}

impl ScanResult {
    pub(crate) fn new(symbols: Vec<Symbol>, finder_region: Option<FinderRegion>) -> Self {
        Self {
            symbols,
            finder_region,
        }
    }

    /// The decoded barcode symbols.
    pub fn symbols(&self) -> &[Symbol] {
        &self.symbols
    }

    /// Consumes the result and returns the decoded symbols.
    pub fn into_symbols(self) -> Vec<Symbol> {
        self.symbols
    }

    /// Region where QR finder patterns were detected but no QR code
    /// was successfully decoded.
    ///
    /// This is the bounding area that likely contains QR codes too small or
    /// distorted to decode at the current resolution. Cropping and upscaling
    /// this region may yield successful decodes.
    ///
    /// Returns `None` when no undecoded regions were found, or when the
    /// `qrcode` feature is disabled.
    pub fn finder_region(&self) -> Option<&FinderRegion> {
        self.finder_region.as_ref()
    }
}

// Backward-compatible: `for symbol in scanner.scan(&mut img)` still works.
impl IntoIterator for ScanResult {
    type Item = Symbol;
    type IntoIter = std::vec::IntoIter<Symbol>;

    fn into_iter(self) -> Self::IntoIter {
        self.symbols.into_iter()
    }
}

impl<'a> IntoIterator for &'a ScanResult {
    type Item = &'a Symbol;
    type IntoIter = std::slice::Iter<'a, Symbol>;

    fn into_iter(self) -> Self::IntoIter {
        self.symbols.iter()
    }
}

// Backward-compatible: `result.is_empty()`, `result.len()`, indexing all work.
impl std::ops::Deref for ScanResult {
    type Target = [Symbol];
    fn deref(&self) -> &[Symbol] {
        &self.symbols
    }
}

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
///     .position_tracking(true)
///     .scan_density(1, 1);
///
/// let mut scanner = Scanner::with_config(config);
///
/// // Scan an image
/// let data = vec![0u8; 640 * 480];
/// let mut image = Image::from_gray(&data, 640, 480).unwrap();
/// let result = scanner.scan(&mut image);
/// ```
pub struct Scanner {
    scanner: ImageScanner,
    retry_undecoded_regions: bool,
}

impl Scanner {
    /// Create a new image scanner with default configuration
    ///
    /// For more control over the configuration, use [`Scanner::with_config()`].
    pub fn new() -> Self {
        Self {
            scanner: ImageScanner::default(),
            retry_undecoded_regions: false,
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
        let retry = config.retry_undecoded_regions;
        Self {
            scanner: ImageScanner::with_config(config),
            retry_undecoded_regions: retry,
        }
    }

    /// Scan an image for barcodes
    ///
    /// Returns a [`ScanResult`] containing decoded symbols and any
    /// undecoded QR finder regions. Check [`ScanResult::finder_region()`]
    /// to find areas that may contain QR codes too small to decode at
    /// the current resolution.
    ///
    /// When [`DecoderConfig::retry_undecoded_regions`] is enabled, each
    /// undecoded finder region is automatically cropped, upscaled 4x, and
    /// re-scanned. Successfully decoded symbols have their coordinates
    /// mapped back to the original image frame.
    pub fn scan(&mut self, image: &mut Image) -> ScanResult {
        let (symbols, raw_region) = self.scanner.scan_image(image.as_mut_image());
        let finder_region = raw_region.map(|(x, y, w, h)| FinderRegion {
            x,
            y,
            width: w,
            height: h,
        });

        if !self.retry_undecoded_regions {
            return ScanResult::new(symbols, finder_region);
        }

        let Some(finder_region) = finder_region else {
            return ScanResult::new(symbols, finder_region);
        };

        // Try multiple scale factors: the adaptive binarization window size
        // is chosen in power-of-2 steps based on image size, and certain
        // intermediate sizes land in a range where the window extends too
        // far into the white quiet zone, causing halo artifacts that break
        // data extraction. Trying 2x, 4x, and 6x covers the common cases.
        const SCALES: &[u32] = &[2, 4, 6];

        // Skip retry for regions that cover more than 10% of the image
        // area — but only for sufficiently large images. On small images
        // (< 200px) the QR code legitimately fills most of the frame.
        // On larger images, a big region is almost always a false positive
        // from a 1D barcode.
        let apply_area_filter = image.width() >= 200 && image.height() >= 200;
        if apply_area_filter {
            let image_area = image.width() as u64 * image.height() as u64;
            let region_area = finder_region.width as u64 * finder_region.height as u64;
            if region_area > image_area / 10 {
                return ScanResult::new(symbols, Some(finder_region));
            }
        }

        // Pad by 50% of the region size on each side for quiet zone
        let pad_x = finder_region.width / 2;
        let pad_y = finder_region.height / 2;
        let cx = finder_region.x.saturating_sub(pad_x);
        let cy = finder_region.y.saturating_sub(pad_y);
        let cw = (finder_region.width + 2 * pad_x).min(image.width().saturating_sub(cx));
        let ch = (finder_region.height + 2 * pad_y).min(image.height().saturating_sub(cy));

        let cropped = match image.crop(cx, cy, cw, ch) {
            Some(c) => c,
            None => return ScanResult::new(symbols, Some(finder_region)),
        };

        for &scale in SCALES {
            if let Some(mut upscaled) = cropped.upscale(scale) {
                let (mut retry_symbols, _) = self.scanner.scan_image(upscaled.as_mut_image());
                if !retry_symbols.is_empty() {
                    // Remap coordinates back to original image frame
                    let half = scale as i32 / 2;
                    for sym in &mut retry_symbols {
                        for pt in &mut sym.pts {
                            pt[0] = (pt[0] + half) / scale as i32 + cx as i32;
                            pt[1] = (pt[1] + half) / scale as i32 + cy as i32;
                        }
                    }
                    // Merge with existing symbols, skip duplicates
                    let mut all = symbols;
                    for sym in retry_symbols {
                        if !all
                            .iter()
                            .any(|s| s.symbol_type() == sym.symbol_type() && s.data == sym.data)
                        {
                            all.push(sym);
                        }
                    }
                    return ScanResult::new(all, None);
                }
            }
        }

        ScanResult::new(symbols, Some(finder_region))
    }
}

impl Default for Scanner {
    fn default() -> Self {
        Self::new()
    }
}
