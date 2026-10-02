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
//! // Create a scanner with every supported symbology enabled.
//! let mut scanner = Scanner::new();
//!
//! // Or build a scanner with only the symbologies you need.
//! use zedbar::{DecoderConfig, config::*};
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
use crate::symbol::{Bounds, Point, Symbol, SymbolType};

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
/// for region in result.finder_regions() {
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
    finder_regions: Vec<FinderRegion>,
}

impl ScanResult {
    pub(crate) fn new(symbols: Vec<Symbol>, finder_regions: Vec<FinderRegion>) -> Self {
        Self {
            symbols,
            finder_regions,
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

    /// Regions where QR finder patterns were detected but no QR code
    /// was successfully decoded.
    ///
    /// Each entry is a separate cluster of finder patterns — typically
    /// one per undecoded QR code in the image. Cropping and upscaling
    /// each region may yield successful decodes.
    ///
    /// The count is capped: an image dense with finder-like patterns can
    /// produce candidates faster than any caller could usefully act on them,
    /// so only the first several dozen are reported.
    ///
    /// Empty when no undecoded regions were found, or when the `qrcode`
    /// feature is disabled.
    pub fn finder_regions(&self) -> &[FinderRegion] {
        &self.finder_regions
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
    retries: Vec<Retry>,
}

/// A way of deriving an image to re-scan when the full-resolution pass
/// leaves a QR code undecoded. Listed cheapest first, which is the order
/// they run in.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Retry {
    /// Crop and upscale each undecoded finder region.
    UndecodedRegions,
    /// Box-filter a large image to half and then quarter size.
    Downscaled,
    /// Smooth the whole image with a 3x3 Gaussian.
    Smoothed,
}

impl Retry {
    fn enabled(config: &DecoderConfig) -> Vec<Retry> {
        let mut retries = Vec::new();
        if config.retry_undecoded_regions {
            retries.push(Retry::UndecodedRegions);
        }
        if config.retry_downscaled {
            retries.push(Retry::Downscaled);
        }
        if config.retry_smoothed {
            retries.push(Retry::Smoothed);
        }
        retries
    }
}

/// Carries a point in a derived image back to the original frame as
/// `p * num / den + offset`, rounded to nearest.
#[derive(Debug, Clone, Copy)]
struct Mapping {
    num: i32,
    den: i32,
    dx: i32,
    dy: i32,
}

impl Mapping {
    const IDENTITY: Self = Self {
        num: 1,
        den: 1,
        dx: 0,
        dy: 0,
    };

    /// A crop whose top-left corner is at `(x, y)`, upscaled by `scale`.
    fn upscaled_crop(scale: u32, x: u32, y: u32) -> Self {
        Self {
            num: 1,
            den: scale as i32,
            dx: x as i32,
            dy: y as i32,
        }
    }

    /// The image downscaled by `scale`, with output pixel `(0, 0)` centered
    /// on source pixel `(offset, offset)`.
    fn downscaled(scale: i32, offset: i32) -> Self {
        Self {
            num: scale,
            den: 1,
            dx: offset,
            dy: offset,
        }
    }

    fn apply(&self, pt: &mut Point) {
        pt.x = (pt.x * self.num + self.den / 2) / self.den + self.dx;
        pt.y = (pt.y * self.num + self.den / 2) / self.den + self.dy;
    }
}

impl Scanner {
    /// Create a new image scanner with every supported symbology enabled.
    ///
    /// Equivalent to `Scanner::with_config(DecoderConfig::all())`. Convenient
    /// for exploratory use; for production use, prefer
    /// [`Scanner::with_config()`] with a
    /// [`DecoderConfig::new()`](DecoderConfig::new) that opts into only the
    /// symbologies you actually need.
    pub fn new() -> Self {
        Self::with_config(DecoderConfig::all())
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
        let retries = Retry::enabled(&config);
        Self {
            scanner: ImageScanner::with_config(config),
            retries,
        }
    }

    /// Scan an image for barcodes
    ///
    /// Returns a [`ScanResult`] containing decoded symbols and any
    /// undecoded QR finder regions. Check [`ScanResult::finder_regions()`]
    /// to find areas that may contain QR codes too small to decode at
    /// the current resolution.
    ///
    /// When [`DecoderConfig::retry_undecoded_regions`] is enabled, each
    /// undecoded finder region is automatically cropped, upscaled, and
    /// re-scanned. Only QR and SQ codes are taken from those re-scans — the
    /// regions come from QR finder patterns, and a linear symbology gains
    /// nothing from upscaling a fragment it already saw at full resolution.
    /// Recovered symbols have their coordinates mapped back to the original
    /// image frame.
    ///
    /// When [`DecoderConfig::retry_downscaled`] is enabled and no QR or SQ
    /// code has been found by then, a large image is box-filtered to half
    /// and then quarter size and re-scanned, which recovers codes in photos
    /// of screens whose pixel grid defeats finder detection at full
    /// resolution. The same symbol filter and coordinate mapping apply.
    ///
    /// When [`DecoderConfig::retry_smoothed`] is enabled, still no QR or SQ
    /// code has been found, and an undecoded finder region remains, the
    /// image is smoothed with a 3x3 Gaussian and re-scanned once, which
    /// closes the halftone gaps that keep a scanned print from decoding. The
    /// same symbol filter applies.
    pub fn scan(&mut self, image: &mut Image) -> ScanResult {
        let (mut symbols, raw_regions) = self.scanner.scan_image(image.as_mut_image());
        let mut finder_regions: Vec<FinderRegion> = raw_regions
            .into_iter()
            .map(|(x, y, w, h)| FinderRegion {
                x,
                y,
                width: w,
                height: h,
            })
            .collect();

        for i in 0..self.retries.len() {
            match self.retries[i] {
                // A leftover region may point at a second, smaller code next
                // to one already decoded, so this runs regardless.
                Retry::UndecodedRegions => {
                    if !finder_regions.is_empty() {
                        let regions = std::mem::take(&mut finder_regions);
                        finder_regions = self.retry_regions(image, &mut symbols, regions);
                    }
                }
                Retry::Downscaled => self.retry_whole_image(
                    image,
                    &mut symbols,
                    &mut finder_regions,
                    Self::retry_downscaled,
                ),
                // Smoothing helps a code that was located but not read, so
                // an image with no finder patterns at all skips the pass.
                Retry::Smoothed if !finder_regions.is_empty() => self.retry_whole_image(
                    image,
                    &mut symbols,
                    &mut finder_regions,
                    Self::retry_smoothed,
                ),
                Retry::Smoothed => {}
            }
        }

        ScanResult::new(symbols, finder_regions)
    }

    /// Scan a derived image and keep the QR and SQ codes it yields, with
    /// their coordinates carried back to the original frame. Returns whether
    /// it yielded any.
    ///
    /// Only 2D symbols are kept: every retry is driven by QR finder patterns,
    /// and a derived image carries no detail the linear decoders did not
    /// already have at full resolution — only more scan lines across a
    /// fragment of the image, which is how a short read happens. Interleaved
    /// 2 of 5 is the clearest case: any even-length substring of one is
    /// itself a valid symbol.
    fn rescan(&mut self, derived: &mut Image, mapping: Mapping, symbols: &mut Vec<Symbol>) -> bool {
        let (found, _) = self.scanner.scan_image(derived.as_mut_image());
        let mut found: Vec<Symbol> = found
            .into_iter()
            .filter(|s| is_2d(s.symbol_type()))
            .collect();
        if found.is_empty() {
            return false;
        }
        for sym in &mut found {
            for pt in &mut sym.pts {
                mapping.apply(pt);
            }
        }
        merge_symbols(symbols, found);
        true
    }

    /// Run a whole-image retry, which exists to find a code the earlier
    /// passes missed entirely, so one already in hand makes it redundant.
    /// A region the full-resolution pass could not decode is resolved once
    /// a recovered code covers it; the rest stay undecoded.
    fn retry_whole_image(
        &mut self,
        image: &Image,
        symbols: &mut Vec<Symbol>,
        finder_regions: &mut Vec<FinderRegion>,
        retry: fn(&mut Self, &Image, &mut Vec<Symbol>),
    ) {
        if symbols.iter().any(|s| is_2d(s.symbol_type())) {
            return;
        }
        let before = symbols.len();
        retry(self, image, symbols);
        let recovered: Vec<Bounds> = symbols[before..]
            .iter()
            .filter_map(Symbol::bounds)
            .collect();
        finder_regions.retain(|region| !recovered.iter().any(|b| intersects(region, b)));
    }

    /// Crop, upscale and re-scan each undecoded finder region. Returns the
    /// regions that still did not decode.
    fn retry_regions(
        &mut self,
        image: &Image,
        symbols: &mut Vec<Symbol>,
        finder_regions: Vec<FinderRegion>,
    ) -> Vec<FinderRegion> {
        // Try multiple scale factors: the adaptive binarization window size
        // is chosen in power-of-2 steps based on image size, and certain
        // intermediate sizes land in a range where the window extends too
        // far into the white quiet zone, causing halo artifacts that break
        // data extraction. Trying 2x, 4x, and 6x covers the common cases.
        const SCALES: &[u32] = &[2, 4, 6];

        // Skip retry for regions that cover more than 10% of the image
        // area — on large images a big region is almost always a false
        // positive from a 1D barcode. Small images (< 200px on either
        // side) are exempt because a legitimate single QR often fills
        // most of the frame.
        let apply_area_filter = image.width() >= 200 && image.height() >= 200;
        let image_area = image.width() as u64 * image.height() as u64;
        let area_limit = image_area / 10;

        // Each retried region costs up to one full rescan per scale, on an
        // upscaled crop. A cluttered image can report dozens of candidates, so
        // cap how many are actually retried; the rest are handed back
        // unresolved for the caller to deal with as it sees fit.
        //
        // Taking the first N is not arbitrary. Finder centers are sorted by
        // (bucketed) edge-point count before the triplet search, so the
        // per-triplet candidates — which are the ones reported whenever any
        // triplet looked like a QR — arrive in descending order of confidence.
        // The cluster-derived fallback is only spatially ordered, but it is
        // used solely when no triplet survived at all.
        const MAX_RETRIED_REGIONS: usize = 16;

        let mut unresolved: Vec<FinderRegion> = Vec::new();
        let mut retried = 0usize;
        for region in &finder_regions {
            if retried >= MAX_RETRIED_REGIONS {
                unresolved.push(*region);
                continue;
            }
            if apply_area_filter {
                let region_area = region.width as u64 * region.height as u64;
                if region_area > area_limit {
                    unresolved.push(*region);
                    continue;
                }
            }
            // Pad by 50% of the region size on each side for quiet zone
            let pad_x = region.width / 2;
            let pad_y = region.height / 2;
            let cx = region.x.saturating_sub(pad_x);
            let cy = region.y.saturating_sub(pad_y);
            let cw = (region.width + 2 * pad_x).min(image.width().saturating_sub(cx));
            let ch = (region.height + 2 * pad_y).min(image.height().saturating_sub(cy));

            let Some(cropped) = image.crop(cx, cy, cw, ch) else {
                unresolved.push(*region);
                continue;
            };

            retried += 1;
            let decoded = SCALES.iter().any(|&scale| {
                cropped.upscale(scale).is_some_and(|mut upscaled| {
                    self.rescan(
                        &mut upscaled,
                        Mapping::upscaled_crop(scale, cx, cy),
                        symbols,
                    )
                })
            });
            if !decoded {
                unresolved.push(*region);
            }
        }

        unresolved
    }

    /// Re-scan the image at half and quarter size, keeping any QR or SQ
    /// codes found. Stops at the first size that decodes one.
    fn retry_downscaled(&mut self, image: &Image, symbols: &mut Vec<Symbol>) {
        // A photographed QR module spans several pixels, so halving a large
        // image keeps the code decodable while the screen's pixel grid
        // averages out. Below this size halving starts costing codes the
        // full-resolution pass could read, and the pixel grid of a screen is
        // no longer resolved anyway. Checked before every step, so the
        // quarter-size pass needs an image about twice this size.
        const MIN_SIDE: u32 = 1024;
        const MAX_STEPS: u32 = 2;

        // Output pixel (x, y) of `downscale(2)` is centered on source pixel
        // (2x + 2, 2y + 2). Compose that over the steps taken so far.
        let mut scale = 1i32;
        let mut offset = 0i32;
        let mut current = None;

        for _ in 0..MAX_STEPS {
            let source = current.as_ref().unwrap_or(image);
            if source.width().min(source.height()) < MIN_SIDE {
                return;
            }
            let Some(mut smaller) = source.downscale(2) else {
                return;
            };
            offset += 2 * scale;
            scale *= 2;

            if self.rescan(&mut smaller, Mapping::downscaled(scale, offset), symbols) {
                return;
            }
            current = Some(smaller);
        }
    }

    /// Re-scan a smoothed copy of the image, keeping any QR or SQ codes
    /// found. Smoothing keeps the image size, so coordinates carry over.
    fn retry_smoothed(&mut self, image: &Image, symbols: &mut Vec<Symbol>) {
        self.rescan(&mut image.smooth(), Mapping::IDENTITY, symbols);
    }
}

/// Whether a finder region and a symbol's bounding box share any area.
fn intersects(region: &FinderRegion, bounds: &Bounds) -> bool {
    let (rx0, ry0) = (region.x as i64, region.y as i64);
    let (rx1, ry1) = (rx0 + region.width as i64, ry0 + region.height as i64);
    let (bx0, by0) = (bounds.x as i64, bounds.y as i64);
    let (bx1, by1) = (bx0 + bounds.width as i64, by0 + bounds.height as i64);
    rx0 < bx1 && bx0 < rx1 && ry0 < by1 && by0 < ry1
}

/// The symbologies a QR-driven retry is allowed to contribute.
fn is_2d(sym: SymbolType) -> bool {
    matches!(sym, SymbolType::QrCode | SymbolType::SqCode)
}

/// Append the retry's symbols that the earlier passes did not already find.
fn merge_symbols(symbols: &mut Vec<Symbol>, retry_symbols: Vec<Symbol>) {
    for sym in retry_symbols {
        if !symbols
            .iter()
            .any(|s| s.symbol_type() == sym.symbol_type() && s.data == sym.data)
        {
            symbols.push(sym);
        }
    }
}

impl Default for Scanner {
    fn default() -> Self {
        Self::new()
    }
}
