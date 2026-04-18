//! WebAssembly bindings for the zedbar barcode scanner.

use wasm_bindgen::prelude::*;

use crate::config::DecoderConfig;
use crate::image::Image;
use crate::scanner::Scanner;

#[cfg(feature = "image")]
use image;

/// Options for barcode scanning.
///
/// All fields are optional and default to sensible values when omitted.
///
/// ```js
/// // Use defaults (retry enabled):
/// const results = scanGrayscale(data, width, height);
///
/// // Disable automatic retry of small QR codes:
/// const results = scanGrayscale(data, width, height, { retryUndecodedRegions: false });
/// ```
#[wasm_bindgen]
#[derive(Default)]
pub struct ScanOptions {
    retry_undecoded_regions: Option<bool>,
}

#[wasm_bindgen]
impl ScanOptions {
    /// Create a new options object with defaults.
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        Self {
            retry_undecoded_regions: None,
        }
    }

    /// Whether to automatically retry undecoded QR finder regions by
    /// cropping and upscaling them. Default: `true`.
    #[wasm_bindgen(setter, js_name = "retryUndecodedRegions")]
    pub fn set_retry_undecoded_regions(&mut self, value: bool) {
        self.retry_undecoded_regions = Some(value);
    }
}

/// A decoded barcode/QR code result.
#[wasm_bindgen]
pub struct DecodeResult {
    symbol_type: String,
    data: Vec<u8>,
}

#[wasm_bindgen]
impl DecodeResult {
    /// The barcode format name (e.g. "QR-Code", "EAN-13").
    #[wasm_bindgen(getter, js_name = "symbolType")]
    pub fn symbol_type(&self) -> String {
        self.symbol_type.clone()
    }

    /// Raw decoded bytes.
    #[wasm_bindgen(getter)]
    pub fn data(&self) -> Vec<u8> {
        self.data.clone()
    }

    /// Decoded data as UTF-8 text, or null if not valid UTF-8.
    #[wasm_bindgen(getter)]
    pub fn text(&self) -> Option<String> {
        std::str::from_utf8(&self.data).ok().map(|s| s.to_string())
    }
}

fn build_scanner(options: Option<ScanOptions>) -> Scanner {
    let retry = options
        .and_then(|o| o.retry_undecoded_regions)
        .unwrap_or(true);

    let config = DecoderConfig::new().retry_undecoded_regions(retry);
    Scanner::with_config(config)
}

/// Scan grayscale image data for barcodes and QR codes.
///
/// `data` must be an array of 8-bit grayscale pixel values,
/// row-major, with dimensions `width` x `height`.
///
/// Returns an array of `DecodeResult` objects.
#[wasm_bindgen(js_name = "scanGrayscale")]
pub fn scan_grayscale(
    data: &[u8],
    width: u32,
    height: u32,
    options: Option<ScanOptions>,
) -> Result<Vec<DecodeResult>, JsValue> {
    let mut image =
        Image::from_gray(data, width, height).map_err(|e| JsValue::from_str(&format!("{e:?}")))?;

    let mut scanner = build_scanner(options);
    let symbols = scanner.scan(&mut image);

    Ok(symbols
        .into_iter()
        .map(|s| DecodeResult {
            symbol_type: s.symbol_type().to_string(),
            data: s.data().to_vec(),
        })
        .collect())
}

/// Scan an encoded image (PNG, JPEG, BMP, WebP) for barcodes and QR codes.
///
/// `bytes` should contain the raw bytes of an image file in one of the
/// supported formats: PNG, JPEG, BMP, or WebP.
///
/// The image will be automatically decoded and converted to grayscale
/// before scanning.
///
/// Returns an array of `DecodeResult` objects.
#[wasm_bindgen(js_name = "scanImageBytes")]
pub fn scan_image_bytes(
    bytes: &[u8],
    options: Option<ScanOptions>,
) -> Result<Vec<DecodeResult>, JsValue> {
    // Decode the image using the image crate
    let img = image::load_from_memory(bytes)
        .map_err(|e| JsValue::from_str(&format!("Failed to decode image: {}", e)))?;

    // Convert to grayscale
    let gray = img.to_luma8();
    let (width, height) = gray.dimensions();
    let data = gray.as_raw();

    scan_grayscale(data, width, height, options)
}
