//! WebAssembly bindings for the zedbar barcode scanner.

use wasm_bindgen::prelude::*;

use crate::SymbolType;
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
/// // Use defaults (every supported symbology, retry enabled):
/// const results = scanGrayscale(data, width, height);
///
/// // Disable automatic retry of small QR codes:
/// const results = scanGrayscale(data, width, height, { retryUndecodedRegions: false });
///
/// // Restrict to QR codes only:
/// const results = scanGrayscale(data, width, height, { symbologies: ["QR-Code"] });
/// ```
#[wasm_bindgen]
#[derive(Default)]
pub struct ScanOptions {
    retry_undecoded_regions: Option<bool>,
    symbologies: Option<Vec<String>>,
}

#[wasm_bindgen]
impl ScanOptions {
    /// Create a new options object with defaults.
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        Self::default()
    }

    /// Whether to automatically retry undecoded QR finder regions by
    /// cropping and upscaling them. Default: `true`.
    #[wasm_bindgen(setter, js_name = "retryUndecodedRegions")]
    pub fn set_retry_undecoded_regions(&mut self, value: bool) {
        self.retry_undecoded_regions = Some(value);
    }

    /// Restrict scanning to the listed symbologies. When omitted, every
    /// supported symbology is enabled.
    ///
    /// Names match the `symbolType` field on decode results: `"QR-Code"`,
    /// `"SQ-Code"`, `"EAN-13"`, `"EAN-8"`, `"EAN-2"`, `"EAN-5"`,
    /// `"UPC-A"`, `"UPC-E"`, `"ISBN-10"`, `"ISBN-13"`, `"I2/5"`,
    /// `"DataBar"`, `"DataBar-Exp"`, `"Codabar"`, `"CODE-39"`,
    /// `"CODE-93"`, `"CODE-128"`.
    #[wasm_bindgen(setter, js_name = "symbologies")]
    pub fn set_symbologies(&mut self, value: Vec<String>) {
        self.symbologies = Some(value);
    }
}

// Non-wasm helpers
impl ScanOptions {
    const DEFAULT_RETRY: bool = true;

    fn retry(&self) -> bool {
        self.retry_undecoded_regions.unwrap_or(Self::DEFAULT_RETRY)
    }
}

/// A decoded barcode/QR code result.
#[wasm_bindgen]
pub struct DecodeResult {
    symbol_type: String,
    data: Vec<u8>,
    text: Option<String>,
}

#[wasm_bindgen]
impl DecodeResult {
    /// The barcode format name (e.g. "QR-Code", "EAN-13").
    #[wasm_bindgen(getter, js_name = "symbolType")]
    pub fn symbol_type(&self) -> String {
        self.symbol_type.clone()
    }

    /// Raw decoded bytes.
    ///
    /// For 2D codes (QR, SQ), these are the original bytes as encoded
    /// in the barcode, before any text encoding or base64 conversion.
    #[wasm_bindgen(getter)]
    pub fn data(&self) -> Vec<u8> {
        self.data.clone()
    }

    /// Decoded data as text, or null if not decodable as text.
    ///
    /// For 2D codes, the text is produced by detecting the encoding (UTF-8,
    /// Shift-JIS, Windows-1252, etc.) and converting to a UTF-8 string.
    /// For linear barcodes, the data is always ASCII text.
    #[wasm_bindgen(getter)]
    pub fn text(&self) -> Option<String> {
        self.text.clone()
    }
}

fn build_scanner(options: Option<ScanOptions>) -> Result<Scanner, JsValue> {
    let options = options.unwrap_or_default();
    let mut config = match &options.symbologies {
        Some(syms) => {
            let mut cfg = DecoderConfig::new();
            for name in syms {
                let sym: SymbolType =
                    name.parse()
                        .map_err(|e: crate::symbol::ParseSymbolTypeError| {
                            JsValue::from_str(&e.to_string())
                        })?;
                cfg = cfg.enable_type(sym);
            }
            cfg
        }
        None => DecoderConfig::all(),
    };
    config = config.retry_undecoded_regions(options.retry());
    Ok(Scanner::with_config(config))
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

    let mut scanner = build_scanner(options)?;
    let symbols = scanner.scan(&mut image);

    Ok(symbols
        .into_iter()
        .map(|s| DecodeResult {
            symbol_type: s.symbol_type().to_string(),
            data: s.raw_data().unwrap_or(s.data()).to_vec(),
            text: s.data_string().map(|t| t.to_string()),
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
