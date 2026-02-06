//! WebAssembly bindings for the zedbar barcode scanner.

use wasm_bindgen::prelude::*;

use crate::image::Image;
use crate::scanner::Scanner;

/// A decoded barcode/QR code result.
#[wasm_bindgen]
pub struct DecodeResult {
    symbol_type: String,
    data: Vec<u8>,
}

#[wasm_bindgen]
impl DecodeResult {
    /// The barcode format name (e.g. "QR-Code", "EAN-13").
    #[wasm_bindgen(getter)]
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

/// Scan grayscale image data for barcodes and QR codes.
///
/// `data` must be an array of 8-bit grayscale pixel values,
/// row-major, with dimensions `width` x `height`.
///
/// Returns an array of `DecodeResult` objects.
#[wasm_bindgen]
pub fn scan_grayscale(data: &[u8], width: u32, height: u32) -> Result<Vec<DecodeResult>, JsValue> {
    let mut image =
        Image::from_gray(data, width, height).map_err(|e| JsValue::from_str(&format!("{e:?}")))?;

    let mut scanner = Scanner::new();
    let symbols = scanner.scan(&mut image);

    Ok(symbols
        .into_iter()
        .map(|s| DecodeResult {
            symbol_type: s.symbol_type().to_string(),
            data: s.data().to_vec(),
        })
        .collect())
}
