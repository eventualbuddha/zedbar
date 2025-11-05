//! Property-based tests for QR code decoding
//!
//! These tests generate random QR code data, encode it to an image,
//! and verify that the decoder can correctly decode it back.
//! Also verifies that zbar matches rqrr's decoding results.

use crate::{Image, Scanner};
use image::{GrayImage, Luma};
use proptest::prelude::*;
use proptest::test_runner::TestCaseError;
use qrcode::QrCode;

/// Generate a QR code image from data
fn generate_qr_image(data: &[u8]) -> Option<GrayImage> {
    // Create QR code
    let code = QrCode::new(data).ok()?;

    // Render to image with some padding for better detection
    let image = code.render::<Luma<u8>>().quiet_zone(true).build();

    Some(image)
}

/// Helper to decode a QR code image with rqrr and return the decoded data
fn decode_qr_image_rqrr(img: &GrayImage) -> Result<Vec<Vec<u8>>, String> {
    let width = img.width() as usize;
    let height = img.height() as usize;
    let raw = img.as_raw();

    let mut prepared_img =
        rqrr::PreparedImage::prepare_from_greyscale(width, height, |x, y| raw[y * width + x]);

    let grids = prepared_img.detect_grids();
    if grids.is_empty() {
        return Err("rqrr: No grids found".to_string());
    }

    let mut results = Vec::new();
    for grid in grids {
        let (_meta, content) = grid
            .decode()
            .map_err(|e| format!("rqrr: Failed to decode grid: {e:?}"))?;
        results.push(content.into_bytes());
    }

    Ok(results)
}

/// Helper to decode a QR code image and return the decoded data
/// Also verifies that zbar and rqrr produce the same results
fn decode_qr_image(img: &GrayImage) -> Result<Vec<Vec<u8>>, String> {
    let (width, height) = img.dimensions();
    let data = img.as_raw();

    // Create ZBar image from grayscale data
    let mut zbar_img = Image::from_gray(data, width, height)
        .map_err(|e| format!("Failed to create ZBar image: {e:?}"))?;

    // Create scanner (QR codes enabled by default)
    let mut scanner = Scanner::new();

    // Scan the image
    let symbols = scanner.scan(&mut zbar_img);

    // Get symbols
    if symbols.is_empty() {
        return Err("No symbols found".to_string());
    }

    // Collect all decoded data
    let zbar_results: Vec<Vec<u8>> = symbols.iter().map(|s| s.data().to_vec()).collect();

    // Also decode with rqrr and verify they match
    let rqrr_results = decode_qr_image_rqrr(img)?;

    // Sort both results for comparison (order may differ)
    let mut zbar_sorted = zbar_results.clone();
    let mut rqrr_sorted = rqrr_results.clone();
    zbar_sorted.sort();
    rqrr_sorted.sort();

    if zbar_sorted != rqrr_sorted {
        return Err(format!(
            "zbar and rqrr produced different results!\nzbar: {} symbols\nrqrr: {} symbols",
            zbar_sorted.len(),
            rqrr_sorted.len()
        ));
    }

    Ok(zbar_results)
}

proptest! {
    /// Test that ASCII text QR codes can be roundtripped
    #[test]
    fn prop_qr_roundtrip_ascii(data in "[a-zA-Z0-9 ]{1,100}") {
        let data_bytes = data.as_bytes();

        // Generate QR code image
        let img = generate_qr_image(data_bytes)
            .ok_or_else(|| TestCaseError::fail("Failed to generate QR code"))?;

        // Decode it
        let decoded = decode_qr_image(&img)
            .map_err(TestCaseError::fail)?;

        // Verify we got exactly one symbol
        prop_assert_eq!(decoded.len(), 1, "Expected exactly one decoded symbol");

        // Verify the data matches
        prop_assert_eq!(&decoded[0], data_bytes, "Decoded data doesn't match original");
    }

    /// Test that numeric strings can be roundtripped
    #[test]
    fn prop_qr_roundtrip_numeric(data in "[0-9]{1,100}") {
        let data_bytes = data.as_bytes();

        let img = generate_qr_image(data_bytes)
            .ok_or_else(|| TestCaseError::fail("Failed to generate QR code"))?;

        let decoded = decode_qr_image(&img)
            .map_err(TestCaseError::fail)?;

        prop_assert_eq!(decoded.len(), 1);
        prop_assert_eq!(&decoded[0], data_bytes);
    }

    /// Test that ASCII-range binary data can be roundtripped
    /// Note: QR codes typically encode data as text, so we limit to ASCII range
    #[test]
    fn prop_qr_roundtrip_binary(data in prop::collection::vec(0u8..128, 1..50)) {
        let img = match generate_qr_image(&data) {
            Some(img) => img,
            None => return Ok(()), // Skip if QR generation fails (data too large, etc)
        };

        let decoded = decode_qr_image(&img)
            .map_err(TestCaseError::fail)?;

        prop_assert_eq!(decoded.len(), 1);
        prop_assert_eq!(&decoded[0], &data);
    }

    /// Test URL patterns
    #[test]
    fn prop_qr_roundtrip_urls(
        protocol in "(https?|ftp)",
        domain in "[a-z]{3,20}",
        tld in "(com|org|net|edu)",
        path in prop::option::of("[a-z0-9/]{0,30}")
    ) {
        let url = match path {
            Some(p) => format!("{protocol}://{domain}.{tld}/{p}"),
            None => format!("{protocol}://{domain}.{tld}"),
        };
        let data_bytes = url.as_bytes();

        let img = generate_qr_image(data_bytes)
            .ok_or_else(|| TestCaseError::fail("Failed to generate QR code"))?;

        let decoded = decode_qr_image(&img)
            .map_err(TestCaseError::fail)?;

        prop_assert_eq!(decoded.len(), 1);
        prop_assert_eq!(&decoded[0], data_bytes);
    }

    /// Test that empty data is handled correctly
    #[test]
    fn prop_qr_empty_or_single_byte(data in prop::collection::vec(0u8..128, 0..2)) {
        if data.is_empty() {
            // Empty QR codes are valid (version 1, minimal)
            if let Some(img) = generate_qr_image(&data) {
                if let Ok(decoded) = decode_qr_image(&img) {
                    prop_assert_eq!(decoded.len(), 1);
                    prop_assert_eq!(&decoded[0], &data);
                }
            }
        } else {
            let img = generate_qr_image(&data)
                .ok_or_else(|| TestCaseError::fail("Failed to generate QR code"))?;

            let decoded = decode_qr_image(&img)
                .map_err(TestCaseError::fail)?;

            prop_assert_eq!(decoded.len(), 1);
            prop_assert_eq!(&decoded[0], &data);
        }
    }
}

#[cfg(test)]
mod unit_tests {
    use super::*;

    #[test]
    fn test_simple_qr_roundtrip() {
        let data = b"Hello, World!";
        let img = generate_qr_image(data).expect("Failed to generate QR code");
        let decoded = decode_qr_image(&img).expect("Failed to decode QR code");

        assert_eq!(decoded.len(), 1);
        assert_eq!(&decoded[0], data);
    }

    #[test]
    fn test_numeric_qr_roundtrip() {
        let data = b"123456789";
        let img = generate_qr_image(data).expect("Failed to generate QR code");
        let decoded = decode_qr_image(&img).expect("Failed to decode QR code");

        assert_eq!(decoded.len(), 1);
        assert_eq!(&decoded[0], data);
    }

    #[test]
    fn test_url_qr_roundtrip() {
        let data = b"https://example.com/path/to/resource";
        let img = generate_qr_image(data).expect("Failed to generate QR code");
        let decoded = decode_qr_image(&img).expect("Failed to decode QR code");

        assert_eq!(decoded.len(), 1);
        assert_eq!(&decoded[0], data);
    }

    #[test]
    fn test_binary_qr_roundtrip() {
        // Use ASCII-range bytes for compatibility with QR code text encoding
        let data = vec![0u8, 1, 2, 3, 65, 66, 67]; // 65, 66, 67 = 'A', 'B', 'C'
        let img = generate_qr_image(&data).expect("Failed to generate QR code");
        let decoded = decode_qr_image(&img).expect("Failed to decode QR code");

        assert_eq!(decoded.len(), 1);
        assert_eq!(&decoded[0], &data);
    }

    #[test]
    fn test_special_chars_qr_roundtrip() {
        let data = b"!@#$%^&*()_+-=[]{}|;:',.<>?/~`";
        let img = generate_qr_image(data).expect("Failed to generate QR code");
        let decoded = decode_qr_image(&img).expect("Failed to decode QR code");

        assert_eq!(decoded.len(), 1);
        assert_eq!(&decoded[0], data);
    }
}
