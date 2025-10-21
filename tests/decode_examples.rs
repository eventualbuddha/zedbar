//! Integration tests for decoding example barcode images

use std::path::Path;
use zbar::{Image, Scanner};

/// Helper function to decode an image file and return the first symbol found
fn decode_image(path: &str) -> Option<(String, String)> {
    let path = Path::new(path);
    if !path.exists() {
        return None;
    }

    let img = image::open(path).ok()?.to_luma8();
    let mut scanner = Scanner::new();

    // Create a zbar Image from the image buffer
    let mut zbar_image = Image::from_gray(img.as_raw(), img.width(), img.height()).ok()?;

    // Scan the image
    scanner.scan(&mut zbar_image).ok()?;

    // Get the symbols
    let symbols = zbar_image.symbols();

    symbols.iter().next().map(|symbol| {
        let symbol_type = format!("{:?}", symbol.symbol_type());
        let data = String::from_utf8_lossy(symbol.data())
            .trim_end_matches('\0') // Remove null terminators
            .to_string();
        (symbol_type, data)
    })
}

#[test]
fn test_qr_simple() {
    let result = decode_image("examples/test-qr.png");
    assert_eq!(
        result,
        Some((
            "QrCode".to_string(),
            "Hello, simplified zbar!\n".to_string()
        ))
    );
}

#[test]
fn test_qr_jpg() {
    let result = decode_image("examples/test-qr.jpg");
    assert_eq!(
        result,
        Some((
            "QrCode".to_string(),
            "Hello, simplified zbar!\n".to_string()
        ))
    );
}

#[test]
fn test_qr_wifi_sharing() {
    let result = decode_image("examples/pixel-wifi-sharing-qr-code.png");
    assert_eq!(
        result,
        Some((
            "QrCode".to_string(),
            "WIFI:S:Not a real network;T:SAE;P:password;H:false;;".to_string()
        ))
    );
}

#[test]
fn test_qr_capstone_interference() {
    let result = decode_image("examples/qr-code-capstone-interference.png");
    assert_eq!(
        result,
        Some((
            "QrCode".to_string(),
            "http://txz.qq.com/p?k=T8sZMvS*JxhU0kQFseMOMQZAKuE7An3u&f=716027609".to_string()
        ))
    );
}

#[test]
fn test_qr_color_bands() {
    let result = decode_image("examples/qr-code-color-bands.png");
    assert!(result.is_some());
    let (symbol_type, data) = result.unwrap();
    assert_eq!(symbol_type, "QrCode");
    // This QR code contains Chinese text, just verify it starts correctly
    assert!(data.starts_with("二维码生成器"));
    assert!(data.contains("https://zh.qr-code-generator.com"));
}

#[test]
fn test_qr_low_contrast() {
    let result = decode_image("examples/qr-code-low-contrast.png");
    assert_eq!(
        result,
        Some((
            "QrCode".to_string(),
            "72f7f23bf1e7428ccaba5366a938f420".to_string()
        ))
    );
}

#[test]
fn test_qr_pacman() {
    let result = decode_image("examples/qr-code-pacman.png");
    assert_eq!(
        result,
        Some((
            "QrCode".to_string(),
            "http://weixin.qq.com/r/gkgQCDrEc2cMrX5r9x2Q".to_string()
        ))
    );
}

#[test]
fn test_codabar() {
    let result = decode_image("examples/test-codabar.png");
    assert_eq!(result, Some(("Codabar".to_string(), "A40156B".to_string())));
}

#[test]
fn test_code128() {
    let result = decode_image("examples/test-code128.png");
    assert_eq!(
        result,
        Some(("Code128".to_string(), "HELLO123".to_string())) // Null terminator trimmed
    );
}

#[test]
fn test_code39() {
    let result = decode_image("examples/test-code39.png");
    assert_eq!(result, Some(("Code39".to_string(), "TEST123".to_string())));
}

#[test]
fn test_code93() {
    let result = decode_image("examples/test-code93.png");
    assert_eq!(result, Some(("Code93".to_string(), "CODE93".to_string())));
}

#[test]
fn test_ean13() {
    let result = decode_image("examples/test-ean13.png");
    assert_eq!(
        result,
        Some(("Ean13".to_string(), "5901234123457".to_string()))
    );
}

#[test]
fn test_ean8_decoded_as_ean13() {
    // Note: This image is decoded as EAN13, not EAN8 or UPCA
    let result = decode_image("examples/test-ean8.png");
    assert_eq!(
        result,
        Some(("Ean13".to_string(), "0000963850742".to_string()))
    );
}

#[test]
fn test_i25() {
    let result = decode_image("examples/test-i25.png");
    assert_eq!(result, Some(("I25".to_string(), "1234567890".to_string())));
}

#[test]
fn test_upca_decoded_as_ean13() {
    // Note: This image is decoded as EAN13, not UPCA
    let result = decode_image("examples/test-upca.png");
    assert_eq!(
        result,
        Some(("Ean13".to_string(), "0012345678905".to_string()))
    );
}

#[test]
fn test_rqrr_crash_2() {
    // Regression test for integer overflow in databar decoder
    // This image previously caused panic due to epoch overflow
    let result = decode_image("examples/rqrr-crash-2.jpeg");
    assert!(result.is_some());
    let (symbol_type, data) = result.unwrap();
    assert_eq!(symbol_type, "QrCode");
    assert!(data.contains("欢迎访问太平洋IT百科栏目"));
}

#[test]
fn test_all_examples_decode() {
    // Verify that all example images can be decoded without errors
    let images = vec![
        "examples/test-qr.png",
        "examples/test-qr.jpg",
        "examples/pixel-wifi-sharing-qr-code.png",
        "examples/qr-code-capstone-interference.png",
        "examples/qr-code-color-bands.png",
        "examples/qr-code-low-contrast.png",
        "examples/qr-code-pacman.png",
        "examples/test-codabar.png",
        "examples/test-code128.png",
        "examples/test-code39.png",
        "examples/test-code93.png",
        "examples/test-ean13.png",
        "examples/test-ean8.png",
        "examples/test-i25.png",
        "examples/test-upca.png",
        "examples/rqrr-crash-2.jpeg",
    ];

    for image_path in images {
        let result = decode_image(image_path);
        assert!(result.is_some(), "Failed to decode {}", image_path);
    }
}
