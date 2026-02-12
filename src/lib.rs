//! Zedbar Barcode Scanning Library
//!
//! A pure Rust barcode scanning library supporting multiple barcode formats including
//! QR codes, EAN, UPC, Code128, Code39, and more. Based on the zbar C library.
//!
//! # Quick Start
//!
//! ```no_run
//! use zedbar::{Image, Scanner};
//!
//! // Load and convert image to grayscale
//! let img = image::open("barcode.png").unwrap();
//! let gray = img.to_luma8();
//! let (width, height) = gray.dimensions();
//!
//! // Create image and scanner
//! let mut img = Image::from_gray(gray.as_raw(), width, height).unwrap();
//! let mut scanner = Scanner::new();
//!
//! // Scan for barcodes
//! let symbols = scanner.scan(&mut img);
//! for symbol in symbols {
//!     println!("{:?}: {}", symbol.symbol_type(), symbol.data_string().unwrap_or(""));
//! }
//! ```
//!
//! # Configuration
//!
//! Use the type-safe [`config`] module to customize decoder behavior:
//!
//! ```
//! use zedbar::config::*;
//! use zedbar::{DecoderConfig, Scanner};
//!
//! let config = DecoderConfig::new()
//!     .enable(QrCode)
//!     .enable(Ean13)
//!     .set_binary(QrCode, true)          // Preserve binary data in QR codes
//!     .set_length_limits(Code39, 4, 20)  // Code39 must be 4-20 chars
//!     .test_inverted(true)               // Try inverted image if no symbols found
//!     .scan_density(2, 2);               // Scan every 2nd line (faster)
//!
//! let mut scanner = Scanner::with_config(config);
//! ```
//!
//! # Supported Formats
//!
//! - **2D Codes**: QR Code, SQCode
//! - **Linear Codes**: EAN-13, EAN-8, UPC-A, UPC-E, ISBN-10, ISBN-13
//! - **Code Family**: Code 128, Code 93, Code 39, Codabar
//! - **Industrial**: Interleaved 2 of 5, DataBar (RSS)
//!
//! # Modules
//!
//! - [`scanner`] - Main scanner API
//! - [`image`] - Image handling
//! - [`symbol`] - Decoded barcode symbols
//! - [`config`] - Type-safe decoder configuration
//! - [`error`] - Error types

#![allow(clippy::missing_safety_doc)]
#![allow(non_camel_case_types)]

// Public modules
pub mod config;
pub mod error;
pub mod image;
pub mod scanner;
pub mod symbol;

// Internal modules
pub(crate) mod color;
pub(crate) mod decoder;
pub(crate) mod decoders;
#[cfg(feature = "qrcode")]
pub(crate) mod finder;
pub(crate) mod image_data;
pub(crate) mod img_scanner;
pub(crate) mod img_scanner_config;
#[cfg(feature = "qrcode")]
pub(crate) mod qrcode;
#[cfg(feature = "sqcode")]
pub(crate) mod sqcode;
#[cfg(feature = "wasm")]
pub mod wasm;

// Re-export main types
pub use config::DecoderConfig;
pub use error::{Error, Result};
pub use image::Image;
pub use scanner::Scanner;
pub use symbol::{Orientation, SymbolType};

#[cfg(all(test, feature = "qrcode"))]
mod proptest_qr;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_qr_decode_png() {
        // Load the test QR code image (using ::image to specify the external crate)
        let img = ::image::ImageReader::open("examples/test-qr.png")
            .expect("Failed to open image")
            .decode()
            .expect("Failed to decode image");

        // Convert to grayscale
        let gray = img.to_luma8();
        let (width, height) = gray.dimensions();
        let data = gray.as_raw();

        // Create image from grayscale data
        let mut img =
            Image::from_gray(data, width, height).expect("Failed to create zedbar image");

        // Create scanner (QR codes enabled by default)
        let mut scanner = Scanner::new();

        // Scan the image
        let symbols = scanner.scan(&mut img);
        assert!(!symbols.is_empty(), "Expected to find at least one QR code");

        let mut zedbar_results = Vec::new();
        for symbol in symbols {
            let symbol_type = symbol.symbol_type();
            let data = symbol.data();
            println!("Decoded {symbol_type:?}: {} bytes", data.len());

            assert_eq!(symbol_type, SymbolType::QrCode);
            assert!(!data.is_empty(), "QR code data should not be empty");
            zedbar_results.push(data.to_vec());
        }

        // Also decode with rqrr and verify they match
        let width_usize = width as usize;
        let height_usize = height as usize;
        let raw = gray.as_raw();
        let mut prepared_img =
            rqrr::PreparedImage::prepare_from_greyscale(width_usize, height_usize, |x, y| {
                raw[y * width_usize + x]
            });
        let grids = prepared_img.detect_grids();
        assert!(
            !grids.is_empty(),
            "rqrr: Expected to find at least one grid"
        );

        let mut rqrr_results = Vec::new();
        for grid in grids {
            let (_meta, content) = grid.decode().expect("rqrr: Failed to decode grid");
            rqrr_results.push(content.into_bytes());
        }

        // Sort both results for comparison
        zedbar_results.sort();
        rqrr_results.sort();

        assert_eq!(
            zedbar_results, rqrr_results,
            "zedbar and rqrr produced different results"
        );
        println!("✓ zbar and rqrr agree on {} symbols", zedbar_results.len());
    }

    #[test]
    fn test_qr_decode_jpg() {
        // Load the test QR code image (using ::image to specify the external crate)
        let img = ::image::ImageReader::open("examples/test-qr.jpg")
            .expect("Failed to open image")
            .decode()
            .expect("Failed to decode image");

        // Convert to grayscale
        let gray = img.to_luma8();
        let (width, height) = gray.dimensions();
        let data = gray.as_raw();

        // Create image from grayscale data
        let mut img =
            Image::from_gray(data, width, height).expect("Failed to create zedbar image");

        // Create scanner (QR codes enabled by default)
        let mut scanner = Scanner::new();

        // Scan the image
        let symbols = scanner.scan(&mut img);
        assert!(!symbols.is_empty(), "Expected to find at least one QR code");

        let mut zedbar_results = Vec::new();
        for symbol in symbols {
            let symbol_type = symbol.symbol_type();
            let data = symbol.data();
            println!("Decoded {symbol_type:?}: {} bytes", data.len());

            assert_eq!(symbol_type, SymbolType::QrCode);
            assert!(!data.is_empty(), "QR code data should not be empty");
            zedbar_results.push(data.to_vec());
        }

        // Also decode with rqrr and verify they match
        let width_usize = width as usize;
        let height_usize = height as usize;
        let raw = gray.as_raw();
        let mut prepared_img =
            rqrr::PreparedImage::prepare_from_greyscale(width_usize, height_usize, |x, y| {
                raw[y * width_usize + x]
            });
        let grids = prepared_img.detect_grids();
        assert!(
            !grids.is_empty(),
            "rqrr: Expected to find at least one grid"
        );

        let mut rqrr_results = Vec::new();
        for grid in grids {
            let (_meta, content) = grid.decode().expect("rqrr: Failed to decode grid");
            rqrr_results.push(content.into_bytes());
        }

        // Sort both results for comparison
        zedbar_results.sort();
        rqrr_results.sort();

        assert_eq!(
            zedbar_results, rqrr_results,
            "zedbar and rqrr produced different results"
        );
        println!("✓ zedbar and rqrr agree on {} symbols", zedbar_results.len());
    }

    #[test]
    fn test_qr_decode_inverted() {
        // Load the test QR code image and convert to grayscale
        let img = ::image::ImageReader::open("examples/test-qr.png")
            .expect("Failed to open image")
            .decode()
            .expect("Failed to decode image");
        let gray = img.to_luma8();
        let (width, height) = gray.dimensions();
        let original_raw = gray.as_raw();

        // Invert the grayscale pixels (swap black/white)
        let mut inverted = Vec::with_capacity(original_raw.len());
        inverted.extend(original_raw.iter().map(|&v| 255u8.saturating_sub(v)));

        // Build a zedbar image from the inverted data
        let mut img =
            Image::from_gray(&inverted, width, height).expect("Failed to create zedbar image");

        // Configure scanner for QR codes with inverted testing
        use crate::config::*;
        let config = DecoderConfig::new().enable(QrCode).test_inverted(true);
        let mut scanner = Scanner::with_config(config);

        // Scan the inverted image
        let symbols = scanner.scan(&mut img);
        assert!(
            !symbols.is_empty(),
            "Expected to find at least one QR code in inverted image"
        );

        let mut zedbar_results = Vec::new();
        for symbol in symbols {
            let symbol_type = symbol.symbol_type();
            let data = symbol.data();
            println!("Decoded {symbol_type:?} (inverted): {} bytes", data.len());

            assert_eq!(symbol_type, SymbolType::QrCode);
            assert!(!data.is_empty(), "QR code data should not be empty");
            zedbar_results.push(data.to_vec());
        }

        // Use rqrr on the original (non-inverted) grayscale image for the ground truth
        let width_usize = width as usize;
        let height_usize = height as usize;
        let mut prepared_img =
            rqrr::PreparedImage::prepare_from_greyscale(width_usize, height_usize, |x, y| {
                original_raw[y * width_usize + x]
            });
        let grids = prepared_img.detect_grids();
        assert!(
            !grids.is_empty(),
            "rqrr: Expected to find at least one grid in original image"
        );

        let mut rqrr_results = Vec::new();
        for grid in grids {
            let (_meta, content) = grid.decode().expect("rqrr: Failed to decode grid");
            rqrr_results.push(content.into_bytes());
        }

        // Sort and compare results
        zedbar_results.sort();
        rqrr_results.sort();

        assert_eq!(
            zedbar_results, rqrr_results,
            "zedbar and rqrr produced different results for inverted image"
        );
        println!(
            "✓ zedbar and rqrr agree on {} symbols for inverted image",
            zedbar_results.len()
        );
    }

    #[test]
    fn test_ean13_decode() {
        // EAN-13 test (13 digits, common retail barcode)
        let img = ::image::ImageReader::open("examples/test-ean13.png")
            .expect("Failed to open test-ean13.png")
            .decode()
            .expect("Failed to decode image");

        let gray = img.to_luma8();
        let (width, height) = gray.dimensions();
        let data = gray.as_raw();

        let mut img =
            Image::from_gray(data, width, height).expect("Failed to create zedbar image");

        // Create scanner (EAN-13 enabled by default)
        let mut scanner = Scanner::new();

        let symbols = scanner.scan(&mut img);
        assert!(
            !symbols.is_empty(),
            "Expected to find at least one EAN-13 barcode"
        );

        for symbol in symbols {
            assert_eq!(symbol.symbol_type(), SymbolType::Ean13);
            println!("Decoded EAN-13: {}", symbol.data_string().unwrap_or(""));
        }
    }

    #[test]
    fn test_ean8_decode() {
        // EAN-8 test (8 digits, compact retail barcode)
        // Note: zedbar may report EAN-8 as EAN-13 with zero padding
        let img = ::image::ImageReader::open("examples/test-ean8.png")
            .expect("Failed to open test-ean8.png")
            .decode()
            .expect("Failed to decode image");

        let gray = img.to_luma8();
        let (width, height) = gray.dimensions();
        let data = gray.as_raw();

        let mut img =
            Image::from_gray(data, width, height).expect("Failed to create zedbar image");

        // Create scanner (EAN-8 enabled by default)
        let mut scanner = Scanner::new();

        let symbols = scanner.scan(&mut img);
        assert!(
            !symbols.is_empty(),
            "Expected to find at least one EAN-8 barcode"
        );

        for symbol in symbols {
            let decoded_data = symbol.data_string().unwrap_or("");
            println!("Decoded as {:?}: {}", symbol.symbol_type(), decoded_data);

            // EAN-8 can be reported as EAN-13 with zero-padding
            assert!(
                symbol.symbol_type() == SymbolType::Ean8
                    || symbol.symbol_type() == SymbolType::Ean13,
                "Expected EAN-8 or EAN-13, got {:?}",
                symbol.symbol_type()
            );

            // Verify the data contains our EAN-8 digits
            assert!(decoded_data.contains("96385074") || decoded_data == "96385074");
        }
    }

    #[test]
    fn test_upca_decode() {
        // UPC-A test (12 digits, common in North America)
        let img = ::image::ImageReader::open("examples/test-upca.png")
            .expect("Failed to open test-upca.png")
            .decode()
            .expect("Failed to decode image");

        let gray = img.to_luma8();
        let (width, height) = gray.dimensions();
        let data = gray.as_raw();

        let mut img =
            Image::from_gray(data, width, height).expect("Failed to create zedbar image");

        // Create scanner with UPC-A explicitly enabled
        use crate::config::*;
        let config = DecoderConfig::new().enable(Upca);
        let mut scanner = Scanner::with_config(config);

        let symbols = scanner.scan(&mut img);
        assert!(
            !symbols.is_empty(),
            "Expected to find at least one UPC-A barcode"
        );

        for symbol in symbols {
            assert_eq!(symbol.symbol_type(), SymbolType::Upca);
            println!("Decoded UPC-A: {}", symbol.data_string().unwrap_or(""));
        }
    }

    #[test]
    fn test_code128_decode() {
        // Code128 test (high-density alphanumeric barcode)
        let img = ::image::ImageReader::open("examples/test-code128.png")
            .expect("Failed to open test-code128.png")
            .decode()
            .expect("Failed to decode image");

        let gray = img.to_luma8();
        let (width, height) = gray.dimensions();
        let data = gray.as_raw();

        let mut img =
            Image::from_gray(data, width, height).expect("Failed to create zedbar image");

        // Create scanner (Code128 enabled by default)
        let mut scanner = Scanner::new();

        let symbols = scanner.scan(&mut img);
        assert!(
            !symbols.is_empty(),
            "Expected to find at least one Code128 barcode"
        );

        for symbol in symbols {
            assert_eq!(symbol.symbol_type(), SymbolType::Code128);
            println!("Decoded Code128: {}", symbol.data_string().unwrap_or(""));
        }
    }

    #[test]
    fn test_code39_decode() {
        // Code39 test with uppercase alphanumeric data
        let img = ::image::ImageReader::open("examples/test-code39.png")
            .expect("Failed to open test-code39.png")
            .decode()
            .expect("Failed to decode image");

        let gray = img.to_luma8();
        let (width, height) = gray.dimensions();
        let data = gray.as_raw();

        let mut img =
            Image::from_gray(data, width, height).expect("Failed to create zedbar image");

        // Create scanner (Code39 enabled by default)
        let mut scanner = Scanner::new();

        let symbols = scanner.scan(&mut img);
        assert!(
            !symbols.is_empty(),
            "Expected to find at least one Code39 barcode"
        );

        for symbol in symbols {
            assert_eq!(symbol.symbol_type(), SymbolType::Code39);
            println!("Decoded Code39: {}", symbol.data_string().unwrap_or(""));
        }
    }

    #[test]
    fn test_code93_decode() {
        // Code93 test
        let img = ::image::ImageReader::open("examples/test-code93.png")
            .expect("Failed to open test-code93.png")
            .decode()
            .expect("Failed to decode image");

        let gray = img.to_luma8();
        let (width, height) = gray.dimensions();
        let data = gray.as_raw();

        let mut img =
            Image::from_gray(data, width, height).expect("Failed to create zedbar image");

        // Create scanner (Code93 enabled by default)
        let mut scanner = Scanner::new();

        let symbols = scanner.scan(&mut img);
        assert!(
            !symbols.is_empty(),
            "Expected to find at least one Code93 barcode"
        );

        for symbol in symbols {
            assert_eq!(symbol.symbol_type(), SymbolType::Code93);
            println!("Decoded Code93: {}", symbol.data_string().unwrap_or(""));
        }
    }

    #[test]
    fn test_codabar_decode() {
        // Codabar test (uses start/stop characters like A/B/C/D)
        let img = ::image::ImageReader::open("examples/test-codabar.png")
            .expect("Failed to open test-codabar.png")
            .decode()
            .expect("Failed to decode image");

        let gray = img.to_luma8();
        let (width, height) = gray.dimensions();
        let data = gray.as_raw();

        let mut img =
            Image::from_gray(data, width, height).expect("Failed to create zedbar image");

        // Create scanner (Codabar enabled by default)
        let mut scanner = Scanner::new();

        let symbols = scanner.scan(&mut img);
        assert!(
            !symbols.is_empty(),
            "Expected to find at least one Codabar barcode"
        );

        for symbol in symbols {
            assert_eq!(symbol.symbol_type(), SymbolType::Codabar);
            println!("Decoded Codabar: {}", symbol.data_string().unwrap_or(""));
        }
    }

    #[test]
    fn test_interleaved2of5_decode() {
        // Interleaved 2 of 5 test (numeric only, even number of digits)
        let img = ::image::ImageReader::open("examples/test-i25.png")
            .expect("Failed to open test-i25.png")
            .decode()
            .expect("Failed to decode image");

        let gray = img.to_luma8();
        let (width, height) = gray.dimensions();
        let data = gray.as_raw();

        let mut img =
            Image::from_gray(data, width, height).expect("Failed to create zedbar image");

        // Create scanner (I25 enabled by default)
        let mut scanner = Scanner::new();

        let symbols = scanner.scan(&mut img);
        assert!(
            !symbols.is_empty(),
            "Expected to find at least one I25 barcode"
        );

        for symbol in symbols {
            assert_eq!(symbol.symbol_type(), SymbolType::I25);
            println!("Decoded I25: {}", symbol.data_string().unwrap_or(""));
        }
    }

    #[test]
    fn test_pixel_wifi_qr_decode() {
        // Test for pixel-wifi-sharing-qr-code.png
        // This is a QR code that previously caused an overflow in the sqcode decoder
        let img = ::image::ImageReader::open("examples/pixel-wifi-sharing-qr-code.png")
            .expect("Failed to open pixel-wifi-sharing-qr-code.png")
            .decode()
            .expect("Failed to decode image");

        let gray = img.to_luma8();
        let (width, height) = gray.dimensions();
        let data = gray.as_raw();

        let mut img =
            Image::from_gray(data, width, height).expect("Failed to create zedbar image");

        // Create scanner (QR codes enabled by default)
        let mut scanner = Scanner::new();

        let symbols = scanner.scan(&mut img);
        assert!(
            !symbols.is_empty(),
            "Expected to find at least one QR code in pixel-wifi-sharing-qr-code.png"
        );

        let symbol = symbols.into_iter().next().unwrap();
        assert_eq!(symbol.symbol_type(), SymbolType::QrCode);
        let data = symbol.data_string().unwrap_or("");
        assert_eq!(data, "WIFI:S:Not a real network;T:SAE;P:password;H:false;;");
    }

    #[test]
    fn test_qr_code_capstone_interference() {
        let img = ::image::ImageReader::open("examples/qr-code-capstone-interference.png")
            .expect("Failed to open image")
            .decode()
            .expect("Failed to decode image");

        let gray = img.to_luma8();
        let (width, height) = gray.dimensions();
        let data = gray.as_raw();

        let mut img =
            Image::from_gray(data, width, height).expect("Failed to create zedbar image");

        // Create scanner (QR codes enabled by default)
        let mut scanner = Scanner::new();

        let symbols = scanner.scan(&mut img);
        assert!(!symbols.is_empty(), "Expected to find at least one QR code");

        let symbol = symbols.into_iter().next().unwrap();
        assert_eq!(symbol.symbol_type(), SymbolType::QrCode);
        let data = symbol.data_string().unwrap_or("");
        assert_eq!(
            data,
            "http://txz.qq.com/p?k=T8sZMvS*JxhU0kQFseMOMQZAKuE7An3u&f=716027609"
        );
    }

    #[test]
    #[ignore = "QR decoder cannot decode low-contrast color QR codes - finder patterns detected but data extraction fails"]
    fn test_qr_code_color_bands() {
        // This QR code has colored bands that result in very low contrast when converted to grayscale
        // (grayscale range 145-255 instead of 0-255). The QR decoder fails during data extraction
        // even though it successfully:
        // - Detects finder patterns (180 horizontal, 191 vertical lines)
        // - Locates the 3 finder centers
        // - Applies histogram stretching to improve contrast
        // The issue appears to be in qr_reader_match_centers or subsequent data decoding functions.
        let img = ::image::ImageReader::open("examples/qr-code-color-bands.png")
            .expect("Failed to open image")
            .decode()
            .expect("Failed to decode image");

        let gray = img.to_luma8();
        let (width, height) = gray.dimensions();
        let data = gray.as_raw();

        // Apply histogram stretching for low-contrast images
        // Find min/max values
        let min_val = *data.iter().min().unwrap_or(&0);
        let max_val = *data.iter().max().unwrap_or(&255);

        let normalized = if max_val > min_val {
            // Stretch the histogram to full 0-255 range
            data.iter()
                .map(|&pixel| {
                    let stretched =
                        ((pixel as u16 - min_val as u16) * 255) / (max_val as u16 - min_val as u16);
                    stretched as u8
                })
                .collect::<Vec<u8>>()
        } else {
            data.to_vec()
        };

        let mut img =
            Image::from_gray(&normalized, width, height).expect("Failed to create zedbar image");

        // Create scanner (QR codes enabled by default)
        let mut scanner = Scanner::new();

        let symbols = scanner.scan(&mut img);
        assert!(!symbols.is_empty(), "Expected to find at least one QR code");

        let symbol = symbols.into_iter().next().unwrap();
        assert_eq!(symbol.symbol_type(), SymbolType::QrCode);
        let data = symbol.data_string().unwrap_or("");
        assert_eq!(
            data,
            "二维码生成器
https://zh.qr-code-generator.com
打印出来的二维码至少要2厘米宽，确保用任何设备或应用都可以成功扫描。 ... 通过三个简单步骤，就能使用二维码生成器在几秒钟内创建一个二维码。首先，选择二维码的 ..."
        );
    }

    #[test]
    fn test_small_qr_code_upscaling() {
        // Test for small QR code (109x89 pixels) that requires automatic upscaling.
        // This QR code has ~3 pixels per module which is too small for reliable
        // detection without upscaling. The scanner should automatically upscale
        // small images to improve detection.
        let img = ::image::ImageReader::open("examples/github-issue-qr.png")
            .expect("Failed to open github-issue-qr.png")
            .decode()
            .expect("Failed to decode image");

        let gray = img.to_luma8();
        let (width, height) = gray.dimensions();

        // Verify the image is small (< 200px)
        assert!(
            width < 200 || height < 200,
            "Test image should be small to test upscaling"
        );

        let mut zedbar_img =
            Image::from_gray(gray.as_raw(), width, height).expect("Failed to create zedbar image");

        // Should decode successfully with automatic upscaling
        let mut scanner = Scanner::new();
        let symbols = scanner.scan(&mut zedbar_img);

        assert!(
            !symbols.is_empty(),
            "Expected to decode small QR code with automatic upscaling"
        );

        let symbol = symbols.into_iter().next().unwrap();
        assert_eq!(symbol.symbol_type(), SymbolType::QrCode);
        assert_eq!(
            symbol.data_string().unwrap_or(""),
            "S;1;019be05c-54ba-7b32-a6f5-c06b288a46e8;A"
        );
    }

    #[test]
    fn test_small_qr_code_upscaling_disabled() {
        // Verify that disabling upscaling prevents detection of small QR codes
        use crate::config::*;

        let img = ::image::ImageReader::open("examples/github-issue-qr.png")
            .expect("Failed to open github-issue-qr.png")
            .decode()
            .expect("Failed to decode image");

        let gray = img.to_luma8();
        let (width, height) = gray.dimensions();

        let mut zedbar_img =
            Image::from_gray(gray.as_raw(), width, height).expect("Failed to create zedbar image");

        // Disable upscaling
        let config = DecoderConfig::new().upscale_small_images(false);
        let mut scanner = Scanner::with_config(config);
        let symbols = scanner.scan(&mut zedbar_img);

        // Without upscaling, the small QR code should NOT be detected
        assert!(
            symbols.is_empty(),
            "Small QR code should not be detected when upscaling is disabled"
        );
    }
}
