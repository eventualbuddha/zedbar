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
//!     .set_length_limits(Code39, 4, 20)  // Code39 must be 4-20 chars
//!     .test_inverted(true)               // Try inverted image if no symbols found
//!     .retry_undecoded_regions(true)     // Crop+upscale small QR codes automatically
//!     .scan_density(2, 2);               // Scan every 2nd line (faster)
//!
//! let mut scanner = Scanner::with_config(config);
//! ```
//!
//! # Small QR Codes
//!
//! When a QR code is too small to decode (e.g. on a scanned page),
//! [`ScanResult::finder_regions()`](scanner::ScanResult::finder_regions) reports
//! where finder patterns were detected. Use [`Image::crop`] and [`Image::upscale`]
//! to retry those regions, or enable
//! [`retry_undecoded_regions`](DecoderConfig::retry_undecoded_regions)
//! to do this automatically.
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
pub use scanner::{FinderRegion, ScanResult, Scanner};
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
        let mut img = Image::from_gray(data, width, height).expect("Failed to create zedbar image");

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
        println!(
            "✓ zedbar and rqrr agree on {} symbols",
            zedbar_results.len()
        );
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
        let mut img = Image::from_gray(data, width, height).expect("Failed to create zedbar image");

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
        println!(
            "✓ zedbar and rqrr agree on {} symbols",
            zedbar_results.len()
        );
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

        let mut img = Image::from_gray(data, width, height).expect("Failed to create zedbar image");

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

        let mut img = Image::from_gray(data, width, height).expect("Failed to create zedbar image");

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

        let mut img = Image::from_gray(data, width, height).expect("Failed to create zedbar image");

        // UPC-A surfaces as a label on the EAN-13 decoder; both must be enabled.
        use crate::config::*;
        let config = DecoderConfig::new().enable(Ean13).enable(Upca);
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

        let mut img = Image::from_gray(data, width, height).expect("Failed to create zedbar image");

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

        let mut img = Image::from_gray(data, width, height).expect("Failed to create zedbar image");

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

        let mut img = Image::from_gray(data, width, height).expect("Failed to create zedbar image");

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

        let mut img = Image::from_gray(data, width, height).expect("Failed to create zedbar image");

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

        let mut img = Image::from_gray(data, width, height).expect("Failed to create zedbar image");

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

        let mut img = Image::from_gray(data, width, height).expect("Failed to create zedbar image");

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

        let mut img = Image::from_gray(data, width, height).expect("Failed to create zedbar image");

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
    fn test_finder_region_manual_workflow() {
        // When a QR code is too small to decode, the scanner reports finder
        // regions. The caller can crop and upscale these to recover the QR.
        let img = ::image::ImageReader::open("examples/synthetic-small-qr-page.png")
            .expect("Failed to open synthetic-small-qr-page.png")
            .decode()
            .expect("Failed to decode image");

        let gray = img.to_luma8();
        let (width, height) = gray.dimensions();

        let mut zedbar_img =
            Image::from_gray(gray.as_raw(), width, height).expect("Failed to create zedbar image");

        let mut scanner = Scanner::new();
        let result = scanner.scan(&mut zedbar_img);

        assert!(result.symbols().is_empty());
        assert!(
            !result.finder_regions().is_empty(),
            "Expected finder region for the small QR code"
        );

        // Crop and upscale finder regions to recover the QR code
        let mut recovered = Vec::new();
        for region in result.finder_regions() {
            let pad = region.width.max(region.height) / 2;
            let x = region.x.saturating_sub(pad);
            let y = region.y.saturating_sub(pad);
            let w = (region.width + 2 * pad).min(width - x);
            let h = (region.height + 2 * pad).min(height - y);

            if let Some(mut upscaled) = zedbar_img.crop(x, y, w, h).and_then(|c| c.upscale(4)) {
                let retry = scanner.scan(&mut upscaled);
                recovered.extend(retry);
            }
        }

        assert!(!recovered.is_empty());
        assert_eq!(
            recovered[0].data_string().unwrap_or(""),
            "https://zedbar.invalid/fixture/small-qr"
        );
    }

    #[test]
    fn test_retry_undecoded_regions() {
        use crate::config::*;

        // Synthetic fixtures generated by examples/generate_test_fixtures.rs.
        // - small standalone: QR fills a tiny frame. Decodes at the top
        //   level here but is still a useful sanity check.
        // - large page:       small QR in a big image — cannot decode at
        //   the top level, must recover via crop+upscale retry.
        let images = [
            "examples/synthetic-small-qr.png",
            "examples/synthetic-small-qr-page.png",
        ];
        for path in images {
            let img = ::image::ImageReader::open(path)
                .unwrap_or_else(|_| panic!("Failed to open {path}"))
                .decode()
                .unwrap_or_else(|_| panic!("Failed to decode {path}"));

            let gray = img.to_luma8();
            let (width, height) = gray.dimensions();
            let mut zedbar_img = Image::from_gray(gray.as_raw(), width, height)
                .expect("Failed to create zedbar image");

            let config = DecoderConfig::new()
                .enable(QrCode)
                .retry_undecoded_regions(true);
            let mut scanner = Scanner::with_config(config);
            let result = scanner.scan(&mut zedbar_img);

            assert!(
                !result.symbols().is_empty(),
                "{path}: expected retry to decode the QR code"
            );
            assert_eq!(
                result.symbols()[0].data_string().unwrap_or(""),
                "https://zedbar.invalid/fixture/small-qr",
                "{path}: wrong decoded data"
            );
        }
    }

    /// Regression test for the "small QR buried in page-wide redaction noise"
    /// failure mode. We build a synthetic redacted document on top of the
    /// `synthetic-small-qr-page.png` fixture by overlaying:
    ///
    /// 1. Large black "redaction" capsules (smooth edges, few finder lines).
    /// 2. Thousands of small specks of salt-and-pepper noise (many false
    ///    1:1:3:1:1 finder-line candidates on every scanline).
    /// 3. A couple of rows of tightly-spaced short bars that trigger real
    ///    finder-pattern ratios in places that aren't actually a QR.
    ///
    /// The union of all raw finder-line positions then spans most of the
    /// page. A single-bbox retry path can't use that — it's rejected by the
    /// per-image area filter. This test ensures the multi-region retry path
    /// still finds the QR because it bboxes each cluster group separately.
    #[test]
    fn test_retry_survives_redaction_noise() {
        use crate::config::*;

        let gray = synthesize_redacted_page();
        let (width, height) = gray.dimensions();
        let mut zedbar_img =
            Image::from_gray(gray.as_raw(), width, height).expect("Failed to create zedbar image");

        let config = DecoderConfig::new()
            .enable(QrCode)
            .retry_undecoded_regions(true);
        let mut scanner = Scanner::with_config(config);
        let result = scanner.scan(&mut zedbar_img);

        assert!(
            !result.symbols().is_empty(),
            "expected retry to recover the QR despite redaction noise"
        );
        assert_eq!(
            result.symbols()[0].data_string().unwrap_or(""),
            "https://zedbar.invalid/fixture/small-qr",
        );
    }

    fn synthesize_redacted_page() -> ::image::GrayImage {
        use ::image::Luma;

        let mut page = ::image::ImageReader::open("examples/synthetic-small-qr-page.png")
            .expect("Failed to open synthetic-small-qr-page.png")
            .decode()
            .expect("Failed to decode image")
            .to_luma8();
        let (page_w, page_h) = page.dimensions();

        // Rectangle occupied by the real QR code (plus quiet zone). We avoid
        // drawing over it so the QR stays decodable post-crop.
        let mut qr_min_x = page_w as i32;
        let mut qr_max_x = -1;
        let mut qr_min_y = page_h as i32;
        let mut qr_max_y = -1;
        for (_, row) in page.enumerate_rows() {
            for (x, y, luma) in row {
                if luma.0[0] < 0x7f {
                    qr_min_x = qr_min_x.min(x as i32);
                    qr_max_x = qr_max_x.max(x as i32);
                    qr_min_y = qr_min_y.min(y as i32);
                    qr_max_y = qr_max_y.max(y as i32);
                }
            }
        }
        let qr_rect = (qr_min_x - 30, qr_min_y - 30, qr_max_x + 30, qr_max_y + 30);

        // Deterministic noise — don't pull in rand as a dev-dep just for this.
        let mut rng = Lcg(0x5A17_E51D_DEAD_BEEF);

        let capsules = [
            (10i32, 20, 380, 110),
            (400, 20, 330, 90),
            (10, 170, 230, 70),
            (10, 300, 620, 40),
            (10, 350, 720, 40),
            (10, 430, 620, 45),
            (10, 490, 510, 45),
            (10, 550, 480, 50),
            (10, 640, 680, 50),
            (10, 710, 680, 50),
            (10, 790, 730, 55),
            (10, 880, 770, 200),
        ];
        for (bx, by, bw, bh) in capsules {
            fill_capsule(&mut page, bx, by, bw, bh);
        }

        for _ in 0..4500 {
            let x = rng.range(0, page_w as i32);
            let y = rng.range(0, page_h as i32);
            if inside_rect(x, y, qr_rect) {
                continue;
            }
            let size = rng.range(1, 4);
            for dy in 0..size {
                for dx in 0..size {
                    let px = x + dx;
                    let py = y + dy;
                    if px >= 0 && py >= 0 && (px as u32) < page_w && (py as u32) < page_h {
                        page.put_pixel(px as u32, py as u32, Luma([0]));
                    }
                }
            }
        }

        for row in 0..20i32 {
            let y = 260 + row * 30;
            if y + 6 > page_h as i32 || (175..=310).contains(&y) {
                continue;
            }
            let mut x = 20i32;
            while x < page_w as i32 - 30 {
                if !(x + 30 < qr_rect.0 || x > qr_rect.2) {
                    x += 120;
                    continue;
                }
                for dy in 0..6 {
                    for dx in 0..3 {
                        let px = x + dx;
                        let py = y + dy;
                        if px >= 0 && py >= 0 && (px as u32) < page_w && (py as u32) < page_h {
                            page.put_pixel(px as u32, py as u32, Luma([0]));
                        }
                    }
                }
                x += rng.range(4, 10);
            }
        }

        page
    }

    fn inside_rect(x: i32, y: i32, r: (i32, i32, i32, i32)) -> bool {
        x >= r.0 && x <= r.2 && y >= r.1 && y <= r.3
    }

    fn fill_capsule(page: &mut ::image::GrayImage, bx: i32, by: i32, bw: i32, bh: i32) {
        use ::image::Luma;
        let (page_w, page_h) = page.dimensions();
        let radius = bh.min(bw) / 2;
        let ry_max = (bh - 1).max(0);
        let rx_max = (bw - 1).max(0);
        for dy in 0..=ry_max {
            for dx in 0..=rx_max {
                let inside_body = dx >= radius && dx <= (rx_max - radius);
                let is_inside = if inside_body {
                    true
                } else {
                    let cx = if dx < radius { radius } else { rx_max - radius };
                    let cy = ry_max / 2;
                    let dxc = dx - cx;
                    let dyc = dy - cy;
                    dxc * dxc + dyc * dyc <= radius * radius
                };
                if !is_inside {
                    continue;
                }
                let px = bx + dx;
                let py = by + dy;
                if px >= 0 && py >= 0 && (px as u32) < page_w && (py as u32) < page_h {
                    page.put_pixel(px as u32, py as u32, Luma([0]));
                }
            }
        }
    }

    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> u32 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (self.0 >> 33) as u32
        }
        fn range(&mut self, lo: i32, hi: i32) -> i32 {
            let span = (hi - lo).max(1) as u32;
            lo + (self.next() % span) as i32
        }
    }
}
