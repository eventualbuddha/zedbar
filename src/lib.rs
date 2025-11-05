//! ZBar Barcode Scanning Library
//!
//! This crate provides barcode scanning functionality for various barcode formats including
//! QR codes, EAN, UPC, Code128, Code39, and more.

#![allow(clippy::missing_safety_doc)]
#![allow(non_camel_case_types)]

pub mod color;
pub mod config;
pub mod decoder;
pub mod decoders;
pub mod error;
pub mod finder;
pub mod image;
pub mod image_ffi;
pub mod img_scanner;
pub mod img_scanner_config;
pub mod qrcode;
pub mod scanner;
pub mod sqcode;
pub mod symbol;

// Re-export main types
pub use config::DecoderConfig;
pub use error::{Error, Result};
pub use image::Image;
pub use scanner::Scanner;
pub use symbol::{Orientation, SymbolType};

#[cfg(test)]
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

        // Create ZBar image from grayscale data
        let mut zbar_img =
            Image::from_gray(data, width, height).expect("Failed to create ZBar image");

        // Create scanner (QR codes enabled by default)
        let mut scanner = Scanner::new();

        // Scan the image
        let num_symbols = scanner.scan(&mut zbar_img).expect("Failed to scan image");

        println!("Found {num_symbols} symbols with zbar");

        // Get and verify symbols
        let symbols = zbar_img.symbols();
        assert!(!symbols.is_empty(), "Expected to find at least one QR code");

        let mut zbar_results = Vec::new();
        for symbol in symbols {
            let symbol_type = symbol.symbol_type();
            let data = symbol.data();
            println!("Decoded {symbol_type:?}: {} bytes", data.len());

            assert_eq!(symbol_type, SymbolType::QrCode);
            assert!(!data.is_empty(), "QR code data should not be empty");
            zbar_results.push(data.to_vec());
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
        zbar_results.sort();
        rqrr_results.sort();

        assert_eq!(
            zbar_results, rqrr_results,
            "zbar and rqrr produced different results"
        );
        println!("✓ zbar and rqrr agree on {} symbols", zbar_results.len());
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

        // Create ZBar image from grayscale data
        let mut zbar_img =
            Image::from_gray(data, width, height).expect("Failed to create ZBar image");

        // Create scanner (QR codes enabled by default)
        let mut scanner = Scanner::new();

        // Scan the image
        let num_symbols = scanner.scan(&mut zbar_img).expect("Failed to scan image");

        println!("Found {num_symbols} symbols with zbar");

        // Get and verify symbols
        let symbols = zbar_img.symbols();
        assert!(!symbols.is_empty(), "Expected to find at least one QR code");

        let mut zbar_results = Vec::new();
        for symbol in symbols {
            let symbol_type = symbol.symbol_type();
            let data = symbol.data();
            println!("Decoded {symbol_type:?}: {} bytes", data.len());

            assert_eq!(symbol_type, SymbolType::QrCode);
            assert!(!data.is_empty(), "QR code data should not be empty");
            zbar_results.push(data.to_vec());
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
        zbar_results.sort();
        rqrr_results.sort();

        assert_eq!(
            zbar_results, rqrr_results,
            "zbar and rqrr produced different results"
        );
        println!("✓ zbar and rqrr agree on {} symbols", zbar_results.len());
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

        // Build a ZBar image from the inverted data
        let mut zbar_img =
            Image::from_gray(&inverted, width, height).expect("Failed to create ZBar image");

        // Configure scanner for QR codes with inverted testing
        use crate::config::*;
        let config = DecoderConfig::new().enable(QrCode).test_inverted(true);
        let mut scanner = Scanner::with_config(config);

        // Scan the inverted image
        let num_symbols = scanner
            .scan(&mut zbar_img)
            .expect("Failed to scan inverted image");
        println!("Found {num_symbols} symbols in inverted image with zbar");

        // Gather decoded symbols
        let symbols = zbar_img.symbols();
        assert!(
            !symbols.is_empty(),
            "Expected to find at least one QR code in inverted image"
        );

        let mut zbar_results = Vec::new();
        for symbol in symbols {
            let symbol_type = symbol.symbol_type();
            let data = symbol.data();
            println!("Decoded {symbol_type:?} (inverted): {} bytes", data.len());

            assert_eq!(symbol_type, SymbolType::QrCode);
            assert!(!data.is_empty(), "QR code data should not be empty");
            zbar_results.push(data.to_vec());
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
        zbar_results.sort();
        rqrr_results.sort();

        assert_eq!(
            zbar_results, rqrr_results,
            "zbar and rqrr produced different results for inverted image"
        );
        println!(
            "✓ zbar and rqrr agree on {} symbols for inverted image",
            zbar_results.len()
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

        let mut zbar_img =
            Image::from_gray(data, width, height).expect("Failed to create ZBar image");

        // Create scanner (EAN-13 enabled by default)
        let mut scanner = Scanner::new();

        let num_symbols = scanner.scan(&mut zbar_img).expect("Failed to scan image");
        assert!(
            num_symbols > 0,
            "Expected to find at least one EAN-13 barcode"
        );

        let symbols = zbar_img.symbols();
        for symbol in symbols {
            assert_eq!(symbol.symbol_type(), SymbolType::Ean13);
            println!("Decoded EAN-13: {}", symbol.data_string().unwrap_or(""));
        }
    }

    #[test]
    fn test_ean8_decode() {
        // EAN-8 test (8 digits, compact retail barcode)
        // Note: ZBar may report EAN-8 as EAN-13 with zero padding
        let img = ::image::ImageReader::open("examples/test-ean8.png")
            .expect("Failed to open test-ean8.png")
            .decode()
            .expect("Failed to decode image");

        let gray = img.to_luma8();
        let (width, height) = gray.dimensions();
        let data = gray.as_raw();

        let mut zbar_img =
            Image::from_gray(data, width, height).expect("Failed to create ZBar image");

        // Create scanner (EAN-8 enabled by default)
        let mut scanner = Scanner::new();

        let num_symbols = scanner.scan(&mut zbar_img).expect("Failed to scan image");
        assert!(
            num_symbols > 0,
            "Expected to find at least one EAN-8 barcode"
        );

        let symbols = zbar_img.symbols();
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

        let mut zbar_img =
            Image::from_gray(data, width, height).expect("Failed to create ZBar image");

        // Create scanner with UPC-A explicitly enabled
        use crate::config::*;
        let config = DecoderConfig::new().enable(Upca);
        let mut scanner = Scanner::with_config(config);

        let num_symbols = scanner.scan(&mut zbar_img).expect("Failed to scan image");
        assert!(
            num_symbols > 0,
            "Expected to find at least one UPC-A barcode"
        );

        let symbols = zbar_img.symbols();
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

        let mut zbar_img =
            Image::from_gray(data, width, height).expect("Failed to create ZBar image");

        // Create scanner (Code128 enabled by default)
        let mut scanner = Scanner::new();

        let num_symbols = scanner.scan(&mut zbar_img).expect("Failed to scan image");
        assert!(
            num_symbols > 0,
            "Expected to find at least one Code128 barcode"
        );

        let symbols = zbar_img.symbols();
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

        let mut zbar_img =
            Image::from_gray(data, width, height).expect("Failed to create ZBar image");

        // Create scanner (Code39 enabled by default)
        let mut scanner = Scanner::new();

        let num_symbols = scanner.scan(&mut zbar_img).expect("Failed to scan image");
        assert!(
            num_symbols > 0,
            "Expected to find at least one Code39 barcode"
        );

        let symbols = zbar_img.symbols();
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

        let mut zbar_img =
            Image::from_gray(data, width, height).expect("Failed to create ZBar image");

        // Create scanner (Code93 enabled by default)
        let mut scanner = Scanner::new();

        let num_symbols = scanner.scan(&mut zbar_img).expect("Failed to scan image");
        assert!(
            num_symbols > 0,
            "Expected to find at least one Code93 barcode"
        );

        let symbols = zbar_img.symbols();
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

        let mut zbar_img =
            Image::from_gray(data, width, height).expect("Failed to create ZBar image");

        // Create scanner (Codabar enabled by default)
        let mut scanner = Scanner::new();

        let num_symbols = scanner.scan(&mut zbar_img).expect("Failed to scan image");
        assert!(
            num_symbols > 0,
            "Expected to find at least one Codabar barcode"
        );

        let symbols = zbar_img.symbols();
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

        let mut zbar_img =
            Image::from_gray(data, width, height).expect("Failed to create ZBar image");

        // Create scanner (I25 enabled by default)
        let mut scanner = Scanner::new();

        let num_symbols = scanner.scan(&mut zbar_img).expect("Failed to scan image");
        assert!(num_symbols > 0, "Expected to find at least one I25 barcode");

        let symbols = zbar_img.symbols();
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

        let mut zbar_img =
            Image::from_gray(data, width, height).expect("Failed to create ZBar image");

        // Create scanner (QR codes enabled by default)
        let mut scanner = Scanner::new();

        let num_symbols = scanner.scan(&mut zbar_img).expect("Failed to scan image");
        assert!(
            num_symbols > 0,
            "Expected to find at least one QR code in pixel-wifi-sharing-qr-code.png"
        );

        let symbol = zbar_img.symbols().into_iter().next().unwrap();
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

        let mut zbar_img =
            Image::from_gray(data, width, height).expect("Failed to create ZBar image");

        // Create scanner (QR codes enabled by default)
        let mut scanner = Scanner::new();

        let num_symbols = scanner.scan(&mut zbar_img).expect("Failed to scan image");
        assert!(num_symbols > 0, "Expected to find at least one QR code");

        let symbol = zbar_img.symbols().into_iter().next().unwrap();
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

        let mut zbar_img =
            Image::from_gray(&normalized, width, height).expect("Failed to create ZBar image");

        // Create scanner (QR codes enabled by default)
        let mut scanner = Scanner::new();

        let num_symbols = scanner.scan(&mut zbar_img).expect("Failed to scan image");
        assert!(num_symbols > 0, "Expected to find at least one QR code");

        let symbol = zbar_img.symbols().into_iter().next().unwrap();
        assert_eq!(symbol.symbol_type(), SymbolType::QrCode);
        let data = symbol.data_string().unwrap_or("");
        assert_eq!(
            data,
            "二维码生成器
https://zh.qr-code-generator.com
打印出来的二维码至少要2厘米宽，确保用任何设备或应用都可以成功扫描。 ... 通过三个简单步骤，就能使用二维码生成器在几秒钟内创建一个二维码。首先，选择二维码的 ..."
        );
    }
}
