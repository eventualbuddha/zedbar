//! ZBar Barcode Scanning Library (Rust Port)
//!
//! This crate provides barcode scanning functionality, originally based on the C ZBar library.
//! The conversion to Rust is being done incrementally.

#![allow(clippy::missing_safety_doc)]
#![allow(non_camel_case_types)]

pub mod decoder;
pub mod decoder_types;
pub mod decoders;
pub mod error;
pub mod finder;
pub mod image;
pub mod image_ffi;
pub mod img_scanner;
pub mod line_scanner;
pub mod qrcode;
pub mod scanner;
pub mod sqcode;
pub mod symbol;

// Re-export main types
pub use decoder::Decoder;
pub use error::{Error, Result};
pub use image::Image;
pub use scanner::Scanner;
pub use symbol::{Symbol, SymbolSet, SymbolType};

mod ffi;

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

        // Create scanner and configure for QR codes
        let mut scanner = Scanner::new();
        scanner
            .set_config(SymbolType::QrCode, scanner::Config::Enable, 1)
            .expect("Failed to configure scanner");

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

        // Create scanner and configure for QR codes
        let mut scanner = Scanner::new();
        scanner
            .set_config(SymbolType::QrCode, scanner::Config::Enable, 1)
            .expect("Failed to configure scanner");

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

        let mut scanner = Scanner::new();
        scanner
            .set_config(SymbolType::Ean13, scanner::Config::Enable, 1)
            .expect("Failed to configure scanner");

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

        let mut scanner = Scanner::new();
        scanner
            .set_config(SymbolType::Ean8, scanner::Config::Enable, 1)
            .expect("Failed to configure scanner");

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

        let mut scanner = Scanner::new();
        scanner
            .set_config(SymbolType::Upca, scanner::Config::Enable, 1)
            .expect("Failed to configure scanner");

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

        let mut scanner = Scanner::new();
        scanner
            .set_config(SymbolType::Code128, scanner::Config::Enable, 1)
            .expect("Failed to configure scanner");

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

        let mut scanner = Scanner::new();
        scanner
            .set_config(SymbolType::Code39, scanner::Config::Enable, 1)
            .expect("Failed to configure scanner");

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

        let mut scanner = Scanner::new();
        scanner
            .set_config(SymbolType::Code93, scanner::Config::Enable, 1)
            .expect("Failed to configure scanner");

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

        let mut scanner = Scanner::new();
        scanner
            .set_config(SymbolType::Codabar, scanner::Config::Enable, 1)
            .expect("Failed to configure scanner");

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

        let mut scanner = Scanner::new();
        scanner
            .set_config(SymbolType::I25, scanner::Config::Enable, 1)
            .expect("Failed to configure scanner");

        let num_symbols = scanner.scan(&mut zbar_img).expect("Failed to scan image");
        assert!(num_symbols > 0, "Expected to find at least one I25 barcode");

        let symbols = zbar_img.symbols();
        for symbol in symbols {
            assert_eq!(symbol.symbol_type(), SymbolType::I25);
            println!("Decoded I25: {}", symbol.data_string().unwrap_or(""));
        }
    }
}
