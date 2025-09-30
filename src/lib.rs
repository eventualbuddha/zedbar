//! ZBar Barcode Scanning Library (Rust Port)
//!
//! This crate provides barcode scanning functionality, originally based on the C ZBar library.
//! The conversion to Rust is being done incrementally using c2rust as a starting point.

pub mod decoder;
pub mod error;
pub mod image;
pub mod scanner;
pub mod symbol;

// Re-export main types
pub use decoder::Decoder;
pub use error::{Error, Result};
pub use image::Image;
pub use scanner::Scanner;
pub use symbol::{Symbol, SymbolSet, SymbolType};

// For now, we'll use FFI bindings to the C library
// These will be replaced as modules are converted
mod ffi;
pub use ffi::*;

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
            .set_config(SymbolType::QrCode as i32, scanner::Config::Enable, 1)
            .expect("Failed to configure scanner");

        // Scan the image
        let num_symbols = scanner.scan(&mut zbar_img).expect("Failed to scan image");

        println!("Found {} symbols", num_symbols);

        // Get and verify symbols
        let symbols = zbar_img.symbols();
        assert!(!symbols.is_empty(), "Expected to find at least one QR code");

        for symbol in symbols {
            let symbol_type = symbol.symbol_type();
            let data = symbol.data_string().unwrap_or("<invalid UTF-8>");
            println!("Decoded {:?}: {}", symbol_type, data);

            assert_eq!(symbol_type, SymbolType::QrCode);
            assert!(!data.is_empty(), "QR code data should not be empty");
        }
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
            .set_config(SymbolType::QrCode as i32, scanner::Config::Enable, 1)
            .expect("Failed to configure scanner");

        // Scan the image
        let num_symbols = scanner.scan(&mut zbar_img).expect("Failed to scan image");

        println!("Found {} symbols", num_symbols);

        // Get and verify symbols
        let symbols = zbar_img.symbols();
        assert!(!symbols.is_empty(), "Expected to find at least one QR code");

        for symbol in symbols {
            let symbol_type = symbol.symbol_type();
            let data = symbol.data_string().unwrap_or("<invalid UTF-8>");
            println!("Decoded {:?}: {}", symbol_type, data);

            assert_eq!(symbol_type, SymbolType::QrCode);
            assert!(!data.is_empty(), "QR code data should not be empty");
        }
    }
}

