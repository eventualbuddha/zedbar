//! Command-line barcode scanner (Rust port of zbarimg)

use std::env;
use std::process;
use zbar::SymbolType;
use zbar::{Image, Scanner};

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() != 2 {
        eprintln!("Usage: {} <image_file>", args[0]);
        process::exit(1);
    }

    let filename = &args[1];

    // Load the image
    let img = match ::image::ImageReader::open(filename) {
        Ok(reader) => match reader.decode() {
            Ok(img) => img,
            Err(e) => {
                eprintln!("Failed to decode image '{}': {}", filename, e);
                process::exit(1);
            }
        },
        Err(e) => {
            eprintln!("Failed to open image '{}': {}", filename, e);
            process::exit(1);
        }
    };

    // Convert to grayscale
    let gray = img.to_luma8();
    let (width, height) = gray.dimensions();
    let data = gray.as_raw();

    // Create ZBar image from grayscale data
    let mut zbar_img = match Image::from_gray(data, width, height) {
        Ok(img) => img,
        Err(e) => {
            eprintln!("Failed to create ZBar image: {}", e);
            process::exit(1);
        }
    };

    // Create scanner and configure for all symbologies
    let mut scanner = Scanner::new();

    // Enable QR codes
    if let Err(e) = scanner.set_config(SymbolType::None, zbar::scanner::Config::Enable, 1) {
        eprintln!("Failed to configure scanner: {}", e);
        process::exit(1);
    }

    // Scan the image
    let num_symbols = match scanner.scan(&mut zbar_img) {
        Ok(n) => n,
        Err(e) => {
            eprintln!("Failed to scan image: {}", e);
            process::exit(1);
        }
    };

    if num_symbols == 0 {
        println!("No barcodes found in {}", filename);
        process::exit(0);
    }

    // Print results
    let symbols = zbar_img.symbols();
    for symbol in symbols {
        let symbol_type = symbol.symbol_type();
        let data_bytes = symbol.data();
        // Try UTF-8 first, otherwise print raw bytes
        if let Ok(data_str) = std::str::from_utf8(data_bytes) {
            println!("{:?}:{}", symbol_type, data_str);
        } else {
            // For binary data, print as raw bytes
            print!("{:?}:", symbol_type);
            std::io::Write::write_all(&mut std::io::stdout(), data_bytes).ok();
            println!();
        }
    }
}
