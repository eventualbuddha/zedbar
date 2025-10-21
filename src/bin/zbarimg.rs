//! Command-line barcode scanner (Rust port of zbarimg)

use clap::Parser;
use std::io::Write;
use std::process;
use zbar::SymbolType;
use zbar::{Image, Scanner};

/// Scan and decode bar codes from one or more image files
#[derive(Parser)]
#[command(name = "zbarimg")]
#[command(version)]
#[command(about = "Scan and decode bar codes from one or more image files", long_about = None)]
struct Args {
    /// Minimal output, only print decoded symbol data
    #[arg(short, long)]
    quiet: bool,

    /// Output decoded symbol data without converting charsets
    #[arg(long)]
    raw: bool,

    /// Image files to scan
    #[arg(required = true)]
    files: Vec<String>,
}

fn main() {
    let args = Args::parse();

    let mut total_symbols = 0;

    for filename in &args.files {
        // Load the image
        let img = match ::image::ImageReader::open(filename) {
            Ok(reader) => match reader.decode() {
                Ok(img) => img,
                Err(e) => {
                    if !args.quiet {
                        eprintln!("Failed to decode image '{}': {}", filename, e);
                    }
                    process::exit(1);
                }
            },
            Err(e) => {
                if !args.quiet {
                    eprintln!("Failed to open image '{}': {}", filename, e);
                }
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
                if !args.quiet {
                    eprintln!("Failed to create ZBar image: {}", e);
                }
                process::exit(1);
            }
        };

        // Create scanner and configure for all symbologies
        let mut scanner = Scanner::new();

        // Enable all barcode types
        if let Err(e) = scanner.set_config(SymbolType::None, zbar::scanner::Config::Enable, 1) {
            if !args.quiet {
                eprintln!("Failed to configure scanner: {}", e);
            }
            process::exit(1);
        }

        // Scan the image
        let num_symbols = match scanner.scan(&mut zbar_img) {
            Ok(n) => n,
            Err(e) => {
                if !args.quiet {
                    eprintln!("Failed to scan image: {}", e);
                }
                process::exit(1);
            }
        };

        total_symbols += num_symbols;

        // Print results
        let symbols = zbar_img.symbols();
        for symbol in symbols {
            let symbol_type = symbol.symbol_type();
            let data_bytes = symbol.data();

            if args.raw {
                // Raw mode: just output the data bytes with newline
                std::io::stdout().write_all(data_bytes).ok();
                println!();
            } else {
                // Normal mode: include symbol type
                // Try UTF-8 first, otherwise print raw bytes
                if let Ok(data_str) = std::str::from_utf8(data_bytes) {
                    println!("{:?}:{}", symbol_type, data_str);
                } else {
                    // For binary data, print as raw bytes
                    print!("{:?}:", symbol_type);
                    std::io::stdout().write_all(data_bytes).ok();
                    println!();
                }
            }
        }
    }

    // Print statistics unless quiet mode
    if !args.quiet {
        if total_symbols == 0 {
            eprintln!("No barcodes found");
            process::exit(1);
        } else {
            eprintln!("scanned {} barcode symbols from {} image(s)", 
                     total_symbols, args.files.len());
        }
    } else if total_symbols == 0 {
        // In quiet mode, still exit with error if no barcodes found
        process::exit(1);
    }
}
