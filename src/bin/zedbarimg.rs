//! Command-line barcode scanner
//!
//! This binary provides a command-line interface to the zedbar library.
//! It requires all symbology features to be enabled to provide full CLI control.
//! 
//! Build with: `cargo build --bin zedbarimg` (uses default features)

// The zedbarimg binary provides CLI flags for all symbologies,
// so it requires all features to be enabled. Users who want a minimal
// library should build with `--lib` and select specific features.
#[cfg(not(all(
    feature = "qrcode",
    feature = "sqcode",
    feature = "ean",
    feature = "code128",
    feature = "code39",
    feature = "code93",
    feature = "codabar",
    feature = "databar",
    feature = "i25"
)))]
compile_error!(
    "zedbarimg binary requires all symbology features enabled. \
     Build with default features: `cargo build --bin zedbarimg` \
     For a minimal library, use: `cargo build --lib --no-default-features --features qrcode`"
);

use clap::Parser;
use std::io::Write;
use std::process;
use zedbar::{config::*, DecoderConfig, Image, Scanner};

/// Scan and decode bar codes from one or more image files
#[derive(Parser)]
#[command(name = "zedbarimg")]
#[command(version)]
#[command(about = "Scan and decode bar codes from one or more image files", long_about = None)]
struct Args {
    /// Minimal output, only print decoded symbol data
    #[arg(short, long)]
    quiet: bool,

    /// Output decoded symbol data without converting charsets
    #[arg(long)]
    raw: bool,

    /// Disable all symbologies (use with --enable-* to select specific ones)
    #[arg(long)]
    disable_all: bool,

    /// Enable EAN-13 barcode
    #[arg(long)]
    enable_ean13: bool,

    /// Enable EAN-8 barcode
    #[arg(long)]
    enable_ean8: bool,

    /// Enable EAN-2 add-on
    #[arg(long)]
    enable_ean2: bool,

    /// Enable EAN-5 add-on
    #[arg(long)]
    enable_ean5: bool,

    /// Enable UPC-A barcode
    #[arg(long)]
    enable_upca: bool,

    /// Enable UPC-E barcode
    #[arg(long)]
    enable_upce: bool,

    /// Enable ISBN-10
    #[arg(long)]
    enable_isbn10: bool,

    /// Enable ISBN-13
    #[arg(long)]
    enable_isbn13: bool,

    /// Enable Code 39
    #[arg(long)]
    enable_code39: bool,

    /// Enable Code 93
    #[arg(long)]
    enable_code93: bool,

    /// Enable Code 128
    #[arg(long)]
    enable_code128: bool,

    /// Enable Interleaved 2 of 5
    #[arg(long)]
    enable_i25: bool,

    /// Enable Codabar
    #[arg(long)]
    enable_codabar: bool,

    /// Enable DataBar (RSS-14)
    #[arg(long)]
    enable_databar: bool,

    /// Enable DataBar Expanded
    #[arg(long)]
    enable_databar_exp: bool,

    /// Enable QR Code
    #[arg(long)]
    enable_qr: bool,

    /// Enable SQ Code
    #[arg(long)]
    enable_sq: bool,

    /// Disable EAN-13 barcode
    #[arg(long)]
    disable_ean13: bool,

    /// Disable EAN-8 barcode
    #[arg(long)]
    disable_ean8: bool,

    /// Disable Code 39
    #[arg(long)]
    disable_code39: bool,

    /// Disable Code 93
    #[arg(long)]
    disable_code93: bool,

    /// Disable Code 128
    #[arg(long)]
    disable_code128: bool,

    /// Disable QR Code
    #[arg(long)]
    disable_qr: bool,

    /// Image files to scan
    #[arg(required = true)]
    files: Vec<String>,
}

fn main() {
    let args = Args::parse();

    // Build configuration based on command-line arguments
    let config = if args.disable_all {
        // Start with all symbologies disabled, then enable selected ones
        let mut cfg = DecoderConfig::new().disable_all();

        // Apply enable flags
        if args.enable_ean13 {
            cfg = cfg.enable(Ean13);
        }
        if args.enable_ean8 {
            cfg = cfg.enable(Ean8);
        }
        if args.enable_ean2 {
            cfg = cfg.enable(Ean2);
        }
        if args.enable_ean5 {
            cfg = cfg.enable(Ean5);
        }
        if args.enable_upca {
            cfg = cfg.enable(Upca);
        }
        if args.enable_upce {
            cfg = cfg.enable(Upce);
        }
        if args.enable_isbn10 {
            cfg = cfg.enable(Isbn10);
        }
        if args.enable_isbn13 {
            cfg = cfg.enable(Isbn13);
        }
        if args.enable_code39 {
            cfg = cfg.enable(Code39);
        }
        if args.enable_code93 {
            cfg = cfg.enable(Code93);
        }
        if args.enable_code128 {
            cfg = cfg.enable(Code128);
        }
        if args.enable_i25 {
            cfg = cfg.enable(I25);
        }
        if args.enable_codabar {
            cfg = cfg.enable(Codabar);
        }
        if args.enable_databar {
            cfg = cfg.enable(Databar);
        }
        if args.enable_databar_exp {
            cfg = cfg.enable(DatabarExp);
        }
        if args.enable_qr {
            cfg = cfg.enable(QrCode);
        }
        if args.enable_sq {
            cfg = cfg.enable(SqCode);
        }
        cfg
    } else {
        // Start with default configuration (all enabled), then apply changes
        let mut cfg = DecoderConfig::new();

        // Apply enable flags (redundant but harmless)
        if args.enable_ean13 {
            cfg = cfg.enable(Ean13);
        }
        if args.enable_ean8 {
            cfg = cfg.enable(Ean8);
        }
        if args.enable_ean2 {
            cfg = cfg.enable(Ean2);
        }
        if args.enable_ean5 {
            cfg = cfg.enable(Ean5);
        }
        if args.enable_upca {
            cfg = cfg.enable(Upca);
        }
        if args.enable_upce {
            cfg = cfg.enable(Upce);
        }
        if args.enable_isbn10 {
            cfg = cfg.enable(Isbn10);
        }
        if args.enable_isbn13 {
            cfg = cfg.enable(Isbn13);
        }
        if args.enable_code39 {
            cfg = cfg.enable(Code39);
        }
        if args.enable_code93 {
            cfg = cfg.enable(Code93);
        }
        if args.enable_code128 {
            cfg = cfg.enable(Code128);
        }
        if args.enable_i25 {
            cfg = cfg.enable(I25);
        }
        if args.enable_codabar {
            cfg = cfg.enable(Codabar);
        }
        if args.enable_databar {
            cfg = cfg.enable(Databar);
        }
        if args.enable_databar_exp {
            cfg = cfg.enable(DatabarExp);
        }
        if args.enable_qr {
            cfg = cfg.enable(QrCode);
        }
        if args.enable_sq {
            cfg = cfg.enable(SqCode);
        }

        // Apply disable flags
        if args.disable_ean13 {
            cfg = cfg.disable(Ean13);
        }
        if args.disable_ean8 {
            cfg = cfg.disable(Ean8);
        }
        if args.disable_code39 {
            cfg = cfg.disable(Code39);
        }
        if args.disable_code93 {
            cfg = cfg.disable(Code93);
        }
        if args.disable_code128 {
            cfg = cfg.disable(Code128);
        }
        if args.disable_qr {
            cfg = cfg.disable(QrCode);
        }
        cfg
    };

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

        // Create zedbar image from grayscale data
        let mut zedbar_img = match Image::from_gray(data, width, height) {
            Ok(img) => img,
            Err(e) => {
                if !args.quiet {
                    eprintln!("Failed to create zedbar image: {}", e);
                }
                process::exit(1);
            }
        };

        // Create scanner with configured symbologies
        let mut scanner = Scanner::with_config(config.clone());

        // Scan the image
        let symbols = scanner.scan(&mut zedbar_img);

        total_symbols += symbols.len();

        // Print results
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
            eprintln!(
                "scanned {} barcode symbols from {} image(s)",
                total_symbols,
                args.files.len()
            );
        }
    } else if total_symbols == 0 {
        // In quiet mode, still exit with error if no barcodes found
        process::exit(1);
    }
}
