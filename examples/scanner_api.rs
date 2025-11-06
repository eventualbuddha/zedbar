//! Example demonstrating the modern Scanner API with type-safe configuration
//!
//! This example shows how to use the Scanner with the type-safe DecoderConfig
//! builder API to configure barcode scanning.

use zedbar::config::*;
use zedbar::{DecoderConfig, Scanner};

fn main() {
    println!("=== Zedbar Modern Scanner API Examples ===\n");

    // Example 1: Basic usage with default configuration
    println!("1. Basic scanner with defaults:");
    let _scanner = Scanner::new();
    println!("   Created scanner with default configuration\n");

    // Example 2: Custom configuration with specific symbologies
    println!("2. Scanner with custom symbology selection:");
    let config = DecoderConfig::new()
        .enable(Ean13)
        .enable(Ean8)
        .enable(QrCode)
        .disable(Code39)
        .disable(Code128);

    let _scanner = Scanner::with_config(config);
    println!("   Enabled: EAN-13, EAN-8, QR Code");
    println!("   Disabled: Code39, Code128\n");

    // Example 3: Configure symbology-specific settings
    println!("3. Scanner with symbology-specific settings:");
    let config = DecoderConfig::new()
        .enable(Code39)
        .set_length_limits(Code39, 4, 20) // Code39 must be 4-20 characters
        .set_checksum(Code39, true, false) // Validate checksum, don't emit
        .set_uncertainty(Code39, 2); // More tolerant edge detection

    let _scanner = Scanner::with_config(config);
    println!("   Code39: length 4-20, checksum validation, uncertainty=2\n");

    // Example 4: Configure QR code for binary data
    println!("4. Scanner optimized for QR codes with binary data:");
    let config = DecoderConfig::new()
        .enable(QrCode)
        .set_binary(QrCode, true) // Preserve binary data
        .disable(Ean13) // Only scan QR codes
        .disable(Code39);

    let _scanner = Scanner::with_config(config);
    println!("   QR Code only, binary mode enabled\n");

    // Example 5: Scanner-level configuration
    println!("5. Scanner with global scan settings:");
    let config = DecoderConfig::new()
        .position_tracking(true) // Track symbol positions
        .test_inverted(true) // Try inverted image if no symbols found
        .scan_density(2, 2); // Scan every 2nd line (faster, less accurate)

    let _scanner = Scanner::with_config(config);
    println!("   Position tracking: ON");
    println!("   Test inverted: ON");
    println!("   Scan density: 2x2\n");

    // Example 6: High-speed scanning configuration
    println!("6. High-speed scanner (optimized for speed):");
    let config = DecoderConfig::new()
        .enable(Ean13)
        .enable(Ean8)
        .scan_density(3, 3) // Scan every 3rd line
        .position_tracking(false); // Skip position calculation

    let _scanner = Scanner::with_config(config);
    println!("   Reduced scan density and position tracking for speed\n");

    // Example 7: Maximum quality scanning configuration
    println!("7. Maximum quality scanner (optimized for accuracy):");
    let config = DecoderConfig::new()
        .scan_density(1, 1) // Scan every line
        .position_tracking(true)
        .test_inverted(true) // Try both normal and inverted
        .set_uncertainty(Code39, 3); // More tolerant for poor quality images

    let _scanner = Scanner::with_config(config);
    println!("   Maximum scan density and inverted testing for accuracy\n");

    // Example 8: Retail scanning configuration
    println!("8. Retail scanner configuration:");
    let config = DecoderConfig::new()
        .enable(Ean13)
        .enable(Ean8)
        .enable(Upca)
        .enable(Upce)
        .set_checksum(Ean13, false, true) // Emit checksum digit
        .position_tracking(true);

    let _scanner = Scanner::with_config(config);
    println!("   EAN and UPC codes only, checksums included in output\n");

    // Example 9: Library/ISBN scanning configuration
    println!("9. Library/ISBN scanner configuration:");
    let config = DecoderConfig::new()
        .enable(Isbn10)
        .enable(Isbn13)
        .enable(Ean13) // ISBN is based on EAN-13
        .set_checksum(Isbn13, true, true); // Validate and emit

    let _scanner = Scanner::with_config(config);
    println!("   ISBN-10, ISBN-13 support with checksum validation\n");

    // Example 10: Industrial/logistics configuration
    println!("10. Industrial/logistics scanner configuration:");
    let config = DecoderConfig::new()
        .enable(Code128)
        .enable(Code39)
        .enable(QrCode)
        .enable(Databar)
        .set_length_limits(Code128, 1, 256) // Accept any length
        .set_binary(QrCode, true); // Binary data in QR codes

    let _scanner = Scanner::with_config(config);
    println!("    Code128, Code39, QR, DataBar with flexible settings\n");

    println!("=== Type Safety Demonstrations ===\n");

    println!("The following configurations are INVALID and won't compile:\n");
    println!("  // config.set_length_limits(Ean13, 1, 20);");
    println!("  // ❌ Compile error: EAN-13 has fixed length!\n");

    println!("  // config.set_binary(Code39, true);");
    println!("  // ❌ Compile error: Code39 is not a 2D code!\n");

    println!("  // config.set_length_limits(QrCode, 1, 100);");
    println!("  // ❌ Compile error: QR codes don't use length limits!\n");

    println!("=== Example Complete ===");
}
