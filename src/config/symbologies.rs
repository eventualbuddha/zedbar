//! Symbology type definitions and capability implementations
//!
//! This module defines zero-sized types for each supported symbology
//! and implements the appropriate capability traits.

use super::*;

// ============================================================================
// Symbology Type Definitions
// ============================================================================

/// EAN-2 add-on (2-digit supplement)
#[derive(Debug, Clone, Copy)]
pub struct Ean2;

/// EAN-5 add-on (5-digit supplement)
#[derive(Debug, Clone, Copy)]
pub struct Ean5;

/// EAN-8 barcode
#[derive(Debug, Clone, Copy)]
pub struct Ean8;

/// UPC-E barcode
#[derive(Debug, Clone, Copy)]
pub struct Upce;

/// ISBN-10 (legacy book code)
#[derive(Debug, Clone, Copy)]
pub struct Isbn10;

/// UPC-A barcode
#[derive(Debug, Clone, Copy)]
pub struct Upca;

/// EAN-13 barcode
#[derive(Debug, Clone, Copy)]
pub struct Ean13;

/// ISBN-13 (book code)
#[derive(Debug, Clone, Copy)]
pub struct Isbn13;

/// Composite symbology (EAN/UPC + add-on)
#[derive(Debug, Clone, Copy)]
pub struct Composite;

/// Interleaved 2 of 5
#[derive(Debug, Clone, Copy)]
pub struct I25;

/// GS1 DataBar (RSS-14)
#[derive(Debug, Clone, Copy)]
pub struct Databar;

/// GS1 DataBar Expanded
#[derive(Debug, Clone, Copy)]
pub struct DatabarExp;

/// Codabar (used in libraries, blood banks)
#[derive(Debug, Clone, Copy)]
pub struct Codabar;

/// Code 39
#[derive(Debug, Clone, Copy)]
pub struct Code39;

/// QR Code
#[derive(Debug, Clone, Copy)]
pub struct QrCode;

/// SQ Code (experimental)
#[derive(Debug, Clone, Copy)]
pub struct SqCode;

/// Code 93
#[derive(Debug, Clone, Copy)]
pub struct Code93;

/// Code 128
#[derive(Debug, Clone, Copy)]
pub struct Code128;

// ============================================================================
// Symbology Trait Implementations
// ============================================================================

impl Symbology for Ean2 {
    const TYPE: SymbolType = SymbolType::Ean2;
    const NAME: &'static str = "EAN-2";
}

impl Symbology for Ean5 {
    const TYPE: SymbolType = SymbolType::Ean5;
    const NAME: &'static str = "EAN-5";
}

impl Symbology for Ean8 {
    const TYPE: SymbolType = SymbolType::Ean8;
    const NAME: &'static str = "EAN-8";
}

impl Symbology for Upce {
    const TYPE: SymbolType = SymbolType::Upce;
    const NAME: &'static str = "UPC-E";
}

impl Symbology for Isbn10 {
    const TYPE: SymbolType = SymbolType::Isbn10;
    const NAME: &'static str = "ISBN-10";
}

impl Symbology for Upca {
    const TYPE: SymbolType = SymbolType::Upca;
    const NAME: &'static str = "UPC-A";
}

impl Symbology for Ean13 {
    const TYPE: SymbolType = SymbolType::Ean13;
    const NAME: &'static str = "EAN-13";
}

impl Symbology for Isbn13 {
    const TYPE: SymbolType = SymbolType::Isbn13;
    const NAME: &'static str = "ISBN-13";
}

impl Symbology for Composite {
    const TYPE: SymbolType = SymbolType::Composite;
    const NAME: &'static str = "Composite";
}

impl Symbology for I25 {
    const TYPE: SymbolType = SymbolType::I25;
    const NAME: &'static str = "Interleaved 2 of 5";
}

impl Symbology for Databar {
    const TYPE: SymbolType = SymbolType::Databar;
    const NAME: &'static str = "DataBar";
}

impl Symbology for DatabarExp {
    const TYPE: SymbolType = SymbolType::DatabarExp;
    const NAME: &'static str = "DataBar Expanded";
}

impl Symbology for Codabar {
    const TYPE: SymbolType = SymbolType::Codabar;
    const NAME: &'static str = "Codabar";
}

impl Symbology for Code39 {
    const TYPE: SymbolType = SymbolType::Code39;
    const NAME: &'static str = "Code 39";
}

impl Symbology for QrCode {
    const TYPE: SymbolType = SymbolType::QrCode;
    const NAME: &'static str = "QR Code";
}

impl Symbology for SqCode {
    const TYPE: SymbolType = SymbolType::SqCode;
    const NAME: &'static str = "SQ Code";
}

impl Symbology for Code93 {
    const TYPE: SymbolType = SymbolType::Code93;
    const NAME: &'static str = "Code 93";
}

impl Symbology for Code128 {
    const TYPE: SymbolType = SymbolType::Code128;
    const NAME: &'static str = "Code 128";
}

// ============================================================================
// Capability Implementations
// ============================================================================

// SupportsEnable - all symbologies can be enabled/disabled
impl SupportsEnable for Ean2 {}
impl SupportsEnable for Ean5 {}
impl SupportsEnable for Ean8 {}
impl SupportsEnable for Upce {}
impl SupportsEnable for Isbn10 {}
impl SupportsEnable for Upca {}
impl SupportsEnable for Ean13 {}
impl SupportsEnable for Isbn13 {}
impl SupportsEnable for I25 {}
impl SupportsEnable for Databar {}
impl SupportsEnable for DatabarExp {}
impl SupportsEnable for Codabar {}
impl SupportsEnable for Code39 {}
impl SupportsEnable for QrCode {}
impl SupportsEnable for SqCode {}
impl SupportsEnable for Code93 {}
impl SupportsEnable for Code128 {}

// SupportsChecksum - symbologies with checksum configuration
impl SupportsChecksum for Ean2 {}
impl SupportsChecksum for Ean5 {}
impl SupportsChecksum for Ean8 {}
impl SupportsChecksum for Upce {}
impl SupportsChecksum for Isbn10 {}
impl SupportsChecksum for Upca {}
impl SupportsChecksum for Ean13 {}
impl SupportsChecksum for Isbn13 {}
impl SupportsChecksum for Databar {}
impl SupportsChecksum for DatabarExp {}
impl SupportsChecksum for Codabar {}
impl SupportsChecksum for Code39 {}
impl SupportsChecksum for Code93 {}
impl SupportsChecksum for Code128 {}

// SupportsLengthLimits - variable-length symbologies only
impl SupportsLengthLimits for I25 {}
impl SupportsLengthLimits for Codabar {}
impl SupportsLengthLimits for Code39 {}
impl SupportsLengthLimits for Code93 {}
impl SupportsLengthLimits for Code128 {}

// SupportsBinary - 2D codes only
impl SupportsBinary for QrCode {}
impl SupportsBinary for SqCode {}

// SupportsUncertainty - most symbologies
impl SupportsUncertainty for I25 {}
impl SupportsUncertainty for Databar {}
impl SupportsUncertainty for DatabarExp {}
impl SupportsUncertainty for Codabar {}
impl SupportsUncertainty for Code39 {}
impl SupportsUncertainty for QrCode {}
impl SupportsUncertainty for SqCode {}
impl SupportsUncertainty for Code93 {}
impl SupportsUncertainty for Code128 {}
impl SupportsUncertainty for Composite {}
