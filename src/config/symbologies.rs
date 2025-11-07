//! Symbology type definitions and capability implementations
//!
//! This module defines zero-sized types for each supported symbology
//! and implements the appropriate capability traits.

use super::*;

// ============================================================================
// Symbology Type Definitions
// ============================================================================

// EAN family (feature = "ean")
#[cfg(feature = "ean")]
/// EAN-2 add-on (2-digit supplement)
#[derive(Debug, Clone, Copy)]
pub struct Ean2;

#[cfg(feature = "ean")]
/// EAN-5 add-on (5-digit supplement)
#[derive(Debug, Clone, Copy)]
pub struct Ean5;

#[cfg(feature = "ean")]
/// EAN-8 barcode
#[derive(Debug, Clone, Copy)]
pub struct Ean8;

#[cfg(feature = "ean")]
/// UPC-E barcode
#[derive(Debug, Clone, Copy)]
pub struct Upce;

#[cfg(feature = "ean")]
/// ISBN-10 (legacy book code)
#[derive(Debug, Clone, Copy)]
pub struct Isbn10;

#[cfg(feature = "ean")]
/// UPC-A barcode
#[derive(Debug, Clone, Copy)]
pub struct Upca;

#[cfg(feature = "ean")]
/// EAN-13 barcode
#[derive(Debug, Clone, Copy)]
pub struct Ean13;

#[cfg(feature = "ean")]
/// ISBN-13 (book code)
#[derive(Debug, Clone, Copy)]
pub struct Isbn13;

#[cfg(feature = "ean")]
/// Composite symbology (EAN/UPC + add-on)
#[derive(Debug, Clone, Copy)]
pub struct Composite;

// Other symbologies
#[cfg(feature = "i25")]
/// Interleaved 2 of 5
#[derive(Debug, Clone, Copy)]
pub struct I25;

#[cfg(feature = "databar")]
/// GS1 DataBar (RSS-14)
#[derive(Debug, Clone, Copy)]
pub struct Databar;

#[cfg(feature = "databar")]
/// GS1 DataBar Expanded
#[derive(Debug, Clone, Copy)]
pub struct DatabarExp;

#[cfg(feature = "codabar")]
/// Codabar (used in libraries, blood banks)
#[derive(Debug, Clone, Copy)]
pub struct Codabar;

#[cfg(feature = "code39")]
/// Code 39
#[derive(Debug, Clone, Copy)]
pub struct Code39;

#[cfg(feature = "qrcode")]
/// QR Code
#[derive(Debug, Clone, Copy)]
pub struct QrCode;

#[cfg(feature = "sqcode")]
/// SQ Code (experimental)
#[derive(Debug, Clone, Copy)]
pub struct SqCode;

#[cfg(feature = "code93")]
/// Code 93
#[derive(Debug, Clone, Copy)]
pub struct Code93;

#[cfg(feature = "code128")]
/// Code 128
#[derive(Debug, Clone, Copy)]
pub struct Code128;

// ============================================================================
// Symbology Trait Implementations
// ============================================================================

#[cfg(feature = "ean")]
impl Symbology for Ean2 {
    const TYPE: SymbolType = SymbolType::Ean2;
    const NAME: &'static str = "EAN-2";
}

#[cfg(feature = "ean")]
impl Symbology for Ean5 {
    const TYPE: SymbolType = SymbolType::Ean5;
    const NAME: &'static str = "EAN-5";
}

#[cfg(feature = "ean")]
impl Symbology for Ean8 {
    const TYPE: SymbolType = SymbolType::Ean8;
    const NAME: &'static str = "EAN-8";
}

#[cfg(feature = "ean")]
impl Symbology for Upce {
    const TYPE: SymbolType = SymbolType::Upce;
    const NAME: &'static str = "UPC-E";
}

#[cfg(feature = "ean")]
impl Symbology for Isbn10 {
    const TYPE: SymbolType = SymbolType::Isbn10;
    const NAME: &'static str = "ISBN-10";
}

#[cfg(feature = "ean")]
impl Symbology for Upca {
    const TYPE: SymbolType = SymbolType::Upca;
    const NAME: &'static str = "UPC-A";
}

#[cfg(feature = "ean")]
impl Symbology for Ean13 {
    const TYPE: SymbolType = SymbolType::Ean13;
    const NAME: &'static str = "EAN-13";
}

#[cfg(feature = "ean")]
impl Symbology for Isbn13 {
    const TYPE: SymbolType = SymbolType::Isbn13;
    const NAME: &'static str = "ISBN-13";
}

#[cfg(feature = "ean")]
impl Symbology for Composite {
    const TYPE: SymbolType = SymbolType::Composite;
    const NAME: &'static str = "Composite";
}

#[cfg(feature = "i25")]
impl Symbology for I25 {
    const TYPE: SymbolType = SymbolType::I25;
    const NAME: &'static str = "Interleaved 2 of 5";
}

#[cfg(feature = "databar")]
impl Symbology for Databar {
    const TYPE: SymbolType = SymbolType::Databar;
    const NAME: &'static str = "DataBar";
}

#[cfg(feature = "databar")]
impl Symbology for DatabarExp {
    const TYPE: SymbolType = SymbolType::DatabarExp;
    const NAME: &'static str = "DataBar Expanded";
}

#[cfg(feature = "codabar")]
impl Symbology for Codabar {
    const TYPE: SymbolType = SymbolType::Codabar;
    const NAME: &'static str = "Codabar";
}

#[cfg(feature = "code39")]
impl Symbology for Code39 {
    const TYPE: SymbolType = SymbolType::Code39;
    const NAME: &'static str = "Code 39";
}

#[cfg(feature = "qrcode")]
impl Symbology for QrCode {
    const TYPE: SymbolType = SymbolType::QrCode;
    const NAME: &'static str = "QR Code";
}

#[cfg(feature = "sqcode")]
impl Symbology for SqCode {
    const TYPE: SymbolType = SymbolType::SqCode;
    const NAME: &'static str = "SQ Code";
}

#[cfg(feature = "code93")]
impl Symbology for Code93 {
    const TYPE: SymbolType = SymbolType::Code93;
    const NAME: &'static str = "Code 93";
}

#[cfg(feature = "code128")]
impl Symbology for Code128 {
    const TYPE: SymbolType = SymbolType::Code128;
    const NAME: &'static str = "Code 128";
}

// ============================================================================
// Capability Implementations
// ============================================================================

// SupportsEnable - all symbologies can be enabled/disabled
#[cfg(feature = "ean")]
impl SupportsEnable for Ean2 {}
#[cfg(feature = "ean")]
impl SupportsEnable for Ean5 {}
#[cfg(feature = "ean")]
impl SupportsEnable for Ean8 {}
#[cfg(feature = "ean")]
impl SupportsEnable for Upce {}
#[cfg(feature = "ean")]
impl SupportsEnable for Isbn10 {}
#[cfg(feature = "ean")]
impl SupportsEnable for Upca {}
#[cfg(feature = "ean")]
impl SupportsEnable for Ean13 {}
#[cfg(feature = "ean")]
impl SupportsEnable for Isbn13 {}
#[cfg(feature = "i25")]
impl SupportsEnable for I25 {}
#[cfg(feature = "databar")]
impl SupportsEnable for Databar {}
#[cfg(feature = "databar")]
impl SupportsEnable for DatabarExp {}
#[cfg(feature = "codabar")]
impl SupportsEnable for Codabar {}
#[cfg(feature = "code39")]
impl SupportsEnable for Code39 {}
#[cfg(feature = "qrcode")]
impl SupportsEnable for QrCode {}
#[cfg(feature = "sqcode")]
impl SupportsEnable for SqCode {}
#[cfg(feature = "code93")]
impl SupportsEnable for Code93 {}
#[cfg(feature = "code128")]
impl SupportsEnable for Code128 {}

// SupportsChecksum - symbologies with checksum configuration
#[cfg(feature = "ean")]
impl SupportsChecksum for Ean2 {}
#[cfg(feature = "ean")]
impl SupportsChecksum for Ean5 {}
#[cfg(feature = "ean")]
impl SupportsChecksum for Ean8 {}
#[cfg(feature = "ean")]
impl SupportsChecksum for Upce {}
#[cfg(feature = "ean")]
impl SupportsChecksum for Isbn10 {}
#[cfg(feature = "ean")]
impl SupportsChecksum for Upca {}
#[cfg(feature = "ean")]
impl SupportsChecksum for Ean13 {}
#[cfg(feature = "ean")]
impl SupportsChecksum for Isbn13 {}
#[cfg(feature = "databar")]
impl SupportsChecksum for Databar {}
#[cfg(feature = "databar")]
impl SupportsChecksum for DatabarExp {}
#[cfg(feature = "codabar")]
impl SupportsChecksum for Codabar {}
#[cfg(feature = "code39")]
impl SupportsChecksum for Code39 {}
#[cfg(feature = "code93")]
impl SupportsChecksum for Code93 {}
#[cfg(feature = "code128")]
impl SupportsChecksum for Code128 {}

// SupportsLengthLimits - variable-length symbologies only
#[cfg(feature = "i25")]
impl SupportsLengthLimits for I25 {}
#[cfg(feature = "codabar")]
impl SupportsLengthLimits for Codabar {}
#[cfg(feature = "code39")]
impl SupportsLengthLimits for Code39 {}
#[cfg(feature = "code93")]
impl SupportsLengthLimits for Code93 {}
#[cfg(feature = "code128")]
impl SupportsLengthLimits for Code128 {}

// SupportsBinary - 2D codes only
#[cfg(feature = "qrcode")]
impl SupportsBinary for QrCode {}
#[cfg(feature = "sqcode")]
impl SupportsBinary for SqCode {}

// SupportsUncertainty - most symbologies
#[cfg(feature = "i25")]
impl SupportsUncertainty for I25 {}
#[cfg(feature = "databar")]
impl SupportsUncertainty for Databar {}
#[cfg(feature = "databar")]
impl SupportsUncertainty for DatabarExp {}
#[cfg(feature = "codabar")]
impl SupportsUncertainty for Codabar {}
#[cfg(feature = "code39")]
impl SupportsUncertainty for Code39 {}
#[cfg(feature = "qrcode")]
impl SupportsUncertainty for QrCode {}
#[cfg(feature = "sqcode")]
impl SupportsUncertainty for SqCode {}
#[cfg(feature = "code93")]
impl SupportsUncertainty for Code93 {}
#[cfg(feature = "code128")]
impl SupportsUncertainty for Code128 {}
#[cfg(feature = "ean")]
impl SupportsUncertainty for Composite {}
