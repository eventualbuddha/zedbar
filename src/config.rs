//! Type-safe configuration system for ZBar decoders
//!
//! This module provides a compile-time verified configuration API that prevents
//! invalid configuration combinations. The type system ensures you can only
//! set configurations that are valid for each symbology.
//!
//! # Examples
//!
//! ## Valid Configuration
//!
//! ```
//! use zbar::config::*;
//! use zbar::DecoderConfig;
//!
//! let config = DecoderConfig::new()
//!     .enable(Ean13)
//!     .enable(Code39)
//!     .set_length_limits(Code39, 4, 20)   // ✓ Code39 supports variable length
//!     .set_binary(QrCode, true)           // ✓ QR codes support binary mode
//!     .position_tracking(true);
//! ```
//!
//! ## Type-Safe Compile Errors
//!
//! The following configurations will NOT compile:
//!
//! ```compile_fail
//! # use zbar::config::*;
//! # use zbar::DecoderConfig;
//! # let config = DecoderConfig::new();
//! // ❌ EAN-13 has fixed length, doesn't support length limits
//! config.set_length_limits(Ean13, 1, 20);
//! ```
//!
//! ```compile_fail
//! # use zbar::config::*;
//! # use zbar::DecoderConfig;
//! # let config = DecoderConfig::new();
//! // ❌ Code39 is not a 2D code, doesn't support binary mode
//! config.set_binary(Code39, true);
//! ```
//!
//! ```compile_fail
//! # use zbar::config::*;
//! # use zbar::DecoderConfig;
//! # let config = DecoderConfig::new();
//! // ❌ QR codes don't use length limits
//! config.set_length_limits(QrCode, 1, 100);
//! ```

use crate::SymbolType;
use std::collections::{HashMap, HashSet};

pub(crate) mod internal;
pub mod symbologies;

// Re-export symbology types for convenience
pub use symbologies::*;

// ============================================================================
// Configuration Value Types
// ============================================================================

/// Enable/disable flag for a symbology
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Enable(pub bool);

/// Add checksum validation
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AddChecksum(pub bool);

/// Emit checksum in decoded data
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EmitChecksum(pub bool);

/// Binary mode (for 2D codes like QR)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Binary(pub bool);

/// Minimum symbol length
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MinLength(pub u32);

/// Maximum symbol length
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MaxLength(pub u32);

/// Edge detection uncertainty threshold
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Uncertainty(pub u32);

/// Position tracking enable/disable
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PositionTracking(pub bool);

/// Test inverted images
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TestInverted(pub bool);

/// Horizontal scan density
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct XDensity(pub u32);

/// Vertical scan density
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct YDensity(pub u32);

// ============================================================================
// Capability Traits
// ============================================================================

/// Marker trait for symbologies that can be enabled/disabled
pub trait SupportsEnable: Symbology {}

/// Marker trait for symbologies that support checksum configuration
pub trait SupportsChecksum: Symbology {}

/// Marker trait for symbologies that support variable length limits
pub trait SupportsLengthLimits: Symbology {}

/// Marker trait for symbologies that support binary mode
pub trait SupportsBinary: Symbology {}

/// Marker trait for symbologies that support uncertainty configuration
pub trait SupportsUncertainty: Symbology {}

/// Base trait that all symbology types must implement
pub trait Symbology: Sized {
    /// The corresponding SymbolType enum value
    const TYPE: SymbolType;

    /// Human-readable name
    const NAME: &'static str;
}

// ============================================================================
// User-Facing Configuration Builder
// ============================================================================

/// Type-safe configuration builder for ZBar decoders
///
/// This builder uses the type system to ensure only valid configurations
/// can be set for each symbology type.
///
/// # Example
/// ```no_run
/// use zbar::config::*;
///
/// let config = DecoderConfig::new()
///     .enable(Ean13)
///     .enable(Code39)
///     .set_checksum(Code39, true, false)
///     .set_length_limits(Code39, 4, 20)
///     .position_tracking(true);
/// ```
#[derive(Debug, Clone)]
pub struct DecoderConfig {
    /// Which symbologies are enabled
    pub(crate) enabled: HashSet<SymbolType>,

    /// Checksum configuration: (add_check, emit_check)
    pub(crate) checksum_flags: HashMap<SymbolType, (bool, bool)>,

    /// Length limits: (min, max)
    pub(crate) length_limits: HashMap<SymbolType, (u32, u32)>,

    /// Binary mode enabled for 2D codes
    pub(crate) binary_mode: HashSet<SymbolType>,

    /// Uncertainty threshold per symbology
    pub(crate) uncertainty: HashMap<SymbolType, u32>,

    /// Global scanner configuration
    pub(crate) position_tracking: bool,
    pub(crate) test_inverted: bool,
    pub(crate) x_density: u32,
    pub(crate) y_density: u32,
}

impl Default for DecoderConfig {
    fn default() -> Self {
        Self::new()
    }
}

impl DecoderConfig {
    /// Create a new configuration with sensible defaults
    ///
    /// By default:
    /// - EAN-13, EAN-8, I25, DataBar, Codabar, Code39, Code93, Code128, QR, SQ are enabled
    /// - UPC-A, UPC-E, ISBN-10, ISBN-13 can be emitted as variants (not separately enabled)
    /// - Position tracking is enabled
    /// - Scan density is 1x1
    /// - Inverted image testing is disabled
    pub fn new() -> Self {
        let mut config = Self {
            enabled: HashSet::new(),
            checksum_flags: HashMap::new(),
            length_limits: HashMap::new(),
            binary_mode: HashSet::new(),
            uncertainty: HashMap::new(),
            position_tracking: true,
            test_inverted: false,
            x_density: 1,
            y_density: 1,
        };

        // Enable common symbologies by default
        config.enabled.insert(SymbolType::Ean13);
        config.enabled.insert(SymbolType::Ean8);
        // Note: UPC-A, UPC-E, ISBN-10, ISBN-13 are NOT enabled by default
        // They can be emitted as variants of EAN-13/EAN-8 when detected
        config.enabled.insert(SymbolType::I25);
        config.enabled.insert(SymbolType::Databar);
        config.enabled.insert(SymbolType::DatabarExp);
        config.enabled.insert(SymbolType::Codabar);
        config.enabled.insert(SymbolType::Code39);
        config.enabled.insert(SymbolType::Code93);
        config.enabled.insert(SymbolType::Code128);
        config.enabled.insert(SymbolType::QrCode);
        config.enabled.insert(SymbolType::SqCode);

        // Set default checksum behavior for EAN/UPC
        // EAN-13 and EAN-8 are enabled with emit_check
        for sym in [SymbolType::Ean13, SymbolType::Ean8] {
            config.checksum_flags.insert(sym, (false, true)); // don't add, do emit
        }

        // UPC-A, UPC-E, ISBN variants are not enabled but have emit_check
        // This allows them to be reported as variants when EAN decoding detects them
        for sym in [
            SymbolType::Upca,
            SymbolType::Upce,
            SymbolType::Isbn10,
            SymbolType::Isbn13,
        ] {
            config.checksum_flags.insert(sym, (false, true)); // don't add, do emit
        }

        // Set default length limits for variable-length symbologies
        config.length_limits.insert(SymbolType::I25, (6, 256));
        config.length_limits.insert(SymbolType::Codabar, (4, 256));
        config.length_limits.insert(SymbolType::Code39, (1, 256));

        // Set default uncertainty values
        config.uncertainty.insert(SymbolType::Codabar, 1);

        config
    }

    // ========================================================================
    // Per-Symbology Configuration
    // ========================================================================

    /// Enable a symbology
    pub fn enable<S: Symbology + SupportsEnable>(mut self, _: S) -> Self {
        self.enabled.insert(S::TYPE);
        self
    }

    /// Disable a symbology
    pub fn disable<S: Symbology + SupportsEnable>(mut self, _: S) -> Self {
        self.enabled.remove(&S::TYPE);
        self
    }

    /// Disable all symbologies
    ///
    /// Useful when you want to start with a clean slate and only enable
    /// specific symbologies.
    ///
    /// # Example
    /// ```
    /// use zbar::config::*;
    /// use zbar::DecoderConfig;
    ///
    /// let config = DecoderConfig::new()
    ///     .disable_all()
    ///     .enable(QrCode)
    ///     .enable(Code39);
    /// ```
    pub fn disable_all(mut self) -> Self {
        self.enabled.clear();
        self
    }

    /// Check if a symbology is enabled
    pub fn is_enabled(&self, sym: SymbolType) -> bool {
        self.enabled.contains(&sym)
    }

    /// Configure checksum behavior for a symbology
    ///
    /// # Arguments
    /// * `add_check` - Validate checksum during decoding
    /// * `emit_check` - Include checksum digit in decoded data
    pub fn set_checksum<S: Symbology + SupportsChecksum>(
        mut self,
        _: S,
        add_check: bool,
        emit_check: bool,
    ) -> Self {
        self.checksum_flags.insert(S::TYPE, (add_check, emit_check));
        self
    }

    /// Set minimum and maximum length limits
    ///
    /// Only valid for variable-length symbologies like Code39, Code128, etc.
    pub fn set_length_limits<S: Symbology + SupportsLengthLimits>(
        mut self,
        _: S,
        min: u32,
        max: u32,
    ) -> Self {
        assert!(min <= max, "min length must be <= max length");
        assert!(max <= 256, "max length must be <= 256");
        self.length_limits.insert(S::TYPE, (min, max));
        self
    }

    /// Enable or disable binary mode for 2D codes
    ///
    /// When enabled, binary data is preserved without text conversion.
    pub fn set_binary<S: Symbology + SupportsBinary>(mut self, _: S, enabled: bool) -> Self {
        if enabled {
            self.binary_mode.insert(S::TYPE);
        } else {
            self.binary_mode.remove(&S::TYPE);
        }
        self
    }

    /// Set uncertainty threshold for edge detection
    ///
    /// Higher values are more tolerant of poor quality images but may
    /// produce more false positives.
    pub fn set_uncertainty<S: Symbology + SupportsUncertainty>(
        mut self,
        _: S,
        threshold: u32,
    ) -> Self {
        self.uncertainty.insert(S::TYPE, threshold);
        self
    }

    // ========================================================================
    // Global Scanner Configuration
    // ========================================================================

    /// Enable or disable position tracking
    ///
    /// When enabled, the scanner records the pixel coordinates of each
    /// detected symbol.
    pub fn position_tracking(mut self, enabled: bool) -> Self {
        self.position_tracking = enabled;
        self
    }

    /// Enable or disable inverted image testing
    ///
    /// When enabled, if no symbols are found in the normal image, the
    /// scanner will try again with an inverted (negative) image.
    pub fn test_inverted(mut self, enabled: bool) -> Self {
        self.test_inverted = enabled;
        self
    }

    /// Set scan density for both axes
    ///
    /// Higher density means more scan lines, which improves detection
    /// but increases processing time. A value of 1 means scan every line.
    pub fn scan_density(mut self, x: u32, y: u32) -> Self {
        assert!(x > 0, "x density must be > 0");
        assert!(y > 0, "y density must be > 0");
        self.x_density = x;
        self.y_density = y;
        self
    }

    /// Set horizontal scan density
    pub fn x_density(mut self, density: u32) -> Self {
        assert!(density > 0, "density must be > 0");
        self.x_density = density;
        self
    }

    /// Set vertical scan density
    pub fn y_density(mut self, density: u32) -> Self {
        assert!(density > 0, "density must be > 0");
        self.y_density = density;
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_config() {
        let config = DecoderConfig::new();
        assert!(config.is_enabled(SymbolType::Ean13));
        assert!(config.is_enabled(SymbolType::Code39));
        assert!(config.position_tracking);
        assert!(!config.test_inverted);
        assert_eq!(config.x_density, 1);
        assert_eq!(config.y_density, 1);
    }

    #[test]
    fn test_builder_pattern() {
        let config = DecoderConfig::new()
            .enable(Ean13)
            .disable(Code39)
            .set_checksum(Code39, true, false)
            .position_tracking(false)
            .scan_density(2, 2);

        assert!(config.is_enabled(SymbolType::Ean13));
        assert!(!config.is_enabled(SymbolType::Code39));
        assert!(!config.position_tracking);
        assert_eq!(config.x_density, 2);
        assert_eq!(config.y_density, 2);
    }

    #[test]
    fn test_type_safe_length_limits() {
        // Variable-length symbologies can have length limits
        let config = DecoderConfig::new()
            .set_length_limits(Code39, 5, 20)
            .set_length_limits(Code128, 1, 50)
            .set_length_limits(I25, 6, 30);

        assert_eq!(
            config.length_limits.get(&SymbolType::Code39),
            Some(&(5, 20))
        );
        assert_eq!(
            config.length_limits.get(&SymbolType::Code128),
            Some(&(1, 50))
        );
        assert_eq!(config.length_limits.get(&SymbolType::I25), Some(&(6, 30)));
    }

    #[test]
    fn test_type_safe_binary_mode() {
        // Only 2D codes support binary mode
        let config = DecoderConfig::new()
            .set_binary(QrCode, true)
            .set_binary(SqCode, false);

        assert!(config.binary_mode.contains(&SymbolType::QrCode));
        assert!(!config.binary_mode.contains(&SymbolType::SqCode));
    }

    #[test]
    fn test_checksum_configuration() {
        let config = DecoderConfig::new()
            .set_checksum(Codabar, true, false) // validate but don't emit
            .set_checksum(Ean13, false, true); // don't validate, do emit

        assert_eq!(
            config.checksum_flags.get(&SymbolType::Codabar),
            Some(&(true, false))
        );
        assert_eq!(
            config.checksum_flags.get(&SymbolType::Ean13),
            Some(&(false, true))
        );
    }

    #[test]
    fn test_uncertainty_configuration() {
        let config = DecoderConfig::new()
            .set_uncertainty(Code39, 2)
            .set_uncertainty(QrCode, 0)
            .set_uncertainty(Codabar, 1);

        assert_eq!(config.uncertainty.get(&SymbolType::Code39), Some(&2));
        assert_eq!(config.uncertainty.get(&SymbolType::QrCode), Some(&0));
        assert_eq!(config.uncertainty.get(&SymbolType::Codabar), Some(&1));
    }

    #[test]
    fn test_config_to_state_conversion() {
        let config = DecoderConfig::new()
            .enable(Ean13)
            .set_checksum(Ean13, false, true)
            .position_tracking(true)
            .scan_density(2, 3);

        let state: internal::DecoderState = (&config).into();

        // Verify EAN is enabled
        assert!(state.is_enabled(SymbolType::Ean13));
        assert!(state.ean_enabled());

        // Verify checksum config
        let ean13_config = state.get(SymbolType::Ean13).unwrap();
        assert!(!ean13_config.checksum.add_check);
        assert!(ean13_config.checksum.emit_check);

        // Verify scanner config
        assert_eq!(state.scanner.x_density, 2);
        assert_eq!(state.scanner.y_density, 3);
        assert!(state.scanner.position_tracking);
    }
}
