//! Internal optimized configuration storage
//!
//! This module provides an efficient, idiomatic Rust representation of configuration
//! optimized for fast runtime access without C compatibility constraints.

use super::DecoderConfig;
use crate::SymbolType;
use std::collections::HashMap;

/// Optimized internal storage for decoder configuration
///
/// Uses idiomatic Rust data structures for clarity and maintainability.
/// The configuration is immutable after creation for thread-safety.
#[derive(Debug, Clone)]
pub(crate) struct DecoderState {
    /// Per-symbology configuration indexed by SymbolType
    symbologies: HashMap<SymbolType, SymbologyConfig>,

    /// Scanner-level configuration
    pub(crate) scanner: ScannerConfig,
}

/// Configuration for a single symbology
#[derive(Debug, Clone, Default)]
pub(crate) struct SymbologyConfig {
    /// Is this symbology enabled?
    pub(crate) enabled: bool,

    /// Checksum configuration
    pub(crate) checksum: ChecksumConfig,

    /// Length limits (None if not applicable)
    pub(crate) length_limits: Option<LengthLimits>,

    /// Binary mode (for 2D codes)
    pub(crate) binary_mode: bool,

    /// Uncertainty threshold for edge detection
    pub(crate) uncertainty: u32,
}

/// Checksum configuration options
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct ChecksumConfig {
    /// Validate checksum during decoding
    pub(crate) add_check: bool,

    /// Include checksum digit in decoded output
    pub(crate) emit_check: bool,
}

/// Length limit constraints for variable-length symbologies
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct LengthLimits {
    pub(crate) min: u32,
    pub(crate) max: u32,
}

/// Scanner-level configuration
#[derive(Debug, Clone, Copy)]
pub(crate) struct ScannerConfig {
    /// Enable position tracking
    pub(crate) position_tracking: bool,

    /// Test inverted images
    pub(crate) test_inverted: bool,

    /// Horizontal scan density
    pub(crate) x_density: u32,

    /// Vertical scan density
    pub(crate) y_density: u32,
}

impl Default for ScannerConfig {
    fn default() -> Self {
        Self {
            position_tracking: true,
            test_inverted: false,
            x_density: 1,
            y_density: 1,
        }
    }
}

impl Default for DecoderState {
    fn default() -> Self {
        (&DecoderConfig::new()).into()
    }
}

impl From<&DecoderConfig> for DecoderState {
    fn from(config: &DecoderConfig) -> Self {
        let mut symbologies = HashMap::new();

        // Build configuration for each symbology
        for &sym in &config.enabled {
            let mut sym_config = SymbologyConfig {
                enabled: true,
                ..Default::default()
            };

            // Set checksum config if present
            if let Some(&(add, emit)) = config.checksum_flags.get(&sym) {
                sym_config.checksum = ChecksumConfig {
                    add_check: add,
                    emit_check: emit,
                };
            }

            // Set length limits if present
            if let Some(&(min, max)) = config.length_limits.get(&sym) {
                sym_config.length_limits = Some(LengthLimits { min, max });
            }

            // Set binary mode if present
            if config.binary_mode.contains(&sym) {
                sym_config.binary_mode = true;
            }

            // Set uncertainty if present
            if let Some(&threshold) = config.uncertainty.get(&sym) {
                sym_config.uncertainty = threshold;
            }

            symbologies.insert(sym, sym_config);
        }

        // Also add disabled symbologies that have non-default config
        for &sym in SymbolType::ALL.iter() {
            if !config.enabled.contains(&sym) {
                let mut sym_config = SymbologyConfig {
                    enabled: false,
                    ..Default::default()
                };

                let has_config = config.checksum_flags.contains_key(&sym)
                    || config.length_limits.contains_key(&sym)
                    || config.binary_mode.contains(&sym)
                    || config.uncertainty.contains_key(&sym);

                if has_config {
                    if let Some(&(add, emit)) = config.checksum_flags.get(&sym) {
                        sym_config.checksum = ChecksumConfig {
                            add_check: add,
                            emit_check: emit,
                        };
                    }

                    if let Some(&(min, max)) = config.length_limits.get(&sym) {
                        sym_config.length_limits = Some(LengthLimits { min, max });
                    }

                    if config.binary_mode.contains(&sym) {
                        sym_config.binary_mode = true;
                    }

                    if let Some(&threshold) = config.uncertainty.get(&sym) {
                        sym_config.uncertainty = threshold;
                    }

                    symbologies.insert(sym, sym_config);
                }
            }
        }

        Self {
            symbologies,
            scanner: ScannerConfig {
                position_tracking: config.position_tracking,
                test_inverted: config.test_inverted,
                x_density: config.x_density,
                y_density: config.y_density,
            },
        }
    }
}

impl DecoderState {
    /// Check if a symbology is enabled
    pub(crate) fn is_enabled(&self, sym: SymbolType) -> bool {
        self.symbologies
            .get(&sym)
            .map(|c| c.enabled)
            .unwrap_or(false)
    }

    /// Get configuration for a symbology
    pub(crate) fn get(&self, sym: SymbolType) -> Option<&SymbologyConfig> {
        self.symbologies.get(&sym)
    }

    /// Check if any EAN/UPC symbology is enabled
    pub(crate) fn ean_enabled(&self) -> bool {
        [
            SymbolType::Ean2,
            SymbolType::Ean5,
            SymbolType::Ean8,
            SymbolType::Ean13,
            SymbolType::Upca,
            SymbolType::Upce,
            SymbolType::Isbn10,
            SymbolType::Isbn13,
        ]
        .iter()
        .any(|&sym| self.is_enabled(sym))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_state() {
        let state = DecoderState::default();
        assert!(state.is_enabled(SymbolType::Ean13));
        assert!(state.is_enabled(SymbolType::Code39));
        assert!(state.scanner.position_tracking);
        assert_eq!(state.scanner.x_density, 1);
        assert_eq!(state.scanner.y_density, 1);
    }

    #[test]
    fn test_symbology_config() {
        use crate::config::{Code39, DecoderConfig};

        let config = DecoderConfig::new()
            .enable(Code39)
            .set_length_limits(Code39, 4, 20)
            .set_checksum(Code39, true, false)
            .set_uncertainty(Code39, 2);

        let state: DecoderState = (&config).into();

        let code39_config = state.get(SymbolType::Code39).unwrap();
        assert!(code39_config.enabled);
        assert_eq!(
            code39_config.length_limits,
            Some(LengthLimits { min: 4, max: 20 })
        );
        assert!(code39_config.checksum.add_check);
        assert!(!code39_config.checksum.emit_check);
        assert_eq!(code39_config.uncertainty, 2);
    }

    #[test]
    fn test_ean_enabled() {
        // Default config has EAN enabled
        let state = DecoderState::default();
        assert!(state.ean_enabled());

        // Manually create state with only non-EAN symbologies enabled
        let mut symbologies = std::collections::HashMap::new();
        symbologies.insert(
            SymbolType::Code39,
            SymbologyConfig {
                enabled: true,
                ..Default::default()
            },
        );
        symbologies.insert(
            SymbolType::Code128,
            SymbologyConfig {
                enabled: true,
                ..Default::default()
            },
        );

        let state = DecoderState {
            symbologies,
            scanner: ScannerConfig::default(),
        };
        assert!(!state.ean_enabled());

        // State with one EAN symbology enabled
        let mut symbologies = std::collections::HashMap::new();
        symbologies.insert(
            SymbolType::Ean13,
            SymbologyConfig {
                enabled: true,
                ..Default::default()
            },
        );

        let state = DecoderState {
            symbologies,
            scanner: ScannerConfig::default(),
        };
        assert!(state.ean_enabled());
    }

    #[test]
    fn test_scanner_config() {
        let config = DecoderConfig::new()
            .position_tracking(false)
            .test_inverted(true)
            .scan_density(2, 3);

        let state: DecoderState = (&config).into();

        assert!(!state.scanner.position_tracking);
        assert!(state.scanner.test_inverted);
        assert_eq!(state.scanner.x_density, 2);
        assert_eq!(state.scanner.y_density, 3);
    }
}
