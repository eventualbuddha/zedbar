//! Barcode decoder implementations
//!
//! This module contains individual decoder implementations for various barcode formats.

#[cfg(feature = "codabar")]
pub(crate) mod codabar;
#[cfg(feature = "code128")]
pub(crate) mod code128;
#[cfg(feature = "code39")]
pub(crate) mod code39;
#[cfg(feature = "code93")]
pub(crate) mod code93;
#[cfg(feature = "databar")]
pub(crate) mod databar;
#[cfg(feature = "ean")]
pub(crate) mod ean;
#[cfg(feature = "i25")]
pub(crate) mod i25;
