//! QR Code decoding support
//!
//! This module contains Rust implementations of QR code decoding functionality

pub mod bch15_5;

pub use bch15_5::{bch15_5_correct, bch15_5_encode};
