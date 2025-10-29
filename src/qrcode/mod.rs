//! QR Code decoding support
//!
//! This module contains Rust implementations of QR code decoding functionality

pub mod bch15_5;
pub mod binarize;
pub mod qrdec;
pub mod qrdectxt;
pub mod util;

pub use bch15_5::bch15_5_encode;
pub use util::{qr_ihypot, qr_ilog, qr_isqrt};
