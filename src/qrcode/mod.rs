//! QR Code decoding support
//!
//! This module contains Rust implementations of QR code decoding functionality

pub mod bch15_5;
pub mod binarize;
pub mod isaac;
pub mod qrdec;
pub mod qrdectxt;
pub mod rs;
pub mod util;

pub use bch15_5::{bch15_5_correct, bch15_5_encode};
pub use binarize::qr_binarize;
pub use isaac::{isaac_init, isaac_next_uint, isaac_next_uint32, IsaacCtx};
pub use qrdec::{
    bch18_6_correct, qr_hamming_dist, qr_line_eval, qr_point_ccw, qr_point_distance2,
    qr_point_translate, QrLine, QrPoint,
};
pub use util::{qr_ihypot, qr_ilog, qr_isqrt};
