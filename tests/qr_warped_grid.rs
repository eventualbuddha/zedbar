//! A screenshot resampled to under three pixels per module can carry a
//! smooth warp of half a module that no finder-anchored homography can
//! express, so plain grid sampling reads too many modules wrong for
//! Reed-Solomon to recover (eventualbuddha/zedbar#58). The decoder measures
//! the warp from the timing patterns and re-samples.

#![cfg(all(feature = "qrcode", feature = "image"))]

use zedbar::config::*;
use zedbar::{DecoderConfig, Image, Scanner};

#[test]
fn warped_grid_fixture_decodes() {
    let img = image::open("examples/qr-code-warped-grid.png").expect("fixture");
    let mut image = Image::from_dynamic(&img).unwrap();
    let result = Scanner::with_config(DecoderConfig::new().enable(QrCode)).scan(&mut image);
    let symbols = result.symbols();
    assert_eq!(symbols.len(), 1);
    assert_eq!(
        symbols[0].data_string(),
        Some("S;1;019f1245-392a-70f0-b274-f655a7af1254;A")
    );
}
