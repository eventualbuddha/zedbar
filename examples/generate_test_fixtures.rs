//! Regenerate the binary fixtures used by the retry-path regression tests.
//!
//! Run via: `cargo run --example generate_test_fixtures`
//!
//! The fixtures are synthetic by design — the original sources we replaced
//! could have contained user-sensitive document content. If you need to
//! tweak the fixtures (e.g. change the QR payload or tune image size),
//! edit the constants below and re-run this binary.

use image::{GrayImage, Luma};
use qrcode::QrCode;

/// Harmless, clearly-synthetic QR payload baked into all fixtures.
const QR_PAYLOAD: &str = "https://zedbar.invalid/fixture/small-qr";

/// Small standalone QR: exercises the "region fills most of the image" case
/// in the retry path. Replaces what used to be `github-issue-qr.png`.
const SMALL_STANDALONE: &str = "examples/synthetic-small-qr.png";

/// Small QR embedded in a mostly-white large page: exercises the
/// "small QR in a big image, needs crop+upscale" case. Replaces what used
/// to be `qr-in-large-page.jpg` and `qr-failing-2.jpg`.
const SMALL_IN_LARGE_PAGE: &str = "examples/synthetic-small-qr-page.png";

fn render_qr(module_size: u32) -> GrayImage {
    let code = QrCode::new(QR_PAYLOAD.as_bytes()).expect("QR encoding failed");
    code.render::<Luma<u8>>()
        .module_dimensions(module_size, module_size)
        .quiet_zone(true)
        .build()
}

fn paste(dst: &mut GrayImage, src: &GrayImage, x: u32, y: u32) {
    for sy in 0..src.height() {
        for sx in 0..src.width() {
            dst.put_pixel(x + sx, y + sy, *src.get_pixel(sx, sy));
        }
    }
}

fn save_png(img: &GrayImage, path: &str) {
    img.save(path)
        .unwrap_or_else(|e| panic!("save {path}: {e}"));
    println!("  wrote {} ({}x{})", path, img.width(), img.height());
}

/// A simple separable 3x3 box blur. Simulates the edge softening of a
/// real scan well enough to defeat pixel-perfect direct detection,
/// without pulling in an imaging-filter crate.
fn box_blur(src: &GrayImage) -> GrayImage {
    let (w, h) = src.dimensions();
    let mut tmp = GrayImage::new(w, h);
    let pixel = |img: &GrayImage, x: i32, y: i32| -> i32 {
        let x = x.clamp(0, w as i32 - 1) as u32;
        let y = y.clamp(0, h as i32 - 1) as u32;
        img.get_pixel(x, y)[0] as i32
    };
    for y in 0..h as i32 {
        for x in 0..w as i32 {
            let sum = pixel(src, x - 1, y) + pixel(src, x, y) + pixel(src, x + 1, y);
            tmp.put_pixel(x as u32, y as u32, Luma([(sum / 3) as u8]));
        }
    }
    let mut dst = GrayImage::new(w, h);
    for y in 0..h as i32 {
        for x in 0..w as i32 {
            let sum = pixel(&tmp, x, y - 1) + pixel(&tmp, x, y) + pixel(&tmp, x, y + 1);
            dst.put_pixel(x as u32, y as u32, Luma([(sum / 3) as u8]));
        }
    }
    dst
}

fn main() {
    // Small standalone: a QR centered in a 160x160 canvas (under the
    // 200px area-filter threshold in Scanner::scan) with a light blur.
    // This fixture decodes at the top level — it's a sanity check that
    // scanning a small image with a QR filling most of the frame still
    // works end-to-end.
    let qr_small = render_qr(3);
    let canvas = 160u32;
    let mut small_standalone = GrayImage::from_pixel(canvas, canvas, Luma([255]));
    let pad_x = (canvas - qr_small.width()) / 2;
    let pad_y = (canvas - qr_small.height()) / 2;
    paste(&mut small_standalone, &qr_small, pad_x, pad_y);
    let small_standalone = box_blur(&box_blur(&small_standalone));
    save_png(&small_standalone, SMALL_STANDALONE);

    // Large page: 793x1121 (typical US-letter scan at ~100 DPI), small QR
    // in the upper right. The QR is kept small (~3px per module) and the
    // image is softened with a box blur to simulate scan edge-softening.
    // The combination of small modules + blur in a big page is enough to
    // defeat the top-level adaptive binarization, forcing the retry
    // (crop+upscale) path to succeed.
    let qr = render_qr(3);
    let mut page = GrayImage::from_pixel(793, 1121, Luma([255]));
    let qx = page.width() - qr.width() - 50;
    let qy = 180;
    paste(&mut page, &qr, qx, qy);
    let blurred = box_blur(&page);
    save_png(&blurred, SMALL_IN_LARGE_PAGE);

    println!("payload: {QR_PAYLOAD}");
}
