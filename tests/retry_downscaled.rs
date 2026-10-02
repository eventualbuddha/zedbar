//! A photo of a QR code on a display carries the screen's pixel grid as a
//! stripe a few pixels wide. Along a scan line that crosses the stripe every
//! module is cut into pieces, so no finder pattern is found at full
//! resolution even though the code is large and sharp. Averaging the image
//! down removes the stripe; `retry_downscaled` does that automatically.

#![cfg(all(feature = "qrcode", feature = "image"))]

use image::Luma;
use zedbar::config::*;
use zedbar::{DecoderConfig, Image, Scanner};

const PAYLOAD: &str = "https://zedbar.invalid/photographed-from-a-screen?with=some&more=data";

/// Where the QR sits in the synthetic photo, in pixels.
const QR_X: u32 = 120;
const QR_Y: u32 = 200;

/// Render the QR at 16px per module into a 1200x1600 frame and overlay a
/// four-pixel-period horizontal stripe with the levels measured in a real
/// phone photo of an LCD: black modules swing between about 30 and 110,
/// white ones between about 130 and 240.
fn screen_photo() -> (Image, u32) {
    screen_photo_sized(1200, 1600)
}

fn screen_photo_sized(w: u32, h: u32) -> (Image, u32) {
    let rendered = qrcode::QrCode::new(PAYLOAD.as_bytes())
        .expect("encode QR")
        .render::<Luma<u8>>()
        .module_dimensions(16, 16)
        .quiet_zone(false)
        .build();
    assert!(rendered.width() + QR_X < w && rendered.height() + QR_Y < h);

    let mut data = vec![255u8; (w * h) as usize];
    for (x, y, p) in rendered.enumerate_pixels() {
        data[((QR_Y + y) * w + QR_X + x) as usize] = p.0[0];
    }
    for y in 0..h {
        let (base, gain) = if y % 4 < 2 { (30.0, 0.4) } else { (110.0, 0.5) };
        for x in 0..w {
            let v = &mut data[(y * w + x) as usize];
            *v = (base + gain * *v as f32).round() as u8;
        }
    }
    let image = Image::from_gray(&data, w, h).expect("valid dimensions");
    (image, rendered.width())
}

fn decoded(config: DecoderConfig, image: &mut Image) -> Vec<zedbar::symbol::Symbol> {
    Scanner::with_config(config)
        .scan(image)
        .into_iter()
        .collect()
}

#[test]
fn stripe_defeats_full_resolution_pass() {
    let (mut image, _) = screen_photo();
    let symbols = decoded(DecoderConfig::new().enable(QrCode), &mut image);
    assert!(
        symbols.is_empty(),
        "the stripe no longer breaks finder detection; the retry test below proves nothing"
    );
}

#[test]
fn retry_downscaled_recovers_the_code_with_original_coordinates() {
    let (mut image, qr_size) = screen_photo();
    let result = Scanner::with_config(
        DecoderConfig::new()
            .enable(QrCode)
            .retry_undecoded_regions(true)
            .retry_downscaled(true),
    )
    .scan(&mut image);
    let symbols = result.symbols();
    assert_eq!(symbols.len(), 1);
    let symbol = &symbols[0];
    assert_eq!(symbol.data_string(), Some(PAYLOAD));

    // The bounding box must land on the code in the full-size frame, which
    // the half-size scan reports in its own coordinates.
    let bounds = symbol.bounds().expect("position tracking is on");
    let slack = 32i32;
    assert!(
        (bounds.x - QR_X as i32).abs() <= slack
            && (bounds.y - QR_Y as i32).abs() <= slack
            && (bounds.width as i32 - qr_size as i32).abs() <= slack
            && (bounds.height as i32 - qr_size as i32).abs() <= slack,
        "bounds {bounds:?} vs QR at ({QR_X}, {QR_Y}) size {qr_size}"
    );

    // The full-resolution pass reported the code's own area as an undecoded
    // finder region. Recovering the code resolves it.
    assert!(
        result.finder_regions().is_empty(),
        "stale regions: {:?}",
        result.finder_regions()
    );
}

/// The size gate is checked before each step, so an image whose shorter side
/// is exactly the threshold still gets the half-size pass, and one just under
/// it gets nothing.
#[test]
fn retry_downscaled_runs_at_exactly_the_threshold() {
    let config = || DecoderConfig::new().enable(QrCode).retry_downscaled(true);

    let (mut image, _) = screen_photo_sized(1024, 1400);
    assert_eq!(decoded(config(), &mut image).len(), 1, "1024px side");

    let (mut image, _) = screen_photo_sized(1023, 1400);
    assert!(decoded(config(), &mut image).is_empty(), "1023px side");
}

/// The retry only ever adds a result on large images that decoded nothing,
/// so on the fixture corpus it must leave every result unchanged.
#[test]
fn retry_downscaled_does_not_change_results_on_fixtures() {
    for path in [
        "examples/test-qr.png",
        "examples/test-qr-version40.png",
        "examples/qr-code-140-grid01.jpg",
        "examples/synthetic-small-qr-page.png",
        "examples/test-ean13.png",
    ] {
        let img = image::open(path).unwrap_or_else(|e| panic!("{path}: {e}"));
        let mut results = Vec::new();
        for retry in [false, true] {
            let mut image = Image::from_dynamic(&img).unwrap();
            let config = DecoderConfig::all()
                .retry_undecoded_regions(true)
                .retry_downscaled(retry);
            let mut data: Vec<_> = Scanner::with_config(config)
                .scan(&mut image)
                .iter()
                .map(|s| (s.symbol_type(), s.data_string().map(str::to_owned)))
                .collect();
            data.sort();
            results.push(data);
        }
        assert_eq!(results[0], results[1], "{path}");
    }
}
