//! A scan or photocopy of a printed QR code renders each dark module as a
//! halftone of specks and gaps a pixel or two across. At full resolution the
//! gaps split every finder-pattern run, so the code is located but never
//! decoded. A 3x3 Gaussian closes the gaps; `retry_smoothed` applies it and
//! re-scans automatically.

#![cfg(feature = "qrcode")]

use zedbar::config::*;
use zedbar::{DecoderConfig, Image, Scanner};

/// Two scans from the same batch of printed labels, one with the code about
/// 200px across and one about 100px, and a 112px render whose modules are
/// under three pixels wide (eventualbuddha/zedbar#60). Smoothing closes
/// the halftone gaps in the first two and evens out the aliased module
/// edges in the third.
const FIXTURES: &[(&str, &str)] = &[
    ("examples/qr-code-halftone-print.jpg", "22931-10-407766"),
    (
        "examples/qr-code-halftone-print-small.jpg",
        "22852-10-407707",
    ),
    (
        "examples/qr-code-tiny-modules.png",
        "10103GGNTBPUXTTKQGCTYD9QT4M0D",
    ),
];

fn load(path: &str) -> Image {
    let img = image::open(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    Image::from_dynamic(&img).unwrap()
}

fn scan(config: DecoderConfig, image: &mut Image) -> zedbar::ScanResult {
    Scanner::with_config(config).scan(image)
}

#[test]
fn halftone_defeats_the_other_passes() {
    for (path, _) in FIXTURES {
        let result = scan(
            DecoderConfig::new()
                .enable(QrCode)
                .retry_undecoded_regions(true)
                .retry_downscaled(true),
            &mut load(path),
        );
        assert!(
            result.symbols().is_empty(),
            "{path}: decodes without smoothing; the retry test below proves nothing"
        );
        assert!(
            !result.finder_regions().is_empty(),
            "{path}: the finder patterns are no longer even located"
        );
    }
}

#[test]
fn retry_smoothed_recovers_the_code_and_resolves_its_region() {
    for (path, payload) in FIXTURES {
        let result = scan(
            DecoderConfig::new()
                .enable(QrCode)
                .retry_undecoded_regions(true)
                .retry_downscaled(true)
                .retry_smoothed(true),
            &mut load(path),
        );
        let symbols = result.symbols();
        assert_eq!(symbols.len(), 1, "{path}");
        assert_eq!(symbols[0].data_string(), Some(*payload), "{path}");
        assert!(
            result.finder_regions().is_empty(),
            "{path}: stale regions: {:?}",
            result.finder_regions()
        );
    }
}

#[test]
fn manual_smooth_decodes_too() {
    for (path, payload) in FIXTURES {
        let image = load(path);
        let mut smoothed = image.smooth();
        assert_eq!(
            (smoothed.width(), smoothed.height()),
            (image.width(), image.height())
        );
        let result = scan(DecoderConfig::new().enable(QrCode), &mut smoothed);
        assert_eq!(result.symbols().len(), 1, "{path}");
        assert_eq!(result.symbols()[0].data_string(), Some(*payload), "{path}");
    }
}

/// The retry only ever adds a result to an image that decoded no 2D code,
/// so on the fixture corpus it must leave every result unchanged.
#[test]
fn retry_smoothed_does_not_change_results_on_fixtures() {
    for path in [
        "examples/test-qr.png",
        "examples/test-qr-version40.png",
        "examples/qr-code-140-grid01.jpg",
        "examples/synthetic-small-qr-page.png",
        "examples/test-ean13.png",
        "examples/nine-barcodes.png",
    ] {
        let img = image::open(path).unwrap_or_else(|e| panic!("{path}: {e}"));
        let mut results = Vec::new();
        for retry in [false, true] {
            let mut image = Image::from_dynamic(&img).unwrap();
            let config = DecoderConfig::all()
                .retry_undecoded_regions(true)
                .retry_smoothed(retry);
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
