//! Degenerate and adversarial inputs must not panic.
//!
//! The decoders are a port of C code that relied on unsigned wraparound and
//! on reads past the end of fixed-size arrays being harmless. In Rust the
//! equivalent slips are panics — arithmetic overflow in debug builds, slice
//! bounds anywhere — so these tests push shapes and pixel patterns that the
//! image fixtures never produce through the whole scanner.

#[cfg(all(
    feature = "code39",
    feature = "code93",
    feature = "code128",
    feature = "i25",
    feature = "codabar"
))]
use zedbar::config::*;
use zedbar::{DecoderConfig, Image, Scanner};

/// Deterministic noise. Avoids a dev-dependency on `rand` and keeps failures
/// reproducible.
struct Lcg(u64);

impl Lcg {
    fn next_u32(&mut self) -> u32 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 33) as u32
    }
}

fn scan(data: &[u8], width: u32, height: u32) {
    let mut image = Image::from_gray(data, width, height).expect("valid dimensions");
    let mut scanner = Scanner::new();
    // The result is irrelevant; not panicking is the assertion.
    let _ = scanner.scan(&mut image);
}

#[test]
fn scans_extreme_image_shapes() {
    // 1xN, Nx1 and 1x1: the scan passes compute borders from the image
    // dimensions and step by density, which underflows easily at these sizes.
    for (w, h) in [
        (0, 0),
        (0, 8),
        (8, 0),
        (1, 1),
        (1, 64),
        (64, 1),
        (2, 2),
        (3, 1),
        (1, 3),
        (1, 1000),
        (1000, 1),
    ] {
        let data = vec![0x80u8; (w * h) as usize];
        scan(&data, w, h);
    }
}

#[test]
fn scans_uniform_images() {
    // A uniform image has no edges at all, so every decoder sits at its
    // start state while the scanner accumulates one enormous element width.
    for fill in [0x00u8, 0x7f, 0x80, 0xff] {
        let data = vec![fill; 200 * 200];
        scan(&data, 200, 200);
    }
}

#[test]
fn scans_maximum_contrast_stripes() {
    // Alternating single-pixel stripes produce the narrowest elements the
    // scanner can see, in both orientations.
    let (w, h) = (256u32, 256u32);
    let mut vertical = vec![0u8; (w * h) as usize];
    let mut horizontal = vec![0u8; (w * h) as usize];
    for y in 0..h {
        for x in 0..w {
            let i = (y * w + x) as usize;
            vertical[i] = if x % 2 == 0 { 0 } else { 255 };
            horizontal[i] = if y % 2 == 0 { 0 } else { 255 };
        }
    }
    scan(&vertical, w, h);
    scan(&horizontal, w, h);
}

/// A single bar surrounded by a very wide quiet zone: element widths are
/// stored as fixed-point, so one element spanning a whole wide image is the
/// largest width the decoders ever multiply together.
#[test]
fn scans_very_wide_elements() {
    let (w, h) = (4000u32, 8u32);
    let mut data = vec![255u8; (w * h) as usize];
    for y in 0..h {
        for x in (w / 2)..(w / 2 + 3) {
            data[(y * w + x) as usize] = 0;
        }
    }
    scan(&data, w, h);
}

#[test]
fn scans_random_noise() {
    let mut rng = Lcg(0x1234_5678_9abc_def0);
    let (w, h) = (300u32, 300u32);

    for _ in 0..8 {
        let mut data = vec![0u8; (w * h) as usize];
        for pixel in data.iter_mut() {
            *pixel = rng.next_u32() as u8;
        }
        scan(&data, w, h);

        // Hard-thresholded noise: far more edges per scan line, which is
        // what drives the decoders deepest into their state machines.
        let bilevel: Vec<u8> = data
            .iter()
            .map(|p| if *p < 128 { 0 } else { 255 })
            .collect();
        scan(&bilevel, w, h);
    }
}

/// Length limits of (0, 0) mean "unlimited", which removes the ceiling the
/// default configuration puts on how many characters a symbol may collect.
#[test]
#[cfg(all(
    feature = "code39",
    feature = "code93",
    feature = "code128",
    feature = "i25",
    feature = "codabar"
))]
fn scans_noise_with_unlimited_length_limits() {
    let mut rng = Lcg(0x0fed_cba9_8765_4321);
    let (w, h) = (600u32, 200u32);

    let config = DecoderConfig::all()
        .set_length_limits(Code39, 0, 0)
        .set_length_limits(Code93, 0, 0)
        .set_length_limits(Code128, 0, 0)
        .set_length_limits(I25, 0, 0)
        .set_length_limits(Codabar, 0, 0);

    for _ in 0..8 {
        let data: Vec<u8> = (0..(w * h))
            .map(|_| if rng.next_u32() & 1 == 0 { 0 } else { 255 })
            .collect();
        let mut image = Image::from_gray(&data, w, h).expect("valid dimensions");
        let mut scanner = Scanner::with_config(config.clone());
        let _ = scanner.scan(&mut image);
    }
}

#[test]
fn scans_with_inversion_and_retry_enabled() {
    let mut rng = Lcg(0xdead_beef_cafe_f00d);
    let (w, h) = (400u32, 400u32);
    let data: Vec<u8> = (0..(w * h)).map(|_| rng.next_u32() as u8).collect();

    let config = DecoderConfig::all()
        .test_inverted(true)
        .retry_undecoded_regions(true);
    let mut image = Image::from_gray(&data, w, h).expect("valid dimensions");
    let mut scanner = Scanner::with_config(config);
    let _ = scanner.scan(&mut image);
}

#[test]
fn scans_at_coarse_densities() {
    let mut rng = Lcg(0x5a5a_a5a5_1234_9876);
    let (w, h) = (150u32, 150u32);
    let data: Vec<u8> = (0..(w * h)).map(|_| rng.next_u32() as u8).collect();

    // Densities larger than the image exercise the border calculation, which
    // takes `(dimension - 1) % density` and halves it.
    for density in [1u32, 2, 3, 7, 149, 150, 151, 1000] {
        let mut image = Image::from_gray(&data, w, h).expect("valid dimensions");
        let mut scanner = Scanner::with_config(DecoderConfig::all().scan_density(density, density));
        let _ = scanner.scan(&mut image);
    }
}

#[test]
fn crop_rejects_out_of_bounds_regions_without_overflowing() {
    let image = Image::from_gray(&[0u8; 64], 8, 8).expect("valid dimensions");

    // Empty, past the edge, and offsets whose extent overflows u32.
    assert!(image.crop(0, 0, 0, 4).is_none());
    assert!(image.crop(0, 0, 4, 0).is_none());
    assert!(image.crop(0, 0, 9, 8).is_none());
    assert!(image.crop(0, 0, 8, 9).is_none());
    assert!(image.crop(5, 0, 4, 4).is_none());
    assert!(image.crop(u32::MAX, 0, 4, 4).is_none());
    assert!(image.crop(0, u32::MAX, 4, 4).is_none());
    assert!(image.crop(u32::MAX - 1, u32::MAX - 1, 4, 4).is_none());

    let cropped = image.crop(2, 2, 4, 4).expect("in-bounds crop");
    assert_eq!(cropped.width(), 4);
    assert_eq!(cropped.height(), 4);
    assert_eq!(cropped.data().len(), 16);
}

#[test]
fn upscale_rejects_invalid_factors() {
    let image = Image::from_gray(&[0u8; 64], 8, 8).expect("valid dimensions");

    assert!(image.upscale(0).is_none());
    assert!(image.upscale(1).is_none());
    assert!(image.upscale(u32::MAX).is_none());

    let upscaled = image.upscale(3).expect("valid factor");
    assert_eq!(upscaled.width(), 24);
    assert_eq!(upscaled.height(), 24);
    assert_eq!(upscaled.data().len(), 24 * 24);
}

#[test]
fn from_gray_rejects_mismatched_dimensions() {
    assert!(Image::from_gray(&[0u8; 10], 4, 4).is_err());
    assert!(Image::from_gray(&[0u8; 16], 4, 4).is_ok());
    assert!(Image::from_gray(&[], 0, 0).is_ok());
    // Dimensions whose product overflows u32 must not wrap into a match.
    assert!(Image::from_gray(&[0u8; 4], u32::MAX, u32::MAX).is_err());
}

/// A barcode rendered from SVG is usually drawn in black on a transparent
/// background, and the transparent pixels commonly carry RGB (0, 0, 0).
/// Converting straight to grayscale drops the alpha and turns the whole
/// image black, so the code becomes invisible. `Image::from_dynamic`
/// composites over white first.
#[test]
#[cfg(feature = "qrcode")]
fn decodes_a_qr_on_a_transparent_background() {
    use image::{DynamicImage, Luma, Rgba, RgbaImage};

    const PAYLOAD: &str = "https://zedbar.invalid/transparent";

    let rendered = qrcode::QrCode::new(PAYLOAD.as_bytes())
        .expect("encode QR")
        .render::<Luma<u8>>()
        .module_dimensions(4, 4)
        .quiet_zone(true)
        .build();

    // Light modules become fully transparent black, as an SVG rasteriser
    // would leave them.
    let mut rgba = RgbaImage::new(rendered.width(), rendered.height());
    for (dst, src) in rgba.pixels_mut().zip(rendered.pixels()) {
        *dst = if src.0[0] < 128 {
            Rgba([0, 0, 0, 255])
        } else {
            Rgba([0, 0, 0, 0])
        };
    }
    let img = DynamicImage::ImageRgba8(rgba);

    // The naive conversion loses the code entirely: every pixel is black.
    let naive = img.to_luma8();
    assert!(
        naive.as_raw().iter().all(|&p| p < 128),
        "expected the alpha-dropping conversion to flatten the image to black"
    );

    let mut zedbar_img = Image::from_dynamic(&img).expect("convert image");
    let mut scanner = Scanner::new();
    let result = scanner.scan(&mut zedbar_img);
    let symbol = result
        .symbols()
        .first()
        .expect("transparent-background QR should decode");
    assert_eq!(symbol.data_string(), Some(PAYLOAD));
}

/// Opaque images must convert exactly as they did before.
#[test]
fn from_dynamic_matches_to_luma8_for_opaque_images() {
    use image::{DynamicImage, Rgb, RgbImage};

    let mut rgb = RgbImage::new(16, 16);
    for (i, p) in rgb.pixels_mut().enumerate() {
        *p = Rgb([(i * 7) as u8, (i * 3) as u8, (i * 11) as u8]);
    }
    let img = DynamicImage::ImageRgb8(rgb);

    let expected = img.to_luma8();
    let actual = Image::from_dynamic(&img).expect("convert image");
    assert_eq!(actual.data(), expected.as_raw());
}

/// A dense grid of bare QR finder patterns produces enormous numbers of
/// candidate finder centers that can never form a decodable code. Every
/// triplet the geometric pre-filter lets through becomes an undecoded region,
/// and the triplet search is cubic in the center count, so this used to
/// report hundreds of thousands of regions — which
/// [`DecoderConfig::retry_undecoded_regions`] then cropped, upscaled and
/// re-scanned one by one, turning a small image into an effective hang.
#[cfg(feature = "qrcode")]
#[test]
fn dense_finder_pattern_grid_stays_bounded() {
    use std::time::Instant;

    const MODULE: usize = 4;
    const CELL: usize = 10 * MODULE; // 7-module finder plus a 3-module gap
    const GRID: usize = 9;
    let dim = GRID * CELL;

    let mut data = vec![255u8; dim * dim];
    for gy in 0..GRID {
        for gx in 0..GRID {
            for y in 0..7 * MODULE {
                for x in 0..7 * MODULE {
                    let (mx, my) = (x / MODULE, y / MODULE);
                    let dark = mx == 0
                        || mx == 6
                        || my == 0
                        || my == 6
                        || (2..=4).contains(&mx) && (2..=4).contains(&my);
                    if dark {
                        data[(gy * CELL + y) * dim + gx * CELL + x] = 0;
                    }
                }
            }
        }
    }

    let config = DecoderConfig::new()
        .enable(QrCode)
        .retry_undecoded_regions(true);
    let mut image = Image::from_gray(&data, dim as u32, dim as u32).expect("valid dimensions");

    let started = Instant::now();
    let result = Scanner::with_config(config).scan(&mut image);
    let elapsed = started.elapsed();

    assert!(
        result.finder_regions().len() <= 64,
        "unbounded region count: {}",
        result.finder_regions().len()
    );
    // Generous: the point is that the work is bounded at all, not that it is
    // fast. Before the caps this did not terminate in any practical time.
    assert!(
        elapsed.as_secs() < 120,
        "scan took {elapsed:?}, expected the retry work to be bounded"
    );
}

/// A `Scanner` may be reused across images, so its result for one image must
/// not depend on what it scanned before.
///
/// The per-scan-line reset is deliberately soft — DataBar pairs segments and
/// EAN pairs halves across the scan lines of a single image — so without a
/// full reset per image, a half found in one image can pair with a half from
/// the next and report a symbol present in neither. zbar leaves this open
/// because it scans video frames behind an inter-frame cache; this API takes
/// one image at a time.
#[test]
fn scan_results_do_not_depend_on_scan_history() {
    let fixtures = [
        "examples/test-databar.png",
        "examples/test-databar-exp.png",
        "examples/test-ean13.png",
        "examples/test-ean8-plain.png",
        "examples/test-ean13-addon5.png",
        "examples/test-code128.png",
        "examples/test-qr.png",
    ];

    let loaded: Vec<(&str, Vec<u8>, u32, u32)> = fixtures
        .iter()
        .filter_map(|p| {
            let img = image::open(p).ok()?.to_luma8();
            Some((*p, img.as_raw().clone(), img.width(), img.height()))
        })
        .collect();
    assert!(!loaded.is_empty(), "no fixtures loaded");

    let decode = |scanner: &mut Scanner, (_, data, w, h): &(&str, Vec<u8>, u32, u32)| {
        let mut image = Image::from_gray(data, *w, *h).expect("valid image");
        let mut out: Vec<String> = scanner
            .scan(&mut image)
            .symbols()
            .iter()
            .map(|s| format!("{:?}:{:?}", s.symbol_type(), s.data()))
            .collect();
        out.sort();
        out
    };

    for first in &loaded {
        for second in &loaded {
            let standalone = decode(&mut Scanner::new(), second);

            let mut reused = Scanner::new();
            let _ = decode(&mut reused, first);
            let after = decode(&mut reused, second);

            assert_eq!(
                standalone, after,
                "scanning {} changed the result for {}",
                first.0, second.0
            );
        }
    }
}

/// The four corners of the code area are checked against the image bounds
/// before a homography is fitted to them, but the bottom-right corner is
/// replaced afterwards by one projected from an alignment pattern — and that
/// projection divides by a determinant that can be near-degenerate, so it can
/// land hundreds of image widths away. The cofactors built from such a corner
/// are quadratic in the corner spacing and overflow, which trapped here and
/// silently produced a meaningless transform in the C original.
///
/// This image reaches that projection on its inverted pass, and decodes once
/// the out-of-range corner is discarded in favour of the edge intersection.
#[test]
fn alignment_pattern_corner_stays_in_range() {
    let img = image::open("examples/qr-alignment-corner.png")
        .expect("fixture missing")
        .to_luma8();
    let (w, h) = (img.width(), img.height());
    let mut image = Image::from_gray(img.as_raw(), w, h).unwrap();

    let config = DecoderConfig::all().test_inverted(true);
    let result = Scanner::with_config(config).scan(&mut image);

    let data: Vec<_> = result.iter().filter_map(|s| s.data_string()).collect();
    assert_eq!(data, ["version 20 payload padding padding padding"]);
}

/// The 5x5 alignment-pattern probe offsets its template positions by the
/// distance between the pattern's projected centre and the position being
/// tried. Both come out of a projection that can collapse, so the offsets are
/// unbounded — they are meant to run off the image, where the sampler clamps
/// them, but adding them trapped first.
///
/// Fixture: a mutated version-40 code that reaches the alignment search with a
/// degenerate cell.
#[test]
fn alignment_pattern_probe_offsets_do_not_trap() {
    let img = image::open("examples/qr-alignment-fetch.png")
        .expect("fixture missing")
        .to_luma8();
    let (w, h) = (img.width(), img.height());
    let mut image = Image::from_gray(img.as_raw(), w, h).unwrap();
    // No assertion on the result: the point is that the scan completes.
    let _ = Scanner::with_config(DecoderConfig::all()).scan(&mut image);
}
