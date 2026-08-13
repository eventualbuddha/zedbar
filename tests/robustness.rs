//! Degenerate and adversarial inputs must not panic.
//!
//! The decoders are a port of C code that relied on unsigned wraparound and
//! on reads past the end of fixed-size arrays being harmless. In Rust the
//! equivalent slips are panics — arithmetic overflow in debug builds, slice
//! bounds anywhere — so these tests push shapes and pixel patterns that the
//! image fixtures never produce through the whole scanner.

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
