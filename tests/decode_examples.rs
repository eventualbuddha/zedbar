//! Integration tests for decoding example barcode images
//!
//! These tests compare results from multiple decoders:
//! - This crate (zedbar)
//! - System zbar (via zbarimg command)
//! - rqrr (for QR codes only)

use image::{DynamicImage, GenericImageView};
use std::path::Path;
use zedbar::config::{DecoderConfig, Upca, Upce};
use zedbar::{Image, Scanner, SymbolType};

/// Downscale an image if it exceeds the maximum dimension
fn downscale_if_needed(img: DynamicImage, max_dimension: u32) -> DynamicImage {
    let (width, height) = img.dimensions();

    if width > max_dimension || height > max_dimension {
        let scale = (max_dimension as f32 / width.max(height) as f32).min(1.0);
        let new_width = (width as f32 * scale) as u32;
        let new_height = (height as f32 * scale) as u32;
        img.resize(new_width, new_height, image::imageops::FilterType::Triangle)
    } else {
        img
    }
}

/// Helper function to decode an image file and return the first symbol found
fn decode_image(path: &str) -> Option<(String, String)> {
    let path = Path::new(path);
    if !path.exists() {
        return None;
    }

    let img = image::open(path).ok()?;
    let img = downscale_if_needed(img, 1280).to_luma8();

    let mut scanner = Scanner::new();

    // Create a zbar Image from the image buffer
    let mut zbar_image = Image::from_gray(img.as_raw(), img.width(), img.height()).ok()?;

    // Scan the image
    let symbols = scanner.scan(&mut zbar_image);

    symbols.first().map(|symbol| {
        let symbol_type = symbol.symbol_type().to_string();
        let data = String::from_utf8_lossy(symbol.data())
            .trim_end_matches('\0') // Remove null terminators
            .to_string();
        (symbol_type, data)
    })
}

/// Whether the system `zbarimg` command is available.
///
/// The cross-checks against the reference C implementation only run when it
/// is installed (CI installs `zbar-tools`); otherwise they are skipped, so a
/// developer without it still gets zedbar's own assertions. Install it with
/// e.g. `apt install zbar-tools` or `brew install zbar` to run them locally.
fn zbarimg_available() -> bool {
    use std::process::Command;
    use std::sync::OnceLock;

    static AVAILABLE: OnceLock<bool> = OnceLock::new();
    *AVAILABLE.get_or_init(|| {
        let found = Command::new("zbarimg")
            .arg("--version")
            .output()
            .is_ok_and(|out| out.status.success());
        if !found {
            eprintln!("note: zbarimg not found; skipping reference cross-checks");
        }
        found
    })
}

/// Assert that zedbar and the reference `zbarimg` decode `path` the same
/// way. Does nothing when `zbarimg` is not installed.
fn assert_matches_zbars(path: &str, expected: &Option<(String, String)>) {
    if !zbarimg_available() {
        return;
    }
    assert_eq!(
        &decode_with_zbars(path),
        expected,
        "zbars failed for {path}"
    );
}

/// Helper function to decode using zbar-rs binary (system zbarimg command)
/// This avoids dependency on the zbars crate which may have build issues
fn decode_with_zbars(path: &str) -> Option<(String, String)> {
    use std::process::Command;

    let path = Path::new(path);
    if !path.exists() {
        return None;
    }

    // Try to use the system zbarimg command (without --raw to get symbol type)
    let output = Command::new("zbarimg")
        .arg("--quiet")
        .arg(path)
        .output()
        .ok()?;

    if !output.status.success() {
        return None;
    }

    // Strip the final newline added by zbarimg, but preserve newlines in the data
    let mut stdout = output.stdout;
    if stdout.last() == Some(&b'\n') {
        stdout.pop();
    }

    let stdout_str = String::from_utf8_lossy(&stdout);

    // Parse the output format: "SYMBOLTYPE:data"
    let mut parts = stdout_str.splitn(2, ':');
    let symbol_type = parts.next()?.to_string();
    let data = parts.next().unwrap_or("").to_string();

    Some((symbol_type, data))
}

/// Helper function to decode QR codes using rqrr
fn decode_with_rqrr(path: &str) -> Option<(String, String)> {
    let path = Path::new(path);
    if !path.exists() {
        return None;
    }

    let img = image::open(path).ok()?;
    let img = downscale_if_needed(img, 1280).to_luma8();

    let mut prepared_img = rqrr::PreparedImage::prepare(img);
    let grids = prepared_img.detect_grids();

    grids.first().and_then(|grid| {
        let (_meta, content) = grid.decode().ok()?;
        Some(("QR-Code".to_string(), content))
    })
}

#[test]
fn test_qr_simple() {
    let expected = Some((
        "QR-Code".to_string(),
        "Hello, simplified zbar!\n".to_string(),
    ));

    let result_this = decode_image("examples/test-qr.png");
    let result_rqrr = decode_with_rqrr("examples/test-qr.png");

    assert_eq!(result_this, expected, "zedbar failed");
    assert_matches_zbars("examples/test-qr.png", &expected);

    // rqrr should match if it succeeds (it may not decode all QR codes)
    if let Some(rqrr_result) = result_rqrr {
        assert_eq!(
            Some(rqrr_result),
            expected,
            "rqrr decoded but gave different result"
        );
    }
}

#[test]
fn test_qr_jpg() {
    let expected = Some((
        "QR-Code".to_string(),
        "Hello, simplified zbar!\n".to_string(),
    ));

    let result_this = decode_image("examples/test-qr.jpg");
    let result_rqrr = decode_with_rqrr("examples/test-qr.jpg");

    assert_eq!(result_this, expected, "zedbar failed");
    assert_matches_zbars("examples/test-qr.jpg", &expected);

    if let Some(rqrr_result) = result_rqrr {
        assert_eq!(
            Some(rqrr_result),
            expected,
            "rqrr decoded but gave different result"
        );
    }
}

#[test]
fn test_qr_webp() {
    let expected = Some((
        "QR-Code".to_string(),
        "Hello, simplified zbar!\n".to_string(),
    ));

    let result_this = decode_image("examples/test-qr.webp");
    let result_rqrr = decode_with_rqrr("examples/test-qr.webp");

    assert_eq!(result_this, expected, "zedbar failed");
    assert_matches_zbars("examples/test-qr.webp", &expected);

    if let Some(rqrr_result) = result_rqrr {
        assert_eq!(
            Some(rqrr_result),
            expected,
            "rqrr decoded but gave different result"
        );
    }
}

/// A version-40 (177x177 module) code, the largest QR defines.
///
/// Regression test for the corner bounds check in `qr_hom_fit`: the C source
/// writes the limit as `_width << QR_FINDER_SUBPREC + 1`, which shifts by
/// SUBPREC+1. Reading it as `(width << SUBPREC) + 1` halved the allowance, and
/// codes filling the frame — everything from version 36 up — had their fitted
/// corners rejected and never decoded at all.
#[test]
fn test_qr_version_40() {
    let expected = Some((
        "QR-Code".to_string(),
        "VERSION 40 QR CODE TEST PAYLOAD".to_string(),
    ));

    let result_this = decode_image("examples/test-qr-version40.png");

    assert_eq!(result_this, expected, "zedbar failed");
    assert_matches_zbars("examples/test-qr-version40.png", &expected);
}

#[test]
fn test_qr_wifi_sharing() {
    let expected = Some((
        "QR-Code".to_string(),
        "WIFI:S:Not a real network;T:SAE;P:password;H:false;;".to_string(),
    ));

    let result_this = decode_image("examples/pixel-wifi-sharing-qr-code.png");
    let result_rqrr = decode_with_rqrr("examples/pixel-wifi-sharing-qr-code.png");

    assert_eq!(result_this, expected, "zedbar failed");
    assert_matches_zbars("examples/pixel-wifi-sharing-qr-code.png", &expected);

    if let Some(rqrr_result) = result_rqrr {
        assert_eq!(
            Some(rqrr_result),
            expected,
            "rqrr decoded but gave different result"
        );
    }
}

#[test]
fn test_qr_capstone_interference() {
    let expected = Some((
        "QR-Code".to_string(),
        "http://txz.qq.com/p?k=T8sZMvS*JxhU0kQFseMOMQZAKuE7An3u&f=716027609".to_string(),
    ));

    let result_this = decode_image("examples/qr-code-capstone-interference.png");
    let result_rqrr = decode_with_rqrr("examples/qr-code-capstone-interference.png");

    assert_eq!(result_this, expected, "zedbar failed");
    assert_matches_zbars("examples/qr-code-capstone-interference.png", &expected);

    if let Some(rqrr_result) = result_rqrr {
        assert_eq!(
            Some(rqrr_result),
            expected,
            "rqrr decoded but gave different result"
        );
    }
}

#[test]
fn test_qr_color_bands() {
    let result_this = decode_image("examples/qr-code-color-bands.png");
    let result_rqrr = decode_with_rqrr("examples/qr-code-color-bands.png");

    assert!(result_this.is_some(), "zedbar failed");

    // Check zedbar result
    let (symbol_type, data) = result_this.as_ref().unwrap();
    assert_eq!(symbol_type, "QR-Code");
    assert!(data.starts_with("二维码生成器"));
    assert!(data.contains("https://zh.qr-code-generator.com"));

    // zedbar and zbars should agree
    assert_matches_zbars("examples/qr-code-color-bands.png", &result_this);

    // rqrr should match if it succeeds
    if let Some(rqrr_result) = result_rqrr {
        assert_eq!(
            Some(rqrr_result),
            result_this,
            "rqrr decoded but gave different result"
        );
    }
}

#[test]
fn test_qr_low_contrast() {
    let expected = Some((
        "QR-Code".to_string(),
        "72f7f23bf1e7428ccaba5366a938f420".to_string(),
    ));

    let result_this = decode_image("examples/qr-code-low-contrast.png");
    let result_rqrr = decode_with_rqrr("examples/qr-code-low-contrast.png");

    assert_eq!(result_this, expected, "zedbar failed");
    assert_matches_zbars("examples/qr-code-low-contrast.png", &expected);

    if let Some(rqrr_result) = result_rqrr {
        assert_eq!(
            Some(rqrr_result),
            expected,
            "rqrr decoded but gave different result"
        );
    }
}

#[test]
fn test_qr_pacman() {
    let expected = Some((
        "QR-Code".to_string(),
        "http://weixin.qq.com/r/gkgQCDrEc2cMrX5r9x2Q".to_string(),
    ));

    let result_this = decode_image("examples/qr-code-pacman.png");
    let result_rqrr = decode_with_rqrr("examples/qr-code-pacman.png");

    assert_eq!(result_this, expected, "zedbar failed");
    assert_matches_zbars("examples/qr-code-pacman.png", &expected);

    if let Some(rqrr_result) = result_rqrr {
        assert_eq!(
            Some(rqrr_result),
            expected,
            "rqrr decoded but gave different result"
        );
    }
}

#[test]
fn test_codabar() {
    let expected = Some(("Codabar".to_string(), "A40156B".to_string()));

    let result_this = decode_image("examples/test-codabar.png");

    assert_eq!(result_this, expected, "zedbar failed");
    assert_matches_zbars("examples/test-codabar.png", &expected);
}

#[test]
fn test_code128() {
    let expected = Some(("CODE-128".to_string(), "HELLO123".to_string()));

    let result_this = decode_image("examples/test-code128.png");

    assert_eq!(result_this, expected, "zedbar failed");
    assert_matches_zbars("examples/test-code128.png", &expected);
}

/// Code 128 packs two digits per character in character set C, so a payload
/// that switches between numeric and alphanumeric runs makes the decoder
/// expand set-C characters mid-buffer. Regression: everything after the
/// first expansion used to be dropped ("1234AB" decoded as "1234").
#[test]
fn test_code128_character_set_switches() {
    for (path, payload) in [
        ("examples/test-code128-setc.png", "12345678"),
        ("examples/test-code128-mixed.png", "1234AB"),
        ("examples/test-code128-mixed2.png", "1234AB5678"),
    ] {
        let expected = Some(("CODE-128".to_string(), payload.to_string()));

        let result_this = decode_image(path);

        assert_eq!(result_this, expected, "zedbar failed for {path}");
        assert_matches_zbars(path, &expected);
    }
}

/// Decoded data must not carry a trailing NUL. Code 128 used to include the
/// C-style string terminator in the symbol data.
#[test]
fn test_decoded_data_has_no_trailing_nul() {
    for path in [
        "examples/test-code128.png",
        "examples/test-code128-mixed.png",
        "examples/test-code39.png",
        "examples/test-code93.png",
        "examples/test-codabar.png",
        "examples/test-i25.png",
        "examples/test-ean13.png",
    ] {
        let img = image::open(path)
            .unwrap_or_else(|e| panic!("open {path}: {e}"))
            .to_luma8();
        let mut zbar_image =
            Image::from_gray(img.as_raw(), img.width(), img.height()).expect("create image");
        let mut scanner = Scanner::new();
        let result = scanner.scan(&mut zbar_image);
        assert!(!result.symbols().is_empty(), "{path}: nothing decoded");
        for symbol in result.symbols() {
            assert!(
                !symbol.data().contains(&0),
                "{path}: {:?} data contains a NUL byte: {:?}",
                symbol.symbol_type(),
                symbol.data()
            );
        }
    }
}

#[test]
fn test_code39() {
    let expected = Some(("CODE-39".to_string(), "TEST123".to_string()));

    let result_this = decode_image("examples/test-code39.png");

    assert_eq!(result_this, expected, "zedbar failed");
    assert_matches_zbars("examples/test-code39.png", &expected);
}

#[test]
fn test_code93() {
    let expected = Some(("CODE-93".to_string(), "CODE93".to_string()));

    let result_this = decode_image("examples/test-code93.png");

    assert_eq!(result_this, expected, "zedbar failed");
    assert_matches_zbars("examples/test-code93.png", &expected);
}

#[test]
fn test_ean13() {
    let expected = Some(("EAN-13".to_string(), "5901234123457".to_string()));

    let result_this = decode_image("examples/test-ean13.png");

    assert_eq!(result_this, expected, "zedbar failed");
    assert_matches_zbars("examples/test-ean13.png", &expected);
}

/// A 978-prefixed EAN-13 is a Bookland ISBN. zbar relabels it ISBN-13 when
/// that symbology is opted in, and reduces it to the 10-digit form (prefix
/// dropped, check digit recomputed mod 11) for ISBN-10.
///
/// ISBN-13 wins when both are enabled, which is the ordering in
/// `integrate_partial`. Cross-checked against
/// `zbarimg -S isbn13.enable=1` / `-S isbn10.enable=1`.
#[test]
fn test_isbn13_and_isbn10() {
    use zedbar::config::{Isbn10, Isbn13};

    let path = "examples/test-isbn13.png";
    let img = image::open(path).expect("open").to_luma8();
    let raw = img.as_raw();
    let (w, h) = (img.width(), img.height());

    let decode = |config: DecoderConfig| {
        let mut image = Image::from_gray(raw, w, h).expect("image");
        let result = Scanner::with_config(config).scan(&mut image);
        let s = result.first().expect("no symbol decoded");
        (s.symbol_type(), s.data_string().unwrap_or("").to_string())
    };

    // Neither opted in: it stays an ordinary EAN-13.
    assert_eq!(
        decode(DecoderConfig::all()),
        (SymbolType::Ean13, "9780306406157".into())
    );

    assert_eq!(
        decode(DecoderConfig::all().enable(Isbn13)),
        (SymbolType::Isbn13, "9780306406157".into())
    );

    assert_eq!(
        decode(DecoderConfig::all().enable(Isbn10)),
        (SymbolType::Isbn10, "0306406152".into())
    );

    // ISBN-13 takes precedence over ISBN-10.
    assert_eq!(
        decode(DecoderConfig::all().enable(Isbn10).enable(Isbn13)),
        (SymbolType::Isbn13, "9780306406157".into())
    );
}

/// EAN-13 with a 2- or 5-digit add-on, which zbar reports as a composite of
/// the main symbol and the add-on concatenated.
///
/// Regression test for the add-on state machine in `decode_pass`: stepping
/// over an add-on character boundary advances both `pass->state` and the local
/// index, and the port advanced only the state, leaving the
/// character-boundary test reading a stale index.
///
/// Add-ons stay opt-in — zbar leaves EAN-2/EAN-5 disabled by default, and
/// [`test_ean13_addon_ignored_by_default`] pins that.
#[test]
fn test_ean13_with_addons() {
    use zedbar::config::{Ean2, Ean5};

    for (path, expected) in [
        ("examples/test-ean13-addon2.png", "590123412345712"),
        ("examples/test-ean13-addon5.png", "590123412345712345"),
    ] {
        let img = image::open(path).expect("open").to_luma8();
        let mut image = Image::from_gray(img.as_raw(), img.width(), img.height()).expect("image");

        let config = DecoderConfig::all()
            .enable(Ean2)
            .enable(Ean5)
            .enable_type(SymbolType::Composite);
        let result = Scanner::with_config(config).scan(&mut image);

        let symbol = result
            .first()
            .unwrap_or_else(|| panic!("no symbol in {path}"));
        assert_eq!(symbol.symbol_type(), SymbolType::Composite, "{path}");
        assert_eq!(symbol.data_string(), Some(expected), "{path}");
    }
}

/// Without `Composite` the add-on is reported as its own EAN-2 / EAN-5
/// symbol alongside the main EAN-13, which is how zbar behaves with
/// `-S ean2.enable=1 -S ean5.enable=1` and composite left off.
#[test]
fn test_ean2_and_ean5_as_standalone_symbols() {
    use zedbar::config::{Ean2, Ean5};

    for (path, addon_type, addon_data) in [
        ("examples/test-ean13-addon2.png", SymbolType::Ean2, "12"),
        ("examples/test-ean13-addon5.png", SymbolType::Ean5, "12345"),
    ] {
        let img = image::open(path).expect("open").to_luma8();
        let mut image = Image::from_gray(img.as_raw(), img.width(), img.height()).expect("image");

        let config = DecoderConfig::all().enable(Ean2).enable(Ean5);
        let result = Scanner::with_config(config).scan(&mut image);

        let found: Vec<(SymbolType, &str)> = result
            .symbols()
            .iter()
            .map(|s| (s.symbol_type(), s.data_string().unwrap_or("")))
            .collect();

        assert!(
            found.contains(&(SymbolType::Ean13, "5901234123457")),
            "{path}: missing the main symbol, got {found:?}"
        );
        assert!(
            found.contains(&(addon_type, addon_data)),
            "{path}: missing the add-on, got {found:?}"
        );
    }
}

/// zbar ships with EAN-2/EAN-5 disabled, so an add-on is not read and only the
/// main symbol comes back.
#[test]
fn test_ean13_addon_ignored_by_default() {
    let expected = Some(("EAN-13".to_string(), "5901234123457".to_string()));

    for path in [
        "examples/test-ean13-addon2.png",
        "examples/test-ean13-addon5.png",
    ] {
        assert_eq!(decode_image(path), expected, "zedbar failed for {path}");
        assert_matches_zbars(path, &expected);
    }
}

/// SQ Code, reported base64-encoded like zbar does.
///
/// No encoder for this symbology exists in any tool I could find — zint,
/// ZXing and zbar all only read it, and zbar ships no sample — so the
/// fixtures are generated by `examples/generate_sqcode_fixture.rs`, which
/// derives the layout from zbar's decoder. `zbarimg` reads them, which is
/// what makes them trustworthy: an independent implementation agrees on
/// both the geometry and the payload.
///
/// The two sizes cover an 8x8 and a 12x12 data area, so the border walk runs
/// over different numbers of alignment dots.
///
/// zedbar reports the sampled bits directly where zbar reports them
/// base64-encoded, so the cross-check encodes before comparing — which also
/// documents the relationship between the two.
#[test]
fn test_sqcode() {
    for (path, payload) in [
        ("examples/test-sqcode.png", &b"zedbar!!"[..]),
        ("examples/test-sqcode-large.png", &b"zedbar SQ Code fix"[..]),
    ] {
        let img = image::open(path).expect("open").to_luma8();
        let mut image = Image::from_gray(img.as_raw(), img.width(), img.height()).expect("image");

        let result = Scanner::new().scan(&mut image);
        let symbol = result
            .first()
            .unwrap_or_else(|| panic!("nothing decoded in {path}"));

        assert_eq!(symbol.symbol_type(), SymbolType::SqCode, "{path}");
        assert_eq!(symbol.data(), payload, "{path}");
        // No text conversion happens, so there is nothing for `raw_data` to
        // hand back that `data` does not already carry.
        assert_eq!(symbol.raw_data(), None, "{path}");

        if zbarimg_available() {
            assert_eq!(
                decode_with_zbars(path),
                Some(("SQ-Code".to_string(), base64(payload))),
                "zbars disagreed for {path}"
            );
        }
    }
}

/// Base64 as zbar emits it, for the SQ cross-check above. Small enough not to
/// be worth a dev-dependency.
fn base64(bytes: &[u8]) -> String {
    const TABLE: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::new();
    for chunk in bytes.chunks(3) {
        let b = [
            chunk[0],
            *chunk.get(1).unwrap_or(&0),
            *chunk.get(2).unwrap_or(&0),
        ];
        let n = u32::from_be_bytes([0, b[0], b[1], b[2]]);
        for i in 0..4 {
            if i <= chunk.len() {
                out.push(TABLE[(n >> (18 - 6 * i)) as usize & 0x3f] as char);
            } else {
                out.push('=');
            }
        }
    }
    out
}

#[test]
fn test_ean8_plain() {
    let expected = Some(("EAN-8".to_string(), "96385074".to_string()));

    let result_this = decode_image("examples/test-ean8-plain.png");

    assert_eq!(result_this, expected, "zedbar failed");
    assert_matches_zbars("examples/test-ean8-plain.png", &expected);
}

/// With EMIT_CHECK turned off, zbar drops the trailing check digit and emits
/// the rest. The port derived that shortened length by subtracting one from
/// the symbology enum and converting back to a `SymbolType`; 7 (EAN-8 minus
/// its check digit) names no symbology, so it collapsed to `None` and the
/// decoder emitted an empty string. UPC-A and UPC-E hit the same hole.
#[test]
fn test_ean8_without_emitted_checksum() {
    use zedbar::config::Ean8;

    let img = image::open("examples/test-ean8-plain.png")
        .expect("open")
        .to_luma8();
    let mut image = Image::from_gray(img.as_raw(), img.width(), img.height()).expect("image");

    let config = DecoderConfig::new().set_checksum(Ean8, false, false);
    let result = Scanner::with_config(config).scan(&mut image);

    let symbol = result.first().expect("no symbol decoded");
    assert_eq!(symbol.symbol_type(), SymbolType::Ean8);
    assert_eq!(symbol.data_string(), Some("9638507"));
}

#[test]
fn test_ean8_decoded_as_ean13() {
    // Note: This image is decoded as EAN13, not EAN8 or UPCA
    let expected = Some(("EAN-13".to_string(), "0000963850742".to_string()));

    let result_this = decode_image("examples/test-ean8.png");

    assert_eq!(result_this, expected, "zedbar failed");
    assert_matches_zbars("examples/test-ean8.png", &expected);
}

#[test]
fn test_i25() {
    let expected = Some(("I2/5".to_string(), "1234567890".to_string()));

    let result_this = decode_image("examples/test-i25.png");

    assert_eq!(result_this, expected, "zedbar failed");
    assert_matches_zbars("examples/test-i25.png", &expected);
}

#[test]
fn test_code128_single_character() {
    // zbar sets no minimum length for Code 128, so a one-character payload
    // (3 characters counting start and check) must still decode.
    let expected = Some(("CODE-128".to_string(), "A".to_string()));

    let result_this = decode_image("examples/test-code128-1char.png");

    assert_eq!(result_this, expected, "zedbar failed");
    assert_matches_zbars("examples/test-code128-1char.png", &expected);
}

#[test]
fn test_code93_single_character() {
    // Likewise for Code 93: one data character plus the C and K check
    // characters is three, which is above zbar's hard floor of two.
    let expected = Some(("CODE-93".to_string(), "A".to_string()));

    let result_this = decode_image("examples/test-code93-1char.png");

    assert_eq!(result_this, expected, "zedbar failed");
    assert_matches_zbars("examples/test-code93-1char.png", &expected);
}

#[test]
fn test_databar() {
    // GTIN-13 2001234567890 plus the emitted check digit, prefixed with AI 01.
    let expected = Some(("DataBar".to_string(), "0120012345678909".to_string()));

    let result_this = decode_image("examples/test-databar.png");

    assert_eq!(result_this, expected, "zedbar failed");
    assert_matches_zbars("examples/test-databar.png", &expected);
}

#[test]
fn test_databar_expanded() {
    // AI 01 (GTIN-14) + AI 3202 (net weight) + AI 15 (best before date).
    let expected = Some((
        "DataBar-Exp".to_string(),
        "0198898765432106320201234515991231".to_string(),
    ));

    let result_this = decode_image("examples/test-databar-exp.png");

    assert_eq!(result_this, expected, "zedbar failed");
    assert_matches_zbars("examples/test-databar-exp.png", &expected);
}

#[test]
fn test_databar_expanded_alphanumeric() {
    // Exercises the general-purpose data compaction tail, which needs the bit
    // feeder to actually advance through the character stream.
    let expected = Some(("DataBar-Exp".to_string(), "10ABC123".to_string()));

    let result_this = decode_image("examples/test-databar-exp-alnum.png");

    assert_eq!(result_this, expected, "zedbar failed");
    assert_matches_zbars("examples/test-databar-exp-alnum.png", &expected);
}

#[test]
fn test_upca_decoded_as_ean13() {
    // Note: This image is decoded as EAN13, not UPCA
    let expected = Some(("EAN-13".to_string(), "0012345678905".to_string()));

    let result_this = decode_image("examples/test-upca.png");

    assert_eq!(result_this, expected, "zedbar failed");
    assert_matches_zbars("examples/test-upca.png", &expected);
}

#[test]
fn test_rqrr_crash_2() {
    // Regression test for integer overflow in databar decoder
    // This image previously caused panic due to epoch overflow
    let result_this = decode_image("examples/rqrr-crash-2.jpeg");
    let result_rqrr = decode_with_rqrr("examples/rqrr-crash-2.jpeg");

    assert!(result_this.is_some(), "zedbar failed");

    let (symbol_type, data) = result_this.as_ref().unwrap();
    assert_eq!(symbol_type, "QR-Code");
    assert!(data.contains("欢迎访问太平洋IT百科栏目"));

    // zedbar and zbars should agree
    assert_matches_zbars("examples/rqrr-crash-2.jpeg", &result_this);

    // rqrr should match if it succeeds
    if let Some(rqrr_result) = result_rqrr {
        assert_eq!(
            Some(rqrr_result),
            result_this,
            "rqrr decoded but gave different result"
        );
    }
}

#[test]
fn test_rqrr_crash_3() {
    let expected_data = "shc:/567629095243206034602924374044603122295953265460346029254077280433602870286471674522280928613331456437653141590640220306450459085643550341424541364037063665417137241236380304375622046737407532323925433443326057360106452931531270742428395038692212766728666731266342087422573776302062041022437658685343255820002167287607585708105505622752282407670809680507692361773323356634342439664440596761410443377667202663224433674530596175400038397052612140292974753658337372662132066669047253044469405210524536242721550377673434280323045475690310233670562227414567090555653507636250537239522776211205312561442568282012726838630039087127042463716936535535602928393065580072763158437500341209546904210458383257586630101033123422114008776058732325243477645920113037325929083272452732223707055550412927584543582550667760036577724025621136525340592771740903663844771261692077697211447057562509437029626707254539002011763240720310114260256672645965627243654061066553770056003044082967606162724306592273682223412466107335331229606157521057357572327529693965670332063208596309543400076452696835713027450728663529345234666377297208583525543653527774072234735706452828641140633528387577054371703966706421520708254156041170353656054471407636552612616834377244090406554327122559623453686207006139712936404138601156656945315611255669116044703333731263580306106975715411702932060511012768634011703371553353213365032550756476005853005224547339310064671161682376335069647622323339523133724171327531702738363650063527592633763908656123314363227707566731311074";
    let expected = Some(("QR-Code".to_string(), expected_data.to_string()));

    let result_this = decode_image("examples/rqrr-crash-3.jpeg");
    let result_rqrr = decode_with_rqrr("examples/rqrr-crash-3.jpeg");

    assert_eq!(result_this, expected, "zedbar failed");
    assert_matches_zbars("examples/rqrr-crash-3.jpeg", &expected);

    if let Some(rqrr_result) = result_rqrr {
        assert_eq!(
            Some(rqrr_result),
            expected,
            "rqrr decoded but gave different result"
        );
    }
}

#[test]
fn test_all_examples_decode() {
    // Verify that all example images can be decoded without errors
    let images = vec![
        "examples/test-qr.png",
        "examples/test-qr.jpg",
        "examples/test-qr.webp",
        "examples/pixel-wifi-sharing-qr-code.png",
        "examples/test-qr-version40.png",
        // Explicitly DON'T pass these to `decode_image`.
        // "examples/qr-code-140-grid01.jpg",
        // "examples/qr-code-140-grid02.jpg",
        // "examples/qr-code-140-grid03.jpg",
        "examples/qr-code-capstone-interference.png",
        "examples/qr-code-color-bands.png",
        "examples/qr-code-low-contrast.png",
        "examples/qr-code-pacman.png",
        "examples/test-codabar.png",
        "examples/test-databar.png",
        "examples/test-databar-exp.png",
        "examples/test-databar-exp-alnum.png",
        "examples/test-code128.png",
        "examples/test-code128-2.png",
        "examples/test-code128-1char.png",
        "examples/test-code93-1char.png",
        "examples/test-code39.png",
        "examples/test-code93.png",
        "examples/test-ean13.png",
        "examples/test-ean8.png",
        "examples/test-ean8-plain.png",
        "examples/test-ean13-addon2.png",
        "examples/test-ean13-addon5.png",
        "examples/test-isbn13.png",
        "examples/test-sqcode.png",
        "examples/test-sqcode-large.png",
        "examples/test-i25.png",
        "examples/test-upca.png",
        "examples/nine-barcodes.png",
    ];

    for image_path in images {
        let result = decode_image(image_path);
        assert!(result.is_some(), "Failed to decode {}", image_path);
    }
}

/// Helper function to decode an image file and return raw binary data
fn decode_image_binary(path: &str) -> Option<(String, Vec<u8>)> {
    let path = Path::new(path);
    if !path.exists() {
        return None;
    }

    let img = image::open(path).ok()?;
    let img = downscale_if_needed(img, 1280).to_luma8();

    let mut scanner = Scanner::new();

    // Create a zbar Image from the image buffer
    let mut zbar_image = Image::from_gray(img.as_raw(), img.width(), img.height()).ok()?;

    // Scan the image
    let symbols = scanner.scan(&mut zbar_image);

    symbols.first().map(|symbol| {
        let symbol_type = symbol.symbol_type().to_string();
        let data = symbol.data().to_vec();
        (symbol_type, data)
    })
}

/// Helper function to decode binary data using zbar binary (system zbarimg command)
fn decode_with_zbars_binary(path: &str) -> Option<(String, Vec<u8>)> {
    use std::process::Command;

    let path = Path::new(path);
    if !path.exists() {
        return None;
    }

    // Get symbol type without --raw
    let type_output = Command::new("zbarimg")
        .arg("--quiet")
        .arg(path)
        .output()
        .ok()?;

    if !type_output.status.success() {
        return None;
    }

    let type_line = String::from_utf8_lossy(&type_output.stdout);
    let symbol_type = type_line.lines().next()?.split(':').next()?.to_string();

    // Get raw binary data with --raw
    let data_output = Command::new("zbarimg")
        .arg("--raw")
        .arg("--quiet")
        .arg(path)
        .output()
        .ok()?;

    if !data_output.status.success() {
        return None;
    }

    // Strip trailing newline if present
    let mut data = data_output.stdout;
    if data.last() == Some(&b'\n') {
        data.pop();
    }

    Some((symbol_type, data))
}

/// Helper function to decode QR codes using rqrr and return binary data
fn decode_with_rqrr_binary(path: &str) -> Option<(String, Vec<u8>)> {
    let path = Path::new(path);
    if !path.exists() {
        return None;
    }

    let img = image::open(path).ok()?.to_luma8();
    let mut prepared_img = rqrr::PreparedImage::prepare(img);
    let grids = prepared_img.detect_grids();

    grids.first().and_then(|grid| {
        let (_meta, content) = grid.decode().ok()?;
        Some(("QR-Code".to_string(), content.into_bytes()))
    })
}

#[test]
fn test_rqrr_crash_4_binary() {
    // This is a binary QR code that previously failed to decode
    // The bug was in qrdectxt.rs where decoding failures were treated as errors
    let result_this = decode_image_binary("examples/rqrr-crash-4.png");
    let result_zbars = decode_with_zbars_binary("examples/rqrr-crash-4.png");
    let result_rqrr = decode_with_rqrr_binary("examples/rqrr-crash-4.png");

    // zedbar should succeed
    assert!(result_this.is_some(), "zedbar should decode the QR code");

    let (symbol_type, data) = result_this.as_ref().unwrap();
    assert_eq!(symbol_type, "QR-Code");

    // data() returns text-decoded output (encoding-detected, converted to UTF-8).
    // The raw QR bytes go through Windows-1252 → UTF-8 conversion, expanding
    // the 1376 raw bytes to 2146 UTF-8 bytes.
    assert_eq!(data.len(), 2146, "Text-decoded data should be 2146 bytes");
    assert_eq!(&data[..4], &[0x07, 0xC3, 0x84, 0x18]);

    // Verify it contains binary data (bytes outside printable ASCII)
    let has_binary = data.iter().any(|b| !(0x20..0x80).contains(b));
    assert!(has_binary, "Should contain binary data");

    // zbars should produce similar text-decoded data (first bytes should match)
    if let Some((zbars_symbol_type, zbars_data)) = result_zbars {
        assert_eq!(zbars_symbol_type, "QR-Code");
        assert_eq!(
            &zbars_data[..4],
            &[0x07, 0xC3, 0x84, 0x18],
            "zbars should decode to similar data"
        );
        assert!(zbars_data.len() >= 2000, "zbars data should be substantial");
    }

    // rqrr must match exactly if it succeeds (strict check)
    if let Some(rqrr_result) = result_rqrr {
        assert_eq!(
            Some(rqrr_result),
            result_this.clone(),
            "rqrr decoded but gave different result"
        );
    }
}

#[test]
fn test_raw_data_binary_qr() {
    // raw_data() should return the original bytes before encoding conversion
    let path = Path::new("examples/rqrr-crash-4.png");
    if !path.exists() {
        return;
    }

    let img = image::open(path).unwrap();
    let img = downscale_if_needed(img, 1280).to_luma8();
    let mut scanner = Scanner::new();
    let mut zbar_image = Image::from_gray(img.as_raw(), img.width(), img.height()).unwrap();
    let symbols = scanner.scan(&mut zbar_image);
    let symbol = symbols.first().expect("should decode the QR code");

    let raw = symbol.raw_data().expect("QR codes should have raw_data");
    let text_decoded = symbol.data();

    // Raw bytes are the pre-encoding-conversion data (1376 bytes)
    assert_eq!(raw.len(), 1376, "Raw data should be 1376 bytes");
    assert_eq!(raw[0], 0x07);
    assert_eq!(
        raw[1], 0xC4,
        "Second raw byte should be 0xC4, not the UTF-8 lead byte 0xC3"
    );

    // Text-decoded data is longer due to Windows-1252 → UTF-8 expansion
    assert_eq!(text_decoded.len(), 2146);
    assert!(raw.len() < text_decoded.len());

    // This particular file's raw bytes are detected as Windows-1252 by the
    // encoding heuristic. If the heuristic changes, the exact text output may
    // differ — but raw bytes should always be shorter than the text-decoded
    // output for this binary payload.
    let (decoded, _enc, _had_errors) = encoding_rs::WINDOWS_1252.decode(raw);
    assert_eq!(
        text_decoded,
        decoded.as_bytes(),
        "Expected text output to match Windows-1252 decode of raw bytes (if this fails, \
         the encoding heuristic may have changed)"
    );
}

#[test]
fn test_raw_data_text_qr() {
    // For ASCII/UTF-8 QR codes, raw_data() and data() should contain the same bytes
    let path = Path::new("examples/test-qr.png");
    if !path.exists() {
        return;
    }

    let img = image::open(path).unwrap();
    let img = downscale_if_needed(img, 1280).to_luma8();
    let mut scanner = Scanner::new();
    let mut zbar_image = Image::from_gray(img.as_raw(), img.width(), img.height()).unwrap();
    let symbols = scanner.scan(&mut zbar_image);
    let symbol = symbols.first().expect("should decode the QR code");

    let raw = symbol.raw_data().expect("QR codes should have raw_data");
    let text_decoded = symbol.data();

    // For ASCII text, raw bytes and text-decoded bytes are identical
    assert_eq!(raw, text_decoded);
    assert_eq!(
        std::str::from_utf8(raw).unwrap(),
        "Hello, simplified zbar!\n"
    );
}

#[test]
fn test_raw_data_linear_barcode() {
    // Linear barcodes don't have raw_data (their data is always ASCII)
    let path = Path::new("examples/test-ean13.png");
    if !path.exists() {
        return;
    }

    let img = image::open(path).unwrap();
    let img = downscale_if_needed(img, 1280).to_luma8();
    let mut scanner = Scanner::new();
    let mut zbar_image = Image::from_gray(img.as_raw(), img.width(), img.height()).unwrap();
    let symbols = scanner.scan(&mut zbar_image);

    if let Some(symbol) = symbols.first() {
        assert!(
            symbol.raw_data().is_none(),
            "Linear barcodes should not have raw_data"
        );
    }
}

#[test]
fn test_qr_140_grids() {
    for path in [
        "examples/qr-code-140-grid01.jpg",
        "examples/qr-code-140-grid02.jpg",
        "examples/qr-code-140-grid03.jpg",
    ] {
        let img = image::open(path).unwrap();
        let img = img.to_luma8();
        let mut scanner = Scanner::new();

        let mut image = Image::from_gray(img.as_raw(), img.width(), img.height()).unwrap();
        let result = scanner.scan(&mut image);

        assert_eq!(result.len(), 140);
        for (i, symbol) in result.iter().enumerate() {
            assert_eq!(
                symbol.symbol_type(),
                SymbolType::QrCode,
                "symbol {i} is not a QR code"
            );
        }
    }
}

#[test]
fn test_nine_barcodes() {
    // Image contains 9 different barcode types. Verify all are decoded correctly.
    let path = "examples/nine-barcodes.png";
    let img = image::open(path).expect("failed to open nine-barcodes.png");
    let img = downscale_if_needed(img, 1280).to_luma8();

    // Kitchen sink + opt in to UPC-A/UPC-E variant labeling on top of EAN.
    let config = DecoderConfig::all().enable(Upca).enable(Upce);
    let mut scanner = Scanner::with_config(config);
    let mut zbar_image = Image::from_gray(img.as_raw(), img.width(), img.height()).unwrap();
    let symbols = scanner.scan(&mut zbar_image);

    // Collect decoded symbols as (type, data) pairs
    let mut results: Vec<(SymbolType, String)> = symbols
        .iter()
        .map(|s| {
            let data = String::from_utf8_lossy(s.data())
                .trim_end_matches('\0')
                .to_string();
            (s.symbol_type(), data)
        })
        .collect();
    results.sort_by_key(|r| r.0.to_string());

    let expected: Vec<(SymbolType, &str)> = {
        let mut v = vec![
            (SymbolType::Code39, "CODE39"),
            (SymbolType::Code93, "CODE93"),
            (SymbolType::Code128, "CODE128"),
            (SymbolType::Codabar, "C012345D"),
            (SymbolType::I25, "00123456"),
            (SymbolType::Ean8, "01234565"),
            (SymbolType::Upca, "012345678905"),
            (SymbolType::Ean13, "1234567890128"),
            (SymbolType::Upce, "01234565"),
        ];
        v.sort_by_key(|r| r.0.to_string());
        v
    };

    assert_eq!(
        results.len(),
        expected.len(),
        "expected {} barcodes, found {}.\nDecoded: {:?}",
        expected.len(),
        results.len(),
        results
    );

    for (actual, expected) in results.iter().zip(expected.iter()) {
        assert_eq!(actual.0, expected.0, "symbol type mismatch");
        assert_eq!(actual.1, expected.1, "data mismatch for {:?}", expected.0);
    }
}
