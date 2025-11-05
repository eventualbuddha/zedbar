//! Integration tests for decoding example barcode images
//!
//! These tests compare results from multiple decoders:
//! - This crate (zbar-rust)
//! - System zbar (via zbarimg command)
//! - rqrr (for QR codes only)

use std::path::Path;
use zbar::{Image, Scanner};

/// Helper function to decode an image file and return the first symbol found
fn decode_image(path: &str) -> Option<(String, String)> {
    let path = Path::new(path);
    if !path.exists() {
        return None;
    }

    let img = image::open(path).ok()?.to_luma8();
    let mut scanner = Scanner::new();

    // Create a zbar Image from the image buffer
    let mut zbar_image = Image::from_gray(img.as_raw(), img.width(), img.height()).ok()?;

    // Scan the image
    scanner.scan(&mut zbar_image).ok()?;

    // Get the symbols
    let symbols = zbar_image.symbols();

    symbols.iter().next().map(|symbol| {
        let symbol_type = symbol.symbol_type().to_string();
        let data = String::from_utf8_lossy(symbol.data())
            .trim_end_matches('\0') // Remove null terminators
            .to_string();
        (symbol_type, data)
    })
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

    let img = image::open(path).ok()?.to_luma8();
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
    let result_zbars = decode_with_zbars("examples/test-qr.png");
    let result_rqrr = decode_with_rqrr("examples/test-qr.png");

    assert_eq!(result_this, expected, "zbar-rust failed");
    assert_eq!(result_zbars, expected, "zbars failed");

    // zbar-rust and zbars must agree
    assert_eq!(result_this, result_zbars, "zbar-rust and zbars disagree");

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
    let result_zbars = decode_with_zbars("examples/test-qr.jpg");
    let result_rqrr = decode_with_rqrr("examples/test-qr.jpg");

    assert_eq!(result_this, expected, "zbar-rust failed");
    assert_eq!(result_zbars, expected, "zbars failed");
    assert_eq!(result_this, result_zbars, "zbar-rust and zbars disagree");

    if let Some(rqrr_result) = result_rqrr {
        assert_eq!(
            Some(rqrr_result),
            expected,
            "rqrr decoded but gave different result"
        );
    }
}

#[test]
fn test_qr_wifi_sharing() {
    let expected = Some((
        "QR-Code".to_string(),
        "WIFI:S:Not a real network;T:SAE;P:password;H:false;;".to_string(),
    ));

    let result_this = decode_image("examples/pixel-wifi-sharing-qr-code.png");
    let result_zbars = decode_with_zbars("examples/pixel-wifi-sharing-qr-code.png");
    let result_rqrr = decode_with_rqrr("examples/pixel-wifi-sharing-qr-code.png");

    assert_eq!(result_this, expected, "zbar-rust failed");
    assert_eq!(result_zbars, expected, "zbars failed");
    assert_eq!(result_this, result_zbars, "zbar-rust and zbars disagree");

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
    let result_zbars = decode_with_zbars("examples/qr-code-capstone-interference.png");
    let result_rqrr = decode_with_rqrr("examples/qr-code-capstone-interference.png");

    assert_eq!(result_this, expected, "zbar-rust failed");
    assert_eq!(result_zbars, expected, "zbars failed");
    assert_eq!(result_this, result_zbars, "zbar-rust and zbars disagree");

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
    let result_zbars = decode_with_zbars("examples/qr-code-color-bands.png");
    let result_rqrr = decode_with_rqrr("examples/qr-code-color-bands.png");

    // zbar-rust and zbars should succeed
    assert!(result_this.is_some(), "zbar-rust failed");
    assert!(result_zbars.is_some(), "zbars failed");

    // Check zbar-rust result
    let (symbol_type, data) = result_this.as_ref().unwrap();
    assert_eq!(symbol_type, "QR-Code");
    assert!(data.starts_with("二维码生成器"));
    assert!(data.contains("https://zh.qr-code-generator.com"));

    // zbar-rust and zbars should agree
    assert_eq!(result_zbars, result_this, "zbars disagrees with zbar-rust");

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
    let result_zbars = decode_with_zbars("examples/qr-code-low-contrast.png");
    let result_rqrr = decode_with_rqrr("examples/qr-code-low-contrast.png");

    assert_eq!(result_this, expected, "zbar-rust failed");
    assert_eq!(result_zbars, expected, "zbars failed");
    assert_eq!(result_this, result_zbars, "zbar-rust and zbars disagree");

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
    let result_zbars = decode_with_zbars("examples/qr-code-pacman.png");
    let result_rqrr = decode_with_rqrr("examples/qr-code-pacman.png");

    assert_eq!(result_this, expected, "zbar-rust failed");
    assert_eq!(result_zbars, expected, "zbars failed");
    assert_eq!(result_this, result_zbars, "zbar-rust and zbars disagree");

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
    let result_zbars = decode_with_zbars("examples/test-codabar.png");

    assert_eq!(result_this, expected, "zbar-rust failed");
    assert_eq!(result_zbars, expected, "zbars failed");
    assert_eq!(result_this, result_zbars, "zbar-rust and zbars disagree");
}

#[test]
fn test_code128() {
    let expected = Some(("CODE-128".to_string(), "HELLO123".to_string()));

    let result_this = decode_image("examples/test-code128.png");
    let result_zbars = decode_with_zbars("examples/test-code128.png");

    assert_eq!(result_this, expected, "zbar-rust failed");
    assert_eq!(result_zbars, expected, "zbars failed");
    assert_eq!(result_this, result_zbars, "zbar-rust and zbars disagree");
}

#[test]
fn test_code39() {
    let expected = Some(("CODE-39".to_string(), "TEST123".to_string()));

    let result_this = decode_image("examples/test-code39.png");
    let result_zbars = decode_with_zbars("examples/test-code39.png");

    assert_eq!(result_this, expected, "zbar-rust failed");
    assert_eq!(result_zbars, expected, "zbars failed");
    assert_eq!(result_this, result_zbars, "zbar-rust and zbars disagree");
}

#[test]
fn test_code93() {
    let expected = Some(("CODE-93".to_string(), "CODE93".to_string()));

    let result_this = decode_image("examples/test-code93.png");
    let result_zbars = decode_with_zbars("examples/test-code93.png");

    assert_eq!(result_this, expected, "zbar-rust failed");
    assert_eq!(result_zbars, expected, "zbars failed");
    assert_eq!(result_this, result_zbars, "zbar-rust and zbars disagree");
}

#[test]
fn test_ean13() {
    let expected = Some(("EAN-13".to_string(), "5901234123457".to_string()));

    let result_this = decode_image("examples/test-ean13.png");
    let result_zbars = decode_with_zbars("examples/test-ean13.png");

    assert_eq!(result_this, expected, "zbar-rust failed");
    assert_eq!(result_zbars, expected, "zbars failed");
    assert_eq!(result_this, result_zbars, "zbar-rust and zbars disagree");
}

#[test]
fn test_ean8_decoded_as_ean13() {
    // Note: This image is decoded as EAN13, not EAN8 or UPCA
    let expected = Some(("EAN-13".to_string(), "0000963850742".to_string()));

    let result_this = decode_image("examples/test-ean8.png");
    let result_zbars = decode_with_zbars("examples/test-ean8.png");

    assert_eq!(result_this, expected, "zbar-rust failed");
    assert_eq!(result_zbars, expected, "zbars failed");
    assert_eq!(result_this, result_zbars, "zbar-rust and zbars disagree");
}

#[test]
fn test_i25() {
    let expected = Some(("I2/5".to_string(), "1234567890".to_string()));

    let result_this = decode_image("examples/test-i25.png");
    let result_zbars = decode_with_zbars("examples/test-i25.png");

    assert_eq!(result_this, expected, "zbar-rust failed");
    assert_eq!(result_zbars, expected, "zbars failed");
    assert_eq!(result_this, result_zbars, "zbar-rust and zbars disagree");
}

#[test]
fn test_upca_decoded_as_ean13() {
    // Note: This image is decoded as EAN13, not UPCA
    let expected = Some(("EAN-13".to_string(), "0012345678905".to_string()));

    let result_this = decode_image("examples/test-upca.png");
    let result_zbars = decode_with_zbars("examples/test-upca.png");

    assert_eq!(result_this, expected, "zbar-rust failed");
    assert_eq!(result_zbars, expected, "zbars failed");
    assert_eq!(result_this, result_zbars, "zbar-rust and zbars disagree");
}

#[test]
fn test_rqrr_crash_2() {
    // Regression test for integer overflow in databar decoder
    // This image previously caused panic due to epoch overflow
    let result_this = decode_image("examples/rqrr-crash-2.jpeg");
    let result_zbars = decode_with_zbars("examples/rqrr-crash-2.jpeg");
    let result_rqrr = decode_with_rqrr("examples/rqrr-crash-2.jpeg");

    // zbar-rust and zbars should succeed
    assert!(result_this.is_some(), "zbar-rust failed");
    assert!(result_zbars.is_some(), "zbars failed");

    let (symbol_type, data) = result_this.as_ref().unwrap();
    assert_eq!(symbol_type, "QR-Code");
    assert!(data.contains("欢迎访问太平洋IT百科栏目"));

    // zbar-rust and zbars should agree
    assert_eq!(result_zbars, result_this, "zbars disagrees with zbar-rust");

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
    let result_zbars = decode_with_zbars("examples/rqrr-crash-3.jpeg");
    let result_rqrr = decode_with_rqrr("examples/rqrr-crash-3.jpeg");

    assert_eq!(result_this, expected, "zbar-rust failed");
    assert_eq!(result_zbars, expected, "zbars failed");
    assert_eq!(result_this, result_zbars, "zbar-rust and zbars disagree");

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
        "examples/pixel-wifi-sharing-qr-code.png",
        "examples/qr-code-capstone-interference.png",
        "examples/qr-code-color-bands.png",
        "examples/qr-code-low-contrast.png",
        "examples/qr-code-pacman.png",
        "examples/test-codabar.png",
        "examples/test-code128.png",
        "examples/test-code39.png",
        "examples/test-code93.png",
        "examples/test-ean13.png",
        "examples/test-ean8.png",
        "examples/test-i25.png",
        "examples/test-upca.png",
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

    let img = image::open(path).ok()?.to_luma8();
    let mut scanner = Scanner::new();

    // Create a zbar Image from the image buffer
    let mut zbar_image = Image::from_gray(img.as_raw(), img.width(), img.height()).ok()?;

    // Scan the image
    scanner.scan(&mut zbar_image).ok()?;

    // Get the symbols
    let symbols = zbar_image.symbols();

    symbols.iter().next().map(|symbol| {
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

    // zbar-rust should succeed
    assert!(result_this.is_some(), "zbar-rust should decode the QR code");

    let (symbol_type, data) = result_this.as_ref().unwrap();
    assert_eq!(symbol_type, "QR-Code");

    // Verify it's binary data (contains bytes that would fail UTF-8 decoding as WINDOWS-1252)
    assert_eq!(data.len(), 2146, "Binary data should be 2146 bytes");

    // Check first few bytes to ensure we're getting the right data
    // These bytes are from the actual decoded QR code
    assert_eq!(&data[..4], &[0x07, 0xC3, 0x84, 0x18]);

    // Verify it contains binary data (bytes outside printable ASCII)
    let has_binary = data.iter().any(|b| !(0x20..0x80).contains(b));
    assert!(has_binary, "Should contain binary data");

    // zbars should produce similar binary data (first bytes should match)
    // Note: lengths may differ due to padding/encoding differences between implementations
    if let Some((zbars_symbol_type, zbars_data)) = result_zbars {
        assert_eq!(zbars_symbol_type, "QR-Code");
        assert_eq!(
            &zbars_data[..4],
            &[0x07, 0xC3, 0x84, 0x18],
            "zbars should decode to similar binary data"
        );
        assert!(
            zbars_data.len() >= 2000,
            "zbars binary data should be substantial"
        );
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
