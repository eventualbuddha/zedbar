//! Recovery of a QR finder with an erased outer bar must retain the normal
//! format and error-correction checks.
#![cfg(feature = "qrcode")]

use image::{GrayImage, Luma, imageops};
use zedbar::{Image, Scanner, SymbolType};

fn decoded(gray: &GrayImage) -> Vec<(SymbolType, Vec<u8>)> {
    let mut image = Image::from_gray(gray.as_raw(), gray.width(), gray.height()).unwrap();
    Scanner::new()
        .scan(&mut image)
        .iter()
        .map(|symbol| (symbol.symbol_type(), symbol.data().to_vec()))
        .collect()
}

fn assert_rotations(gray: &GrayImage, expected: &[u8]) {
    for mirrored in [false, true] {
        let mut rotated = if mirrored {
            imageops::flip_horizontal(gray)
        } else {
            gray.clone()
        };
        for rotation in 0..4 {
            assert_eq!(
                decoded(&rotated),
                [(SymbolType::QrCode, expected.to_vec())],
                "rotation {rotation}, mirrored {mirrored}"
            );
            rotated = imageops::rotate90(&rotated);
        }
    }
}

/// Original attachment from https://github.com/eventualbuddha/zedbar/issues/57.
/// The lower-left finder's left outer bar is almost entirely white. ZBar and
/// rqrr also miss it, so the expected payload comes from the reporter's iPhone
/// decode rather than a reference cross-check.
#[test]
fn issue57_erased_finder_bar() {
    let gray = image::open("examples/qr-erased-finder-bar.png")
        .unwrap()
        .to_luma8();
    assert_rotations(&gray, b"S;1;019ed9b6-6cde-7e33-a057-95a00735f1f4;A");
}

fn damaged_qr(payload: &[u8], module: u32) -> (GrayImage, u32) {
    let code = qrcode::QrCode::with_error_correction_level(payload, qrcode::EcLevel::M).unwrap();
    let dim = code.width() as u32;
    let mut gray = code
        .render::<Luma<u8>>()
        .module_dimensions(module, module)
        .build();
    // Four-module quiet zone. Erase the left outer bar of the bottom-left
    // finder, keeping its top/bottom edges and inner square intact.
    for y in (4 + dim - 6) * module..(4 + dim - 1) * module {
        for x in 4 * module..5 * module {
            gray.put_pixel(x, y, Luma([255]));
        }
    }
    (gray, dim)
}

#[test]
fn generated_erased_finder_bar() {
    for payload in [
        b"damaged finder".as_slice(),
        b"https://example.com/qr/damaged-finder/regression",
    ] {
        for module in [3, 5, 8] {
            let (gray, _) = damaged_qr(payload, module);
            assert_rotations(&gray, payload);
        }
    }
}

#[test]
fn damaged_finder_does_not_bypass_error_correction() {
    let module = 5;
    let (mut gray, dim) = damaged_qr(b"https://example.com/qr/damaged-finder/regression", module);
    // Preserve all three finder patterns but destroy most payload modules.
    for y in 12 * module..(4 + dim) * module {
        for x in 12 * module..(4 + dim) * module {
            gray.put_pixel(x, y, Luma([255]));
        }
    }
    for _ in 0..4 {
        assert!(decoded(&gray).is_empty());
        gray = imageops::rotate90(&gray);
    }
}

#[test]
fn many_partial_finders_do_not_create_symbols() {
    // Two complete finders enable recovery; the remaining 98 are plausible
    // damaged markers. There is no QR payload anywhere in the image. The
    // work done is bounded by caps in the recovery itself, not by this test.
    const MODULE: u32 = 4;
    const CELL: u32 = 11 * MODULE;
    let mut gray = GrayImage::from_pixel(10 * CELL, 10 * CELL, Luma([255]));
    for gy in 0..10 {
        for gx in 0..10 {
            for y in 0..7 * MODULE {
                for x in 0..7 * MODULE {
                    let (mx, my) = (x / MODULE, y / MODULE);
                    let dark = mx == 0
                        || mx == 6
                        || my == 0
                        || my == 6
                        || (2..=4).contains(&mx) && (2..=4).contains(&my);
                    let erased = (gy != 0 || gx > 1) && mx == 0 && (1..=5).contains(&my);
                    if dark && !erased {
                        gray.put_pixel(gx * CELL + x + MODULE, gy * CELL + y + MODULE, Luma([0]));
                    }
                }
            }
        }
    }
    assert!(decoded(&gray).is_empty());
}
