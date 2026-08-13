//! Regenerate the Code 128 fixtures used by the character-set-switch tests.
//!
//! Run via: `cargo run --example generate_code128_fixtures`
//!
//! Code 128 packs two digits into one symbol character in character set C,
//! so a payload that switches between numeric and alphanumeric runs exercises
//! the decoder's set-C expansion path. These fixtures cover one expansion
//! (`1234AB`) and two (`1234AB5678`), which regressed once because the
//! expansion grows the decode buffer mid-post-process.

use image::{GrayImage, Luma};

/// Module patterns for symbol values 0..=105, from the Code 128 symbology
/// (ISO/IEC 15417). `1` is a bar module, `0` a space module.
static CODES: [&str; 106] = [
    "11011001100",
    "11001101100",
    "11001100110",
    "10010011000",
    "10010001100",
    "10001001100",
    "10011001000",
    "10011000100",
    "10001100100",
    "11001001000",
    "11001000100",
    "11000100100",
    "10110011100",
    "10011011100",
    "10011001110",
    "10111001100",
    "10011101100",
    "10011100110",
    "11001110010",
    "11001011100",
    "11001001110",
    "11011100100",
    "11001110100",
    "11101101110",
    "11101001100",
    "11100101100",
    "11100100110",
    "11101100100",
    "11100110100",
    "11100110010",
    "11011011000",
    "11011000110",
    "11000110110",
    "10100011000",
    "10001011000",
    "10001000110",
    "10110001000",
    "10001101000",
    "10001100010",
    "11010001000",
    "11000101000",
    "11000100010",
    "10110111000",
    "10110001110",
    "10001101110",
    "10111011000",
    "10111000110",
    "10001110110",
    "11101110110",
    "11010001110",
    "11000101110",
    "11011101000",
    "11011100010",
    "11011101110",
    "11101011000",
    "11101000110",
    "11100010110",
    "11101101000",
    "11101100010",
    "11100011010",
    "11101111010",
    "11001000010",
    "11110001010",
    "10100110000",
    "10100001100",
    "10010110000",
    "10010000110",
    "10000101100",
    "10000100110",
    "10110010000",
    "10110000100",
    "10011010000",
    "10011000010",
    "10000110100",
    "10000110010",
    "11000010010",
    "11001010000",
    "11110111010",
    "11000010100",
    "10001111010",
    "10100111100",
    "10010111100",
    "10010011110",
    "10111100100",
    "10011110100",
    "10011110010",
    "11110100100",
    "11110010100",
    "11110010010",
    "11011011110",
    "11011110110",
    "11110110110",
    "10101111000",
    "10100011110",
    "10001011110",
    "10111101000",
    "10111100010",
    "11110101000",
    "11110100010",
    "10111011110",
    "10111101110",
    "11101011110",
    "11110101110",
    "11010000100",
    "11010010000",
    "11010011100",
];

/// Stop pattern (symbol value 106) including its trailing termination bar.
const STOP: &str = "1100011101011";

const CODE_C: usize = 99;
const CODE_B: usize = 100;
const START_C: usize = 105;

/// Encode `parts`, alternating numeric runs (even digit count, character set
/// C) and text runs (character set B), starting with a numeric run.
fn encode(parts: &[&str]) -> Vec<usize> {
    let mut codes = vec![START_C];
    let mut numeric = true;
    for (i, part) in parts.iter().enumerate() {
        if i > 0 {
            codes.push(if numeric { CODE_C } else { CODE_B });
        }
        if numeric {
            let digits = part.as_bytes();
            assert!(digits.len() % 2 == 0, "numeric run must have even length");
            for pair in digits.chunks(2) {
                codes.push(((pair[0] - b'0') * 10 + (pair[1] - b'0')) as usize);
            }
        } else {
            for &c in part.as_bytes() {
                codes.push((c - 32) as usize);
            }
        }
        numeric = !numeric;
    }

    // Modulo-103 check character: start value plus position-weighted data.
    let mut sum = codes[0];
    for (i, &c) in codes.iter().enumerate().skip(1) {
        sum += i * c;
    }
    codes.push(sum % 103);
    codes
}

fn render(codes: &[usize], module: u32, height: u32) -> GrayImage {
    let mut modules = String::new();
    for &c in codes {
        modules.push_str(CODES[c]);
    }
    modules.push_str(STOP);

    let quiet = 10 * module;
    let width = modules.len() as u32 * module + 2 * quiet;
    let mut img = GrayImage::from_pixel(width, height, Luma([255]));
    for (i, m) in modules.chars().enumerate() {
        if m == '1' {
            for dx in 0..module {
                for y in 0..height {
                    img.put_pixel(quiet + i as u32 * module + dx, y, Luma([0]));
                }
            }
        }
    }
    img
}

fn main() {
    for (path, parts, payload) in [
        (
            "examples/test-code128-setc.png",
            vec!["12345678"],
            "12345678",
        ),
        (
            "examples/test-code128-mixed.png",
            vec!["1234", "AB"],
            "1234AB",
        ),
        (
            "examples/test-code128-mixed2.png",
            vec!["1234", "AB", "5678"],
            "1234AB5678",
        ),
    ] {
        let img = render(&encode(&parts), 3, 60);
        img.save(path)
            .unwrap_or_else(|e| panic!("save {path}: {e}"));
        println!("  wrote {path} ({payload})");
    }
}
