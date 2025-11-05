//! Code 39 barcode decoder
//!
//! This module implements decoding for Code 39 barcodes.

use crate::{color::Color, finder::decode_e, img_scanner::zbar_image_scanner_t, SymbolType};
use libc::c_int;

// Number of characters in Code 39
const NUM_CHARS: usize = 0x2c;

// Assertion macro
macro_rules! zassert {
    ($condition:expr, $retval:expr, $($arg:tt)*) => {
        if !$condition {
            return $retval;
        }
    };
}

// ============================================================================
// Code 39 lookup tables
// ============================================================================

static CODE39_HI: [u8; 32] = [
    0x80,        // 2 next
    0x40 | 0x02, // 4
    0x80 | 0x06, // 2 next
    0xc0 | 0x08, // 2 skip
    0x40 | 0x0a, // 4
    0x80 | 0x0e, // 2 next
    0xc0 | 0x10, // 2 skip
    0x12,        // direct
    0x80 | 0x13, // 2 next
    0xc0 | 0x15, // 2 skip
    0x80 | 0x17, // 2 next
    0xff,
    0xc0 | 0x19, // 2 skip
    0x1b,        // direct
    0xff,
    0xff,
    0x40 | 0x1c, // 4
    0x80 | 0x20, // 2 next
    0xc0 | 0x22, // 2 skip
    0x24,        // direct
    0x80 | 0x25, // 2 next
    0xff,
    0x27, // direct
    0xff,
    0xc0 | 0x28, // 2 skip
    0x2a,        // direct
    0xff,
    0xff,
    0x2b, // direct
    0xff,
    0xff,
    0xff,
];

struct Char39 {
    chk: u8,
    rev: u8,
    fwd: u8,
}

static CODE39_ENCODINGS: [Char39; NUM_CHARS] = [
    Char39 {
        chk: 0x07,
        rev: 0x1a,
        fwd: 0x20,
    }, // 00
    Char39 {
        chk: 0x0d,
        rev: 0x10,
        fwd: 0x03,
    }, // 01
    Char39 {
        chk: 0x13,
        rev: 0x17,
        fwd: 0x22,
    }, // 02
    Char39 {
        chk: 0x16,
        rev: 0x1d,
        fwd: 0x23,
    }, // 03
    Char39 {
        chk: 0x19,
        rev: 0x0d,
        fwd: 0x05,
    }, // 04
    Char39 {
        chk: 0x1c,
        rev: 0x13,
        fwd: 0x06,
    }, // 05
    Char39 {
        chk: 0x25,
        rev: 0x07,
        fwd: 0x0c,
    }, // 06
    Char39 {
        chk: 0x2a,
        rev: 0x2a,
        fwd: 0x27,
    }, // 07
    Char39 {
        chk: 0x31,
        rev: 0x04,
        fwd: 0x0e,
    }, // 08
    Char39 {
        chk: 0x34,
        rev: 0x00,
        fwd: 0x0f,
    }, // 09
    Char39 {
        chk: 0x43,
        rev: 0x15,
        fwd: 0x25,
    }, // 0a
    Char39 {
        chk: 0x46,
        rev: 0x1c,
        fwd: 0x26,
    }, // 0b
    Char39 {
        chk: 0x49,
        rev: 0x0b,
        fwd: 0x08,
    }, // 0c
    Char39 {
        chk: 0x4c,
        rev: 0x12,
        fwd: 0x09,
    }, // 0d
    Char39 {
        chk: 0x52,
        rev: 0x19,
        fwd: 0x2b,
    }, // 0e
    Char39 {
        chk: 0x58,
        rev: 0x0f,
        fwd: 0x00,
    }, // 0f
    Char39 {
        chk: 0x61,
        rev: 0x02,
        fwd: 0x11,
    }, // 10
    Char39 {
        chk: 0x64,
        rev: 0x09,
        fwd: 0x12,
    }, // 11
    Char39 {
        chk: 0x70,
        rev: 0x06,
        fwd: 0x13,
    }, // 12
    Char39 {
        chk: 0x85,
        rev: 0x24,
        fwd: 0x16,
    }, // 13
    Char39 {
        chk: 0x8a,
        rev: 0x29,
        fwd: 0x28,
    }, // 14
    Char39 {
        chk: 0x91,
        rev: 0x21,
        fwd: 0x18,
    }, // 15
    Char39 {
        chk: 0x94,
        rev: 0x2b,
        fwd: 0x19,
    }, // 16
    Char39 {
        chk: 0xa2,
        rev: 0x28,
        fwd: 0x29,
    }, // 17
    Char39 {
        chk: 0xa8,
        rev: 0x27,
        fwd: 0x2a,
    }, // 18
    Char39 {
        chk: 0xc1,
        rev: 0x1f,
        fwd: 0x1b,
    }, // 19
    Char39 {
        chk: 0xc4,
        rev: 0x26,
        fwd: 0x1c,
    }, // 1a
    Char39 {
        chk: 0xd0,
        rev: 0x23,
        fwd: 0x1d,
    }, // 1b
    Char39 {
        chk: 0x03,
        rev: 0x14,
        fwd: 0x1e,
    }, // 1c
    Char39 {
        chk: 0x06,
        rev: 0x1b,
        fwd: 0x1f,
    }, // 1d
    Char39 {
        chk: 0x09,
        rev: 0x0a,
        fwd: 0x01,
    }, // 1e
    Char39 {
        chk: 0x0c,
        rev: 0x11,
        fwd: 0x02,
    }, // 1f
    Char39 {
        chk: 0x12,
        rev: 0x18,
        fwd: 0x21,
    }, // 20
    Char39 {
        chk: 0x18,
        rev: 0x0e,
        fwd: 0x04,
    }, // 21
    Char39 {
        chk: 0x21,
        rev: 0x01,
        fwd: 0x0a,
    }, // 22
    Char39 {
        chk: 0x24,
        rev: 0x08,
        fwd: 0x0b,
    }, // 23
    Char39 {
        chk: 0x30,
        rev: 0x05,
        fwd: 0x0d,
    }, // 24
    Char39 {
        chk: 0x42,
        rev: 0x16,
        fwd: 0x24,
    }, // 25
    Char39 {
        chk: 0x48,
        rev: 0x0c,
        fwd: 0x07,
    }, // 26
    Char39 {
        chk: 0x60,
        rev: 0x03,
        fwd: 0x10,
    }, // 27
    Char39 {
        chk: 0x81,
        rev: 0x1e,
        fwd: 0x14,
    }, // 28
    Char39 {
        chk: 0x84,
        rev: 0x25,
        fwd: 0x15,
    }, // 29
    Char39 {
        chk: 0x90,
        rev: 0x22,
        fwd: 0x17,
    }, // 2a
    Char39 {
        chk: 0xc0,
        rev: 0x20,
        fwd: 0x1a,
    }, // 2b
];

static CODE39_CHARACTERS: &[u8; NUM_CHARS] = b"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ-. $/+%*";

// ============================================================================
// Code 39 Decoder functions
// ============================================================================

/// Decode a single element
fn code39_decode1(enc: u8, e: u32, s: u32) -> u8 {
    let e_val = decode_e(e, s, 72);
    if !(0..=18).contains(&e_val) {
        return 0xff;
    }
    let mut enc = enc << 1;
    if e_val > 6 {
        enc |= 1;
    }
    enc
}

/// Decode 9 elements into a character (5 bars + 4 spaces or vice versa)
fn code39_decode9(dcode: &mut zbar_image_scanner_t) -> i8 {
    let s9 = dcode.code39.s9;

    if s9 < 9 {
        return -1;
    }

    let mut enc: u8 = 0;

    // Threshold bar width ratios for first 5 elements
    for i in 0..5 {
        enc = code39_decode1(enc, dcode.get_width(i), s9);
        if enc == 0xff {
            return -1;
        }
    }
    zassert!((enc as usize) < 0x20, -1, " enc={:x} s9={:x}\n", enc, s9);

    // Lookup first 5 encoded widths for coarse decode
    let mut idx = CODE39_HI[enc as usize];
    if idx == 0xff {
        return -1;
    }

    // Encode remaining widths (NB first encoded width is lost)
    for i in 5..9 {
        enc = code39_decode1(enc, dcode.get_width(i), s9);
        if enc == 0xff {
            return -1;
        }
    }

    if (idx & 0xc0) == 0x80 {
        idx = (idx & 0x3f) + ((enc >> 3) & 1);
    } else if (idx & 0xc0) == 0xc0 {
        idx = (idx & 0x3f) + ((enc >> 2) & 1);
    } else if (idx & 0xc0) != 0 {
        idx = (idx & 0x3f) + ((enc >> 2) & 3);
    }
    zassert!(
        (idx as usize) < 0x2c,
        -1,
        " idx={:x} enc={:x} s9={:x}\n",
        idx,
        enc,
        s9
    );

    let c = &CODE39_ENCODINGS[idx as usize];
    if enc != c.chk {
        return -1;
    }

    dcode.code39.width = s9;
    if dcode.code39.direction() {
        c.rev as i8
    } else {
        c.fwd as i8
    }
}

/// Decode start pattern
fn code39_decode_start(dcode: &mut zbar_image_scanner_t) -> SymbolType {
    let c = code39_decode9(dcode);
    if c != 0x19 && c != 0x2b {
        return SymbolType::None;
    }
    dcode
        .code39
        .set_direction(dcode.code39.direction() ^ (c == 0x19));

    // Check leading quiet zone - spec is 10x
    let quiet = dcode.get_width(9);
    if quiet != 0 && quiet < dcode.code39.s9 / 2 {
        return SymbolType::None;
    }

    dcode.code39.set_element(9);
    dcode.code39.set_character(0);
    SymbolType::Partial
}

/// Post-process decoded buffer
fn code39_postprocess(dcode: &mut zbar_image_scanner_t) -> i32 {
    let character = dcode.code39.character() as usize;
    let direction = dcode.code39.direction();

    dcode.direction = 1 - 2 * (direction as c_int);

    let buffer = &mut dcode.code39.buffer[..character];

    if direction {
        // Reverse buffer
        buffer.reverse();
    }

    for c in buffer.iter_mut() {
        *c = *CODE39_CHARACTERS.get(*c as usize).unwrap_or(&b'?');
    }

    let buffer = buffer.to_vec();
    dcode
        .buffer_mut_slice(character)
        .unwrap()
        .copy_from_slice(&buffer);

    dcode.modifiers = 0;
    0
}

/// Check width against reference
fn check_width(ref_width: u32, w: u32) -> bool {
    let dref = ref_width;
    let ref_4 = ref_width * 4;
    let w_4 = w * 4;
    ref_4.wrapping_sub(dref) <= w_4 && w_4 <= ref_4.wrapping_add(dref)
}

/// Main Code 39 decode function
pub(crate) fn _zbar_decode_code39(dcode: &mut zbar_image_scanner_t) -> SymbolType {
    // Update latest character width
    let w9 = dcode.get_width(9);
    let w0 = dcode.get_width(0);
    dcode.code39.s9 = dcode.code39.s9.wrapping_sub(w9).wrapping_add(w0);

    if dcode.code39.character() < 0 {
        if dcode.color() != Color::Bar {
            return SymbolType::None;
        }
        return code39_decode_start(dcode);
    }

    let element = dcode.code39.element();
    dcode.code39.set_element(element + 1);

    if dcode.code39.element() < 9 {
        return SymbolType::None;
    }

    if dcode.code39.element() == 10 {
        let space = dcode.get_width(0);
        let character = dcode.code39.character();

        // Check if STOP character is in the buffer
        let has_stop = if character > 0 {
            let buf_slice = &dcode.code39.buffer;
            buf_slice
                .get((character - 1) as usize)
                .is_some_and(|&b| b == 0x2b)
        } else {
            false
        };

        if has_stop {
            // STOP character found
            let mut sym = SymbolType::None;

            // Trim STOP character
            dcode.code39.set_character(character - 1);
            let character = dcode.code39.character();

            // Trailing quiet zone check
            if space != 0 && space < dcode.code39.width / 2 {
                // Failed quiet zone check
            } else {
                let (min_len, max_len) = dcode
                    .get_length_limits(SymbolType::Code39)
                    .unwrap_or((4, 0)); // Default: min=4, max=0 (unlimited)

                if character < min_len as i16 || (max_len > 0 && character > max_len as i16) {
                    // Failed length check
                } else if code39_postprocess(dcode) == 0 {
                    // FIXME checksum
                    sym = SymbolType::Code39;
                }
            }

            dcode.code39.set_character(-1);
            if sym == SymbolType::None {
                dcode.release_lock(SymbolType::Code39);
            }
            return sym;
        }

        if space > dcode.code39.width / 2 {
            // Inter-character space check failure
            if character > 0 {
                dcode.release_lock(SymbolType::Code39);
            }
            dcode.code39.set_character(-1);
        }
        dcode.code39.set_element(0);
        return SymbolType::None;
    }

    if !check_width(dcode.code39.width, dcode.code39.s9) {
        let character = dcode.code39.character();
        if character > 0 {
            dcode.release_lock(SymbolType::Code39);
        }
        dcode.code39.set_character(-1);
        return SymbolType::None;
    }

    let c = code39_decode9(dcode);

    let character = dcode.code39.character();

    // Lock shared resources
    if character == 0 && !dcode.acquire_lock(SymbolType::Code39) {
        dcode.code39.set_character(-1);
        return SymbolType::Partial;
    }

    if c < 0 {
        dcode.release_lock(SymbolType::Code39);
        dcode.code39.set_character(-1);
        return SymbolType::None;
    }

    zassert!(
        c < 0x2c,
        SymbolType::None,
        "c={:02x} s9={:x}\n",
        c,
        dcode.code39.s9
    );

    dcode.code39.set_byte(character as usize, c as u8);
    dcode.code39.set_character(character + 1);

    SymbolType::None
}
