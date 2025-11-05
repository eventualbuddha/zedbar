//! Interleaved 2 of 5 barcode decoder
//!
//! This module implements decoding for Interleaved 2 of 5 (I25) barcodes.

use crate::{color::Color, finder::decode_e, img_scanner::zbar_image_scanner_t, SymbolType};
use libc::c_int;

// ============================================================================
// I25 Decoder functions
// ============================================================================

/// Decode a single element
fn i25_decode1(enc: u8, e: u32, s: u32) -> u8 {
    let e_val = decode_e(e, s, 45);
    if e_val > 7 {
        return 0xff;
    }
    let mut enc = enc << 1;
    if e_val > 2 {
        enc |= 1;
    }
    enc
}

/// Decode 10 elements into a digit (5 bars + 5 spaces)
fn i25_decode10(dcode: &zbar_image_scanner_t, offset: u8) -> u8 {
    let dcode25 = &dcode.i25;

    if dcode25.s10 < 10 {
        return 0xff;
    }

    let mut enc: u8 = 0;
    let mut par: u8 = 0;

    // Threshold bar width ratios
    let mut i: i8 = 8;
    while i >= 0 {
        let j = offset
            + if dcode25.direction() {
                i as u8
            } else {
                8 - i as u8
            };
        enc = i25_decode1(enc, dcode.get_width(j), dcode25.s10);
        if enc == 0xff {
            return 0xff;
        }
        if enc & 1 != 0 {
            par += 1;
        }
        i -= 2;
    }

    // Parity check
    if par != 2 {
        return 0xff;
    }

    // Decode binary weights
    enc &= 0xf;
    if enc & 8 != 0 {
        if enc == 12 {
            enc = 0;
        } else {
            enc -= 1;
            if enc > 9 {
                return 0xff;
            }
        }
    }

    enc
}

/// Decode start pattern
fn i25_decode_start(dcode: &mut zbar_image_scanner_t) -> SymbolType {
    let s10 = dcode.i25.s10;

    if s10 < 10 {
        return SymbolType::None;
    }

    let mut enc: u8 = 0;
    let mut i: u8 = 10;

    enc = i25_decode1(enc, dcode.get_width(i), s10);
    i += 1;
    enc = i25_decode1(enc, dcode.get_width(i), s10);
    i += 1;
    enc = i25_decode1(enc, dcode.get_width(i), s10);
    i += 1;

    let valid = if dcode.color() == Color::Bar {
        enc == 4
    } else {
        enc = i25_decode1(enc, dcode.get_width(i), s10);
        i += 1;
        enc == 0
    };

    if !valid {
        return SymbolType::None;
    }

    // Check leading quiet zone - spec is 10n(?)
    // we require 5.25n for w=2n to 6.75n for w=3n
    let quiet = dcode.get_width(i);
    if quiet != 0 && quiet < s10 * 3 / 8 {
        return SymbolType::None;
    }

    dcode.i25.set_direction(dcode.color() != Color::Space);
    dcode.i25.set_element(1);
    dcode.i25.set_character(0);
    SymbolType::Partial
}

/// Acquire lock and copy holding buffer
fn i25_acquire_lock(dcode: &mut zbar_image_scanner_t) -> bool {
    // Lock shared resources
    if !dcode.acquire_lock(SymbolType::I25) {
        dcode.i25.set_character(-1);
        return true;
    }

    false
}

/// Decode end pattern and validate
fn i25_decode_end(dcode: &mut zbar_image_scanner_t) -> SymbolType {
    let width = dcode.i25.width;
    let direction = dcode.i25.direction();
    let character = dcode.i25.character();

    // Check trailing quiet zone
    let quiet = dcode.get_width(0);
    if (quiet != 0 && quiet < width * 3 / 8)
        || decode_e(dcode.get_width(1), width, 45) > 2
        || decode_e(dcode.get_width(2), width, 45) > 2
    {
        return SymbolType::None;
    }

    // Check exit condition
    let e = decode_e(dcode.get_width(3), width, 45);
    let valid = if !direction {
        (e as u32).wrapping_sub(3) <= 4
    } else {
        e <= 2 && decode_e(dcode.get_width(4), width, 45) <= 2
    };

    if !valid {
        return SymbolType::None;
    }

    if character <= 4 && i25_acquire_lock(dcode) {
        return SymbolType::Partial;
    }

    dcode.direction = 1 - 2 * (direction as c_int);

    let char_count = character as usize;

    // Get length limits from config
    let (min_len, max_len) = dcode.get_length_limits(SymbolType::I25).unwrap_or((6, 0)); // Default: min=6, max=0 (unlimited)

    let buffer = &mut dcode.i25.buffer;

    if direction {
        // Reverse buffer
        buffer[..char_count].reverse();
    }

    if character < min_len as i16 || (max_len > 0 && character > max_len as i16) {
        dcode.release_lock(SymbolType::I25);
        dcode.i25.set_character(-1);
        return SymbolType::None;
    }

    let buffer = buffer.clone();

    dcode
        .buffer_mut_slice(char_count)
        .unwrap()
        .copy_from_slice(&buffer[..char_count]);

    dcode.modifiers = 0;
    dcode.i25.set_character(-1);
    SymbolType::I25
}

/// Main I25 decode function
pub(crate) fn _zbar_decode_i25(dcode: &mut zbar_image_scanner_t) -> SymbolType {
    // Update latest character width
    let w10 = dcode.get_width(10);
    let w0 = dcode.get_width(0);
    dcode.i25.s10 = dcode.i25.s10.wrapping_sub(w10).wrapping_add(w0);

    if dcode.i25.character() < 0 && i25_decode_start(dcode) == SymbolType::None {
        return SymbolType::None;
    }

    let element = dcode.i25.element();
    dcode.i25.set_element(element.wrapping_sub(1));

    let check_element = 6u8.wrapping_sub(dcode.i25.direction() as u8);
    if dcode.i25.element() == check_element {
        return i25_decode_end(dcode);
    } else if dcode.i25.element() != 0 {
        return SymbolType::None;
    }

    // FIXME check current character width against previous
    dcode.i25.width = dcode.i25.s10;

    let _direction = dcode.i25.direction();
    let _character = dcode.i25.character();
    let _element = dcode.i25.element();

    if dcode.i25.character() == 4 && i25_acquire_lock(dcode) {
        return SymbolType::Partial;
    }

    let c = i25_decode10(dcode, 1);
    if c > 9 {
        // goto reset
        if dcode.i25.character() >= 4 {
            dcode.release_lock(SymbolType::I25);
        }
        dcode.i25.set_character(-1);
        return SymbolType::None;
    }

    let character = dcode.i25.character();

    dcode.i25.set_byte(character as usize, c + b'0');
    dcode.i25.set_character(character + 1);

    let c = i25_decode10(dcode, 0);
    if c > 9 {
        // goto reset
        if dcode.i25.character() >= 4 {
            dcode.release_lock(SymbolType::I25);
        }
        dcode.i25.set_character(-1);
        return SymbolType::None;
    }

    let character = dcode.i25.character();

    dcode.i25.set_byte(character as usize, c + b'0');
    dcode.i25.set_character(character + 1);
    dcode.i25.set_element(10);

    if dcode.i25.character() == 2 {
        SymbolType::Partial
    } else {
        SymbolType::None
    }
}
