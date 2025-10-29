//! Code 128 barcode decoder
//!
//! This module implements decoding for Code 128 barcodes.

use crate::{
    decoder::{
        code128_decoder_t, zbar_decoder_t, ZBAR_CFG_MAX_LEN, ZBAR_CFG_MIN_LEN, ZBAR_MOD_AIM,
        ZBAR_MOD_GS1,
    },
    line_scanner::zbar_color_t,
    SymbolType,
};
use libc::{c_int, c_uint};

// Character count
const NUM_CHARS: usize = 108;

// Code 128 character codes
const FNC3: u8 = 0x60;
const FNC2: u8 = 0x61;
const SHIFT: u8 = 0x62;
const CODE_C: u8 = 0x63;
// CODE_B = 0x64 - unused in decoder (part of spec)
const CODE_A: u8 = 0x65;
const FNC1: u8 = 0x66;
const START_A: u8 = 0x67;
// START_B = 0x68 - unused in decoder (part of spec)
const START_C: u8 = 0x69;
const STOP_FWD: u8 = 0x6a;
const STOP_REV: u8 = 0x6b;
// FNC4 = 0x6c - unused in decoder (FIXME: extended ASCII not implemented)

// Assertion macro
macro_rules! zassert {
    ($condition:expr, $retval:expr, $($arg:tt)*) => {
        if !$condition {
            return $retval;
        }
    };
}

// ============================================================================
// Code 128 lookup tables
// ============================================================================

static CHARACTERS: [u8; NUM_CHARS] = [
    0x5c, 0xbf, 0xa1, // [00] 00
    0x2a, 0xc5, 0x0c, 0xa4, // [03] 01
    0x2d, 0xe3, 0x0f, // [07] 02
    0x5f, 0xe4, // [0a] 03
    0x6b, 0xe8, 0x69, 0xa7, 0xe7, // [0c] 10
    0xc1, 0x51, 0x1e, 0x83, 0xd9, 0x00, 0x84, 0x1f, // [11] 11
    0xc7, 0x0d, 0x33, 0x86, 0xb5, 0x0e, 0x15, 0x87, // [19] 12
    0x10, 0xda, 0x11, // [21] 13
    0x36, 0xe5, 0x18, 0x37, // [24] 20
    0xcc, 0x13, 0x39, 0x89, 0x97, 0x14, 0x1b, 0x8a, 0x3a, 0xbd, // [28] 21
    0xa2, 0x5e, 0x01, 0x85, 0xb0, 0x02, 0xa3, // [32] 22
    0xa5, 0x2c, 0x16, 0x88, 0xbc, 0x12, 0xa6, // [39] 23
    0x61, 0xe6, 0x56, 0x62, // [40] 30
    0x19, 0xdb, 0x1a, // [44] 31
    0xa8, 0x32, 0x1c, 0x8b, 0xcd, 0x1d, 0xa9, // [47] 32
    0xc3, 0x20, 0xc4, // [4e] 33
    0x50, 0x5d, 0xc0, // [51] 0014 0025 0034
    0x2b, 0xc6, // [54] 0134 0143
    0x2e, // [56] 0243
    0x53, 0x60, // [57] 0341 0352
    0x31, // [59] 1024
    0x52, 0xc2, // [5a] 1114 1134
    0x34, 0xc8, // [5c] 1242 1243
    0x55, // [5e] 1441
    0x57, 0x3e, 0xce, // [5f] 4100 5200 4300
    0x3b, 0xc9, // [62] 4310 3410
    0x6a, // [64] 3420
    0x54, 0x4f, // [65] 1430 2530
    0x38, // [67] 4201
    0x58, 0xcb, // [68] 4111 4311
    0x2f, 0xca, // [6a] 2421 3421
];

static LO_BASE: [u8; 8] = [0x00, 0x07, 0x0c, 0x19, 0x24, 0x32, 0x40, 0x47];

static LO_OFFSET: [u8; 0x80] = [
    0xff, 0xf0, 0xff, 0x1f, 0xff, 0xf2, 0xff, 0xff, // 00 [00]
    0xff, 0xff, 0xff, 0x3f, 0xf4, 0xf5, 0xff, 0x6f, // 01
    0xff, 0xff, 0xff, 0xff, 0xf0, 0xf1, 0xff, 0x2f, // 02 [07]
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x3f, 0x4f, // 03
    0xff, 0x0f, 0xf1, 0xf2, 0xff, 0x3f, 0xff, 0xf4, // 10 [0c]
    0xf5, 0xf6, 0xf7, 0x89, 0xff, 0xab, 0xff, 0xfc, // 11
    0xff, 0xff, 0x0f, 0x1f, 0x23, 0x45, 0xf6, 0x7f, // 12 [19]
    0xff, 0xff, 0xff, 0xff, 0xf8, 0xff, 0xf9, 0xaf, // 13
    0xf0, 0xf1, 0xff, 0x2f, 0xff, 0xf3, 0xff, 0xff, // 20 [24]
    0x4f, 0x5f, 0x67, 0x89, 0xfa, 0xbf, 0xff, 0xcd, // 21
    0xf0, 0xf1, 0xf2, 0x3f, 0xf4, 0x56, 0xff, 0xff, // 22 [32]
    0xff, 0xff, 0x7f, 0x8f, 0x9a, 0xff, 0xbc, 0xdf, // 23
    0x0f, 0x1f, 0xf2, 0xff, 0xff, 0x3f, 0xff, 0xff, // 30 [40]
    0xf4, 0xff, 0xf5, 0x6f, 0xff, 0xff, 0xff, 0xff, // 31
    0x0f, 0x1f, 0x23, 0xff, 0x45, 0x6f, 0xff, 0xff, // 32 [47]
    0xf7, 0xff, 0xf8, 0x9f, 0xff, 0xff, 0xff, 0xff, // 33
];

// ============================================================================
// Helper functions from decoder.h
// ============================================================================

/// Fixed character width decode assist
#[inline]
fn decode_e(e: c_uint, s: c_uint, n: c_uint) -> i32 {
    let e_val = ((e * n * 2 + 1) / s).wrapping_sub(3) / 2;
    if e_val >= n - 3 {
        -1
    } else {
        e_val as i32
    }
}

/// Access config value by index
#[inline]
fn cfg(decoder: &code128_decoder_t, cfg: c_int) -> c_int {
    decoder.configs[(cfg - ZBAR_CFG_MIN_LEN) as usize]
}

// ============================================================================
// Code 128 Decoder functions
// ============================================================================

/// Decode low signature
#[inline]
fn decode_lo(sig: i32) -> i8 {
    let offset =
        (((sig >> 1) & 0x01) | ((sig >> 3) & 0x06) | ((sig >> 5) & 0x18) | ((sig >> 7) & 0x60))
            as u8;
    let mut idx = LO_OFFSET[offset as usize];

    if (sig & 1) != 0 {
        idx &= 0xf;
    } else {
        idx >>= 4;
    }
    if idx == 0xf {
        return -1;
    }

    let base = ((sig >> 11) | ((sig >> 9) & 1)) as u8;
    zassert!(
        base < 8,
        -1,
        "sig={:x} offset={:x} idx={:x} base={:x}\n",
        sig,
        offset,
        idx,
        base
    );
    let idx = idx + LO_BASE[base as usize];

    zassert!(
        idx <= 0x50,
        -1,
        "sig={:x} offset={:x} base={:x} idx={:x}\n",
        sig,
        offset,
        base,
        idx
    );
    CHARACTERS[idx as usize] as i8
}

/// Decode high signature
#[inline]
fn decode_hi(mut sig: i32) -> i8 {
    let mut rev = (sig & 0x4400) != 0;
    if rev {
        sig = ((sig >> 12) & 0x000f)
            | ((sig >> 4) & 0x00f0)
            | ((sig << 4) & 0x0f00)
            | ((sig << 12) & 0xf000);
    }

    let mut idx = match sig {
        0x0014 => 0x0,
        0x0025 => 0x1,
        0x0034 => 0x2,
        0x0134 => 0x3,
        0x0143 => 0x4,
        0x0243 => 0x5,
        0x0341 => 0x6,
        0x0352 => 0x7,
        0x1024 => 0x8,
        0x1114 => 0x9,
        0x1134 => 0xa,
        0x1242 => 0xb,
        0x1243 => 0xc,
        0x1441 => {
            rev = false;
            0xd
        }
        _ => return -1,
    };
    if rev {
        idx += 0xe;
    }
    CHARACTERS[0x51 + idx] as i8
}

/// Calculate check value
#[inline]
fn calc_check(c: u8) -> u8 {
    if (c & 0x80) == 0 {
        return 0x18;
    }
    let c = c & 0x7f;
    if c < 0x3d {
        return if c < 0x30 && c != 0x17 { 0x10 } else { 0x20 };
    }
    if c < 0x50 {
        return if c == 0x4d { 0x20 } else { 0x10 };
    }
    if c < 0x67 {
        0x20
    } else {
        0x10
    }
}

/// Decode 6 elements
#[inline]
fn decode6(dcode: &zbar_decoder_t) -> i8 {
    let s = dcode.code128.s6;
    if s < 5 {
        return -1;
    }

    // Build edge signature of character
    let sig = if dcode.color() == zbar_color_t::ZBAR_BAR {
        (decode_e(dcode.get_width(0) + dcode.get_width(1), s, 11) << 12)
            | (decode_e(dcode.get_width(1) + dcode.get_width(2), s, 11) << 8)
            | (decode_e(dcode.get_width(2) + dcode.get_width(3), s, 11) << 4)
            | decode_e(dcode.get_width(3) + dcode.get_width(4), s, 11)
    } else {
        (decode_e(dcode.get_width(5) + dcode.get_width(4), s, 11) << 12)
            | (decode_e(dcode.get_width(4) + dcode.get_width(3), s, 11) << 8)
            | (decode_e(dcode.get_width(3) + dcode.get_width(2), s, 11) << 4)
            | decode_e(dcode.get_width(2) + dcode.get_width(1), s, 11)
    };

    if sig < 0 {
        return -1;
    }

    // Lookup edge signature
    let c = if (sig & 0x4444) != 0 {
        decode_hi(sig)
    } else {
        decode_lo(sig)
    };

    if c == -1 {
        return -1;
    }

    // Character validation
    let bars = if dcode.color() == zbar_color_t::ZBAR_BAR {
        dcode.get_width(0) + dcode.get_width(2) + dcode.get_width(4)
    } else {
        dcode.get_width(1) + dcode.get_width(3) + dcode.get_width(5)
    };
    let bars = bars * 11 * 4 / s;
    let chk = calc_check(c as u8);

    if (chk as i32 - 7) > bars as i32 || bars as i32 > (chk as i32 + 7) {
        return -1;
    }

    (c & 0x7f) as i8
}

/// Validate checksum
#[inline]
fn validate_checksum(dcode: &zbar_decoder_t) -> bool {
    if dcode.code128.character() < 3 {
        return true;
    }

    let buf = dcode.buffer_slice();

    // Add in irregularly weighted start character
    let idx = if dcode.code128.direction() != 0 {
        (dcode.code128.character() - 1) as usize
    } else {
        0
    };
    let mut sum = buf[idx] as c_uint;
    if sum >= 103 {
        sum -= 103;
    }

    // Calculate sum in reverse to avoid multiply operations
    let mut acc: c_uint = 0;
    for i in (1..(dcode.code128.character() - 2) as usize).rev() {
        zassert!(
            sum < 103,
            true,
            "dir={:x} i={:x} sum={:x} acc={:x}\n",
            dcode.code128.direction(),
            i,
            sum,
            acc
        );
        let idx = if dcode.code128.direction() != 0 {
            (dcode.code128.character() as usize) - 1 - i
        } else {
            i
        };
        acc += buf[idx] as c_uint;
        if acc >= 103 {
            acc -= 103;
        }
        zassert!(
            acc < 103,
            true,
            "dir={:x} i={:x} sum={:x} acc={:x}\n",
            dcode.code128.direction(),
            i,
            sum,
            acc
        );
        sum += acc;
        if sum >= 103 {
            sum -= 103;
        }
    }

    // Compare to check character
    let idx = if dcode.code128.direction() != 0 {
        1
    } else {
        (dcode.code128.character() - 2) as usize
    };
    let check = buf[idx] as c_uint;
    sum != check
}

/// Expand and decode character set C
#[inline]
fn postprocess_c(dcode: &mut zbar_decoder_t, start: usize, end: usize, dst: usize) -> c_uint {
    // Expand buffer to accommodate 2x set C characters (2 digits per-char)
    let delta = end - start;
    let old_len = dcode.code128.character() as usize;
    let newlen = old_len + delta;
    if dcode.set_buffer_capacity(newlen as c_uint).is_err() {
        return SymbolType::None as c_uint;
    }

    let buf = match dcode.buffer_mut_slice(newlen) {
        Ok(buf) => buf,
        Err(_) => return SymbolType::None as c_uint,
    };

    // Relocate unprocessed data to end of buffer
    buf.copy_within(start..old_len, start + delta);

    for i in 0..delta {
        let j = dst + i * 2;
        // Convert each set C character into two ASCII digits
        let mut code = buf[start + delta + i];
        buf[j] = b'0';
        if code >= 50 {
            code -= 50;
            buf[j] += 5;
        }
        if code >= 30 {
            code -= 30;
            buf[j] += 3;
        }
        if code >= 20 {
            code -= 20;
            buf[j] += 2;
        }
        if code >= 10 {
            code -= 10;
            buf[j] += 1;
        }
        zassert!(
            buf[j] <= b'9',
            delta as c_uint,
            "start={:x} end={:x} i={:x} j={:x}\n",
            start,
            end,
            i,
            j
        );
        zassert!(
            code <= 9,
            delta as c_uint,
            "start={:x} end={:x} i={:x} j={:x}\n",
            start,
            end,
            i,
            j
        );
        buf[j + 1] = b'0' + code;
    }

    dcode.code128.set_character(newlen as i16);
    delta as c_uint
}

/// Resolve scan direction and convert to ASCII
#[inline]
fn postprocess(dcode: &mut zbar_decoder_t) -> bool {
    dcode.modifiers = 0;
    dcode.direction = 1 - 2 * (dcode.code128.direction() as c_int);
    let character_count = dcode.code128.character() as usize;
    let direction = dcode.code128.direction();

    // First phase: reverse buffer if needed and validate
    {
        let buf = match dcode.buffer_mut_slice(character_count) {
            Ok(buf) => buf,
            Err(_) => return true,
        };

        if direction != 0 {
            // Reverse buffer
            let half = character_count / 2;
            for i in 0..half {
                let j = character_count - 1 - i;
                buf.swap(i, j);
            }
            zassert!(
                buf[character_count - 1] == STOP_REV,
                true,
                "dir={:x}\n",
                direction
            );
        } else {
            zassert!(
                buf[character_count - 1] == STOP_FWD,
                true,
                "dir={:x}\n",
                direction
            );
        }

        let code = buf[0];
        zassert!(
            (START_A..=START_C).contains(&code),
            true,
            "code={:x}\n",
            code
        );
    } // buf borrow ends here

    // Second phase: convert characters
    let start_code = {
        let buf = match dcode.buffer_mut_slice(character_count) {
            Ok(buf) => buf,
            Err(_) => return true,
        };
        buf[0]
    };

    let mut charset = start_code - START_A;
    let mut cexp = if start_code == START_C { 1 } else { 0 };

    let mut i = 1usize;
    let mut j = 0usize;
    while i < (character_count - 2) {
        let code = {
            let buf = match dcode.buffer_mut_slice(character_count.max(j + 1)) {
                Ok(buf) => buf,
                Err(_) => return true,
            };
            buf[i]
        };

        zassert!(
            (code & 0x80) == 0,
            true,
            "i={:x} j={:x} code={:02x} charset={:x} cexp={:x}\n",
            i,
            j,
            code,
            charset,
            cexp
        );

        if (charset & 0x2) != 0 && code < 100 {
            // Defer character set C for expansion
            i += 1;
            continue;
        } else if code < 0x60 {
            // Convert character set B to ASCII
            let mut ascii_code = code + 0x20;
            if (charset == 0 || charset == 0x81) && ascii_code >= 0x60 {
                // Convert character set A to ASCII
                ascii_code -= 0x60;
            }

            let buf = match dcode.buffer_mut_slice(character_count.max(j + 1)) {
                Ok(buf) => buf,
                Err(_) => return true,
            };
            buf[j] = ascii_code;
            j += 1;
            if (charset & 0x80) != 0 {
                charset &= 0x7f;
            }
        } else {
            if (charset & 0x2) != 0 {
                // Expand character set C to ASCII
                zassert!(
                    cexp != 0,
                    true,
                    "i={:x} j={:x} code={:02x} charset={:x} cexp={:x}\n",
                    i,
                    j,
                    code,
                    charset,
                    cexp
                );
                let delta = postprocess_c(dcode, cexp, i, j);
                // Buffer was modified by postprocess_c, no need to re-acquire
                i += delta as usize;
                j += (delta * 2) as usize;
                cexp = 0;
            }
            if code < CODE_C {
                if code == SHIFT {
                    charset |= 0x80;
                } else if code == FNC2 {
                    // FIXME FNC2 - message append
                } else if code == FNC3 {
                    // FIXME FNC3 - initialize
                }
            } else if code == FNC1 {
                // FNC1 - Code 128 subsets or ASCII 0x1d
                if i == 1 {
                    dcode.modifiers |= 1 << ZBAR_MOD_GS1;
                } else if i == 2 {
                    dcode.modifiers |= 1 << ZBAR_MOD_AIM;
                } else if i < (character_count - 3) {
                    let buf = match dcode.buffer_mut_slice(j + 1) {
                        Ok(buf) => buf,
                        Err(_) => return true,
                    };
                    buf[j] = 0x1d;
                    j += 1;
                }
                // else drop trailing FNC1
            } else if code >= START_A {
                return true;
            } else {
                let newset = CODE_A - code;
                zassert!(
                    (CODE_C..=CODE_A).contains(&code),
                    true,
                    "i={:x} j={:x} code={:02x} charset={:x} cexp={:x}\n",
                    i,
                    j,
                    code,
                    charset,
                    cexp
                );
                if newset != charset {
                    charset = newset;
                } else {
                    // FIXME FNC4 - extended ASCII
                }
            }
            if (charset & 0x2) != 0 {
                cexp = i + 1;
            }
        }
        i += 1;
    }

    if (charset & 0x2) != 0 {
        zassert!(
            cexp != 0,
            true,
            "i={:x} j={:x} code={:02x} charset={:x} cexp={:x}\n",
            i,
            j,
            0u8, // code is out of scope here
            charset,
            cexp
        );
        let delta = postprocess_c(dcode, cexp, i, j);
        j += (delta * 2) as usize;
    }

    zassert!(
        (j as c_uint) < dcode.buffer_capacity(),
        true,
        "j={:02x}\n",
        j
    );
    dcode.set_buffer_len(j as c_uint);

    // Write null terminator
    let buf = match dcode.buffer_mut_slice(j + 1) {
        Ok(buf) => buf,
        Err(_) => return true,
    };
    buf[j] = 0;

    dcode.code128.set_character(j as i16);
    false
}

/// Main Code 128 decode function
pub(crate) fn _zbar_decode_code128(dcode: &mut zbar_decoder_t) -> SymbolType {
    // Update latest character width
    dcode.code128.s6 = dcode
        .code128
        .s6
        .wrapping_sub(dcode.get_width(6))
        .wrapping_add(dcode.get_width(0));

    if (dcode.code128.character() < 0 && dcode.color() != zbar_color_t::ZBAR_SPACE)
        || (dcode.code128.character() >= 0
            && (dcode.code128.element() + 1 != 6
                || dcode.color() as u8 != dcode.code128.direction()))
    {
        if dcode.code128.character() >= 0 {
            dcode.code128.set_element(dcode.code128.element() + 1);
        }
        return SymbolType::None;
    }
    dcode.code128.set_element(0);

    let c = decode6(dcode);

    if dcode.code128.character() < 0 {
        let qz = dcode.get_width(6);
        if c < START_A as i8 || c > STOP_REV as i8 || c == STOP_FWD as i8 {
            return SymbolType::None;
        }
        if qz != 0 && qz < (dcode.code128.s6 * 3) / 4 {
            return SymbolType::None;
        }
        // Decoded valid start/stop - initialize state
        dcode.code128.set_character(1);
        if c == STOP_REV as i8 {
            dcode.code128.set_direction(zbar_color_t::ZBAR_BAR as u8);
            dcode.code128.set_element(7);
        } else {
            dcode.code128.set_direction(zbar_color_t::ZBAR_SPACE as u8);
        }
        dcode.code128.set_start(c as u8);
        dcode.code128.width = dcode.code128.s6;
        return SymbolType::None;
    } else if c < 0
        || dcode
            .set_buffer_capacity((dcode.code128.character() + 1) as c_uint)
            .is_err()
    {
        if dcode.code128.character() > 1 {
            dcode.release_lock(SymbolType::Code128);
        }
        dcode.code128.set_character(-1);
        return SymbolType::None;
    } else {
        let dw = dcode.code128.width.abs_diff(dcode.code128.s6);
        let dw = dw * 4;
        if dw > dcode.code128.width {
            if dcode.code128.character() > 1 {
                dcode.release_lock(SymbolType::Code128);
            }
            dcode.code128.set_character(-1);
            return SymbolType::None;
        }
    }
    dcode.code128.width = dcode.code128.s6;

    let capacity = dcode.buffer_capacity();
    zassert!(
        (capacity as i16) > dcode.code128.character(),
        SymbolType::None,
        "alloc={:x} idx={:x} c={:02x}\n",
        capacity,
        dcode.code128.character(),
        c
    );

    if dcode.code128.character() == 1 {
        // Lock shared resources
        if dcode.acquire_lock(SymbolType::Code128) {
            dcode.code128.set_character(-1);
            return SymbolType::None;
        }
        let start = dcode.code128.start();
        if let Ok(buf) = dcode.buffer_mut_slice(1) {
            buf[0] = start;
        }
    }

    let character = dcode.code128.character();
    if let Ok(buf) = dcode.buffer_mut_slice((character + 1) as usize) {
        buf[character as usize] = c as u8;
    }
    dcode.code128.set_character(character + 1);

    if dcode.code128.character() > 2
        && ((dcode.code128.direction() != 0 && c >= START_A as i8 && c <= START_C as i8)
            || (dcode.code128.direction() == 0 && c == STOP_FWD as i8))
    {
        // FIXME STOP_FWD should check extra bar (and QZ!)
        let mut sym = SymbolType::Code128;
        #[allow(clippy::if_same_then_else)]
        if validate_checksum(dcode) || postprocess(dcode) {
            sym = SymbolType::None;
        } else if dcode.code128.character() < cfg(&dcode.code128, ZBAR_CFG_MIN_LEN) as i16
            || (cfg(&dcode.code128, ZBAR_CFG_MAX_LEN) > 0
                && dcode.code128.character() > cfg(&dcode.code128, ZBAR_CFG_MAX_LEN) as i16)
        {
            sym = SymbolType::None;
        }
        dcode.code128.set_character(-1);
        if sym == SymbolType::None {
            dcode.release_lock(SymbolType::Code128);
        }
        return sym;
    }
    SymbolType::None
}
