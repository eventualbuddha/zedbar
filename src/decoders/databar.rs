use std::ptr::null_mut;

use libc::{c_int, c_uint};

use crate::{
    decoder::{_zbar_decoder_decode_e, decoder_realloc_databar_segments},
    decoder_types::{
        databar_decoder_t, databar_segment_t, zbar_decoder_t, zbar_symbol_type_t,
        ZBAR_CFG_EMIT_CHECK, ZBAR_CFG_ENABLE, ZBAR_DATABAR, ZBAR_DATABAR_EXP, ZBAR_MOD_GS1,
        ZBAR_NONE, ZBAR_PARTIAL,
    },
    line_scanner::zbar_color_t,
};

const DATABAR_MAX_SEGMENTS: usize = 32;
const GS: u8 = 0x1d;
const SCH_NUM: i32 = 0;
const SCH_ALNUM: i32 = 1;
const SCH_ISO646: i32 = 2;

#[inline]
fn var_max(len: i32, offset: i32) -> i32 {
    (((len * 12 + offset) * 2) + 6) / 7
}

#[inline]
unsafe fn push_char(buf: &mut *mut u8, ch: u8) {
    **buf = ch;
    *buf = (*buf).add(1);
}

#[inline]
unsafe fn push_char4(buf: &mut *mut u8, c0: u8, c1: u8, c2: u8, c3: u8) {
    push_char(buf, c0);
    push_char(buf, c1);
    push_char(buf, c2);
    push_char(buf, c3);
}

#[inline]
unsafe fn feed_bits(
    d: &mut u64,
    bit_count: &mut i32,
    len: &mut i32,
    data_ptr: &mut *mut i32,
    required: i32,
) {
    while *bit_count < required && *len > 0 {
        let next = (**data_ptr as u32 & 0x0fff) as u64;
        *data_ptr = (*data_ptr).add(1);
        *d = (*d << 12) | next;
        *bit_count += 12;
        *len -= 1;
    }
}

#[inline]
unsafe fn segment_index(seg: *mut crate::decoder_types::databar_segment_t) -> i32 {
    (((*seg).finder() as i32) << 2)
        | (((*seg).color() as i32) << 1)
        | ((((*seg).color() as u8 ^ (*seg).side()) as i32) & 1)
}

/// DataBar finder pattern hash table
static FINDER_HASH: [i8; 0x20] = [
    0x16, 0x1f, 0x02, 0x00, 0x03, 0x00, 0x06, 0x0b, 0x1f, 0x0e, 0x17, 0x0c, 0x0b, 0x14, 0x11, 0x0c,
    0x1f, 0x03, 0x13, 0x08, 0x00, 0x0a, -1, 0x16, 0x0c, 0x09, -1, 0x1a, 0x1f, 0x1c, 0x00, -1,
];

/// DataBar expanded sequences
static EXP_SEQUENCES: [u8; 30] = [
    // sequence Group 1
    0x01, 0x23, 0x25, 0x07, 0x29, 0x47, 0x29, 0x67, 0x0b, 0x29, 0x87, 0xab,
    // sequence Group 2
    0x21, 0x43, 0x65, 0x07, 0x21, 0x43, 0x65, 0x89, 0x21, 0x43, 0x65, 0xa9, 0x0b, 0x21, 0x43, 0x67,
    0x89, 0xab,
];

static EXP_CHECKSUMS: [u8; 12] = [1, 189, 62, 113, 46, 43, 109, 134, 6, 79, 161, 45];

struct GroupS {
    sum: u16,
    wmax: u8,
    todd: u8,
    teven: u8,
}

static GROUPS: [GroupS; 14] = [
    // (17,4) DataBar Expanded character groups
    GroupS {
        sum: 0,
        wmax: 7,
        todd: 87,
        teven: 4,
    },
    GroupS {
        sum: 348,
        wmax: 5,
        todd: 52,
        teven: 20,
    },
    GroupS {
        sum: 1388,
        wmax: 4,
        todd: 30,
        teven: 52,
    },
    GroupS {
        sum: 2948,
        wmax: 3,
        todd: 10,
        teven: 104,
    },
    GroupS {
        sum: 3988,
        wmax: 1,
        todd: 1,
        teven: 204,
    },
    // (16,4) DataBar outer character groups
    GroupS {
        sum: 0,
        wmax: 8,
        todd: 161,
        teven: 1,
    },
    GroupS {
        sum: 161,
        wmax: 6,
        todd: 80,
        teven: 10,
    },
    GroupS {
        sum: 961,
        wmax: 4,
        todd: 31,
        teven: 34,
    },
    GroupS {
        sum: 2015,
        wmax: 3,
        todd: 10,
        teven: 70,
    },
    GroupS {
        sum: 2715,
        wmax: 1,
        todd: 1,
        teven: 126,
    },
    // (15,4) DataBar inner character groups
    GroupS {
        sum: 1516,
        wmax: 8,
        todd: 81,
        teven: 1,
    },
    GroupS {
        sum: 1036,
        wmax: 6,
        todd: 48,
        teven: 10,
    },
    GroupS {
        sum: 336,
        wmax: 4,
        todd: 20,
        teven: 35,
    },
    GroupS {
        sum: 0,
        wmax: 2,
        todd: 4,
        teven: 84,
    },
];

/// Append checksum digit to 13-digit buffer
///
/// Calculates and appends a check digit to a 13-digit numeric string
/// using weighted sum modulo 10 algorithm.
///
/// # Parameters
/// - `buf`: Buffer containing 13 ASCII digits, with space for 14th digit
///
/// # Safety
/// Buffer must contain at least 14 bytes, with first 13 being ASCII digits '0'-'9'
pub unsafe fn _zbar_databar_append_check14(buf: *mut u8) {
    if buf.is_null() {
        return;
    }

    let mut chk: u8 = 0;
    let mut ptr = buf;

    for i in (0..13).rev() {
        let d = *ptr - b'0';
        chk = chk.wrapping_add(d);
        if (i & 1) == 0 {
            chk = chk.wrapping_add(d << 1);
        }
        ptr = ptr.add(1);
    }

    chk %= 10;
    if chk != 0 {
        chk = 10 - chk;
    }
    *ptr = chk + b'0';
}

/// Decode a number into decimal digits
///
/// Converts an unsigned long value into ASCII decimal digits and stores
/// them in the provided buffer.
///
/// # Parameters
/// - `buf`: Buffer to write digits to
/// - `n`: Number to decode
/// - `i`: Number of digits to write
///
/// # Safety
/// Buffer must have at least `i` bytes available
pub unsafe fn _zbar_databar_decode10(buf: *mut u8, mut n: u64, i: c_int) {
    if buf.is_null() || i <= 0 {
        return;
    }

    let mut pos = buf.add(i as usize);
    let mut remaining = i;

    while remaining > 0 {
        let d = (n % 10) as u8;
        n /= 10;
        pos = pos.sub(1);
        *pos = b'0' + d;
        remaining -= 1;
    }
}

unsafe fn databar_postprocess_exp(dcode: *mut zbar_decoder_t, data: *mut i32) -> c_int {
    if dcode.is_null() || data.is_null() {
        return -1;
    }

    let mut data_ptr = data;
    let first = *data_ptr as u64;
    data_ptr = data_ptr.add(1);
    let mut len = (first / 211 + 4) as i32;

    let mut d = *data_ptr as u64;
    data_ptr = data_ptr.add(1);

    let mut i_bits: i32;
    let enc: i32;
    let buflen: i32;

    let mut n = ((d >> 4) & 0x7f) as u32;
    if n >= 0x40 {
        i_bits = 10;
        enc = 1;
        buflen = 2 + 14 + var_max(len, 10 - 2 - 44 + 6) + 2;
    } else if n >= 0x38 {
        i_bits = 4;
        enc = 6 + (n as i32 & 7);
        buflen = 2 + 14 + 4 + 6 + 2 + 6 + 2;
    } else if n >= 0x30 {
        i_bits = 6;
        enc = 2 + (((n >> 2) & 1) as i32);
        buflen = 2 + 14 + 4 + 3 + var_max(len, 6 - 2 - 44 - 2 - 10) + 2;
    } else if n >= 0x20 {
        i_bits = 7;
        enc = 4 + (((n >> 3) & 1) as i32);
        buflen = 2 + 14 + 4 + 6;
    } else {
        i_bits = 9;
        enc = 0;
        buflen = var_max(len, 9 - 2) + 2;
    }

    if buflen <= 2 {
        return -1;
    }

    if enc < 4 {
        i_bits -= 1;
        let parity_bit = ((d >> i_bits) & 1) as i32;
        if ((len ^ parity_bit) & 1) != 0 {
            return -1;
        }

        i_bits -= 1;
        let size_group = ((d >> i_bits) & 1) as i32;
        if size_group != (len > 14) as i32 {
            return -1;
        }
    }

    len -= 2;

    if (*dcode).set_buffer_capacity(buflen as c_uint).is_err() {
        return -1;
    }

    let base_buf = (*dcode).buffer_mut_ptr() as *mut u8;
    let mut buf = base_buf;

    if enc != 0 {
        push_char(&mut buf, b'0');
        push_char(&mut buf, b'1');
    }

    if enc == 1 {
        i_bits -= 4;
        if i_bits >= 10 {
            return -1;
        }
        let digit = ((d >> i_bits) & 0xf) as u8;
        push_char(&mut buf, b'0' + digit);
    } else if enc != 0 {
        push_char(&mut buf, b'9');
    }

    if enc != 0 {
        for _ in 0..4 {
            feed_bits(&mut d, &mut i_bits, &mut len, &mut data_ptr, 10);
            i_bits -= 10;
            if i_bits < 0 {
                return -1;
            }
            n = ((d >> i_bits) & 0x3ff) as u32;
            if n >= 1000 {
                return -1;
            }
            _zbar_databar_decode10(buf, n as u64, 3);
            buf = buf.add(3);
        }
        _zbar_databar_append_check14(buf.sub(13));
        buf = buf.add(1);
    }

    match enc {
        2 => {
            feed_bits(&mut d, &mut i_bits, &mut len, &mut data_ptr, 2);
            i_bits -= 2;
            let val = ((d >> i_bits) & 0x3) as u8;
            push_char4(&mut buf, b'3', b'9', b'2', b'0' + val);
        }
        3 => {
            feed_bits(&mut d, &mut i_bits, &mut len, &mut data_ptr, 12);
            i_bits -= 2;
            let val = ((d >> i_bits) & 0x3) as u8;
            push_char4(&mut buf, b'3', b'9', b'3', b'0' + val);
            i_bits -= 10;
            if i_bits < 0 {
                return -1;
            }
            n = ((d >> i_bits) & 0x3ff) as u32;
            if n >= 1000 {
                return -1;
            }
            _zbar_databar_decode10(buf, n as u64, 3);
            buf = buf.add(3);
        }
        4 => {
            feed_bits(&mut d, &mut i_bits, &mut len, &mut data_ptr, 15);
            i_bits -= 15;
            if i_bits < 0 {
                return -1;
            }
            n = ((d >> i_bits) & 0x7fff) as u32;
            push_char4(&mut buf, b'3', b'1', b'0', b'3');
            _zbar_databar_decode10(buf, n as u64, 6);
            buf = buf.add(6);
        }
        5 => {
            feed_bits(&mut d, &mut i_bits, &mut len, &mut data_ptr, 15);
            i_bits -= 15;
            if i_bits < 0 {
                return -1;
            }
            n = ((d >> i_bits) & 0x7fff) as u32;
            let mut prefix = b'2';
            if n >= 10000 {
                prefix = b'3';
            }
            push_char4(&mut buf, b'3', b'2', b'0', prefix);
            if n >= 10000 {
                n -= 10000;
            }
            _zbar_databar_decode10(buf, n as u64, 6);
            buf = buf.add(6);
        }
        _ => {}
    }

    if enc >= 6 {
        let mut temp_enc = enc & 1;
        push_char4(&mut buf, b'3', b'1' + temp_enc as u8, b'0', b'x');

        feed_bits(&mut d, &mut i_bits, &mut len, &mut data_ptr, 20);
        i_bits -= 20;
        if i_bits < 0 {
            return -1;
        }
        n = ((d >> i_bits) & 0xfffff) as u32;
        if n >= 1_000_000 {
            return -1;
        }
        let buf_start = buf;
        _zbar_databar_decode10(buf_start, n as u64, 6);
        *buf_start.offset(-1) = *buf_start;
        *buf_start = b'0';
        buf = buf_start.add(6);

        feed_bits(&mut d, &mut i_bits, &mut len, &mut data_ptr, 16);
        i_bits -= 16;
        if i_bits < 0 {
            return -1;
        }
        n = ((d >> i_bits) & 0xffff) as u32;
        if n < 38400 {
            let dd = n % 32;
            let rem = n / 32;
            let mm = rem % 12 + 1;
            let yy = rem / 12;

            push_char(&mut buf, b'1');
            temp_enc = enc - 6;
            push_char(&mut buf, b'0' + ((temp_enc | 1) as u8));
            _zbar_databar_decode10(buf, yy as u64, 2);
            buf = buf.add(2);
            _zbar_databar_decode10(buf, mm as u64, 2);
            buf = buf.add(2);
            _zbar_databar_decode10(buf, dd as u64, 2);
            buf = buf.add(2);
        } else if n > 38400 {
            return -1;
        }
    }

    if enc < 4 {
        let mut scheme = SCH_NUM;
        while i_bits > 0 || len > 0 {
            feed_bits(&mut d, &mut i_bits, &mut len, &mut data_ptr, 8);

            if scheme == SCH_NUM {
                i_bits -= 4;
                if i_bits < 0 {
                    break;
                }
                if ((d >> i_bits) & 0xf) == 0 {
                    scheme = SCH_ALNUM;
                    continue;
                }

                if len == 0 && i_bits < 3 {
                    let digit = ((d >> i_bits) & 0xf) as i32 - 1;
                    if digit > 9 {
                        return -1;
                    }
                    push_char(&mut buf, b'0' + digit as u8);
                    break;
                }

                i_bits -= 3;
                if i_bits < 0 {
                    break;
                }
                let mut val = ((d >> i_bits) & 0x7f) as i32 - 8;
                let n1 = val % 11;
                val /= 11;
                push_char(&mut buf, if val < 10 { b'0' + val as u8 } else { GS });
                push_char(&mut buf, if n1 < 10 { b'0' + n1 as u8 } else { GS });
            } else {
                let mut ch: u8 = 0;
                i_bits -= 3;
                if i_bits < 0 {
                    break;
                }
                if ((d >> i_bits) & 0x7) == 0 {
                    scheme = SCH_NUM;
                    continue;
                }

                i_bits -= 2;
                if i_bits < 0 {
                    break;
                }
                let mut val = ((d >> i_bits) & 0x1f) as u32;
                if val == 0x04 {
                    scheme ^= 0x3;
                } else if val == 0x0f {
                    ch = GS;
                } else if val < 0x0f {
                    ch = 43 + val as u8;
                } else if scheme == SCH_ALNUM {
                    i_bits -= 1;
                    if i_bits < 0 {
                        return -1;
                    }
                    val = ((d >> i_bits) & 0x1f) as u32;
                    ch = if val < 0x1a {
                        b'A' + val as u8
                    } else if val == 0x1a {
                        b'*'
                    } else if val < 0x1f {
                        b',' + (val as u8) - 0x1b
                    } else {
                        return -1;
                    };
                } else if scheme == SCH_ISO646 && val < 0x1d {
                    i_bits -= 2;
                    if i_bits < 0 {
                        return -1;
                    }
                    val = ((d >> i_bits) & 0x3f) as u32;
                    ch = if val < 0x1a {
                        b'A' + val as u8
                    } else if val < 0x34 {
                        b'a' + (val as u8) - 0x1a
                    } else {
                        return -1;
                    };
                } else if scheme == SCH_ISO646 {
                    i_bits -= 3;
                    if i_bits < 0 {
                        return -1;
                    }
                    val = ((d >> i_bits) & 0x1f) as u32;
                    ch = match val {
                        0x00..=0x09 => b'!' + val as u8 - 8,
                        0x0a..=0x14 => b'%' + val as u8 - 0x0a,
                        0x15..=0x1a => b':' + val as u8 - 0x15,
                        0x1b => b'_',
                        0x1c => b' ',
                        _ => return -1,
                    };
                } else {
                    return -1;
                }

                if ch != 0 {
                    push_char(&mut buf, ch);
                }
            }
        }
    }

    let total_len = buf.offset_from(base_buf) as c_int;
    if total_len < 0 || (total_len as u32) >= (*dcode).buffer_capacity() {
        return -1;
    }

    *buf = 0;
    (*dcode).set_buffer_len(total_len as c_uint);
    if total_len > 0 {
        let last = buf.sub(1);
        if *last == GS {
            *last = 0;
            let new_len = (*dcode).buffer_len().saturating_sub(1);
            (*dcode).set_buffer_len(new_len);
        }
    }

    0
}

/// Convert DataBar data from heterogeneous base {1597,2841} to base 10 character representation
pub unsafe fn _zbar_databar_postprocess(dcode: *mut zbar_decoder_t, d: *mut c_uint) {
    let dcode = &mut *dcode;
    let mut d_array = [*d.offset(0), *d.offset(1), *d.offset(2), *d.offset(3)];

    // Get config before borrowing buffer
    let emit_check = ((dcode.databar.config >> ZBAR_CFG_EMIT_CHECK) & 1) != 0;

    let buf = match dcode.buffer_mut_slice(16) {
        Ok(buf) => buf,
        Err(_) => return,
    };
    let mut chk = 0u32;

    // Write "01" prefix
    buf[0] = b'0';
    buf[1] = b'1';

    // Start at position 15 and work backwards
    let mut buf_idx = 15;

    // Write two null terminators
    buf[buf_idx] = 0;
    buf_idx -= 1;
    buf[buf_idx] = 0;
    buf_idx -= 1;

    // First conversion
    let mut r = (d_array[0] as u64) * 1597 + (d_array[1] as u64);
    d_array[1] = (r / 10000) as c_uint;
    r %= 10000;
    r = r * 2841 + (d_array[2] as u64);
    d_array[2] = (r / 10000) as c_uint;
    r %= 10000;
    r = r * 1597 + (d_array[3] as u64);
    d_array[3] = (r / 10000) as c_uint;

    // Extract 4 decimal digits
    for i in (0..4).rev() {
        let c = (r % 10) as u32;
        chk += c;
        if (i & 1) != 0 {
            chk += c << 1;
        }
        buf[buf_idx] = c as u8 + b'0';
        buf_idx -= 1;
        if i != 0 {
            r /= 10;
        }
    }

    // Second conversion
    r = (d_array[1] as u64) * 2841 + (d_array[2] as u64);
    d_array[2] = (r / 10000) as c_uint;
    r %= 10000;
    r = r * 1597 + (d_array[3] as u64);
    d_array[3] = (r / 10000) as c_uint;

    // Extract 4 more decimal digits
    for i in (0..4).rev() {
        let c = (r % 10) as u32;
        chk += c;
        if (i & 1) != 0 {
            chk += c << 1;
        }
        buf[buf_idx] = c as u8 + b'0';
        buf_idx -= 1;
        if i != 0 {
            r /= 10;
        }
    }

    // Third conversion
    r = (d_array[2] as u64) * 1597 + (d_array[3] as u64);

    // Extract 5 decimal digits
    for i in (0..5).rev() {
        let c = (r % 10) as u32;
        chk += c;
        if (i & 1) == 0 {
            chk += c << 1;
        }
        buf[buf_idx] = c as u8 + b'0';
        buf_idx -= 1;
        if i != 0 {
            r /= 10;
        }
    }

    // Add check digit if configured
    if emit_check {
        chk %= 10;
        if chk != 0 {
            chk = 10 - chk;
        }
        buf[13] = chk as u8 + b'0';
        dcode.set_buffer_len(14);
    } else {
        dcode.set_buffer_len(13);
    }
}

/// Check if two widths are compatible within tolerance
///
/// Validates that two measured widths (wf and wd) are within an acceptable
/// tolerance range for n modules. Used to match DataBar segments.
///
/// # Parameters
/// - `wf`: First width measurement
/// - `wd`: Second width measurement
/// - `n`: Number of modules
///
/// # Returns
/// 1 if widths match within tolerance, 0 otherwise
pub fn _zbar_databar_check_width(wf: u32, wd: u32, n: u32) -> c_int {
    let dwf = wf * 3;
    let wd = wd * 14;
    let wf = wf * n;

    // Check: wf - dwf <= wd && wd <= wf + dwf
    // In C, this relies on unsigned wraparound if wf < dwf
    // For unsigned subtraction: wf - dwf will wrap if wf < dwf,
    // resulting in a very large number, making the condition false
    if wf.wrapping_sub(dwf) <= wd && wd <= wf + dwf {
        1
    } else {
        0
    }
}

/// Merge or update a DataBar segment with existing segments
pub unsafe fn _zbar_databar_merge_segment(
    db: *mut databar_decoder_t,
    seg: *mut crate::decoder_types::databar_segment_t,
) {
    let csegs = (*db).csegs() as isize;

    for i in 0..csegs {
        let s = ((*db).segs).offset(i);

        // Skip if this is the same segment
        if s == seg {
            continue;
        }

        // Check if this segment matches and should be merged
        if (*s).finder() == (*seg).finder()
            && (*s).exp() == (*seg).exp()
            && (*s).color() == (*seg).color()
            && (*s).side() == (*seg).side()
            && (*s).data == (*seg).data
            && (*s).check() == (*seg).check()
            && _zbar_databar_check_width((*seg).width as u32, (*s).width as u32, 14) != 0
        {
            // Found a matching segment - merge with it
            let mut cnt = (*s).count();
            if cnt < 0x7F {
                cnt += 1;
            }
            (*seg).set_count(cnt);

            // Merge partial flags (bitwise AND)
            let new_partial = (*seg).partial() && (*s).partial();
            (*seg).set_partial(new_partial);

            // Average the widths (weighted average favoring new measurement)
            let new_width = (3 * (*seg).width + (*s).width + 2) / 4;
            (*seg).width = new_width;

            // Mark old segment as unused
            (*s).set_finder(-1);
        } else if (*s).finder() >= 0 {
            // Not a match, check if we should age it out
            let age = (*db).epoch().wrapping_sub((*s).epoch());
            if age >= 248 || (age >= 128 && (*s).count() < 2) {
                (*s).set_finder(-1);
            }
        }
    }
}

/// Match DataBar segment to find complete symbol
pub unsafe fn match_segment(
    dcode: *mut zbar_decoder_t,
    seg: *mut crate::decoder_types::databar_segment_t,
) -> zbar_symbol_type_t {
    let db = &mut (*dcode).databar;
    let csegs = db.csegs();
    let mut maxage = 0xfff;
    let mut maxcnt = 0;
    let mut smax: [*mut crate::decoder_types::databar_segment_t; 3] = [std::ptr::null_mut(); 3];
    let mut d = [0u32; 4];

    if (*seg).partial() && (*seg).count() < 4 {
        return ZBAR_PARTIAL;
    }

    for i0 in 0..(csegs as usize) {
        let s0 = db.segs.add(i0);
        if s0 == seg
            || (*s0).finder() != (*seg).finder()
            || (*s0).exp()
            || (*s0).color() != (*seg).color()
            || (*s0).side() == (*seg).side()
            || ((*s0).partial() && (*s0).count() < 4)
            || _zbar_databar_check_width((*seg).width as c_uint, (*s0).width as c_uint, 14) == 0
        {
            continue;
        }

        for i1 in 0..(csegs as usize) {
            let s1 = db.segs.add(i1);
            if i1 == i0
                || (*s1).finder() < 0
                || (*s1).exp()
                || (*s1).color() == (*seg).color()
                || ((*s1).partial() && (*s1).count() < 4)
                || _zbar_databar_check_width((*seg).width as c_uint, (*s1).width as c_uint, 14) == 0
            {
                continue;
            }

            let mut chkf = if (*seg).color() != zbar_color_t::ZBAR_SPACE {
                (*seg).finder() as i32 + (*s1).finder() as i32 * 9
            } else {
                (*s1).finder() as i32 + (*seg).finder() as i32 * 9
            };
            if chkf > 72 {
                chkf -= 1;
            }
            if chkf > 8 {
                chkf -= 1;
            }

            let chks =
                (((*seg).check() as i32) + ((*s0).check() as i32) + ((*s1).check() as i32)) % 79;

            let chk = if chkf >= chks {
                chkf - chks
            } else {
                79 + chkf - chks
            };

            let age1 = ((db.epoch().wrapping_sub((*s0).epoch())) as u32)
                + ((db.epoch().wrapping_sub((*s1).epoch())) as u32);

            for i2 in (i1 + 1)..(csegs as usize) {
                let s2 = db.segs.add(i2);
                if i2 == i0
                    || (*s2).finder() != (*s1).finder()
                    || (*s2).exp()
                    || (*s2).color() != (*s1).color()
                    || (*s2).side() == (*s1).side()
                    || (*s2).check() as i32 != chk
                    || ((*s2).partial() && (*s2).count() < 4)
                    || _zbar_databar_check_width((*seg).width as c_uint, (*s2).width as c_uint, 14)
                        == 0
                {
                    continue;
                }
                let age2 = db.epoch().wrapping_sub((*s2).epoch()) as u32;
                let age = age1 + age2;
                let cnt = (*s0).count() as u32 + (*s1).count() as u32 + (*s2).count() as u32;
                if maxcnt < cnt as i32 || (maxcnt == cnt as i32 && (maxage as i32) > (age as i32)) {
                    maxcnt = cnt as i32;
                    maxage = age;
                    smax[0] = s0;
                    smax[1] = s1;
                    smax[2] = s2;
                }
            }
        }
    }

    if smax[0].is_null() {
        return ZBAR_PARTIAL;
    }

    d[(((*seg).color() as usize) << 1) | ((*seg).side() as usize)] = (*seg).data as u32;
    for i0 in 0..3 {
        d[(((*smax[i0]).color() as usize) << 1) | ((*smax[i0]).side() as usize)] =
            (*smax[i0]).data as u32;
        let new_count = (*smax[i0]).count().wrapping_sub(1);
        (*smax[i0]).set_count(new_count);
        if new_count == 0 {
            (*smax[i0]).set_finder(-1);
        }
    }
    (*seg).set_finder(-1);

    if (*dcode).set_buffer_capacity(18).is_err() {
        return ZBAR_PARTIAL;
    }

    if (*dcode)._zbar_decoder_acquire_lock(ZBAR_DATABAR) != 0 {
        return ZBAR_PARTIAL;
    }

    _zbar_databar_postprocess(dcode, d.as_mut_ptr());
    (*dcode).modifiers = 1 << ZBAR_MOD_GS1;
    (*dcode).direction = 1 - 2 * (((*seg).side() as i32) ^ ((*seg).color() as i32) ^ 1);
    ZBAR_DATABAR
}

/// Lookup DataBar expanded sequence
/// Returns -1 on error, 0 or 1 on success
pub unsafe fn lookup_sequence(
    seg: *mut crate::decoder_types::databar_segment_t,
    fixed: i32,
    seq: *mut i32,
    maxsize: usize,
) -> i32 {
    let mut n = ((*seg).data as u32 / 211) as usize;
    let mut i = n.div_ceil(2) + 1;
    n += 4;
    i = (i * i) / 4;
    let p = &EXP_SEQUENCES[i..];

    if n >= maxsize - 1 {
        // The loop below checks i<n and increments i by one within the loop
        // when accessing seq[22]. For this to be safe, n needs to be < 21.
        // See CVE-2023-40890.
        return -1;
    }

    let mut fixed = fixed >> 1;
    *seq.offset(0) = 0;
    *seq.offset(1) = 1;
    let mut i = 2;
    let mut p_idx = 0;
    while i < n {
        let mut s = p[p_idx] as i32;
        if (i & 2) == 0 {
            p_idx += 1;
            s >>= 4;
        } else {
            s &= 0xf;
        }
        if s == fixed {
            fixed = -1;
        }
        s <<= 1;
        *seq.add(i) = s;
        i += 1;
        s += 1;
        *seq.add(i) = s;
        i += 1;
    }
    *seq.add(n) = -1;
    if fixed < 1 {
        1
    } else {
        0
    }
}

pub unsafe fn match_segment_exp(
    dcode: *mut zbar_decoder_t,
    seg: *mut databar_segment_t,
    dir: c_int,
) -> zbar_symbol_type_t {
    let db = &mut (*dcode).databar;
    let csegs = db.csegs() as usize;
    if csegs == 0 {
        return ZBAR_PARTIAL;
    }

    let ifixed = seg.offset_from(db.segs) as usize;
    let fixed = segment_index(seg);
    let mut bestsegs = [-1i32; 22];
    let mut segs_idx = [-1i32; 22];
    let mut seq = [-1i32; 22];
    let mut iseg = [-1i32; DATABAR_MAX_SEGMENTS];
    let mut width_stack = [0u32; 22];

    seq[0] = 0;
    seq[1] = -1;
    segs_idx[0] = -1;
    bestsegs[0] = -1;
    width_stack[0] = (*seg).width as u32;

    #[allow(clippy::needless_range_loop)]
    for j in 0..csegs {
        let s = db.segs.add(j);
        iseg[j] =
            if (*s).exp() && (*s).finder() >= 0 && (!(*s).partial() || (*s).count() as i32 >= 4) {
                segment_index(s)
            } else {
                -1
            };
    }

    let mut maxcnt = 0i32;
    let mut maxage = 0x7fff_u32;
    let mut i: i32 = 0;
    let mut seg_ptr = seg;

    loop {
        while i >= 0 && seq[i as usize] >= 0 {
            let idx = i as usize;
            let target = seq[idx];
            let mut found: Option<usize> = None;
            let mut candidate: *mut databar_segment_t = std::ptr::null_mut();
            let current_width = width_stack[idx];

            if target == fixed {
                candidate = db.segs.add(ifixed);
                if segs_idx[idx] < 0
                    && _zbar_databar_check_width(current_width, (*candidate).width as u32, 14) != 0
                {
                    found = Some(ifixed);
                }
            } else {
                let mut start = if segs_idx[idx] >= 0 {
                    (segs_idx[idx] + 1) as usize
                } else {
                    0
                };
                while start < csegs {
                    if iseg[start] == target {
                        let cand = db.segs.add(start);
                        if idx == 0
                            || _zbar_databar_check_width(current_width, (*cand).width as u32, 14)
                                != 0
                        {
                            found = Some(start);
                            candidate = cand;
                            break;
                        }
                    }
                    start += 1;
                }
            }

            if let Some(jidx) = found {
                if candidate.is_null() {
                    candidate = db.segs.add(jidx);
                }

                if idx == 0 {
                    let lu = lookup_sequence(candidate, fixed, seq.as_mut_ptr(), seq.len());
                    if lu == 0 {
                        i -= 1;
                        continue;
                    }
                    if lu < 0 {
                        return ZBAR_NONE;
                    }
                }

                let next_width = if idx == 0 {
                    (*candidate).width as u32
                } else {
                    (current_width + (*candidate).width as u32) / 2
                };
                width_stack[idx + 1] = next_width;
                segs_idx[idx] = jidx as i32;
                segs_idx[idx + 1] = -1;
                seg_ptr = candidate;
                i = idx as i32 + 1;
            } else {
                i -= 1;
            }
        }

        if i < 0 {
            break;
        }

        let mut cnt = 0u32;
        let mut chk = 0u32;
        let mut age: u32;

        let first_idx = segs_idx[0] as usize;
        let mut seg_eval = db.segs.add(first_idx);
        age = ((*db).epoch().wrapping_sub((*seg_eval).epoch())) as u32 & 0xff;

        let mut pos = 1usize;
        while pos < segs_idx.len() && segs_idx[pos] >= 0 {
            let sidx = segs_idx[pos] as usize;
            seg_eval = db.segs.add(sidx);
            chk += (*seg_eval).check() as u32;
            cnt += (*seg_eval).count() as u32;
            age += ((*db).epoch().wrapping_sub((*seg_eval).epoch())) as u32 & 0xff;
            pos += 1;
        }

        let mut chk0 = (*db.segs.add(first_idx)).data as i32 % 211;
        if chk0 < 0 {
            chk0 += 211;
        }
        chk %= 211;
        if chk != chk0 as u32 {
            i -= 1;
            continue;
        }
        if maxcnt > cnt as i32 || (maxcnt == cnt as i32 && maxage <= age) {
            i -= 1;
            continue;
        }

        maxcnt = cnt as i32;
        maxage = age;
        bestsegs[..pos].copy_from_slice(&segs_idx[..pos]);
        bestsegs[pos] = -1;
        i -= 1;
    }

    if bestsegs[0] < 0 {
        return ZBAR_PARTIAL;
    }

    if (*dcode)._zbar_decoder_acquire_lock(ZBAR_DATABAR_EXP) != 0 {
        return ZBAR_PARTIAL;
    }

    let mut data_vals = [0i32; 22];
    let mut count = 0usize;
    while count < bestsegs.len() && bestsegs[count] >= 0 {
        let sidx = bestsegs[count] as usize;
        let s = db.segs.add(sidx);
        data_vals[count] = (*s).data as i32;
        count += 1;
    }

    if databar_postprocess_exp(dcode, data_vals.as_mut_ptr()) != 0 {
        (*dcode)._zbar_decoder_release_lock(ZBAR_DATABAR_EXP);
        return ZBAR_PARTIAL;
    }

    for sidx in bestsegs {
        let s = db.segs.add(sidx as usize);
        seg_ptr = s;
        if sidx as usize != ifixed {
            let mut cnt = (*s).count();
            if cnt > 0 {
                cnt -= 1;
                (*s).set_count(cnt);
                if cnt == 0 {
                    (*s).set_finder(-1);
                }
            }
        }
    }

    (*dcode).direction = (1 - 2 * (((*seg_ptr).side() ^ (*seg_ptr).color() as u8) as i32)) * dir;
    (*dcode).modifiers = 1 << ZBAR_MOD_GS1;
    ZBAR_DATABAR_EXP
}

/// Calculate DataBar checksum
///
/// Computes a checksum value for DataBar symbols based on signature values.
///
/// # Parameters
/// - `sig0`: First signature value (4 nibbles)
/// - `sig1`: Second signature value (4 nibbles)
/// - `side`: Side indicator (0 or 1)
/// - `mod_val`: Modulus value for checksum calculation
///
/// # Returns
/// Calculated checksum value
pub fn _zbar_databar_calc_check(mut sig0: u32, mut sig1: u32, side: u32, mod_val: u32) -> u32 {
    let mut chk: u32 = 0;

    for i in (0..4).rev() {
        chk = (chk * 3 + (sig1 & 0xf) + 1) * 3 + (sig0 & 0xf) + 1;
        sig1 >>= 4;
        sig0 >>= 4;
        if (i & 1) == 0 {
            chk %= mod_val;
        }
    }

    if side != 0 {
        chk = (chk * (6561 % mod_val)) % mod_val;
    }

    chk
}

/// Calculate DataBar character value from 4-element signature
/// Returns -1 on error
pub fn calc_value4(sig: c_uint, mut n: c_uint, wmax: c_uint, mut nonarrow: c_uint) -> c_int {
    let mut v = 0u32;
    n = n.wrapping_sub(1);

    let w0 = (sig >> 12) & 0xF;
    if w0 > 1 {
        if w0 > wmax {
            return -1;
        }
        let n0 = n.wrapping_sub(w0);
        let sk20 = n
            .wrapping_sub(1)
            .wrapping_mul(n)
            .wrapping_mul(n.wrapping_mul(2).wrapping_sub(1));
        let sk21 = n0
            .wrapping_mul(n0.wrapping_add(1))
            .wrapping_mul(n0.wrapping_mul(2).wrapping_add(1));
        v = sk20.wrapping_sub(sk21).wrapping_sub(
            w0.wrapping_sub(1)
                .wrapping_mul(3)
                .wrapping_mul(n.wrapping_mul(2).wrapping_sub(w0)),
        );

        if nonarrow == 0 && w0 > 2 && n > 4 {
            let mut k = n
                .wrapping_sub(2)
                .wrapping_mul(n.wrapping_sub(1))
                .wrapping_mul(n.wrapping_mul(2).wrapping_sub(3))
                .wrapping_sub(sk21);
            k = k.wrapping_sub(
                w0.wrapping_sub(2).wrapping_mul(3).wrapping_mul(
                    n.wrapping_mul(14)
                        .wrapping_sub(w0.wrapping_mul(7))
                        .wrapping_sub(31),
                ),
            );
            v = v.wrapping_sub(k);
        }

        if n.wrapping_sub(2) > wmax {
            let wm20 = wmax.wrapping_mul(2).wrapping_mul(wmax.wrapping_add(1));
            let wm21 = wmax.wrapping_mul(2).wrapping_add(1);
            let mut k = sk20;
            if n0 > wmax {
                k = k.wrapping_sub(sk21);
                k = k.wrapping_add(w0.wrapping_sub(1).wrapping_mul(3).wrapping_mul(
                    wm20.wrapping_sub(wm21.wrapping_mul(n.wrapping_mul(2).wrapping_sub(w0))),
                ));
            } else {
                k = k.wrapping_sub(
                    wmax.wrapping_add(1)
                        .wrapping_mul(wmax.wrapping_add(2))
                        .wrapping_mul(wmax.wrapping_mul(2).wrapping_add(3)),
                );
                k =
                    k.wrapping_add(
                        n.wrapping_sub(wmax)
                            .wrapping_sub(2)
                            .wrapping_mul(3)
                            .wrapping_mul(wm20.wrapping_sub(
                                wm21.wrapping_mul(n.wrapping_add(wmax).wrapping_add(1)),
                            )),
                    );
            }
            k = k.wrapping_mul(3);
            v = v.wrapping_sub(k);
        }
        v /= 12;
    } else {
        nonarrow = 1;
    }
    n = n.wrapping_sub(w0);

    let w1 = (sig >> 8) & 0xF;
    if w1 > 1 {
        if w1 > wmax {
            return -1;
        }
        v = v.wrapping_add(
            n.wrapping_mul(2)
                .wrapping_sub(w1)
                .wrapping_mul(w1.wrapping_sub(1))
                / 2,
        );
        if nonarrow == 0 && w1 > 2 && n > 3 {
            v = v.wrapping_sub(
                n.wrapping_mul(2)
                    .wrapping_sub(w1)
                    .wrapping_sub(5)
                    .wrapping_mul(w1.wrapping_sub(2))
                    / 2,
            );
        }
        if n.wrapping_sub(1) > wmax {
            if n.wrapping_sub(w1) > wmax {
                v = v.wrapping_sub(
                    w1.wrapping_sub(1).wrapping_mul(
                        n.wrapping_mul(2)
                            .wrapping_sub(w1)
                            .wrapping_sub(wmax.wrapping_mul(2)),
                    ),
                );
            } else {
                v = v.wrapping_sub(
                    n.wrapping_sub(wmax)
                        .wrapping_mul(n.wrapping_sub(wmax).wrapping_sub(1)),
                );
            }
        }
    } else {
        nonarrow = 1;
    }
    n = n.wrapping_sub(w1);

    let w2 = (sig >> 4) & 0xF;
    if w2 > 1 {
        if w2 > wmax {
            return -1;
        }
        v = v.wrapping_add(w2.wrapping_sub(1));
        if nonarrow == 0 && w2 > 2 && n > 2 {
            v = v.wrapping_sub(n.wrapping_sub(2));
        }
        if n > wmax {
            v = v.wrapping_sub(n.wrapping_sub(wmax));
        }
    } else {
        nonarrow = 1;
    }

    let w3 = sig & 0xF;
    if w3 == 1 {
        nonarrow = 1;
    } else if w3 > wmax {
        return -1;
    }

    if nonarrow == 0 {
        return -1;
    }

    v as c_int
}

/// Decode a DataBar character from width measurements
pub unsafe fn decode_char(
    dcode: *mut zbar_decoder_t,
    seg: *mut crate::decoder_types::databar_segment_t,
    off: c_int,
    dir: c_int,
) -> zbar_symbol_type_t {
    let db = &mut (*dcode).databar;
    let s = (*dcode).calc_s(if dir > 0 { off } else { off - 6 } as u8, 8);
    let mut emin = [0i32, 0i32];
    let mut sum = 0i32;
    let mut sig0 = 0u32;
    let mut sig1 = 0u32;

    let n = if (*seg).exp() {
        17
    } else if (*seg).side() != 0 {
        15
    } else {
        16
    };
    emin[1] = -(n as i32);

    if s < 13 || _zbar_databar_check_width((*seg).width as c_uint, s, n) == 0 {
        return ZBAR_NONE;
    }

    let mut off = off;
    for i in (0..4).rev() {
        let e = _zbar_decoder_decode_e((*dcode).pair_width(off as u8), s, n);
        if e < 0 {
            return ZBAR_NONE;
        }
        sum = e - sum;
        off += dir;
        sig1 <<= 4;
        if emin[1] < -sum {
            emin[1] = -sum;
        }
        sig1 = sig1.wrapping_add(sum as u32);
        if i == 0 {
            break;
        }

        let e = _zbar_decoder_decode_e((*dcode).pair_width(off as u8), s, n);
        if e < 0 {
            return ZBAR_NONE;
        }
        sum = e - sum;
        off += dir;
        sig0 <<= 4;
        if emin[0] > sum {
            emin[0] = sum;
        }
        sig0 = sig0.wrapping_add(sum as u32);
    }

    let mut diff = emin[(!(n as i32) & 1) as usize];
    diff = diff + (diff << 4);
    diff = diff + (diff << 8);

    sig0 = sig0.wrapping_sub(diff as u32);
    sig1 = sig1.wrapping_add(diff as u32);

    let mut sum0 = sig0.wrapping_add(sig0 >> 8);
    let mut sum1 = sig1.wrapping_add(sig1 >> 8);
    sum0 = sum0.wrapping_add(sum0 >> 4);
    sum1 = sum1.wrapping_add(sum1 >> 4);
    sum0 &= 0xf;
    sum1 &= 0xf;

    if sum0.wrapping_add(sum1).wrapping_add(8) as c_int != n as c_int {
        return ZBAR_NONE;
    }

    if ((sum0 ^ (n >> 1)) | (sum1 ^ (n >> 1) ^ n)) & 1 != 0 {
        return ZBAR_NONE;
    }

    let i = ((n & 0x3) ^ 1) * 5 + (sum1 >> 1);
    if i as usize >= GROUPS.len() {
        return ZBAR_NONE;
    }
    let g = &GROUPS[i as usize];

    let vodd = calc_value4(
        sig0.wrapping_add(0x1111),
        sum0.wrapping_add(4),
        g.wmax as c_uint,
        (!(n as i32) & 1) as c_uint,
    );
    if vodd < 0 || vodd > g.todd as i32 {
        return ZBAR_NONE;
    }

    let veven = calc_value4(
        sig1.wrapping_add(0x1111),
        sum1.wrapping_add(4),
        (9 - g.wmax) as c_uint,
        (n & 1) as c_uint,
    );
    if veven < 0 || veven > g.teven as i32 {
        return ZBAR_NONE;
    }

    let mut v = g.sum as i32;
    if (n & 2) != 0 {
        v = v
            .wrapping_add(vodd)
            .wrapping_add(veven.wrapping_mul(g.todd as i32));
    } else {
        v = v
            .wrapping_add(veven)
            .wrapping_add(vodd.wrapping_mul(g.teven as i32));
    }

    let mut chk;
    if (*seg).exp() {
        let side = (*seg).color() as u8 ^ (*seg).side() ^ 1;
        if v >= 4096 {
            return ZBAR_NONE;
        }
        chk = _zbar_databar_calc_check(sig0, sig1, side as c_uint, 211);
        if (*seg).finder() != 0 || (*seg).color() != zbar_color_t::ZBAR_SPACE || (*seg).side() != 0
        {
            let i = ((*seg).finder() as i32) * 2 - (side as i32) + ((*seg).color() as i32);
            if !(0..12).contains(&i) {
                return ZBAR_NONE;
            }
            chk = (chk * EXP_CHECKSUMS[i as usize] as u32) % 211;
        } else if v >= 4009 {
            return ZBAR_NONE;
        } else {
            chk = 0;
        }
    } else {
        chk = _zbar_databar_calc_check(sig0, sig1, (*seg).side() as c_uint, 79);
        if (*seg).color() != zbar_color_t::ZBAR_SPACE {
            chk = (chk * 16) % 79;
        }
    }

    (*seg).set_check(chk as u8);
    (*seg).data = v as i16;

    _zbar_databar_merge_segment(db, seg);

    if (*seg).exp() {
        return match_segment_exp(dcode, seg, dir);
    } else if dir > 0 {
        return match_segment(dcode, seg);
    }
    ZBAR_PARTIAL
}

/// Allocate a new DataBar segment (or reuse an old one)
/// Returns the index of the allocated segment, or -1 on failure
pub unsafe fn _zbar_databar_alloc_segment(db: *mut databar_decoder_t) -> c_int {
    let mut maxage = 0u32;
    let csegs = (*db).csegs() as usize;
    let mut old: c_int = -1;

    // First pass: look for empty slots or very old segments
    for i in 0..csegs {
        let seg = ((*db).segs).add(i);

        if (*seg).finder() < 0 {
            return i as c_int;
        }

        let age = (*db).epoch().wrapping_sub((*seg).epoch());
        if age >= 128 && (*seg).count() < 2 {
            (*seg).set_finder(-1);
            return i as c_int;
        }

        // Score based on both age and count
        let score = if age > (*seg).count() {
            age - (*seg).count() + 1
        } else {
            1
        };

        if maxage < score as u32 {
            maxage = score as u32;
            old = i as c_int;
        }
    }

    // Try to grow the segment array if not at max
    if csegs < DATABAR_MAX_SEGMENTS {
        let i = csegs;
        let mut new_csegs = csegs * 2;
        if new_csegs > DATABAR_MAX_SEGMENTS {
            new_csegs = DATABAR_MAX_SEGMENTS;
        }

        if new_csegs != csegs {
            // Reallocate segment array
            let new_ptr = decoder_realloc_databar_segments((*db).segs, new_csegs);

            if new_ptr.is_null() {
                // Allocation failed, fall through to reuse old segment
            } else {
                (*db).segs = new_ptr;
                (*db).set_csegs(new_csegs as u8);

                // Initialize new segments
                for j in i..new_csegs {
                    let seg = ((*db).segs).add(j);
                    (*seg).set_finder(-1);
                    (*seg).set_exp(false);
                    (*seg).set_color(zbar_color_t::ZBAR_SPACE);
                    (*seg).set_side(0);
                    (*seg).set_partial(false);
                    (*seg).set_count(0);
                    (*seg).set_epoch(0);
                    (*seg).set_check(0);
                }
                return i as c_int;
            }
        }
    }

    // Reuse oldest segment
    if old >= 0 {
        let seg = ((*db).segs).offset(old as isize);
        (*seg).set_finder(-1);
    }
    old
}

/// Decode DataBar finder pattern
pub unsafe fn decode_finder(dcode: *mut zbar_decoder_t) -> zbar_symbol_type_t {
    let db = &mut (*dcode).databar;
    let e0 = (*dcode).pair_width(1);
    let e2 = (*dcode).pair_width(3);
    let (dir, e2, e3) = if e0 < e2 {
        let e = e2 * 4;
        if e < 15 * e0 || e > 34 * e0 {
            return ZBAR_NONE;
        }
        (0, e2, (*dcode).pair_width(4))
    } else {
        let e = e0 * 4;
        if e < 15 * e2 || e > 34 * e2 {
            return ZBAR_NONE;
        }
        (1, e0, (*dcode).pair_width(0))
    };
    let e1 = (*dcode).pair_width(2);

    let s = e1 + e3;
    if s < 12 {
        return ZBAR_NONE;
    }

    let sig = (_zbar_decoder_decode_e(e3, s, 14) << 8)
        | (_zbar_decoder_decode_e(e2, s, 14) << 4)
        | _zbar_decoder_decode_e(e1, s, 14);
    if sig < 0
        || ((sig >> 4) & 0xf) < 8
        || ((sig >> 4) & 0xf) > 10
        || (sig & 0xf) >= 10
        || ((sig >> 8) & 0xf) >= 10
        || (((sig >> 8) + sig) & 0xf) != 10
    {
        return ZBAR_NONE;
    }

    let finder = (FINDER_HASH[((sig - (sig >> 5)) & 0x1f) as usize]
        + FINDER_HASH[((sig >> 1) & 0x1f) as usize])
        & 0x1f;
    if finder == 0x1f
        || (((if finder < 9 { db.config } else { db.config_exp }) >> ZBAR_CFG_ENABLE) & 1) == 0
    {
        return ZBAR_NONE;
    }

    if finder < 0 {
        return ZBAR_NONE;
    }

    let iseg = _zbar_databar_alloc_segment(db);
    if iseg < 0 {
        return ZBAR_NONE;
    }

    let seg = &mut (*db.segs.offset(iseg as isize));
    seg.set_finder(if finder >= 9 { finder - 9 } else { finder });
    seg.set_exp(finder >= 9);
    seg.set_color((((*dcode).color() as u8) ^ dir as u8 ^ 1).into());
    seg.set_side(dir as u8);
    seg.set_partial(false);
    seg.set_count(1);
    seg.width = s as i16;
    seg.set_epoch(db.epoch());

    let rc = decode_char(dcode, seg, 12 - dir, -1);
    if rc == 0 {
        seg.set_partial(true);
    } else {
        db.set_epoch(db.epoch().wrapping_add(1));
    }

    let i = (((*dcode).idx as c_int + 8 + dir) & 0xf) as usize;
    if db.chars[i] != -1 {
        return ZBAR_NONE;
    }
    db.chars[i] = iseg as i8;
    rc
}

pub unsafe fn _zbar_decode_databar(dcode: *mut zbar_decoder_t) -> zbar_symbol_type_t {
    let db = &mut (*dcode).databar;
    let i = (*dcode).idx & 0xf;

    let mut sym = decode_finder(dcode);

    let iseg = db.chars[i as usize];
    if iseg < 0 {
        return sym;
    }

    db.chars[i as usize] = -1;
    let mut seg = db.segs.offset(iseg as isize);

    let pair;
    if (*seg).partial() {
        pair = null_mut();
        (*seg).set_side(!(*seg).side());
    } else {
        let jseg = _zbar_databar_alloc_segment(db as *mut _);
        pair = db.segs.offset(iseg as isize);
        seg = db.segs.offset(jseg as isize);
        (*seg).set_finder((*pair).finder());
        (*seg).set_exp((*pair).exp());
        (*seg).set_color((*pair).color());
        (*seg).set_side(!(*pair).side());
        (*seg).set_partial(false);
        (*seg).set_count(1);
        (*seg).width = (*pair).width;
        (*seg).set_epoch(db.epoch());
    }

    sym = decode_char(dcode, seg, 1, 1);
    if sym == 0 {
        (*seg).set_finder(-1);
        if !pair.is_null() {
            (*pair).set_partial(true);
        }
    } else {
        db.set_epoch(db.epoch() + 1);
    }

    sym
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_check_width() {
        // Test exact match: wf=100, wd=100, n=14
        // After scaling: dwf=300, wd=1400, wf=1400
        // Check: (1400-300 <= 1400) && (1400 <= 1400+300) => (1100<=1400) && (1400<=1700) => true
        assert_eq!(_zbar_databar_check_width(100, 100, 14), 1);

        // Test within tolerance: wf=100, wd=105, n=14
        // After scaling: dwf=300, wd=1470, wf=1400
        // Check: (1400-300 <= 1470) && (1470 <= 1400+300) => (1100<=1470) && (1470<=1700) => true
        assert_eq!(_zbar_databar_check_width(100, 105, 14), 1);

        // Test within tolerance: wf=100, wd=95, n=14
        // After scaling: dwf=300, wd=1330, wf=1400
        // Check: (1400-300 <= 1330) && (1330 <= 1400+300) => (1100<=1330) && (1330<=1700) => true
        assert_eq!(_zbar_databar_check_width(100, 95, 14), 1);

        // Test outside tolerance: wf=100, wd=130, n=14
        // After scaling: dwf=300, wd=1820, wf=1400
        // Check: (1400-300 <= 1820) && (1820 <= 1400+300) => (1100<=1820) && (1820<=1700) => false
        assert_eq!(_zbar_databar_check_width(100, 130, 14), 0);

        // Test outside tolerance: wf=100, wd=70, n=14
        // After scaling: dwf=300, wd=980, wf=1400
        // Check: (1400-300 <= 980) && (980 <= 1400+300) => (1100<=980) && (980<=1700) => false
        assert_eq!(_zbar_databar_check_width(100, 70, 14), 0);
    }

    #[test]
    fn test_decode10() {
        unsafe {
            let mut buf = [0u8; 10];

            // Test simple number
            _zbar_databar_decode10(buf.as_mut_ptr(), 123, 3);
            assert_eq!(&buf[0..3], b"123");

            // Test with leading zeros
            _zbar_databar_decode10(buf.as_mut_ptr(), 45, 6);
            assert_eq!(&buf[0..6], b"000045");

            // Test zero
            _zbar_databar_decode10(buf.as_mut_ptr(), 0, 3);
            assert_eq!(&buf[0..3], b"000");
        }
    }

    #[test]
    fn test_append_check14() {
        unsafe {
            // Test with 13 ASCII digits
            let mut buf = [
                b'9', b'7', b'8', b'0', b'1', b'4', b'3', b'0', b'0', b'7', b'2', b'3', b'0', 0,
            ];
            _zbar_databar_append_check14(buf.as_mut_ptr());
            // Check digit: (9+7+8+0+1+4+3+0+0+7+2+3) + (9+8+1+3+0+2)*2 = 44 + 46 = 90, 90%10=0, check=0
            assert_eq!(buf[13], b'0');

            // Test another example: 1234567890120
            let mut buf2 = [
                b'1', b'2', b'3', b'4', b'5', b'6', b'7', b'8', b'9', b'0', b'1', b'2', b'0', 0,
            ];
            _zbar_databar_append_check14(buf2.as_mut_ptr());
            // Check: (1+2+3+4+5+6+7+8+9+0+1+2) + (1+3+5+7+9+1)*2 = 48 + 52 = 100, 100%10=0, check=0
            assert_eq!(buf2[13], b'0');
        }
    }

    #[test]
    fn test_calc_check() {
        // Test basic checksum calculation
        // These values are based on understanding the algorithm
        let chk1 = _zbar_databar_calc_check(0x1234, 0x5678, 0, 211);
        assert!(chk1 < 211);

        let chk2 = _zbar_databar_calc_check(0x1234, 0x5678, 1, 211);
        assert!(chk2 < 211);

        // Side should affect the result
        assert_ne!(chk1, chk2);

        // Test with different modulus
        let chk3 = _zbar_databar_calc_check(0x1234, 0x5678, 0, 79);
        assert!(chk3 < 79);

        // Same inputs should give same output
        let chk4 = _zbar_databar_calc_check(0x1234, 0x5678, 0, 211);
        assert_eq!(chk1, chk4);
    }
}
