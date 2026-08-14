use crate::{
    Error, Result, SymbolType,
    color::Color,
    decoder::{DatabarDecoder, DatabarSegment, Modifier, decode_e},
    img_scanner::{BUFFER_MIN, ImageScanner},
};

const DATABAR_MAX_SEGMENTS: usize = 32;
const GS: u8 = 0x1d;

/// Encoding scheme for DataBar Expanded data
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DatabarScheme {
    Numeric = 0,
    Alphanumeric = 1,
    Iso646 = 2,
}

fn var_max(len: i32, offset: i32) -> i32 {
    (((len * 12 + offset) * 2) + 6) / 7
}

/// Sliding-window reader over the 12-bit character stream of a DataBar
/// Expanded symbol.
///
/// Mirrors the C decoder's `FEED_BITS` macro: `window` accumulates characters
/// 12 bits at a time and `bits` counts how many of them are still unconsumed.
/// [`advance`](Self::advance) plus [`peek`](Self::peek) reproduce the C idiom
/// `i -= n; (d >> i) & mask` — note that C often peeks a *wider* field than it
/// just advanced past, deliberately re-reading bits.
struct BitReader<'a> {
    data: &'a [i32],
    pos: usize,
    window: u64,
    bits: i32,
    /// Characters left to feed into the window.
    remaining: i32,
}

impl<'a> BitReader<'a> {
    fn new(data: &'a [i32], window: u64, bits: i32, remaining: i32) -> Self {
        Self {
            data,
            pos: 0,
            window,
            bits,
            remaining,
        }
    }

    /// Pull characters in until at least `required` bits are buffered.
    fn feed(&mut self, required: i32) {
        while self.bits < required && self.remaining > 0 {
            let next = self.data.get(self.pos).copied().unwrap_or(0) as u32 & 0x0fff;
            self.pos += 1;
            self.window = (self.window << 12) | u64::from(next);
            self.bits += 12;
            self.remaining -= 1;
        }
    }

    /// Move the read cursor down by `n` bits.
    ///
    /// Returns `None` once the window has run dry; that is where the C decoder
    /// tests `i < 0` (and, at the sites where it forgets to, invokes undefined
    /// behavior by shifting by a negative amount).
    fn advance(&mut self, n: i32) -> Option<()> {
        self.bits -= n;
        (self.bits >= 0).then_some(())
    }

    /// Read the `n` bits at the current cursor without moving it.
    fn peek(&self, n: i32) -> u64 {
        debug_assert!(self.bits >= 0 && (1..=63).contains(&n));
        (self.window >> self.bits) & ((1u64 << n) - 1)
    }

    /// [`feed`](Self::feed) then [`advance`](Self::advance) then
    /// [`peek`](Self::peek) the same width — the common case.
    fn read(&mut self, n: i32) -> Option<u64> {
        self.feed(n);
        self.advance(n)?;
        Some(self.peek(n))
    }
}

/// Bounds-checked write cursor over the decoder's output buffer.
struct OutBuf<'a> {
    buf: &'a mut [u8],
    len: usize,
}

impl OutBuf<'_> {
    /// Reserve `n` bytes at the cursor and advance past them.
    fn take_mut(&mut self, n: usize) -> Result<&mut [u8]> {
        let end = self.len.checked_add(n).ok_or(Error::Invalid)?;
        let slice = self.buf.get_mut(self.len..end).ok_or(Error::Invalid)?;
        self.len = end;
        Ok(slice)
    }

    fn push(&mut self, byte: u8) -> Result<()> {
        self.take_mut(1)?[0] = byte;
        Ok(())
    }

    fn push_all(&mut self, bytes: &[u8]) -> Result<()> {
        self.take_mut(bytes.len())?.copy_from_slice(bytes);
        Ok(())
    }

    /// Append `digits` zero-padded decimal digits of `n`.
    fn push_decimal(&mut self, n: u64, digits: usize) -> Result<()> {
        decode10(self.take_mut(digits)?, n, digits);
        Ok(())
    }
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
fn append_check14(buf: &mut [u8]) {
    let mut chk: u8 = 0;
    let mut ptr = buf;

    for i in (0..13).rev() {
        let d = ptr[0] - b'0';
        chk = chk.wrapping_add(d);
        if (i & 1) == 0 {
            chk = chk.wrapping_add(d << 1);
        }
        ptr = &mut ptr[1..];
    }

    chk %= 10;
    if chk != 0 {
        chk = 10 - chk;
    }
    ptr[0] = chk + b'0';
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
fn decode10(buf: &mut [u8], mut n: u64, mut i: usize) {
    let mut remaining = i;

    while remaining > 0 {
        let d = (n % 10) as u8;
        n /= 10;
        i -= 1;
        buf[i] = b'0' + d;
        remaining -= 1;
    }
}

/// Decode the general-purpose data compaction stream that follows the
/// compressed AI fields (C: the `enc < 4` tail of `databar_postprocess_exp`).
fn postprocess_exp_tail(r: &mut BitReader<'_>, out: &mut OutBuf<'_>) -> Result<()> {
    let mut scheme = DatabarScheme::Numeric;

    while r.bits > 0 || r.remaining > 0 {
        r.feed(8);

        if scheme == DatabarScheme::Numeric {
            if r.advance(4).is_none() {
                break;
            }
            if r.peek(4) == 0 {
                scheme = DatabarScheme::Alphanumeric;
                continue;
            }

            if r.remaining == 0 && r.bits < 3 {
                // special case last digit
                let digit = r.peek(4) - 1;
                if digit > 9 {
                    return Err(Error::Invalid);
                }
                out.push(b'0' + digit as u8)?;
                break;
            }

            if r.advance(3).is_none() {
                break;
            }
            // Re-reads the 4 bits above plus the 3 just consumed.
            let val = r.peek(7) as i32 - 8;
            let low = val % 11;
            let high = val / 11;
            out.push(if high < 10 { b'0' + high as u8 } else { GS })?;
            out.push(if low < 10 { b'0' + low as u8 } else { GS })?;
            continue;
        }

        // Alphanumeric / ISO-646
        if r.advance(3).is_none() {
            break;
        }
        if r.peek(3) == 0 {
            scheme = DatabarScheme::Numeric;
            continue;
        }

        if r.advance(2).is_none() {
            break;
        }
        let val = r.peek(5);
        let ch: u8 = if val == 0x04 {
            // C toggles with `scheme ^= 0x3`, swapping ALNUM <-> ISO646.
            scheme = if scheme == DatabarScheme::Alphanumeric {
                DatabarScheme::Iso646
            } else {
                DatabarScheme::Alphanumeric
            };
            0
        } else if val == 0x0f {
            GS
        } else if val < 0x0f {
            43 + val as u8
        } else if scheme == DatabarScheme::Alphanumeric {
            r.advance(1).ok_or(Error::Invalid)?;
            match r.peek(5) {
                v if v < 0x1a => b'A' + v as u8,
                0x1a => b'*',
                v if v < 0x1f => b',' + v as u8 - 0x1b,
                _ => return Err(Error::Invalid),
            }
        } else if val < 0x1d {
            r.advance(2).ok_or(Error::Invalid)?;
            match r.peek(6) {
                v if v < 0x1a => b'A' + v as u8,
                v if v < 0x34 => b'a' + v as u8 - 0x1a,
                _ => return Err(Error::Invalid),
            }
        } else {
            r.advance(3).ok_or(Error::Invalid)?;
            match r.peek(5) {
                v @ 0x00..=0x09 => b'!' + v as u8 - 8,
                v @ 0x0a..=0x14 => b'%' + v as u8 - 0x0a,
                v @ 0x15..=0x1a => b':' + v as u8 - 0x15,
                0x1b => b'_',
                0x1c => b' ',
                _ => return Err(Error::Invalid),
            }
        };

        if ch != 0 {
            out.push(ch)?;
        }
    }

    Ok(())
}

fn postprocess_exp(dcode: &mut ImageScanner, data: &[i32]) -> Result<usize> {
    let first = *data.first().ok_or(Error::Invalid)? as u64;
    let mut len = (first / 211 + 4) as i32;

    // Encodation method lives in the second character; the rest of the stream
    // is fed in on demand.
    let d = *data.get(1).ok_or(Error::Invalid)? as u64;
    let mut r = BitReader::new(&data[2..], d, 12, len);

    let n = (d >> 4) & 0x7f;
    let (enc, buflen) = if n >= 0x40 {
        r.bits = 10;
        (1, 2 + 14 + var_max(len, 10 - 2 - 44 + 6) + 2)
    } else if n >= 0x38 {
        r.bits = 4;
        (6 + (n as i32 & 7), 2 + 14 + 4 + 6 + 2 + 6 + 2)
    } else if n >= 0x30 {
        r.bits = 6;
        (2 + ((n >> 2) & 1) as i32, {
            2 + 14 + 4 + 3 + var_max(len, 6 - 2 - 44 - 2 - 10) + 2
        })
    } else if n >= 0x20 {
        r.bits = 7;
        (4 + ((n >> 3) & 1) as i32, 2 + 14 + 4 + 6)
    } else {
        r.bits = 9;
        (0, var_max(len, 9 - 2) + 2)
    };

    if buflen <= 2 {
        return Err(Error::Invalid);
    }

    if enc < 4 {
        // variable length symbol bit field
        r.advance(1).ok_or(Error::Invalid)?;
        if ((len ^ r.peek(1) as i32) & 1) != 0 {
            // even/odd length mismatch
            return Err(Error::Invalid);
        }

        r.advance(1).ok_or(Error::Invalid)?;
        if r.peek(1) as i32 != (len > 14) as i32 {
            // size group mismatch
            return Err(Error::Invalid);
        }
    }

    len -= 2;
    r.remaining = len;

    let buflen = buflen as usize;
    if dcode.set_buffer_capacity(buflen).is_err() {
        return Err(Error::Invalid);
    }

    // C's `size_buf` never shrinks below BUFFER_MIN, so the writable region is
    // at least that large even when `buflen` is smaller.
    let alloc = buflen.max(BUFFER_MIN);
    let buf = dcode.buffer_mut_slice(alloc).map_err(|()| Error::Invalid)?;
    let mut out = OutBuf { buf, len: 0 };

    // handle compressed fields
    if enc != 0 {
        out.push_all(b"01")?;
    }

    if enc == 1 {
        r.advance(4).ok_or(Error::Invalid)?;
        if r.bits >= 10 {
            return Err(Error::Invalid);
        }
        out.push(b'0' + r.peek(4) as u8)?;
    } else if enc != 0 {
        out.push(b'9')?;
    }

    if enc != 0 {
        for _ in 0..4 {
            let n = r.read(10).ok_or(Error::Invalid)?;
            if n >= 1000 {
                return Err(Error::Invalid);
            }
            out.push_decimal(n, 3)?;
        }
        // The 12 digits just written, plus the AI's leading digit, are the
        // first 13 of the GTIN-14; the check digit lands right after them.
        let start = out.len.checked_sub(13).ok_or(Error::Invalid)?;
        let end = out.len + 1;
        append_check14(out.buf.get_mut(start..end).ok_or(Error::Invalid)?);
        out.len = end;
    }

    match enc {
        // 01100: AI 392x
        2 => {
            let val = r.read(2).ok_or(Error::Invalid)?;
            out.push_all(&[b'3', b'9', b'2', b'0' + val as u8])?;
        }
        // 01101: AI 393x
        3 => {
            r.feed(12);
            r.advance(2).ok_or(Error::Invalid)?;
            let val = r.peek(2);
            out.push_all(&[b'3', b'9', b'3', b'0' + val as u8])?;
            r.advance(10).ok_or(Error::Invalid)?;
            let n = r.peek(10);
            if n >= 1000 {
                return Err(Error::Invalid);
            }
            out.push_decimal(n, 3)?;
        }
        // 0100: AI 3103
        4 => {
            let n = r.read(15).ok_or(Error::Invalid)?;
            out.push_all(b"3103")?;
            out.push_decimal(n, 6)?;
        }
        // 0101: AI 3202/3203
        5 => {
            let mut n = r.read(15).ok_or(Error::Invalid)?;
            out.push_all(&[b'3', b'2', b'0', if n >= 10000 { b'3' } else { b'2' }])?;
            if n >= 10000 {
                n -= 10000;
            }
            out.push_decimal(n, 6)?;
        }
        _ => {}
    }

    // 0111000 - 0111111: AI 310x/320x + AI 11/13/15/17
    if enc >= 6 {
        out.push_all(&[b'3', b'1' + (enc & 1) as u8, b'0', b'x'])?;

        let n = r.read(20).ok_or(Error::Invalid)?;
        if n >= 1_000_000 {
            return Err(Error::Invalid);
        }
        // The decimal point position replaces the placeholder 'x', shifting
        // the first digit of the value back one slot.
        let start = out.len;
        out.push_decimal(n, 6)?;
        out.buf[start - 1] = out.buf[start];
        out.buf[start] = b'0';

        let n = r.read(16).ok_or(Error::Invalid)?;
        if n < 38400 {
            let dd = n % 32;
            let rem = n / 32;
            let mm = rem % 12 + 1;
            let yy = rem / 12;

            out.push(b'1')?;
            out.push(b'0' + ((enc - 6) | 1) as u8)?;
            out.push_decimal(yy, 2)?;
            out.push_decimal(mm, 2)?;
            out.push_decimal(dd, 2)?;
        } else if n > 38400 {
            return Err(Error::Invalid);
        }
    }

    if enc < 4 {
        // remainder is general-purpose data compaction
        postprocess_exp_tail(&mut r, &mut out)?;
    }

    // C NUL-terminates and asserts the terminator is still in bounds.
    let total_len = out.len;
    out.push(0)?;

    // Trailing separators are not part of the data.
    let final_len = if total_len > 0 && out.buf[total_len - 1] == GS {
        out.buf[total_len - 1] = 0;
        total_len - 1
    } else {
        total_len
    };

    dcode.set_buffer_len(final_len);
    Ok(final_len)
}

/// Convert DataBar data from heterogeneous base {1597,2841} to base 10 character representation
fn postprocess(dcode: &mut ImageScanner, mut d: [u32; 4]) {
    // Get config before borrowing buffer
    let emit_check = dcode.should_emit_checksum(SymbolType::Databar);

    // Layout: "01" at [0..2], 13 digits at [2..15] written back to front,
    // then the check digit at [15] and a NUL at [16].
    let buf = match dcode.buffer_mut_slice(17) {
        Ok(buf) => buf,
        Err(_) => return,
    };
    let mut chk = 0u32;

    // Write "01" prefix
    buf[0] = b'0';
    buf[1] = b'1';

    // Write two null terminators, then work backwards from just below them
    let mut buf_idx = 17;
    buf_idx -= 1;
    buf[buf_idx] = 0;
    buf_idx -= 1;
    buf[buf_idx] = 0;

    // First conversion
    let mut r = (d[0] as u64) * 1597 + (d[1] as u64);
    d[1] = (r / 10000) as u32;
    r %= 10000;
    r = r * 2841 + (d[2] as u64);
    d[2] = (r / 10000) as u32;
    r %= 10000;
    r = r * 1597 + (d[3] as u64);
    d[3] = (r / 10000) as u32;

    // Extract 4 decimal digits
    for i in (0..4).rev() {
        let c = (r % 10) as u32;
        chk += c;
        if (i & 1) != 0 {
            chk += c << 1;
        }
        buf_idx -= 1;
        buf[buf_idx] = c as u8 + b'0';
        if i != 0 {
            r /= 10;
        }
    }

    // Second conversion
    r = (d[1] as u64) * 2841 + (d[2] as u64);
    d[2] = (r / 10000) as u32;
    r %= 10000;
    r = r * 1597 + (d[3] as u64);
    d[3] = (r / 10000) as u32;

    // Extract 4 more decimal digits
    for i in (0..4).rev() {
        let c = (r % 10) as u32;
        chk += c;
        if (i & 1) != 0 {
            chk += c << 1;
        }
        buf_idx -= 1;
        buf[buf_idx] = c as u8 + b'0';
        if i != 0 {
            r /= 10;
        }
    }

    // Third conversion
    r = (d[2] as u64) * 1597 + (d[3] as u64);

    // Extract 5 decimal digits
    for i in (0..5).rev() {
        let c = (r % 10) as u32;
        chk += c;
        if (i & 1) == 0 {
            chk += c << 1;
        }
        buf_idx -= 1;
        buf[buf_idx] = c as u8 + b'0';
        if i != 0 {
            r /= 10;
        }
    }

    debug_assert_eq!(buf_idx, 2, "13 digits should land in buf[2..15]");

    // Add check digit if configured
    if emit_check {
        chk %= 10;
        if chk != 0 {
            chk = 10 - chk;
        }
        buf[buf_idx + 13] = chk as u8 + b'0';
        dcode.set_buffer_len(buf_idx + 14);
    } else {
        dcode.set_buffer_len(buf_idx + 13);
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
fn check_width(wf: u32, wd: u32, n: u32) -> i32 {
    let dwf = wf.wrapping_mul(3);
    let wd = wd.wrapping_mul(14);
    let wf = wf.wrapping_mul(n);

    // Check: wf - dwf <= wd && wd <= wf + dwf
    // In C, this relies on unsigned wraparound if wf < dwf
    // For unsigned subtraction: wf - dwf will wrap if wf < dwf,
    // resulting in a very large number, making the condition false
    if wf.wrapping_sub(dwf) <= wd && wd <= wf.wrapping_add(dwf) {
        1
    } else {
        0
    }
}

/// Merge or update a DataBar segment with existing segments
fn merge_segment(db: &mut DatabarDecoder, seg_idx: usize) {
    let csegs = db.csegs();
    let epoch = db.epoch();

    for i in 0..csegs {
        // Skip if this is the same segment
        if i == seg_idx {
            continue;
        }

        // Read values we need from the target segment
        let seg_finder = db.seg(seg_idx).finder();
        let seg_exp = db.seg(seg_idx).exp();
        let seg_color = db.seg(seg_idx).color();
        let seg_side = db.seg(seg_idx).side();
        let seg_data = db.seg(seg_idx).data;
        let seg_check = db.seg(seg_idx).check();
        let seg_width = db.seg(seg_idx).width;
        let seg_partial = db.seg(seg_idx).partial();

        let s = db.seg_mut(i);
        // Check if this segment matches and should be merged
        if s.finder() == seg_finder
            && s.exp() == seg_exp
            && s.color() == seg_color
            && s.side() == seg_side
            && s.data == seg_data
            && s.check() == seg_check
            && check_width(seg_width as u32, s.width as u32, 14) != 0
        {
            // Found a matching segment - merge with it
            let mut cnt = s.count();
            if cnt < 0x7F {
                cnt += 1;
            }

            // Merge partial flags (bitwise AND)
            let new_partial = seg_partial && s.partial();

            // Average the widths (weighted average favoring new measurement).
            // C promotes both `unsigned short`s to `int`, so the intermediate
            // `3 * width` must not be computed in 16 bits.
            let new_width =
                ((3 * seg_width as u32 + s.width as u32 + 2) / 4).min(u16::MAX as u32) as u16;

            // Mark old segment as unused
            s.set_finder(-1);

            // Now update the target segment
            let seg = &mut db.seg_mut(seg_idx);
            seg.set_count(cnt);
            seg.set_partial(new_partial);
            seg.width = new_width;
        } else if s.finder() >= 0 {
            // Not a match, check if we should age it out
            let age = epoch.wrapping_sub(s.epoch());
            if age >= 248 || (age >= 128 && s.count() < 2) {
                s.set_finder(-1);
            }
        }
    }
}

/// Match DataBar segment to find complete symbol
fn match_segment(dcode: &mut ImageScanner, seg_idx: usize) -> SymbolType {
    let db = &mut dcode.databar;
    let csegs = db.csegs();
    let mut maxage = 0xfff;
    let mut maxcnt = 0;
    let mut smax: Option<[usize; 3]> = None;
    let mut d = [0u32; 4];

    let seg = db.seg(seg_idx);
    if seg.partial() && seg.count() < 4 {
        return SymbolType::Partial;
    }

    // Cache values we need from seg
    let seg_finder = seg.finder();
    let seg_color = seg.color();
    let seg_side = seg.side();
    let seg_width = seg.width;
    let seg_check = seg.check();

    for i0 in 0..csegs {
        if i0 == seg_idx {
            continue;
        }
        let s0 = db.seg(i0);
        if s0.finder() != seg_finder
            || s0.exp()
            || s0.color() != seg_color
            || s0.side() == seg_side
            || (s0.partial() && s0.count() < 4)
            || check_width(seg_width as u32, s0.width as u32, 14) == 0
        {
            continue;
        }

        for i1 in 0..csegs {
            if i1 == i0 {
                continue;
            }
            let s1 = db.seg(i1);
            if s1.finder() < 0
                || s1.exp()
                || s1.color() == seg_color
                || (s1.partial() && s1.count() < 4)
                || check_width(seg_width as u32, s1.width as u32, 14) == 0
            {
                continue;
            }

            let mut chkf = if seg_color != Color::Space {
                seg_finder as i32 + s1.finder() as i32 * 9
            } else {
                s1.finder() as i32 + seg_finder as i32 * 9
            };
            if chkf > 72 {
                chkf -= 1;
            }
            if chkf > 8 {
                chkf -= 1;
            }

            let chks = ((seg_check as i32) + (s0.check() as i32) + (s1.check() as i32)) % 79;

            let chk = if chkf >= chks {
                chkf - chks
            } else {
                79 + chkf - chks
            };

            let age1 = ((db.epoch().wrapping_sub(s0.epoch())) as u32)
                + ((db.epoch().wrapping_sub(s1.epoch())) as u32);

            for i2 in (i1 + 1)..csegs {
                if i2 == i0 {
                    continue;
                }
                let s2 = db.seg(i2);
                if s2.finder() != s1.finder()
                    || s2.exp()
                    || s2.color() != s1.color()
                    || s2.side() == s1.side()
                    || s2.check() as i32 != chk
                    || (s2.partial() && s2.count() < 4)
                    || check_width(seg_width as u32, s2.width as u32, 14) == 0
                {
                    continue;
                }
                let age2 = db.epoch().wrapping_sub(s2.epoch()) as u32;
                let age = age1 + age2;
                let cnt = s0.count() as u32 + s1.count() as u32 + s2.count() as u32;
                if maxcnt < cnt as i32 || (maxcnt == cnt as i32 && (maxage as i32) > (age as i32)) {
                    maxcnt = cnt as i32;
                    maxage = age;
                    smax = Some([i0, i1, i2]);
                }
            }
        }
    }

    let Some(smax) = smax else {
        return SymbolType::Partial;
    };

    let seg = db.seg(seg_idx);
    d[((seg.color() as usize) << 1) | (seg.side() as usize)] = seg.data as u32;

    for idx in smax {
        let s = db.seg(idx);
        d[((s.color() as usize) << 1) | (s.side() as usize)] = s.data as u32;
        let new_count = s.count().wrapping_sub(1);
        let s = db.seg_mut(idx);
        s.set_count(new_count);
        if new_count == 0 {
            s.set_finder(-1);
        }
    }
    db.seg_mut(seg_idx).set_finder(-1);

    // Extract values before dropping db borrow
    let seg_side = db.seg(seg_idx).side();
    let seg_color = db.seg(seg_idx).color();

    if dcode.set_buffer_capacity(18).is_err() {
        return SymbolType::Partial;
    }

    if !dcode.acquire_lock(SymbolType::Databar) {
        return SymbolType::Partial;
    }

    postprocess(dcode, d);
    dcode.modifiers = Modifier::Gs1.bit();
    dcode.direction = 1 - 2 * ((seg_side as i32) ^ (seg_color as i32) ^ 1);
    SymbolType::Databar
}

/// Lookup DataBar expanded sequence
/// Returns -1 on error, 0 or 1 on success
fn lookup_sequence(seg: &mut DatabarSegment, fixed: i32, seq: &mut [i32], maxsize: usize) -> i32 {
    let mut n = (seg.data as u32 / 211) as usize;
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
    seq[0] = 0;
    seq[1] = 1;
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
        seq[i] = s;
        i += 1;
        s += 1;
        seq[i] = s;
        i += 1;
    }
    seq[n] = -1;
    if fixed < 1 { 1 } else { 0 }
}

fn match_segment_exp(dcode: &mut ImageScanner, seg_idx: usize, dir: i32) -> SymbolType {
    let db = &mut dcode.databar;
    let csegs = db.csegs();
    if csegs == 0 {
        return SymbolType::Partial;
    }

    let ifixed = seg_idx;
    let fixed = db.seg(seg_idx).segment_index();
    let mut bestsegs = [-1i32; 22];
    let mut segs_idx = [-1i32; 22];
    let mut seq = [-1i32; 22];
    let mut iseg = [-1i32; DATABAR_MAX_SEGMENTS];
    let mut width_stack = [0u32; 22];

    seq[0] = 0;
    seq[1] = -1;
    segs_idx[0] = -1;
    bestsegs[0] = -1;
    width_stack[0] = db.seg(seg_idx).width as u32;

    #[allow(clippy::needless_range_loop)]
    for j in 0..csegs {
        let s = db.seg(j);
        iseg[j] = if s.exp() && s.finder() >= 0 && (!s.partial() || s.count() as i32 >= 4) {
            s.segment_index()
        } else {
            -1
        };
    }

    let mut maxcnt = 0i32;
    let mut maxage = 0x7fff_u32;
    let mut i: i32 = 0;

    loop {
        while i >= 0 && seq[i as usize] >= 0 {
            let idx = i as usize;
            let target = seq[idx];
            let mut found: Option<usize> = None;
            let current_width = width_stack[idx];

            if target == fixed {
                let candidate = db.seg(ifixed);
                if segs_idx[idx] < 0 && check_width(current_width, candidate.width as u32, 14) != 0
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
                        let cand = db.seg(start);
                        if idx == 0 || check_width(current_width, cand.width as u32, 14) != 0 {
                            found = Some(start);
                            break;
                        }
                    }
                    start += 1;
                }
            }

            if let Some(jidx) = found {
                if idx == 0 {
                    let maxsize = seq.len();
                    let lu = lookup_sequence(db.seg_mut(jidx), fixed, &mut seq, maxsize);
                    if lu == 0 {
                        i -= 1;
                        continue;
                    }
                    if lu < 0 {
                        return SymbolType::None;
                    }
                }

                let candidate = db.seg(jidx);
                let next_width = if idx == 0 {
                    candidate.width as u32
                } else {
                    (current_width + candidate.width as u32) / 2
                };
                width_stack[idx + 1] = next_width;
                segs_idx[idx] = jidx as i32;
                segs_idx[idx + 1] = -1;
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
        let seg_eval = db.seg(first_idx);
        age = (db.epoch().wrapping_sub(seg_eval.epoch())) as u32 & 0xff;

        let mut pos = 1usize;
        while pos < segs_idx.len() && segs_idx[pos] >= 0 {
            let sidx = segs_idx[pos] as usize;
            let seg_eval = db.seg(sidx);
            chk += seg_eval.check() as u32;
            cnt += seg_eval.count() as u32;
            age += (db.epoch().wrapping_sub(seg_eval.epoch())) as u32 & 0xff;
            pos += 1;
        }

        let mut chk0 = db.seg(first_idx).data as i32 % 211;
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
        return SymbolType::Partial;
    }

    // Extract data values before dropping db borrow
    let mut data_vals = [0i32; 22];
    let mut count = 0usize;
    while count < bestsegs.len() && bestsegs[count] >= 0 {
        let sidx = bestsegs[count] as usize;
        let s = db.seg(sidx);
        data_vals[count] = s.data as i32;
        count += 1;
    }

    // Update segments before dropping db borrow
    let mut last_idx = ifixed;
    for sidx in bestsegs {
        if sidx < 0 {
            break;
        }
        let idx = sidx as usize;
        last_idx = idx;
        if idx != ifixed {
            let s = db.seg_mut(idx);
            let mut cnt = s.count();
            if cnt > 0 {
                cnt -= 1;
                s.set_count(cnt);
                if cnt == 0 {
                    s.set_finder(-1);
                }
            }
        }
    }

    let seg_side = db.seg(last_idx).side();
    let seg_color = db.seg(last_idx).color();

    if !dcode.acquire_lock(SymbolType::DatabarExp) {
        return SymbolType::Partial;
    }

    if postprocess_exp(dcode, &data_vals).is_err() {
        dcode.release_lock(SymbolType::DatabarExp);
        return SymbolType::Partial;
    }

    dcode.direction = (1 - 2 * ((seg_side ^ seg_color as u8) as i32)) * dir;
    dcode.modifiers = Modifier::Gs1.bit();
    SymbolType::DatabarExp
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
fn calc_check(mut sig0: u32, mut sig1: u32, side: u32, mod_val: u32) -> u32 {
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
fn calc_value4(sig: u32, mut n: u32, wmax: u32, mut nonarrow: u32) -> i32 {
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

    v as i32
}

/// Decode a DataBar character from width measurements
fn decode_char(dcode: &mut ImageScanner, seg_idx: usize, off: i32, dir: i32) -> SymbolType {
    // Read segment values we need before taking other borrows
    let seg_exp = dcode.databar.seg(seg_idx).exp();
    let seg_side = dcode.databar.seg(seg_idx).side();
    let seg_width = dcode.databar.seg(seg_idx).width;
    let seg_finder = dcode.databar.seg(seg_idx).finder();
    let seg_color = dcode.databar.seg(seg_idx).color();

    let s = dcode.calc_s(if dir > 0 { off } else { off - 6 } as u8, 8);
    let mut emin = [0i32, 0i32];
    let mut sum = 0i32;
    let mut sig0 = 0u32;
    let mut sig1 = 0u32;

    let n = if seg_exp {
        17
    } else if seg_side != 0 {
        15
    } else {
        16
    };
    emin[1] = -(n as i32);

    if s < 13 || check_width(seg_width as u32, s, n) == 0 {
        return SymbolType::None;
    }

    let mut off = off;
    for i in (0..4).rev() {
        let e = decode_e(dcode.pair_width(off as u8), s, n);
        if e < 0 {
            return SymbolType::None;
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

        let e = decode_e(dcode.pair_width(off as u8), s, n);
        if e < 0 {
            return SymbolType::None;
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

    if sum0.wrapping_add(sum1).wrapping_add(8) as i32 != n as i32 {
        return SymbolType::None;
    }

    if ((sum0 ^ (n >> 1)) | (sum1 ^ (n >> 1) ^ n)) & 1 != 0 {
        return SymbolType::None;
    }

    let i = ((n & 0x3) ^ 1) * 5 + (sum1 >> 1);
    if i as usize >= GROUPS.len() {
        return SymbolType::None;
    }
    let g = &GROUPS[i as usize];

    let vodd = calc_value4(
        sig0.wrapping_add(0x1111),
        sum0.wrapping_add(4),
        g.wmax as u32,
        (!(n as i32) & 1) as u32,
    );
    if vodd < 0 || vodd > g.todd as i32 {
        return SymbolType::None;
    }

    let veven = calc_value4(
        sig1.wrapping_add(0x1111),
        sum1.wrapping_add(4),
        (9 - g.wmax) as u32,
        n & 1,
    );
    if veven < 0 || veven > g.teven as i32 {
        return SymbolType::None;
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
    if seg_exp {
        let side = seg_color as u8 ^ seg_side ^ 1;
        if v >= 4096 {
            return SymbolType::None;
        }
        chk = calc_check(sig0, sig1, side as u32, 211);
        if seg_finder != 0 || seg_color != Color::Space || seg_side != 0 {
            let i = (seg_finder as i32) * 2 - (side as i32) + (seg_color as i32);
            if !(0..12).contains(&i) {
                return SymbolType::None;
            }
            chk = (chk * EXP_CHECKSUMS[i as usize] as u32) % 211;
        } else if v >= 4009 {
            return SymbolType::None;
        } else {
            chk = 0;
        }
    } else {
        chk = calc_check(sig0, sig1, seg_side as u32, 79);
        if seg_color != Color::Space {
            chk = (chk * 16) % 79;
        }
    }

    let seg = dcode.databar.seg_mut(seg_idx);
    seg.set_check(chk as u8);
    seg.data = v as i16;

    merge_segment(&mut dcode.databar, seg_idx);

    let is_exp = dcode.databar.seg(seg_idx).exp();
    if is_exp {
        return match_segment_exp(dcode, seg_idx, dir);
    } else if dir > 0 {
        return match_segment(dcode, seg_idx);
    }
    SymbolType::Partial
}

/// Allocate a new DataBar segment (or reuse an old one)
/// Returns the index of the allocated segment, or -1 on failure
fn alloc_segment(db: &mut DatabarDecoder) -> i32 {
    let mut maxage = 0u32;
    let csegs = db.csegs();
    let epoch = db.epoch();
    let mut old: i32 = -1;

    // First pass: look for empty slots or very old segments
    for i in 0..csegs {
        let seg = db.seg_mut(i);

        if seg.finder() < 0 {
            return i as i32;
        }

        let age = epoch.wrapping_sub(seg.epoch());
        if age >= 128 && seg.count() < 2 {
            seg.set_finder(-1);
            return i as i32;
        }

        // Score based on both age and count
        let score = if age > seg.count() {
            age - seg.count() + 1
        } else {
            1
        };

        if maxage < score as u32 {
            maxage = score as u32;
            old = i as i32;
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
            // Grow the segment array
            db.resize_segs(new_csegs);
            return i as i32;
        }
    }

    // Reuse oldest segment
    if old >= 0 {
        db.seg_mut(old as usize).set_finder(-1);
    }
    old
}

/// Decode DataBar finder pattern
fn decode_finder(dcode: &mut ImageScanner) -> SymbolType {
    let e0 = dcode.pair_width(1);
    let e2 = dcode.pair_width(3);
    let (dir, e2, e3) = if e0 < e2 {
        let e = e2 * 4;
        if e < 15 * e0 || e > 34 * e0 {
            return SymbolType::None;
        }
        (0, e2, dcode.pair_width(4))
    } else {
        let e = e0 * 4;
        if e < 15 * e2 || e > 34 * e2 {
            return SymbolType::None;
        }
        (1, e0, dcode.pair_width(0))
    };
    let e1 = dcode.pair_width(2);

    let s = e1 + e3;
    if s < 12 {
        return SymbolType::None;
    }

    let sig = (decode_e(e3, s, 14) << 8) | (decode_e(e2, s, 14) << 4) | decode_e(e1, s, 14);
    if sig < 0
        || ((sig >> 4) & 0xf) < 8
        || ((sig >> 4) & 0xf) > 10
        || (sig & 0xf) >= 10
        || ((sig >> 8) & 0xf) >= 10
        || (((sig >> 8) + sig) & 0xf) != 10
    {
        return SymbolType::None;
    }

    let finder = (FINDER_HASH[((sig - (sig >> 5)) & 0x1f) as usize]
        + FINDER_HASH[((sig >> 1) & 0x1f) as usize])
        & 0x1f;
    if finder == 0x1f
        || !dcode.is_enabled(if finder < 9 {
            SymbolType::Databar
        } else {
            SymbolType::DatabarExp
        })
    {
        return SymbolType::None;
    }

    if finder < 0 {
        return SymbolType::None;
    }

    let iseg = alloc_segment(&mut dcode.databar);
    if iseg < 0 {
        return SymbolType::None;
    }

    let color = dcode.color();
    let epoch = dcode.databar.epoch();
    let seg = dcode.databar.seg_mut(iseg as usize);
    seg.set_finder(if finder >= 9 { finder - 9 } else { finder });
    seg.set_exp(finder >= 9);
    seg.set_color(((color as u8) ^ dir as u8 ^ 1).into());
    seg.set_side(dir as u8);
    seg.set_partial(false);
    seg.set_count(1);
    seg.width = s as u16;
    seg.set_epoch(epoch);

    let rc = decode_char(dcode, iseg as usize, 12 - dir, -1);
    if rc == SymbolType::None {
        let seg = dcode.databar.seg_mut(iseg as usize);
        seg.set_partial(true);
    } else {
        dcode
            .databar
            .set_epoch(dcode.databar.epoch().wrapping_add(1));
    }

    let i = ((dcode.idx as i32 + 8 + dir) & 0xf) as usize;
    if dcode.databar.char(i) != -1 {
        return SymbolType::None;
    }
    dcode.databar.set_char(i, iseg as i8);
    rc
}

pub(crate) fn decode_databar(dcode: &mut ImageScanner) -> SymbolType {
    let i = dcode.idx & 0xf;

    let mut sym = decode_finder(dcode);

    let iseg = dcode.databar.char(i as usize);
    if iseg < 0 {
        return sym;
    }

    dcode.databar.set_char(i as usize, -1);

    let seg_idx: usize;
    let pair_idx: Option<usize>;

    if dcode.databar.seg(iseg as usize).partial() {
        pair_idx = None;
        seg_idx = iseg as usize;
        let old_side = dcode.databar.seg(seg_idx).side();
        dcode.databar.seg_mut(seg_idx).set_side(!old_side);
    } else {
        let jseg = alloc_segment(&mut dcode.databar);
        pair_idx = Some(iseg as usize);
        seg_idx = jseg as usize;

        // Extract values from pair before borrowing seg mutably
        let pair_finder = dcode.databar.seg(iseg as usize).finder();
        let pair_exp = dcode.databar.seg(iseg as usize).exp();
        let pair_color = dcode.databar.seg(iseg as usize).color();
        let pair_side = dcode.databar.seg(iseg as usize).side();
        let pair_width = dcode.databar.seg(iseg as usize).width;
        let epoch = dcode.databar.epoch();

        let seg = dcode.databar.seg_mut(jseg as usize);
        seg.set_finder(pair_finder);
        seg.set_exp(pair_exp);
        seg.set_color(pair_color);
        seg.set_side(!pair_side);
        seg.set_partial(false);
        seg.set_count(1);
        seg.width = pair_width;
        seg.set_epoch(epoch);
    }

    sym = decode_char(dcode, seg_idx, 1, 1);
    if sym == SymbolType::None {
        dcode.databar.seg_mut(seg_idx).set_finder(-1);
        if let Some(pidx) = pair_idx {
            dcode.databar.seg_mut(pidx).set_partial(true);
        }
    } else {
        dcode
            .databar
            .set_epoch(dcode.databar.epoch().wrapping_add(1));
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
        assert_eq!(check_width(100, 100, 14), 1);

        // Test within tolerance: wf=100, wd=105, n=14
        // After scaling: dwf=300, wd=1470, wf=1400
        // Check: (1400-300 <= 1470) && (1470 <= 1400+300) => (1100<=1470) && (1470<=1700) => true
        assert_eq!(check_width(100, 105, 14), 1);

        // Test within tolerance: wf=100, wd=95, n=14
        // After scaling: dwf=300, wd=1330, wf=1400
        // Check: (1400-300 <= 1330) && (1330 <= 1400+300) => (1100<=1330) && (1330<=1700) => true
        assert_eq!(check_width(100, 95, 14), 1);

        // Test outside tolerance: wf=100, wd=130, n=14
        // After scaling: dwf=300, wd=1820, wf=1400
        // Check: (1400-300 <= 1820) && (1820 <= 1400+300) => (1100<=1820) && (1820<=1700) => false
        assert_eq!(check_width(100, 130, 14), 0);

        // Test outside tolerance: wf=100, wd=70, n=14
        // After scaling: dwf=300, wd=980, wf=1400
        // Check: (1400-300 <= 980) && (980 <= 1400+300) => (1100<=980) && (980<=1700) => false
        assert_eq!(check_width(100, 70, 14), 0);
    }

    #[test]
    fn test_decode10() {
        let mut buf = [0u8; 10];

        // Test simple number
        decode10(&mut buf, 123, 3);
        assert_eq!(&buf[0..3], b"123");

        // Test with leading zeros
        decode10(&mut buf, 45, 6);
        assert_eq!(&buf[0..6], b"000045");

        // Test zero
        decode10(&mut buf, 0, 3);
        assert_eq!(&buf[0..3], b"000");
    }

    #[test]
    fn test_append_check14() {
        // Test with 13 ASCII digits
        let mut buf = [
            b'9', b'7', b'8', b'0', b'1', b'4', b'3', b'0', b'0', b'7', b'2', b'3', b'0', 0,
        ];
        append_check14(&mut buf);
        // Check digit: (9+7+8+0+1+4+3+0+0+7+2+3) + (9+8+1+3+0+2)*2 = 44 + 46 = 90, 90%10=0, check=0
        assert_eq!(buf[13], b'0');

        // Test another example: 1234567890120
        let mut buf2 = [
            b'1', b'2', b'3', b'4', b'5', b'6', b'7', b'8', b'9', b'0', b'1', b'2', b'0', 0,
        ];
        append_check14(&mut buf2);
        // Check: (1+2+3+4+5+6+7+8+9+0+1+2) + (1+3+5+7+9+1)*2 = 48 + 52 = 100, 100%10=0, check=0
        assert_eq!(buf2[13], b'0');
    }

    #[test]
    fn test_calc_check() {
        // Test basic checksum calculation
        // These values are based on understanding the algorithm
        let chk1 = calc_check(0x1234, 0x5678, 0, 211);
        assert!(chk1 < 211);

        let chk2 = calc_check(0x1234, 0x5678, 1, 211);
        assert!(chk2 < 211);

        // Side should affect the result
        assert_ne!(chk1, chk2);

        // Test with different modulus
        let chk3 = calc_check(0x1234, 0x5678, 0, 79);
        assert!(chk3 < 79);

        // Same inputs should give same output
        let chk4 = calc_check(0x1234, 0x5678, 0, 211);
        assert_eq!(chk1, chk4);
    }
}
