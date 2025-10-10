//! Low-level barcode decoder

use crate::{
    databar_utils::{_zbar_databar_calc_check, _zbar_databar_check_width},
    decoder_types::{
        codabar_decoder_t, code128_decoder_t, code39_decoder_t, code93_decoder_t,
        databar_decoder_t, ean_decoder_t, i25_decoder_t, qr_finder_t, zbar_decoder_t,
        zbar_symbol_type_t, BUFFER_INCR, BUFFER_MAX, BUFFER_MIN, DECODE_WINDOW,
        ZBAR_CFG_EMIT_CHECK, ZBAR_CFG_ENABLE, ZBAR_CFG_MAX_LEN, ZBAR_CFG_MIN_LEN, ZBAR_CODABAR,
        ZBAR_CODE128, ZBAR_CODE39, ZBAR_CODE93, ZBAR_COMPOSITE, ZBAR_DATABAR, ZBAR_DATABAR_EXP,
        ZBAR_EAN13, ZBAR_EAN2, ZBAR_EAN5, ZBAR_EAN8, ZBAR_I25, ZBAR_ISBN10, ZBAR_ISBN13,
        ZBAR_MOD_GS1, ZBAR_NONE, ZBAR_PARTIAL, ZBAR_QRCODE, ZBAR_SQCODE, ZBAR_UPCA, ZBAR_UPCE,
    },
    decoders::{
        codabar::_zbar_decode_codabar, code128::_zbar_decode_code128, code39::_zbar_decode_code39,
        code93::_zbar_decode_code93, ean::_zbar_decode_ean, i25::_zbar_decode_i25,
    },
};
use libc::{c_char, c_int, c_uint, c_void};

// Config constant not in decoder_types
const ZBAR_CFG_NUM: c_int = 5;

// DataBar constants
const DATABAR_MAX_SEGMENTS: usize = 32;

#[repr(C)]
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

static EXP_CHECKSUMS: [u8; 12] = [1, 189, 62, 113, 46, 43, 109, 134, 6, 79, 161, 45];

/// DataBar expanded sequences
static EXP_SEQUENCES: [u8; 30] = [
    // sequence Group 1
    0x01, 0x23, 0x25, 0x07, 0x29, 0x47, 0x29, 0x67, 0x0b, 0x29, 0x87, 0xab,
    // sequence Group 2
    0x21, 0x43, 0x65, 0x07, 0x21, 0x43, 0x65, 0x89, 0x21, 0x43, 0x65, 0xa9, 0x0b, 0x21, 0x43, 0x67,
    0x89, 0xab,
];

/// DataBar finder pattern hash table
static FINDER_HASH: [i8; 0x20] = [
    0x16, 0x1f, 0x02, 0x00, 0x03, 0x00, 0x06, 0x0b, 0x1f, 0x0e, 0x17, 0x0c, 0x0b, 0x14, 0x11, 0x0c,
    0x1f, 0x03, 0x13, 0x08, 0x00, 0x0a, -1, 0x16, 0x0c, 0x09, -1, 0x1a, 0x1f, 0x1c, 0x00, -1,
];

/// Calculate DataBar character value from 4-element signature
/// Returns -1 on error
#[no_mangle]
pub unsafe extern "C" fn calc_value4(
    sig: c_uint,
    mut n: c_uint,
    wmax: c_uint,
    mut nonarrow: c_uint,
) -> c_int {
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

// Forward declarations for complex C functions not yet ported
extern "C" {
    fn match_segment_exp(
        dcode: *mut zbar_decoder_t,
        seg: *mut crate::decoder_types::databar_segment_t,
        dir: c_int,
    ) -> zbar_symbol_type_t;
}

/// Match DataBar segment to find complete symbol
unsafe fn match_segment(
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

            let mut chkf = if (*seg).color() != 0 {
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

    if _zbar_decoder_size_buf(dcode, 18) != 0 {
        return ZBAR_PARTIAL;
    }

    if _zbar_decoder_acquire_lock(dcode, ZBAR_DATABAR) != 0 {
        return ZBAR_PARTIAL;
    }

    _zbar_databar_postprocess(dcode, d.as_mut_ptr());
    (*dcode).modifiers = 1 << ZBAR_MOD_GS1;
    (*dcode).direction = 1 - 2 * (((*seg).side() as i32) ^ ((*seg).color() as i32) ^ 1);
    ZBAR_DATABAR
}

/// Lookup DataBar expanded sequence
/// Returns -1 on error, 0 or 1 on success
#[no_mangle]
pub unsafe extern "C" fn lookup_sequence(
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

/// Decode a DataBar character from width measurements
unsafe fn decode_char(
    dcode: *mut zbar_decoder_t,
    seg: *mut crate::decoder_types::databar_segment_t,
    off: c_int,
    dir: c_int,
) -> zbar_symbol_type_t {
    let db = &mut (*dcode).databar;
    let s = _zbar_decoder_calc_s(
        dcode as *const zbar_decoder_t,
        if dir > 0 { off } else { off - 6 } as u8,
        8,
    );
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
        let e = _zbar_decoder_decode_e(_zbar_decoder_pair_width(dcode, off as u8), s, n);
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

        let e = _zbar_decoder_decode_e(_zbar_decoder_pair_width(dcode, off as u8), s, n);
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
        let side = (*seg).color() ^ (*seg).side() ^ 1;
        if v >= 4096 {
            return ZBAR_NONE;
        }
        chk = _zbar_databar_calc_check(sig0, sig1, side as c_uint, 211);
        if (*seg).finder() != 0 || (*seg).color() != 0 || (*seg).side() != 0 {
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
        if (*seg).color() != 0 {
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

/// Decode DataBar finder pattern
#[no_mangle]
pub unsafe extern "C" fn decode_finder(dcode: *mut zbar_decoder_t) -> zbar_symbol_type_t {
    let db = &mut (*dcode).databar;
    let e0 = _zbar_decoder_pair_width(dcode, 1);
    let e2 = _zbar_decoder_pair_width(dcode, 3);
    let (dir, e2, e3) = if e0 < e2 {
        let e = e2 * 4;
        if e < 15 * e0 || e > 34 * e0 {
            return ZBAR_NONE;
        }
        (0, e2, _zbar_decoder_pair_width(dcode, 4))
    } else {
        let e = e0 * 4;
        if e < 15 * e2 || e > 34 * e2 {
            return ZBAR_NONE;
        }
        (1, e0, _zbar_decoder_pair_width(dcode, 0))
    };
    let e1 = _zbar_decoder_pair_width(dcode, 2);

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
    seg.set_color(((_zbar_decoder_get_color(dcode) as c_int) ^ dir ^ 1) as u8);
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

// External C functions for decoders not yet converted
extern "C" {
    fn _zbar_find_qr(dcode: *mut zbar_decoder_t) -> zbar_symbol_type_t;
    fn _zbar_decode_databar(dcode: *mut zbar_decoder_t) -> zbar_symbol_type_t;
}

// Macro equivalents
#[inline]
unsafe fn cfg_set(configs: &mut [c_int; 2], cfg: c_int, val: c_int) {
    configs[(cfg - ZBAR_CFG_MIN_LEN) as usize] = val;
}

#[inline]
fn test_cfg(config: c_uint, cfg: c_int) -> bool {
    ((config >> cfg) & 1) != 0
}

/// Get the color (bar/space) of the current element
#[no_mangle]
pub unsafe extern "C" fn _zbar_decoder_get_color(dcode: *const zbar_decoder_t) -> c_char {
    ((*dcode).idx & 1) as c_char
}

/// Get width of a specific element from the decoder's history window
#[no_mangle]
pub unsafe extern "C" fn _zbar_decoder_get_width(
    dcode: *const zbar_decoder_t,
    offset: u8,
) -> c_uint {
    (*dcode).w[(((*dcode).idx as usize).wrapping_sub(offset as usize)) & (DECODE_WINDOW - 1)]
}

/// Get the combined width of two consecutive elements
#[no_mangle]
pub unsafe extern "C" fn _zbar_decoder_pair_width(
    dcode: *const zbar_decoder_t,
    offset: u8,
) -> c_uint {
    _zbar_decoder_get_width(dcode, offset) + _zbar_decoder_get_width(dcode, offset + 1)
}

/// Calculate sum of n consecutive element widths
#[no_mangle]
pub unsafe extern "C" fn _zbar_decoder_calc_s(
    dcode: *const zbar_decoder_t,
    mut offset: u8,
    mut n: u8,
) -> c_uint {
    let mut s = 0;
    while n > 0 {
        s += _zbar_decoder_get_width(dcode, offset);
        offset += 1;
        n -= 1;
    }
    s
}

/// Decode element width into a discrete value
/// Returns -1 if the element width is invalid
#[no_mangle]
pub unsafe extern "C" fn _zbar_decoder_decode_e(e: c_uint, s: c_uint, n: c_uint) -> c_int {
    let big_e = ((e * n * 2 + 1) / s).wrapping_sub(3) / 2;
    if big_e >= n - 3 {
        -1
    } else {
        big_e as c_int
    }
}

/// Sort 3 elements by width and return their indices
#[no_mangle]
pub unsafe extern "C" fn _zbar_decoder_decode_sort3(
    dcode: *mut zbar_decoder_t,
    i0: c_int,
) -> c_uint {
    let w0 = _zbar_decoder_get_width(dcode, i0 as u8);
    let w2 = _zbar_decoder_get_width(dcode, (i0 + 2) as u8);
    let w4 = _zbar_decoder_get_width(dcode, (i0 + 4) as u8);

    if w0 < w2 {
        if w2 < w4 {
            ((i0 << 8) | ((i0 + 2) << 4) | (i0 + 4)) as c_uint
        } else if w0 < w4 {
            ((i0 << 8) | ((i0 + 4) << 4) | (i0 + 2)) as c_uint
        } else {
            (((i0 + 4) << 8) | (i0 << 4) | (i0 + 2)) as c_uint
        }
    } else if w4 < w2 {
        (((i0 + 4) << 8) | ((i0 + 2) << 4) | i0) as c_uint
    } else if w0 < w4 {
        (((i0 + 2) << 8) | (i0 << 4) | (i0 + 4)) as c_uint
    } else {
        (((i0 + 2) << 8) | ((i0 + 4) << 4) | i0) as c_uint
    }
}

/// Sort n elements by width and return their indices packed into a value
#[no_mangle]
pub unsafe extern "C" fn _zbar_decoder_decode_sortn(
    dcode: *mut zbar_decoder_t,
    n: c_int,
    i0: c_int,
) -> c_uint {
    let mut mask: c_uint = 0;
    let mut sort: c_uint = 0;

    for _ in (0..n).rev() {
        let mut wmin = c_uint::MAX;
        let mut jmin: c_int = -1;

        for j in (0..n).rev() {
            if (mask >> j) & 1 != 0 {
                continue;
            }
            let w = _zbar_decoder_get_width(dcode, (i0 + j * 2) as u8);
            if wmin >= w {
                wmin = w;
                jmin = j;
            }
        }

        debug_assert!(jmin >= 0, "sortn({}, {}) jmin={}", n, i0, jmin);
        sort <<= 4;
        mask |= 1 << jmin;
        sort |= (i0 + jmin * 2) as c_uint;
    }

    sort
}

/// Acquire a decoder lock for a specific symbology type
/// Returns 1 if already locked, 0 if lock acquired
#[no_mangle]
pub unsafe extern "C" fn _zbar_decoder_acquire_lock(
    dcode: *mut zbar_decoder_t,
    req: zbar_symbol_type_t,
) -> c_char {
    if (*dcode).lock != 0 {
        return 1;
    }
    (*dcode).lock = req;
    0
}

/// Release a decoder lock
/// Returns 0 on success
#[no_mangle]
pub unsafe extern "C" fn _zbar_decoder_release_lock(
    dcode: *mut zbar_decoder_t,
    req: zbar_symbol_type_t,
) -> c_char {
    debug_assert_eq!((*dcode).lock, req, "lock={} req={}", (*dcode).lock, req);
    (*dcode).lock = 0;
    0
}

/// Resize the decoder's data buffer if needed
/// Returns 1 on allocation failure or if max size exceeded, 0 on success
#[no_mangle]
pub unsafe extern "C" fn _zbar_decoder_size_buf(dcode: *mut zbar_decoder_t, len: c_uint) -> c_char {
    if len <= BUFFER_MIN {
        return 0;
    }
    if len < (*dcode).buf_alloc {
        // FIXME: size reduction heuristic?
        return 0;
    }
    if len > BUFFER_MAX {
        return 1;
    }

    let mut new_len = len;
    if len < (*dcode).buf_alloc + BUFFER_INCR {
        new_len = (*dcode).buf_alloc + BUFFER_INCR;
        if new_len > BUFFER_MAX {
            new_len = BUFFER_MAX;
        }
    }

    let buf = libc::realloc((*dcode).buf as *mut libc::c_void, new_len as usize) as *mut c_char;
    if buf.is_null() {
        return 1;
    }

    (*dcode).buf = buf;
    (*dcode).buf_alloc = new_len;
    0
}

// ============================================================================
// Decoder-specific reset functions
// ============================================================================

/// Reset codabar decoder state
#[no_mangle]
pub unsafe extern "C" fn _zbar_codabar_reset(codabar: *mut codabar_decoder_t) {
    (*codabar).set_direction(false);
    (*codabar).set_element(0);
    (*codabar).set_character(-1);
    (*codabar).s7 = 0;
}

/// Reset code128 decoder state
#[no_mangle]
pub unsafe extern "C" fn _zbar_code128_reset(dcode128: *mut code128_decoder_t) {
    (*dcode128).set_direction(0);
    (*dcode128).set_element(0);
    (*dcode128).set_character(-1);
    (*dcode128).s6 = 0;
}

/// Reset code39 decoder state
#[no_mangle]
pub unsafe extern "C" fn _zbar_code39_reset(dcode39: *mut code39_decoder_t) {
    (*dcode39).set_direction(false);
    (*dcode39).set_element(0);
    (*dcode39).set_character(-1);
    (*dcode39).s9 = 0;
}

/// Reset code93 decoder state
#[no_mangle]
pub unsafe extern "C" fn _zbar_code93_reset(dcode93: *mut code93_decoder_t) {
    (*dcode93).set_direction(false);
    (*dcode93).set_element(0);
    (*dcode93).set_character(-1);
}

/// Reset i25 decoder state
#[no_mangle]
pub unsafe extern "C" fn _zbar_i25_reset(i25: *mut i25_decoder_t) {
    (*i25).set_direction(false);
    (*i25).set_element(0);
    (*i25).set_character(-1);
    (*i25).s10 = 0;
}

/// Reset QR finder state
#[no_mangle]
pub unsafe extern "C" fn _zbar_qr_finder_reset(qrf: *mut qr_finder_t) {
    (*qrf).s5 = 0;
}

/// Prepare DataBar decoder for new scan
#[no_mangle]
pub unsafe extern "C" fn _zbar_databar_new_scan(db: *mut databar_decoder_t) {
    for i in 0..16 {
        if (*db).chars[i] >= 0 {
            let seg = ((*db).segs).offset((*db).chars[i] as isize);
            if (*seg).partial() {
                (*seg).set_finder(-1);
            }
            (*db).chars[i] = -1;
        }
    }
}

/// Reset DataBar decoder state
#[no_mangle]
pub unsafe extern "C" fn _zbar_databar_reset(db: *mut databar_decoder_t) {
    let n = (*db).csegs() as isize;
    _zbar_databar_new_scan(db);
    for i in 0..n {
        let seg = ((*db).segs).offset(i);
        (*seg).set_finder(-1);
    }
}

/// Allocate a new DataBar segment (or reuse an old one)
/// Returns the index of the allocated segment, or -1 on failure
#[no_mangle]
pub unsafe extern "C" fn _zbar_databar_alloc_segment(db: *mut databar_decoder_t) -> c_int {
    use crate::decoder_types::databar_segment_t;
    use std::mem::size_of;

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
            let new_ptr = libc::realloc(
                (*db).segs as *mut c_void,
                new_csegs * size_of::<databar_segment_t>(),
            ) as *mut databar_segment_t;

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
                    (*seg).set_color(0);
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

/// Convert DataBar data from heterogeneous base {1597,2841} to base 10 character representation
#[no_mangle]
pub unsafe extern "C" fn _zbar_databar_postprocess(dcode: *mut zbar_decoder_t, d: *mut c_uint) {
    let db = &mut (*dcode).databar;
    let mut d_array = [*d.offset(0), *d.offset(1), *d.offset(2), *d.offset(3)];

    let buf = (*dcode).buf;
    let mut chk = 0u32;

    // Write "01" prefix
    *buf.offset(0) = b'0' as c_char;
    *buf.offset(1) = b'1' as c_char;

    // Start at position 15 and work backwards
    let mut buf_idx = 15;

    // Write two null terminators
    *buf.offset(buf_idx) = 0;
    buf_idx -= 1;
    *buf.offset(buf_idx) = 0;
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
        *buf.offset(buf_idx) = (c as u8 + b'0') as c_char;
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
        *buf.offset(buf_idx) = (c as u8 + b'0') as c_char;
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
        *buf.offset(buf_idx) = (c as u8 + b'0') as c_char;
        buf_idx -= 1;
        if i != 0 {
            r /= 10;
        }
    }

    // Add check digit if configured
    if ((db.config >> ZBAR_CFG_EMIT_CHECK) & 1) != 0 {
        chk %= 10;
        if chk != 0 {
            chk = 10 - chk;
        }
        *buf.offset(13) = (chk as u8 + b'0') as c_char;
        (*dcode).buflen = 14;
    } else {
        (*dcode).buflen = 13;
    }
}

/// Merge or update a DataBar segment with existing segments
#[no_mangle]
pub unsafe extern "C" fn _zbar_databar_merge_segment(
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

/// Prepare EAN decoder for new scan
pub unsafe fn _zbar_ean_new_scan(ean: *mut ean_decoder_t) {
    (*ean).pass[0].state = -1;
    (*ean).pass[1].state = -1;
    (*ean).pass[2].state = -1;
    (*ean).pass[3].state = -1;
    (*ean).s4 = 0;
}

/// Reset EAN decoder state
pub unsafe fn _zbar_ean_reset(ean: *mut ean_decoder_t) {
    _zbar_ean_new_scan(ean);
    (*ean).left = 0; // ZBAR_NONE
    (*ean).right = 0; // ZBAR_NONE
}

// ============================================================================
// Decoder lifecycle functions
// ============================================================================

/// Create a new decoder instance
#[no_mangle]
pub unsafe extern "C" fn zbar_decoder_create() -> *mut zbar_decoder_t {
    let dcode = libc::calloc(1, std::mem::size_of::<zbar_decoder_t>()) as *mut zbar_decoder_t;

    (*dcode).buf_alloc = BUFFER_MIN;
    (*dcode).buf = libc::malloc((*dcode).buf_alloc as usize) as *mut c_char;

    // Initialize default configs
    (*dcode).ean.enable = 1;
    (*dcode).ean.ean13_config = (1 << ZBAR_CFG_ENABLE) | (1 << ZBAR_CFG_EMIT_CHECK);
    (*dcode).ean.ean8_config = (1 << ZBAR_CFG_ENABLE) | (1 << ZBAR_CFG_EMIT_CHECK);
    (*dcode).ean.upca_config = 1 << ZBAR_CFG_EMIT_CHECK;
    (*dcode).ean.upce_config = 1 << ZBAR_CFG_EMIT_CHECK;
    (*dcode).ean.isbn10_config = 1 << ZBAR_CFG_EMIT_CHECK;
    (*dcode).ean.isbn13_config = 1 << ZBAR_CFG_EMIT_CHECK;
    // FIXME_ADDON_SYNC not defined, skip ean2/ean5 config

    (*dcode).i25.config = 1 << ZBAR_CFG_ENABLE;
    cfg_set(&mut (*dcode).i25.configs, ZBAR_CFG_MIN_LEN, 6);

    (*dcode).databar.config = (1 << ZBAR_CFG_ENABLE) | (1 << ZBAR_CFG_EMIT_CHECK);
    (*dcode).databar.config_exp = (1 << ZBAR_CFG_ENABLE) | (1 << ZBAR_CFG_EMIT_CHECK);
    (*dcode).databar.set_csegs(4);
    (*dcode).databar.segs = libc::calloc(
        4,
        std::mem::size_of::<crate::decoder_types::databar_segment_t>(),
    ) as *mut crate::decoder_types::databar_segment_t;

    (*dcode).codabar.config = 1 << ZBAR_CFG_ENABLE;
    cfg_set(&mut (*dcode).codabar.configs, ZBAR_CFG_MIN_LEN, 4);

    (*dcode).code39.config = 1 << ZBAR_CFG_ENABLE;
    cfg_set(&mut (*dcode).code39.configs, ZBAR_CFG_MIN_LEN, 1);

    (*dcode).code93.config = 1 << ZBAR_CFG_ENABLE;
    (*dcode).code128.config = 1 << ZBAR_CFG_ENABLE;
    (*dcode).qrf.config = 1 << ZBAR_CFG_ENABLE;
    (*dcode).sqf.config = 1 << ZBAR_CFG_ENABLE;

    zbar_decoder_reset(dcode);
    dcode
}

/// Destroy a decoder instance
#[no_mangle]
pub unsafe extern "C" fn zbar_decoder_destroy(dcode: *mut zbar_decoder_t) {
    if !(*dcode).databar.segs.is_null() {
        libc::free((*dcode).databar.segs as *mut c_void);
    }
    if !(*dcode).buf.is_null() {
        libc::free((*dcode).buf as *mut c_void);
    }
    libc::free(dcode as *mut c_void);
}

/// Reset decoder to initial state
#[no_mangle]
pub unsafe extern "C" fn zbar_decoder_reset(dcode: *mut zbar_decoder_t) {
    // Calculate the offset to buf_alloc field
    let reset_size = &(*dcode).buf_alloc as *const _ as usize - dcode as usize;
    std::ptr::write_bytes(dcode as *mut u8, 0, reset_size);

    _zbar_ean_reset(&mut (*dcode).ean);
    _zbar_i25_reset(&mut (*dcode).i25);
    _zbar_databar_reset(&mut (*dcode).databar);
    _zbar_codabar_reset(&mut (*dcode).codabar);
    _zbar_code39_reset(&mut (*dcode).code39);
    _zbar_code93_reset(&mut (*dcode).code93);
    _zbar_code128_reset(&mut (*dcode).code128);
    _zbar_qr_finder_reset(&mut (*dcode).qrf);
}

/// Mark start of a new scan pass
///
/// Clears any intra-symbol state and resets color to ZBAR_SPACE.
/// Any partially decoded symbol state is retained.
pub unsafe fn zbar_decoder_new_scan(dcode: *mut zbar_decoder_t) {
    // Soft reset decoder
    std::ptr::write_bytes((*dcode).w.as_mut_ptr(), 0, (*dcode).w.len());
    (*dcode).lock = 0;
    (*dcode).idx = 0;
    (*dcode).s6 = 0;

    _zbar_ean_new_scan(&mut (*dcode).ean);
    _zbar_i25_reset(&mut (*dcode).i25);
    _zbar_databar_new_scan(&mut (*dcode).databar);
    _zbar_codabar_reset(&mut (*dcode).codabar);
    _zbar_code39_reset(&mut (*dcode).code39);
    _zbar_code93_reset(&mut (*dcode).code93);
    _zbar_code128_reset(&mut (*dcode).code128);
    _zbar_qr_finder_reset(&mut (*dcode).qrf);
}

// ============================================================================
// Decoder accessor functions
// ============================================================================

/// Get current decoder color (bar/space)
#[no_mangle]
pub unsafe extern "C" fn zbar_decoder_get_color(dcode: *const zbar_decoder_t) -> c_int {
    _zbar_decoder_get_color(dcode) as c_int
}

/// Get decoded data buffer
#[no_mangle]
pub unsafe extern "C" fn zbar_decoder_get_data(dcode: *const zbar_decoder_t) -> *const c_char {
    (*dcode).buf as *const c_char
}

/// Get length of decoded data
#[no_mangle]
pub unsafe extern "C" fn zbar_decoder_get_data_length(dcode: *const zbar_decoder_t) -> c_uint {
    (*dcode).buflen
}

/// Get decode direction
#[no_mangle]
pub unsafe extern "C" fn zbar_decoder_get_direction(dcode: *const zbar_decoder_t) -> c_int {
    (*dcode).direction
}

/// Set decoder callback handler
#[no_mangle]
pub unsafe extern "C" fn zbar_decoder_set_handler(
    dcode: *mut zbar_decoder_t,
    handler: Option<crate::decoder_types::zbar_decoder_handler_t>,
) -> Option<crate::decoder_types::zbar_decoder_handler_t> {
    let result = (*dcode).handler;
    (*dcode).handler = handler;
    result
}

/// Set user data pointer
#[no_mangle]
pub unsafe extern "C" fn zbar_decoder_set_userdata(
    dcode: *mut zbar_decoder_t,
    userdata: *mut c_void,
) {
    (*dcode).userdata = userdata;
}

/// Get user data pointer
#[no_mangle]
pub unsafe extern "C" fn zbar_decoder_get_userdata(dcode: *const zbar_decoder_t) -> *mut c_void {
    (*dcode).userdata
}

/// Get decoded symbol type
#[no_mangle]
pub unsafe extern "C" fn zbar_decoder_get_type(dcode: *const zbar_decoder_t) -> zbar_symbol_type_t {
    (*dcode).type_
}

/// Get decoded symbol modifiers
#[no_mangle]
pub unsafe extern "C" fn zbar_decoder_get_modifiers(dcode: *const zbar_decoder_t) -> c_uint {
    (*dcode).modifiers
}

// ============================================================================
// Main decode function
// ============================================================================

/// Process next bar/space width from input stream
///
/// The width is in arbitrary relative units. First value of a scan
/// is ZBAR_SPACE width, alternating from there.
///
/// # Returns
/// - Appropriate symbol type if width completes decode of a symbol (data is available for retrieval)
/// - ZBAR_PARTIAL as a hint if part of a symbol was decoded
/// - ZBAR_NONE (0) if no new symbol data is available
pub unsafe fn zbar_decode_width(dcode: *mut zbar_decoder_t, w: c_uint) -> zbar_symbol_type_t {
    let mut sym = 0; // ZBAR_NONE

    // Store width in circular buffer
    (*dcode).w[((*dcode).idx & (DECODE_WINDOW - 1) as u8) as usize] = w;

    // Update shared character width
    (*dcode).s6 = (*dcode).s6.wrapping_sub(_zbar_decoder_get_width(dcode, 7));
    (*dcode).s6 = (*dcode).s6.wrapping_add(_zbar_decoder_get_width(dcode, 1));

    // Each decoder processes width stream in parallel
    if test_cfg((*dcode).qrf.config, ZBAR_CFG_ENABLE) {
        let tmp = _zbar_find_qr(dcode);
        if tmp > ZBAR_PARTIAL {
            sym = tmp;
        }
    }

    if (*dcode).ean.enable != 0 {
        let tmp = _zbar_decode_ean(dcode);
        if tmp != 0 {
            sym = tmp;
        }
    }

    if test_cfg((*dcode).code39.config, ZBAR_CFG_ENABLE) {
        let tmp = _zbar_decode_code39(dcode);
        if tmp > ZBAR_PARTIAL {
            sym = tmp;
        }
    }

    if test_cfg((*dcode).code93.config, ZBAR_CFG_ENABLE) {
        let tmp = _zbar_decode_code93(dcode);
        if tmp > ZBAR_PARTIAL {
            sym = tmp;
        }
    }

    if test_cfg((*dcode).code128.config, ZBAR_CFG_ENABLE) {
        let tmp = _zbar_decode_code128(dcode);
        if tmp > ZBAR_PARTIAL {
            sym = tmp;
        }
    }

    if test_cfg(
        (*dcode).databar.config | (*dcode).databar.config_exp,
        ZBAR_CFG_ENABLE,
    ) {
        let tmp = _zbar_decode_databar(dcode);
        if tmp > ZBAR_PARTIAL {
            sym = tmp;
        }
    }

    if test_cfg((*dcode).codabar.config, ZBAR_CFG_ENABLE) {
        let tmp = _zbar_decode_codabar(dcode);
        if tmp > ZBAR_PARTIAL {
            sym = tmp;
        }
    }

    if test_cfg((*dcode).i25.config, ZBAR_CFG_ENABLE) {
        let tmp = _zbar_decode_i25(dcode);
        if tmp > ZBAR_PARTIAL {
            sym = tmp;
        }
    }

    (*dcode).idx = (*dcode).idx.wrapping_add(1);
    (*dcode).type_ = sym;

    if sym != 0 {
        if (*dcode).lock != 0 && sym > ZBAR_PARTIAL && sym != ZBAR_QRCODE {
            _zbar_decoder_release_lock(dcode, sym);
        }
        if let Some(handler) = (*dcode).handler {
            handler(dcode);
        }
    }

    sym
}

// ============================================================================
// Configuration functions
// ============================================================================

/// Get configuration pointer for a symbology (internal helper)
unsafe fn decoder_get_configp(
    dcode: *const zbar_decoder_t,
    sym: zbar_symbol_type_t,
) -> *const c_uint {
    match sym {
        ZBAR_EAN13 => &(*dcode).ean.ean13_config as *const c_uint,
        ZBAR_EAN2 => &(*dcode).ean.ean2_config as *const c_uint,
        ZBAR_EAN5 => &(*dcode).ean.ean5_config as *const c_uint,
        ZBAR_EAN8 => &(*dcode).ean.ean8_config as *const c_uint,
        ZBAR_UPCA => &(*dcode).ean.upca_config as *const c_uint,
        ZBAR_UPCE => &(*dcode).ean.upce_config as *const c_uint,
        ZBAR_ISBN10 => &(*dcode).ean.isbn10_config as *const c_uint,
        ZBAR_ISBN13 => &(*dcode).ean.isbn13_config as *const c_uint,
        ZBAR_I25 => &(*dcode).i25.config as *const c_uint,
        ZBAR_DATABAR => &(*dcode).databar.config as *const c_uint,
        ZBAR_DATABAR_EXP => &(*dcode).databar.config_exp as *const c_uint,
        ZBAR_CODABAR => &(*dcode).codabar.config as *const c_uint,
        ZBAR_CODE39 => &(*dcode).code39.config as *const c_uint,
        ZBAR_CODE93 => &(*dcode).code93.config as *const c_uint,
        ZBAR_CODE128 => &(*dcode).code128.config as *const c_uint,
        ZBAR_QRCODE => &(*dcode).qrf.config as *const c_uint,
        ZBAR_SQCODE => &(*dcode).sqf.config as *const c_uint,
        _ => std::ptr::null(),
    }
}

/// Get all configurations for a symbology
#[no_mangle]
pub unsafe extern "C" fn zbar_decoder_get_configs(
    dcode: *const zbar_decoder_t,
    sym: zbar_symbol_type_t,
) -> c_uint {
    let config = decoder_get_configp(dcode, sym);
    if config.is_null() {
        0
    } else {
        *config
    }
}

/// Set boolean configuration (internal helper)
unsafe fn decoder_set_config_bool(
    dcode: *mut zbar_decoder_t,
    sym: zbar_symbol_type_t,
    cfg: c_int,
    val: c_int,
) -> c_int {
    let config = decoder_get_configp(dcode, sym) as *mut c_uint;
    if config.is_null() || cfg >= ZBAR_CFG_NUM {
        return 1;
    }

    if val == 0 {
        *config &= !(1 << cfg);
    } else if val == 1 {
        *config |= 1 << cfg;
    } else {
        return 1;
    }

    // Update EAN enable flag
    (*dcode).ean.enable = if test_cfg(
        (*dcode).ean.ean13_config
            | (*dcode).ean.ean2_config
            | (*dcode).ean.ean5_config
            | (*dcode).ean.ean8_config
            | (*dcode).ean.upca_config
            | (*dcode).ean.upce_config
            | (*dcode).ean.isbn10_config
            | (*dcode).ean.isbn13_config,
        ZBAR_CFG_ENABLE,
    ) {
        1
    } else {
        0
    };

    0
}

/// Set integer configuration (internal helper)
unsafe fn decoder_set_config_int(
    dcode: *mut zbar_decoder_t,
    sym: zbar_symbol_type_t,
    cfg: c_int,
    val: c_int,
) -> c_int {
    match sym {
        ZBAR_I25 => {
            (*dcode).i25.configs[(cfg - ZBAR_CFG_MIN_LEN) as usize] = val;
            0
        }
        ZBAR_CODABAR => {
            (*dcode).codabar.configs[(cfg - ZBAR_CFG_MIN_LEN) as usize] = val;
            0
        }
        ZBAR_CODE39 => {
            (*dcode).code39.configs[(cfg - ZBAR_CFG_MIN_LEN) as usize] = val;
            0
        }
        ZBAR_CODE93 => {
            (*dcode).code93.configs[(cfg - ZBAR_CFG_MIN_LEN) as usize] = val;
            0
        }
        ZBAR_CODE128 => {
            (*dcode).code128.configs[(cfg - ZBAR_CFG_MIN_LEN) as usize] = val;
            0
        }
        _ => 1,
    }
}

/// Get decoder configuration value
#[no_mangle]
pub unsafe extern "C" fn zbar_decoder_get_config(
    dcode: *mut zbar_decoder_t,
    sym: zbar_symbol_type_t,
    cfg: c_int,
    val: *mut c_int,
) -> c_int {
    let config = decoder_get_configp(dcode, sym);

    // Return error if symbol doesn't have config
    if sym <= ZBAR_PARTIAL || sym > ZBAR_CODE128 || sym == ZBAR_COMPOSITE {
        return 1;
    }

    // Return decoder boolean configs
    if cfg < ZBAR_CFG_NUM {
        *val = if (*config & (1 << cfg)) != 0 { 1 } else { 0 };
        return 0;
    }

    // Return decoder integer configs
    if (ZBAR_CFG_MIN_LEN..=ZBAR_CFG_MAX_LEN).contains(&cfg) {
        match sym {
            ZBAR_I25 => {
                *val = (*dcode).i25.configs[(cfg - ZBAR_CFG_MIN_LEN) as usize];
                0
            }
            ZBAR_CODABAR => {
                *val = (*dcode).codabar.configs[(cfg - ZBAR_CFG_MIN_LEN) as usize];
                0
            }
            ZBAR_CODE39 => {
                *val = (*dcode).code39.configs[(cfg - ZBAR_CFG_MIN_LEN) as usize];
                0
            }
            ZBAR_CODE93 => {
                *val = (*dcode).code93.configs[(cfg - ZBAR_CFG_MIN_LEN) as usize];
                0
            }
            ZBAR_CODE128 => {
                *val = (*dcode).code128.configs[(cfg - ZBAR_CFG_MIN_LEN) as usize];
                0
            }
            _ => 1,
        }
    } else {
        1
    }
}

/// Set decoder configuration
#[no_mangle]
pub unsafe extern "C" fn zbar_decoder_set_config(
    dcode: *mut zbar_decoder_t,
    sym: zbar_symbol_type_t,
    cfg: c_int,
    val: c_int,
) -> c_int {
    // If ZBAR_NONE, set config for all symbologies
    if sym == ZBAR_NONE {
        const ALL: [zbar_symbol_type_t; 17] = [
            ZBAR_EAN13,
            ZBAR_EAN2,
            ZBAR_EAN5,
            ZBAR_EAN8,
            ZBAR_UPCA,
            ZBAR_UPCE,
            ZBAR_ISBN10,
            ZBAR_ISBN13,
            ZBAR_I25,
            ZBAR_DATABAR,
            ZBAR_DATABAR_EXP,
            ZBAR_CODABAR,
            ZBAR_CODE39,
            ZBAR_CODE93,
            ZBAR_CODE128,
            ZBAR_QRCODE,
            ZBAR_SQCODE,
        ];
        for &s in &ALL {
            zbar_decoder_set_config(dcode, s, cfg, val);
        }
        return 0;
    }

    if (0..ZBAR_CFG_NUM).contains(&cfg) {
        decoder_set_config_bool(dcode, sym, cfg, val)
    } else if (ZBAR_CFG_MIN_LEN..=ZBAR_CFG_MAX_LEN).contains(&cfg) {
        decoder_set_config_int(dcode, sym, cfg, val)
    } else {
        1
    }
}

// ============================================================================
// Debug helper
// ============================================================================

use std::sync::Mutex;

static DECODER_DUMP: Mutex<Option<Vec<u8>>> = Mutex::new(None);

/// Format decoder buffer as hex string (for debugging)
#[no_mangle]
pub unsafe extern "C" fn _zbar_decoder_buf_dump(buf: *mut u8, buflen: c_uint) -> *const c_char {
    let dumplen = (buflen * 3) + 12;
    let mut dump = DECODER_DUMP.lock().unwrap();

    // Allocate or reallocate buffer
    if dump.is_none() || dump.as_ref().unwrap().len() < dumplen as usize {
        *dump = Some(Vec::with_capacity(dumplen as usize));
    }

    let dump_vec = dump.as_mut().unwrap();
    dump_vec.clear();

    // Format header
    let len_display = if buflen > 0xffff { 0xffff } else { buflen };
    let header = format!("buf[{:04x}]=", len_display);
    dump_vec.extend_from_slice(header.as_bytes());

    // Format buffer contents as hex
    let slice = std::slice::from_raw_parts(buf, buflen as usize);
    for (i, &byte) in slice.iter().enumerate() {
        if i > 0 {
            dump_vec.push(b' ');
        }
        let hex = format!("{:02x}", byte);
        dump_vec.extend_from_slice(hex.as_bytes());
    }

    dump_vec.push(0); // Null terminator
    dump_vec.as_ptr() as *const c_char
}

/// Low-level decoder for processing bar/space width streams
pub struct Decoder {
    _ptr: *mut std::ffi::c_void,
}

impl Decoder {
    /// Create a new decoder
    pub fn new() -> Self {
        // Note: This is a placeholder - the C library doesn't expose the decoder directly
        // In a full conversion, we'd implement the decoder logic in Rust
        todo!("Decoder needs to be implemented as part of the Rust conversion")
    }
}

impl Default for Decoder {
    fn default() -> Self {
        Self::new()
    }
}
