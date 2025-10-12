/*
    This is a Rust port of zbar/qrcode/qrdectxt.c
    Copyright (C) 2008-2010  Timothy B. Terriberry (tterribe@xiph.org)
*/

use std::collections::VecDeque;
use std::os::raw::{c_char, c_int, c_uchar, c_uint};
use std::ptr::null_mut;
use std::slice;

use encoding_rs::{Encoding, BIG5, SHIFT_JIS, UTF_8, WINDOWS_1252};

use crate::decoder_types::{
    ZBAR_CFG_BINARY, ZBAR_MOD_AIM, ZBAR_MOD_GS1, ZBAR_PARTIAL, ZBAR_QRCODE,
};
use crate::image_ffi::{zbar_image_t, zbar_symbol_t};
use crate::img_scanner::{
    _zbar_image_scanner_add_sym, _zbar_image_scanner_alloc_sym, _zbar_image_scanner_recycle_syms,
    zbar_image_scanner_get_config, zbar_image_scanner_t,
};
use crate::symbol::{_zbar_symbol_add_point, _zbar_symbol_set_create};

use super::qr_point;
use super::qrdec::qr_mode;

fn qr_mode_has_data(_mode: qr_mode) -> bool {
    matches!(
        _mode,
        qr_mode::QR_MODE_NUM
            | qr_mode::QR_MODE_ALNUM
            | qr_mode::QR_MODE_BYTE
            | qr_mode::QR_MODE_KANJI
    )
}

pub union qr_code_data_entry_payload {
    pub data: qr_code_data_entry_data,
    pub eci: c_uint,
    pub ai: c_int,
    pub sa: qr_code_data_entry_sa,
}

#[derive(Clone, Copy)]
pub struct qr_code_data_entry_data {
    pub buf: *mut c_uchar,
    pub len: c_int,
}

#[derive(Clone, Copy)]
pub struct qr_code_data_entry_sa {
    pub sa_index: c_uchar,
    pub sa_size: c_uchar,
    pub sa_parity: c_uchar,
}

pub struct qr_code_data_entry {
    pub mode: qr_mode,
    pub payload: qr_code_data_entry_payload,
}

pub struct qr_code_data {
    pub entries: *mut qr_code_data_entry,
    pub nentries: c_int,
    pub version: c_uchar,
    pub ecc_level: c_uchar,
    pub sa_index: c_uchar,
    pub sa_size: c_uchar,
    pub sa_parity: c_uchar,
    pub self_parity: c_uchar,
    pub bbox: [qr_point; 4],
}

pub struct qr_code_data_list {
    pub qrdata: *mut qr_code_data,
    pub nqrdata: c_int,
    pub cqrdata: c_int,
}

fn text_is_ascii(text: &[u8]) -> bool {
    text.iter().all(|&c| c < 0x80)
}

fn text_is_latin1(text: &[u8]) -> bool {
    text.iter().all(|&c| !(0x80..0xA0).contains(&c))
}

fn text_is_big5(text: &[u8]) -> bool {
    let mut i = 0;
    while i < text.len() {
        if text[i] == 0xFF {
            return false;
        } else if text[i] >= 0x80 {
            i += 1;
            if i >= text.len() {
                return false;
            }
            if text[i] < 0x40 || (text[i] > 0x7E && text[i] < 0xA1) || text[i] > 0xFE {
                return false;
            }
        }
        i += 1;
    }
    true
}

fn enc_list_mtf(enc_list: &mut VecDeque<&'static Encoding>, enc: &'static Encoding) {
    if let Some(pos) = enc_list.iter().position(|&e| e == enc) {
        let e = enc_list.remove(pos).unwrap();
        enc_list.push_front(e);
    }
}

unsafe fn sym_add_point(sym: *mut zbar_symbol_t, x: c_int, y: c_int) {
    _zbar_symbol_add_point(sym, x, y);
}

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers passed from C.
/// The caller must ensure that the pointers are valid and that the data they
/// point to has the expected layout.
pub unsafe fn qr_code_data_list_extract_text(
    _qrlist: *const qr_code_data_list,
    iscn: *mut zbar_image_scanner_t,
    _img: *mut zbar_image_t,
) -> c_int {
    let mut raw_binary: c_int = 0;
    zbar_image_scanner_get_config(iscn, ZBAR_QRCODE, ZBAR_CFG_BINARY, &mut raw_binary);

    let qrlist = &*_qrlist;
    let qrdata = slice::from_raw_parts(qrlist.qrdata, qrlist.nqrdata as usize);
    let mut mark = vec![0u8; qrlist.nqrdata as usize];
    let ntext = 0;

    for i in 0..qrlist.nqrdata as usize {
        if mark[i] == 0 {
            let mut sa: [c_int; 16] = [-1; 16];
            let sa_size;

            if qrdata[i].sa_size > 0 {
                sa_size = qrdata[i].sa_size as usize;
                let sa_parity = qrdata[i].sa_parity;
                for j in i..qrlist.nqrdata as usize {
                    if mark[j] == 0
                        && qrdata[j].sa_size as usize == sa_size
                        && qrdata[j].sa_parity == sa_parity
                        && sa[qrdata[j].sa_index as usize] < 0
                    {
                        sa[qrdata[j].sa_index as usize] = j as c_int;
                        mark[j] = 1;
                    }
                }
            } else {
                sa[0] = i as c_int;
                sa_size = 1;
            }

            let mut sa_ctext = 0;
            let mut fnc1 = 0;
            let mut fnc1_2ai = 0;
            let mut has_kanji = false;

            for j in 0..sa_size {
                if sa[j] >= 0 {
                    let qrdataj = &qrdata[sa[j] as usize];
                    let entries = slice::from_raw_parts(qrdataj.entries, qrdataj.nentries as usize);
                    for entry in entries {
                        let shift = match entry.mode {
                            qr_mode::QR_MODE_FNC1_1ST => {
                                if fnc1 == 0 {
                                    fnc1 = 1 << ZBAR_MOD_GS1;
                                }
                                0
                            }
                            qr_mode::QR_MODE_FNC1_2ND => {
                                if fnc1 == 0 {
                                    fnc1 = 1 << ZBAR_MOD_AIM;
                                    fnc1_2ai = entry.payload.ai;
                                    sa_ctext += 2;
                                }
                                0
                            }
                            qr_mode::QR_MODE_KANJI => {
                                has_kanji = true;
                                2
                            }
                            qr_mode::QR_MODE_BYTE => 2,
                            _ => {
                                if qr_mode_has_data(entry.mode) {
                                    sa_ctext += entry.payload.data.len as usize;
                                }
                                0
                            }
                        };
                        if shift > 0 {
                            sa_ctext += (entry.payload.data.len as usize) << shift;
                        }
                    }
                }
            }

            let mut sa_text: Vec<u8> = Vec::with_capacity(sa_ctext + 1);
            if fnc1 == (1 << ZBAR_MOD_AIM) {
                if fnc1_2ai < 100 {
                    sa_text.push(b'0' + (fnc1_2ai / 10) as u8);
                    sa_text.push(b'0' + (fnc1_2ai % 10) as u8);
                } else {
                    sa_text.push((fnc1_2ai - 100) as u8);
                }
            }

            let mut eci: Option<&'static Encoding> = None;
            let mut enc_list: VecDeque<&'static Encoding> =
                VecDeque::from(vec![SHIFT_JIS, WINDOWS_1252, BIG5, UTF_8]);
            let mut err = false;
            let mut bytebuf: Vec<u8> = Vec::with_capacity(sa_ctext + 1);
            let mut syms_head: *mut zbar_symbol_t = null_mut();
            let mut sym_cur = &mut syms_head;

            let mut j = 0;
            while j < sa_size {
                if err {
                    break;
                }
                let sym = _zbar_image_scanner_alloc_sym(iscn, ZBAR_QRCODE, 0);
                *sym_cur = sym;
                (*sym).datalen = sa_text.len() as u32;

                if sa[j] < 0 {
                    (*sym).symbol_type = ZBAR_PARTIAL;
                    let mut k = j + 1;
                    while k < sa_size && sa[k] < 0 {
                        k += 1;
                    }
                    j = k;
                    if j >= sa_size {
                        break;
                    }
                    sa_text.push(0);
                    (*sym).datalen = sa_text.len() as u32;
                    sym_cur = &mut (*sym).next;
                    let next_sym = _zbar_image_scanner_alloc_sym(iscn, ZBAR_QRCODE, 0);
                    *sym_cur = next_sym;
                }

                let qrdataj = &qrdata[sa[j] as usize];
                sym_add_point(sym, qrdataj.bbox[0][0], qrdataj.bbox[0][1]);
                sym_add_point(sym, qrdataj.bbox[2][0], qrdataj.bbox[2][1]);
                sym_add_point(sym, qrdataj.bbox[3][0], qrdataj.bbox[3][1]);
                sym_add_point(sym, qrdataj.bbox[1][0], qrdataj.bbox[1][1]);

                let entries = slice::from_raw_parts(qrdataj.entries, qrdataj.nentries as usize);
                for entry in entries.iter() {
                    // Process byte buffer if needed
                    if !bytebuf.is_empty()
                        && (entry.mode != qr_mode::QR_MODE_BYTE
                            && entry.mode != qr_mode::QR_MODE_KANJI)
                    {
                        // convert bytes to text
                        if let Some(enc) = eci {
                            let (res, _enc, had_errors) = enc.decode(&bytebuf);
                            if had_errors {
                                err = true;
                                break;
                            }
                            sa_text.extend_from_slice(res.as_bytes());
                        } else if raw_binary != 0 {
                            sa_text.extend_from_slice(&bytebuf);
                        } else {
                            if has_kanji {
                                enc_list_mtf(&mut enc_list, SHIFT_JIS);
                            } else if bytebuf.starts_with(&[0xEF, 0xBB, 0xBF]) {
                                let (res, _enc, had_errors) = UTF_8.decode(&bytebuf[3..]);
                                if !had_errors {
                                    sa_text.extend_from_slice(res.as_bytes());
                                    enc_list_mtf(&mut enc_list, UTF_8);
                                } else {
                                    err = true;
                                }
                            } else if text_is_ascii(&bytebuf) {
                                enc_list_mtf(&mut enc_list, UTF_8);
                            } else if text_is_big5(&bytebuf) {
                                enc_list_mtf(&mut enc_list, BIG5);
                            }

                            if !err {
                                let mut decoded = false;
                                for &enc in &enc_list {
                                    if enc == WINDOWS_1252 && !text_is_latin1(&bytebuf) {
                                        continue;
                                    }
                                    let (res, _enc, had_errors) = enc.decode(&bytebuf);
                                    if !had_errors {
                                        sa_text.extend_from_slice(res.as_bytes());
                                        enc_list_mtf(&mut enc_list, enc);
                                        decoded = true;
                                        break;
                                    }
                                }
                                if !decoded {
                                    err = true;
                                }
                            }
                        }
                        bytebuf.clear();
                    }
                    if err {
                        break;
                    }

                    match entry.mode {
                        qr_mode::QR_MODE_NUM | qr_mode::QR_MODE_ALNUM => {
                            let data = slice::from_raw_parts(
                                entry.payload.data.buf,
                                entry.payload.data.len as usize,
                            );
                            sa_text.extend_from_slice(data);
                        }
                        qr_mode::QR_MODE_BYTE | qr_mode::QR_MODE_KANJI => {
                            let data = slice::from_raw_parts(
                                entry.payload.data.buf,
                                entry.payload.data.len as usize,
                            );
                            bytebuf.extend_from_slice(data);
                        }
                        qr_mode::QR_MODE_ECI => {
                            // Simplified ECI handling
                            eci = match entry.payload.eci {
                                3..=13 | 15..=18 => Some(WINDOWS_1252), // approx
                                20 => Some(SHIFT_JIS),
                                26 => Some(UTF_8),
                                _ => None,
                            };
                        }
                        _ => {}
                    }
                }
                sym_cur = &mut (*sym).next;
                j += 1;
            }

            if !err && !bytebuf.is_empty() {
                // convert bytes to text
                if let Some(enc) = eci {
                    let (res, _enc, had_errors) = enc.decode(&bytebuf);
                    if had_errors {
                        err = true;
                    } else {
                        sa_text.extend_from_slice(res.as_bytes());
                    }
                } else if raw_binary != 0 {
                    sa_text.extend_from_slice(&bytebuf);
                } else {
                    if has_kanji {
                        enc_list_mtf(&mut enc_list, SHIFT_JIS);
                    } else if bytebuf.starts_with(&[0xEF, 0xBB, 0xBF]) {
                        let (res, _enc, had_errors) = UTF_8.decode(&bytebuf[3..]);
                        if !had_errors {
                            sa_text.extend_from_slice(res.as_bytes());
                            enc_list_mtf(&mut enc_list, UTF_8);
                        } else {
                            err = true;
                        }
                    } else if text_is_ascii(&bytebuf) {
                        enc_list_mtf(&mut enc_list, UTF_8);
                    } else if text_is_big5(&bytebuf) {
                        enc_list_mtf(&mut enc_list, BIG5);
                    }

                    if !err {
                        let mut decoded = false;
                        for &enc in &enc_list {
                            if enc == WINDOWS_1252 && !text_is_latin1(&bytebuf) {
                                continue;
                            }
                            let (res, _enc, had_errors) = enc.decode(&bytebuf);
                            if !had_errors {
                                sa_text.extend_from_slice(res.as_bytes());
                                enc_list_mtf(&mut enc_list, enc);
                                decoded = true;
                                break;
                            }
                        }
                        if !decoded {
                            err = true;
                        }
                    }
                }
                bytebuf.clear();
            }

            if !err {
                sa_text.push(0);
                sa_text.shrink_to_fit();
                let len = sa_text.len();
                let ptr = sa_text.as_mut_ptr();
                std::mem::forget(sa_text);

                if sa_size == 1 {
                    let sym = syms_head;
                    (*sym).data = ptr as *mut c_char;
                    (*sym).datalen = (len - 1) as u32;
                    (*sym).data_alloc = len as u32;
                    (*sym).modifiers = fnc1 as u32;
                    _zbar_image_scanner_add_sym(iscn, sym);
                } else {
                    let sa_sym = _zbar_image_scanner_alloc_sym(iscn, ZBAR_QRCODE, 0);
                    (*sa_sym).syms = _zbar_symbol_set_create();
                    // This part is complex, involving restructuring the symbol set.
                    // For now, we just add the combined text to a single symbol.
                    // A full port would need to replicate the logic from the C code.
                    (*sa_sym).data = ptr as *mut c_char;
                    (*sa_sym).datalen = (len - 1) as u32;
                    (*sa_sym).data_alloc = len as u32;
                    (*sa_sym).modifiers = fnc1 as u32;
                    _zbar_image_scanner_add_sym(iscn, sa_sym);
                    _zbar_image_scanner_recycle_syms(iscn, syms_head);
                }
            } else {
                _zbar_image_scanner_recycle_syms(iscn, syms_head);
            }
        }
    }

    ntext
}
