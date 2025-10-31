/*
    This is a Rust port of zbar/qrcode/qrdectxt.c
    Copyright (C) 2008-2010  Timothy B. Terriberry (tterribe@xiph.org)
*/

use std::collections::VecDeque;
use std::os::raw::{c_int, c_uchar};

use encoding_rs::{Encoding, BIG5, SHIFT_JIS, UTF_8, WINDOWS_1252};

use crate::decoder::{ZBAR_CFG_BINARY, ZBAR_MOD_AIM, ZBAR_MOD_GS1};
use crate::img_scanner::{zbar_image_scanner_get_config, zbar_image_scanner_t, zbar_symbol_set_t};
use crate::qrcode::qrdec::{qr_code_data_list, qr_code_data_payload};
use crate::symbol::zbar_symbol_t;
use crate::SymbolType;

#[derive(Clone, Copy)]
pub struct qr_code_data_entry_sa {
    pub sa_index: c_uchar,
    pub sa_size: c_uchar,
    pub sa_parity: c_uchar,
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

unsafe fn sym_add_point(sym: &mut zbar_symbol_t, x: c_int, y: c_int) {
    sym.pts.push([x, y]);
}

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers passed from C.
/// The caller must ensure that the pointers are valid and that the data they
/// point to has the expected layout.
pub(crate) unsafe fn qr_code_data_list_extract_text(
    qrlist: &qr_code_data_list,
    iscn: &mut zbar_image_scanner_t,
) -> c_int {
    let raw_binary: c_int =
        zbar_image_scanner_get_config(iscn, SymbolType::QrCode, ZBAR_CFG_BINARY).unwrap_or(0);

    let qrdata = &qrlist.qrdata;
    let mut mark = vec![0u8; qrdata.len()];
    let ntext = 0;

    for i in 0..qrdata.len() {
        if mark[i] == 0 {
            let mut sa: [c_int; 16] = [-1; 16];
            let sa_size;

            if qrdata[i].sa_size > 0 {
                sa_size = qrdata[i].sa_size as usize;
                let sa_parity = qrdata[i].sa_parity;
                for j in i..qrdata.len() {
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
                    for entry in &qrdataj.entries {
                        match entry.payload {
                            qr_code_data_payload::Fnc1FirstPositionMarker => {
                                if fnc1 == 0 {
                                    fnc1 = 1 << ZBAR_MOD_GS1;
                                }
                            }
                            qr_code_data_payload::ApplicationIndicator(ai) => {
                                if fnc1 == 0 {
                                    fnc1 = 1 << ZBAR_MOD_AIM;
                                    fnc1_2ai = ai;
                                    sa_ctext += 2;
                                }
                            }
                            qr_code_data_payload::Kanji(ref data) => {
                                has_kanji = true;
                                sa_ctext += data.len() << 2;
                            }
                            qr_code_data_payload::Bytes(ref data) => {
                                sa_ctext += data.len() << 2;
                            }
                            qr_code_data_payload::Numeric(ref data)
                            | qr_code_data_payload::Alphanumeric(ref data) => {
                                sa_ctext += data.len();
                            }
                            qr_code_data_payload::StructuredAppendedHeaderData(_)
                            | qr_code_data_payload::ExtendedChannelInterpretation(_) => {
                                // does not count toward `sa_ctext`
                            }
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
            // Note: encoding_rs treats ISO-8859-1 as WINDOWS-1252 per WHATWG spec,
            // but for QR codes we want the actual ISO-8859-1 behavior (bytes 0x80-0x9F as-is)
            // So we use WINDOWS_1252 which is close enough for most cases
            // Order: Try UTF-8 first (strict), then SHIFT_JIS, then others, WINDOWS_1252 last (accepts anything)
            let mut enc_list: VecDeque<&'static Encoding> =
                VecDeque::from(vec![UTF_8, SHIFT_JIS, BIG5, WINDOWS_1252]);
            let mut err = false;
            let mut bytebuf: Vec<u8> = Vec::with_capacity(sa_ctext + 1);
            let mut component_syms = Vec::new();

            let mut j = 0;
            while j < sa_size {
                if err {
                    break;
                }
                let mut sym = iscn.alloc_sym(SymbolType::QrCode);

                if sa[j] < 0 {
                    sym.symbol_type = SymbolType::Partial;
                    component_syms.push(sym);
                    let mut k = j + 1;
                    while k < sa_size && sa[k] < 0 {
                        k += 1;
                    }
                    j = k;
                    if j >= sa_size {
                        break;
                    }
                    sym = iscn.alloc_sym(SymbolType::QrCode);
                }

                let qrdataj = &qrdata[sa[j] as usize];
                sym_add_point(&mut sym, qrdataj.bbox[0][0], qrdataj.bbox[0][1]);
                sym_add_point(&mut sym, qrdataj.bbox[2][0], qrdataj.bbox[2][1]);
                sym_add_point(&mut sym, qrdataj.bbox[3][0], qrdataj.bbox[3][1]);
                sym_add_point(&mut sym, qrdataj.bbox[1][0], qrdataj.bbox[1][1]);

                let entries = &qrdataj.entries;
                for entry in entries.iter() {
                    // Process byte buffer if needed
                    if !bytebuf.is_empty()
                        && !matches!(
                            entry.payload,
                            qr_code_data_payload::Bytes(_) | qr_code_data_payload::Kanji(_)
                        )
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

                    match &entry.payload {
                        qr_code_data_payload::Numeric(ref data)
                        | qr_code_data_payload::Alphanumeric(ref data) => {
                            sa_text.extend_from_slice(data);
                        }
                        qr_code_data_payload::Bytes(ref data)
                        | qr_code_data_payload::Kanji(ref data) => {
                            bytebuf.extend_from_slice(data);
                        }
                        qr_code_data_payload::ExtendedChannelInterpretation(val) => {
                            // Simplified ECI handling
                            eci = match val {
                                3..=13 | 15..=18 => Some(WINDOWS_1252), // approx
                                20 => Some(SHIFT_JIS),
                                26 => Some(UTF_8),
                                _ => None,
                            };
                        }
                        _ => {}
                    }
                }

                // Add the symbol to our collection
                component_syms.push(sym);
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
                        // Move WINDOWS_1252 to end if it has C1 control chars
                        if enc_list.front() == Some(&WINDOWS_1252) && !text_is_latin1(&bytebuf) {
                            enc_list.pop_front();
                            enc_list.push_back(WINDOWS_1252);
                        }

                        for &enc in &enc_list {
                            let (res, _enc, had_errors) = enc.decode(&bytebuf);
                            if !had_errors {
                                sa_text.extend_from_slice(res.as_bytes());
                                enc_list_mtf(&mut enc_list, enc);
                                decoded = true;
                                break;
                            }
                        }
                        // Note: Unlike some encoding errors, C does not set err=1 if decoding fails
                        // If no encoding worked, just copy the raw bytes
                        if !decoded {
                            sa_text.extend_from_slice(&bytebuf);
                        }
                    }
                }
                bytebuf.clear();
            }

            if !err {
                sa_text.shrink_to_fit();

                if sa_size == 1 {
                    // Single QR code - add it directly
                    if let Some(mut sym) = component_syms.into_iter().next() {
                        sym.data = sa_text;
                        sym.modifiers = fnc1 as u32;
                        iscn.add_symbol(sym);
                    }
                } else {
                    // Multiple QR codes - create structured append symbol
                    let mut sa_sym = iscn.alloc_sym(SymbolType::QrCode);
                    let component_set = zbar_symbol_set_t {
                        symbols: component_syms,
                    };
                    sa_sym.components = Some(component_set);
                    sa_sym.data = sa_text;
                    sa_sym.modifiers = fnc1 as u32;
                    iscn.add_symbol(sa_sym);
                }
            }
        }
    }

    ntext
}
