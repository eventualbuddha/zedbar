//! Low-level barcode decoder

use crate::{
    decoder_types::{
        databar_segment_t, zbar_decoder_t, zbar_symbol_type_t, DECODE_WINDOW, ZBAR_CFG_ENABLE,
        ZBAR_CFG_MAX_LEN, ZBAR_CFG_MIN_LEN, ZBAR_CODABAR, ZBAR_CODE128, ZBAR_CODE39, ZBAR_CODE93,
        ZBAR_COMPOSITE, ZBAR_DATABAR, ZBAR_DATABAR_EXP, ZBAR_EAN13, ZBAR_EAN2, ZBAR_EAN5,
        ZBAR_EAN8, ZBAR_I25, ZBAR_ISBN10, ZBAR_ISBN13, ZBAR_NONE, ZBAR_PARTIAL, ZBAR_QRCODE,
        ZBAR_SQCODE, ZBAR_UPCA, ZBAR_UPCE,
    },
    decoders::{
        codabar::_zbar_decode_codabar, code128::_zbar_decode_code128, code39::_zbar_decode_code39,
        code93::_zbar_decode_code93, databar::_zbar_decode_databar, ean::_zbar_decode_ean,
        i25::_zbar_decode_i25,
    },
    finder::find_qr,
};
use libc::{c_int, c_uint, c_void};
use std::mem::size_of;

// Config constant not in decoder_types
const ZBAR_CFG_NUM: c_int = 5;

#[inline]
fn test_cfg(config: c_uint, cfg: c_int) -> bool {
    ((config >> cfg) & 1) != 0
}

#[inline]
pub(crate) unsafe fn decoder_alloc_databar_segments(count: usize) -> *mut databar_segment_t {
    libc::calloc(count, size_of::<databar_segment_t>()) as *mut databar_segment_t
}

#[inline]
pub(crate) unsafe fn decoder_free_databar_segments(ptr: *mut databar_segment_t) {
    libc::free(ptr as *mut c_void);
}

#[inline]
pub(crate) unsafe fn decoder_realloc_databar_segments(
    ptr: *mut databar_segment_t,
    count: usize,
) -> *mut databar_segment_t {
    libc::realloc(ptr as *mut c_void, count * size_of::<databar_segment_t>())
        as *mut databar_segment_t
}

/// Decode element width into a discrete value
/// Returns -1 if the element width is invalid
pub unsafe fn _zbar_decoder_decode_e(e: c_uint, s: c_uint, n: c_uint) -> c_int {
    let big_e = ((e * n * 2 + 1) / s).wrapping_sub(3) / 2;
    if big_e >= n - 3 {
        -1
    } else {
        big_e as c_int
    }
}

// ============================================================================
// Decoder accessor functions
// ============================================================================

/// Set decoder callback handler
pub unsafe fn zbar_decoder_set_handler(
    dcode: *mut zbar_decoder_t,
    handler: Option<crate::decoder_types::zbar_decoder_handler_t>,
) -> Option<crate::decoder_types::zbar_decoder_handler_t> {
    if dcode.is_null() {
        return None;
    }
    let dcode = &mut *dcode;
    let result = dcode.handler;
    dcode.handler = handler;
    result
}

/// Set user data pointer
pub unsafe fn zbar_decoder_set_userdata(dcode: *mut zbar_decoder_t, userdata: *mut c_void) {
    if dcode.is_null() {
        return;
    }
    let dcode = &mut *dcode;
    dcode.userdata = userdata;
}

/// Get user data pointer
pub unsafe fn zbar_decoder_get_userdata(dcode: *const zbar_decoder_t) -> *mut c_void {
    if dcode.is_null() {
        return std::ptr::null_mut();
    }
    let dcode = &*dcode;
    dcode.userdata
}

/// Get decoded symbol type
pub unsafe fn zbar_decoder_get_type(dcode: *const zbar_decoder_t) -> zbar_symbol_type_t {
    if dcode.is_null() {
        return ZBAR_NONE;
    }
    let dcode = &*dcode;
    dcode.type_
}

/// Get decoded symbol modifiers
pub unsafe fn zbar_decoder_get_modifiers(dcode: *const zbar_decoder_t) -> c_uint {
    if dcode.is_null() {
        return 0;
    }
    let dcode = &*dcode;
    dcode.modifiers
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
    if dcode.is_null() {
        return ZBAR_NONE;
    }

    let dcode_ref = &mut *dcode;
    let mut sym = 0; // ZBAR_NONE

    // Store width in circular buffer
    dcode_ref.w[(dcode_ref.idx & (DECODE_WINDOW - 1) as u8) as usize] = w;

    // Update shared character width
    dcode_ref.s6 = dcode_ref.s6.wrapping_sub(dcode_ref.get_width(7));
    dcode_ref.s6 = dcode_ref.s6.wrapping_add(dcode_ref.get_width(1));

    // Each decoder processes width stream in parallel
    if test_cfg(dcode_ref.qrf.config, ZBAR_CFG_ENABLE) {
        let tmp = find_qr(dcode_ref);
        if tmp > ZBAR_PARTIAL {
            sym = tmp;
        }
    }

    if dcode_ref.ean.enable != 0 {
        let tmp = _zbar_decode_ean(&mut *dcode);
        if tmp != 0 {
            sym = tmp;
        }
    }

    if test_cfg(dcode_ref.code39.config, ZBAR_CFG_ENABLE) {
        let tmp = _zbar_decode_code39(&mut *dcode);
        if tmp > ZBAR_PARTIAL {
            sym = tmp;
        }
    }

    if test_cfg(dcode_ref.code93.config, ZBAR_CFG_ENABLE) {
        let tmp = _zbar_decode_code93(&mut *dcode);
        if tmp > ZBAR_PARTIAL {
            sym = tmp;
        }
    }

    if test_cfg(dcode_ref.code128.config, ZBAR_CFG_ENABLE) {
        let tmp = _zbar_decode_code128(&mut *dcode);
        if tmp > ZBAR_PARTIAL {
            sym = tmp;
        }
    }

    if test_cfg(
        dcode_ref.databar.config | dcode_ref.databar.config_exp,
        ZBAR_CFG_ENABLE,
    ) {
        let tmp = _zbar_decode_databar(&mut *dcode);
        if tmp > ZBAR_PARTIAL {
            sym = tmp;
        }
    }

    if test_cfg(dcode_ref.codabar.config, ZBAR_CFG_ENABLE) {
        let tmp = _zbar_decode_codabar(&mut *dcode);
        if tmp > ZBAR_PARTIAL {
            sym = tmp;
        }
    }

    if test_cfg(dcode_ref.i25.config, ZBAR_CFG_ENABLE) {
        let tmp = _zbar_decode_i25(&mut *dcode);
        if tmp > ZBAR_PARTIAL {
            sym = tmp;
        }
    }

    dcode_ref.idx = dcode_ref.idx.wrapping_add(1);
    dcode_ref.type_ = sym;

    if sym != 0 {
        if dcode_ref.lock != 0 && sym > ZBAR_PARTIAL && sym != ZBAR_QRCODE {
            dcode_ref._zbar_decoder_release_lock(sym);
        }
        if let Some(handler) = dcode_ref.handler {
            handler(dcode);
        }
    }

    sym
}

// ============================================================================
// Configuration functions
// ============================================================================

/// Get configuration reference for a symbology (internal helper)
fn decoder_get_config(dcode: &zbar_decoder_t, sym: zbar_symbol_type_t) -> Option<&c_uint> {
    match sym {
        ZBAR_EAN13 => Some(&dcode.ean.ean13_config),
        ZBAR_EAN2 => Some(&dcode.ean.ean2_config),
        ZBAR_EAN5 => Some(&dcode.ean.ean5_config),
        ZBAR_EAN8 => Some(&dcode.ean.ean8_config),
        ZBAR_UPCA => Some(&dcode.ean.upca_config),
        ZBAR_UPCE => Some(&dcode.ean.upce_config),
        ZBAR_ISBN10 => Some(&dcode.ean.isbn10_config),
        ZBAR_ISBN13 => Some(&dcode.ean.isbn13_config),
        ZBAR_I25 => Some(&dcode.i25.config),
        ZBAR_DATABAR => Some(&dcode.databar.config),
        ZBAR_DATABAR_EXP => Some(&dcode.databar.config_exp),
        ZBAR_CODABAR => Some(&dcode.codabar.config),
        ZBAR_CODE39 => Some(&dcode.code39.config),
        ZBAR_CODE93 => Some(&dcode.code93.config),
        ZBAR_CODE128 => Some(&dcode.code128.config),
        ZBAR_QRCODE => Some(&dcode.qrf.config),
        ZBAR_SQCODE => Some(&dcode.sqf.config),
        _ => None,
    }
}

/// Get mutable configuration reference for a symbology (internal helper)
fn decoder_get_config_mut(
    dcode: &mut zbar_decoder_t,
    sym: zbar_symbol_type_t,
) -> Option<&mut c_uint> {
    match sym {
        ZBAR_EAN13 => Some(&mut dcode.ean.ean13_config),
        ZBAR_EAN2 => Some(&mut dcode.ean.ean2_config),
        ZBAR_EAN5 => Some(&mut dcode.ean.ean5_config),
        ZBAR_EAN8 => Some(&mut dcode.ean.ean8_config),
        ZBAR_UPCA => Some(&mut dcode.ean.upca_config),
        ZBAR_UPCE => Some(&mut dcode.ean.upce_config),
        ZBAR_ISBN10 => Some(&mut dcode.ean.isbn10_config),
        ZBAR_ISBN13 => Some(&mut dcode.ean.isbn13_config),
        ZBAR_I25 => Some(&mut dcode.i25.config),
        ZBAR_DATABAR => Some(&mut dcode.databar.config),
        ZBAR_DATABAR_EXP => Some(&mut dcode.databar.config_exp),
        ZBAR_CODABAR => Some(&mut dcode.codabar.config),
        ZBAR_CODE39 => Some(&mut dcode.code39.config),
        ZBAR_CODE93 => Some(&mut dcode.code93.config),
        ZBAR_CODE128 => Some(&mut dcode.code128.config),
        ZBAR_QRCODE => Some(&mut dcode.qrf.config),
        ZBAR_SQCODE => Some(&mut dcode.sqf.config),
        _ => None,
    }
}

/// Get all configurations for a symbology
pub fn zbar_decoder_get_configs(dcode: &zbar_decoder_t, sym: zbar_symbol_type_t) -> c_uint {
    decoder_get_config(dcode, sym).copied().unwrap_or(0)
}

/// Set boolean configuration (internal helper)
fn decoder_set_config_bool(
    dcode: &mut zbar_decoder_t,
    sym: zbar_symbol_type_t,
    cfg: c_int,
    val: c_int,
) -> c_int {
    let config = match decoder_get_config_mut(dcode, sym) {
        Some(c) => c,
        None => return 1,
    };

    if cfg >= ZBAR_CFG_NUM {
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
    dcode.ean.enable = if test_cfg(
        dcode.ean.ean13_config
            | dcode.ean.ean2_config
            | dcode.ean.ean5_config
            | dcode.ean.ean8_config
            | dcode.ean.upca_config
            | dcode.ean.upce_config
            | dcode.ean.isbn10_config
            | dcode.ean.isbn13_config,
        ZBAR_CFG_ENABLE,
    ) {
        1
    } else {
        0
    };

    0
}

/// Set integer configuration (internal helper)
fn decoder_set_config_int(
    dcode: &mut zbar_decoder_t,
    sym: zbar_symbol_type_t,
    cfg: c_int,
    val: c_int,
) -> c_int {
    match sym {
        ZBAR_I25 => {
            dcode.i25.configs[(cfg - ZBAR_CFG_MIN_LEN) as usize] = val;
            0
        }
        ZBAR_CODABAR => {
            dcode.codabar.configs[(cfg - ZBAR_CFG_MIN_LEN) as usize] = val;
            0
        }
        ZBAR_CODE39 => {
            dcode.code39.configs[(cfg - ZBAR_CFG_MIN_LEN) as usize] = val;
            0
        }
        ZBAR_CODE93 => {
            dcode.code93.configs[(cfg - ZBAR_CFG_MIN_LEN) as usize] = val;
            0
        }
        ZBAR_CODE128 => {
            dcode.code128.configs[(cfg - ZBAR_CFG_MIN_LEN) as usize] = val;
            0
        }
        _ => 1,
    }
}

/// Get decoder configuration value
pub fn zbar_decoder_get_config(
    dcode: &mut zbar_decoder_t,
    sym: zbar_symbol_type_t,
    cfg: c_int,
) -> Result<c_int, c_int> {
    let config = match decoder_get_config(dcode, sym) {
        Some(c) => c,
        None => return Err(1),
    };

    // Return error if symbol doesn't have config
    if sym <= ZBAR_PARTIAL || sym > ZBAR_CODE128 || sym == ZBAR_COMPOSITE {
        return Err(1);
    }

    // Return decoder boolean configs
    if cfg < ZBAR_CFG_NUM {
        return Ok(if (*config & (1 << cfg)) != 0 { 1 } else { 0 });
    }

    // Return decoder integer configs
    if (ZBAR_CFG_MIN_LEN..=ZBAR_CFG_MAX_LEN).contains(&cfg) {
        Ok(match sym {
            ZBAR_I25 => dcode.i25.configs[(cfg - ZBAR_CFG_MIN_LEN) as usize],
            ZBAR_CODABAR => dcode.codabar.configs[(cfg - ZBAR_CFG_MIN_LEN) as usize],
            ZBAR_CODE39 => dcode.code39.configs[(cfg - ZBAR_CFG_MIN_LEN) as usize],
            ZBAR_CODE93 => dcode.code93.configs[(cfg - ZBAR_CFG_MIN_LEN) as usize],
            ZBAR_CODE128 => dcode.code128.configs[(cfg - ZBAR_CFG_MIN_LEN) as usize],
            _ => return Err(1),
        })
    } else {
        Err(1)
    }
}

/// Set decoder configuration
pub unsafe fn zbar_decoder_set_config(
    dcode: *mut zbar_decoder_t,
    sym: zbar_symbol_type_t,
    cfg: c_int,
    val: c_int,
) -> c_int {
    if dcode.is_null() {
        return 1;
    }

    let dcode = &mut *dcode;

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
            // We need to recursively call this, but we can't reborrow dcode
            // So we need to use the raw pointer API
            unsafe {
                zbar_decoder_set_config(dcode as *mut _, s, cfg, val);
            }
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
