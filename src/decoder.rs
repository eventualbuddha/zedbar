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
    finder::_zbar_find_qr,
};
use libc::{c_char, c_int, c_uint, c_void};
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

/// Get decoded data buffer
pub unsafe fn zbar_decoder_get_data(dcode: *const zbar_decoder_t) -> *const c_char {
    (*dcode).buf as *const c_char
}

/// Get length of decoded data
pub unsafe fn zbar_decoder_get_data_length(dcode: *const zbar_decoder_t) -> c_uint {
    (*dcode).buflen
}

/// Get decode direction
pub unsafe fn zbar_decoder_get_direction(dcode: *const zbar_decoder_t) -> c_int {
    (*dcode).direction
}

/// Set decoder callback handler
pub unsafe fn zbar_decoder_set_handler(
    dcode: *mut zbar_decoder_t,
    handler: Option<crate::decoder_types::zbar_decoder_handler_t>,
) -> Option<crate::decoder_types::zbar_decoder_handler_t> {
    let result = (*dcode).handler;
    (*dcode).handler = handler;
    result
}

/// Set user data pointer
pub unsafe fn zbar_decoder_set_userdata(dcode: *mut zbar_decoder_t, userdata: *mut c_void) {
    (*dcode).userdata = userdata;
}

/// Get user data pointer
pub unsafe fn zbar_decoder_get_userdata(dcode: *const zbar_decoder_t) -> *mut c_void {
    (*dcode).userdata
}

/// Get decoded symbol type
pub unsafe fn zbar_decoder_get_type(dcode: *const zbar_decoder_t) -> zbar_symbol_type_t {
    (*dcode).type_
}

/// Get decoded symbol modifiers
pub unsafe fn zbar_decoder_get_modifiers(dcode: *const zbar_decoder_t) -> c_uint {
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
    (*dcode).s6 = (*dcode).s6.wrapping_sub((*dcode).get_width(7));
    (*dcode).s6 = (*dcode).s6.wrapping_add((*dcode).get_width(1));

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
            (*dcode)._zbar_decoder_release_lock(sym);
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
pub unsafe fn zbar_decoder_get_configs(
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
pub unsafe fn zbar_decoder_get_config(
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
pub unsafe fn zbar_decoder_set_config(
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
