//! Low-level barcode decoder

use crate::{
    decoder_types::{
        codabar_decoder_t, code128_decoder_t, code39_decoder_t, code93_decoder_t,
        databar_decoder_t, databar_segment_t, ean_decoder_t, i25_decoder_t, qr_finder_t,
        zbar_decoder_t, zbar_symbol_type_t, BUFFER_INCR, BUFFER_MAX, BUFFER_MIN, DECODE_WINDOW,
        ZBAR_CFG_EMIT_CHECK, ZBAR_CFG_ENABLE, ZBAR_CFG_MAX_LEN, ZBAR_CFG_MIN_LEN, ZBAR_CODABAR,
        ZBAR_CODE128, ZBAR_CODE39, ZBAR_CODE93, ZBAR_COMPOSITE, ZBAR_DATABAR, ZBAR_DATABAR_EXP,
        ZBAR_EAN13, ZBAR_EAN2, ZBAR_EAN5, ZBAR_EAN8, ZBAR_I25, ZBAR_ISBN10, ZBAR_ISBN13, ZBAR_NONE,
        ZBAR_PARTIAL, ZBAR_QRCODE, ZBAR_SQCODE, ZBAR_UPCA, ZBAR_UPCE,
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

// Macro equivalents
#[inline]
unsafe fn cfg_set(configs: &mut [c_int; 2], cfg: c_int, val: c_int) {
    configs[(cfg - ZBAR_CFG_MIN_LEN) as usize] = val;
}

#[inline]
fn test_cfg(config: c_uint, cfg: c_int) -> bool {
    ((config >> cfg) & 1) != 0
}

#[inline]
pub(crate) unsafe fn decoder_alloc_zeroed() -> *mut zbar_decoder_t {
    let decoder = Box::new(zbar_decoder_t::default());
    Box::into_raw(decoder)
}

#[inline]
pub(crate) unsafe fn decoder_free_struct(dcode: *mut zbar_decoder_t) {
    drop(Box::from_raw(dcode))
}

#[inline]
pub(crate) unsafe fn decoder_alloc_buffer(size: usize) -> *mut c_char {
    libc::malloc(size) as *mut c_char
}

#[inline]
pub(crate) unsafe fn decoder_free_buffer(buf: *mut c_char) {
    libc::free(buf as *mut c_void);
}

#[inline]
pub(crate) unsafe fn decoder_realloc_buffer(buf: *mut c_char, new_len: usize) -> *mut c_char {
    libc::realloc(buf as *mut c_void, new_len) as *mut c_char
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

/// Get the color (bar/space) of the current element
pub unsafe fn _zbar_decoder_get_color(dcode: *const zbar_decoder_t) -> c_char {
    ((*dcode).idx & 1) as c_char
}

/// Get width of a specific element from the decoder's history window
pub unsafe fn _zbar_decoder_get_width(dcode: *const zbar_decoder_t, offset: u8) -> c_uint {
    (*dcode).w[(((*dcode).idx as usize).wrapping_sub(offset as usize)) & (DECODE_WINDOW - 1)]
}

/// Get the combined width of two consecutive elements
pub unsafe fn _zbar_decoder_pair_width(dcode: *const zbar_decoder_t, offset: u8) -> c_uint {
    _zbar_decoder_get_width(dcode, offset) + _zbar_decoder_get_width(dcode, offset + 1)
}

/// Calculate sum of n consecutive element widths
pub unsafe fn _zbar_decoder_calc_s(
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
pub unsafe fn _zbar_decoder_decode_e(e: c_uint, s: c_uint, n: c_uint) -> c_int {
    let big_e = ((e * n * 2 + 1) / s).wrapping_sub(3) / 2;
    if big_e >= n - 3 {
        -1
    } else {
        big_e as c_int
    }
}

/// Sort 3 elements by width and return their indices
pub unsafe fn _zbar_decoder_decode_sort3(dcode: *mut zbar_decoder_t, i0: c_int) -> c_uint {
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
pub unsafe fn _zbar_decoder_decode_sortn(
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
pub unsafe fn _zbar_decoder_acquire_lock(
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
pub unsafe fn _zbar_decoder_release_lock(
    dcode: *mut zbar_decoder_t,
    req: zbar_symbol_type_t,
) -> c_char {
    debug_assert_eq!((*dcode).lock, req, "lock={} req={}", (*dcode).lock, req);
    (*dcode).lock = 0;
    0
}

/// Resize the decoder's data buffer if needed
/// Returns 1 on allocation failure or if max size exceeded, 0 on success
pub unsafe fn _zbar_decoder_size_buf(dcode: *mut zbar_decoder_t, len: c_uint) -> c_char {
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

    let buf = decoder_realloc_buffer((*dcode).buf, new_len as usize);
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
pub unsafe fn _zbar_codabar_reset(codabar: *mut codabar_decoder_t) {
    (*codabar).set_direction(false);
    (*codabar).set_element(0);
    (*codabar).set_character(-1);
    (*codabar).s7 = 0;
}

/// Reset code128 decoder state
pub unsafe fn _zbar_code128_reset(dcode128: *mut code128_decoder_t) {
    (*dcode128).set_direction(0);
    (*dcode128).set_element(0);
    (*dcode128).set_character(-1);
    (*dcode128).s6 = 0;
}

/// Reset code39 decoder state
pub unsafe fn _zbar_code39_reset(dcode39: *mut code39_decoder_t) {
    (*dcode39).set_direction(false);
    (*dcode39).set_element(0);
    (*dcode39).set_character(-1);
    (*dcode39).s9 = 0;
}

/// Reset code93 decoder state
pub unsafe fn _zbar_code93_reset(dcode93: *mut code93_decoder_t) {
    (*dcode93).set_direction(false);
    (*dcode93).set_element(0);
    (*dcode93).set_character(-1);
}

/// Reset i25 decoder state
pub unsafe fn _zbar_i25_reset(i25: *mut i25_decoder_t) {
    (*i25).set_direction(false);
    (*i25).set_element(0);
    (*i25).set_character(-1);
    (*i25).s10 = 0;
}

/// Reset QR finder state
pub unsafe fn _zbar_qr_finder_reset(qrf: *mut qr_finder_t) {
    (*qrf).s5 = 0;
}

/// Prepare DataBar decoder for new scan
pub unsafe fn _zbar_databar_new_scan(db: *mut databar_decoder_t) {
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
pub unsafe fn _zbar_databar_reset(db: *mut databar_decoder_t) {
    let n = (*db).csegs() as isize;
    _zbar_databar_new_scan(db);
    for i in 0..n {
        let seg = ((*db).segs).offset(i);
        (*seg).set_finder(-1);
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
pub unsafe fn zbar_decoder_create() -> *mut zbar_decoder_t {
    let dcode = decoder_alloc_zeroed();
    if dcode.is_null() {
        return std::ptr::null_mut();
    }

    (*dcode).buf_alloc = BUFFER_MIN;
    (*dcode).buf = decoder_alloc_buffer((*dcode).buf_alloc as usize);
    if (*dcode).buf.is_null() {
        decoder_free_struct(dcode);
        return std::ptr::null_mut();
    }

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
    (*dcode).databar.segs = decoder_alloc_databar_segments(4);
    if (*dcode).databar.segs.is_null() {
        decoder_free_buffer((*dcode).buf);
        decoder_free_struct(dcode);
        return std::ptr::null_mut();
    }

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
pub unsafe fn zbar_decoder_destroy(dcode: *mut zbar_decoder_t) {
    if !(*dcode).databar.segs.is_null() {
        decoder_free_databar_segments((*dcode).databar.segs);
    }
    if !(*dcode).buf.is_null() {
        decoder_free_buffer((*dcode).buf);
    }
    decoder_free_struct(dcode);
}

/// Reset decoder to initial state
pub unsafe fn zbar_decoder_reset(dcode: *mut zbar_decoder_t) {
    (*dcode).idx = 0;
    (*dcode).w.fill(0);
    (*dcode).type_ = ZBAR_NONE;
    (*dcode).lock = ZBAR_NONE;
    (*dcode).modifiers = 0;
    (*dcode).direction = 0;
    (*dcode).s6 = 0;
    (*dcode).buflen = 0;

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
    (*dcode).w.fill(0);
    (*dcode).lock = ZBAR_NONE;
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
pub unsafe fn zbar_decoder_get_color(dcode: *const zbar_decoder_t) -> c_int {
    _zbar_decoder_get_color(dcode) as c_int
}

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

// ============================================================================
// Debug helper
// ============================================================================

use std::sync::Mutex;

static DECODER_DUMP: Mutex<Option<Vec<u8>>> = Mutex::new(None);

/// Format decoder buffer as hex string (for debugging)
pub unsafe fn _zbar_decoder_buf_dump(buf: *mut u8, buflen: c_uint) -> *const c_char {
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
