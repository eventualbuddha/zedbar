//! ZBar Barcode Scanning Library (Rust Port)
//!
//! This crate provides barcode scanning functionality, originally based on the C ZBar library.
//! The conversion to Rust is being done incrementally.

use std::{
    mem::swap,
    ptr::{null, NonNull},
};

use libc::{c_int, c_uint};

use crate::{
    img_scanner::zbar_symbol_set_t,
    symbol::{zbar_symbol_set_ref, zbar_symbol_t},
};

#[derive(Default)]
#[allow(non_camel_case_types)]
pub struct zbar_image_t {
    pub format: u32,
    pub width: c_uint,
    pub height: c_uint,
    pub data: Vec<u8>,
    pub refcnt: c_int,
    pub seq: c_uint,
    syms: Option<NonNull<zbar_symbol_set_t>>,
}

impl zbar_image_t {
    #[inline]
    pub(crate) fn set_syms_ptr(&mut self, ptr: *mut zbar_symbol_set_t) {
        self.syms = NonNull::new(ptr);
    }

    pub(crate) fn swap_symbols_with(&mut self, other: &mut Self) {
        swap(&mut self.syms, &mut other.syms);
    }
}

impl Drop for zbar_image_t {
    fn drop(&mut self) {
        unsafe {
            if let Some(syms) = self.syms.take() {
                zbar_symbol_set_ref(syms.as_ptr(), -1);
            }
        }
    }
}

pub(crate) unsafe fn zbar_image_first_symbol(img: &zbar_image_t) -> *const zbar_symbol_t {
    img.syms
        .map_or(null(), |syms| unsafe { (*syms.as_ptr()).head })
}

pub(crate) fn _zbar_image_copy(src: &zbar_image_t, inverted: c_int) -> Option<zbar_image_t> {
    const FOURCC_Y800: u32 = 0x30303859; // fourcc('Y', '8', '0', '0')
    const FOURCC_GREY: u32 = 0x59455247; // fourcc('G', 'R', 'E', 'Y')

    if inverted != 0 && src.format != FOURCC_Y800 && src.format != FOURCC_GREY {
        return None;
    }

    let mut dst = zbar_image_t::default();
    dst.format = src.format;
    dst.width = src.width;
    dst.height = src.height;
    dst.data = vec![0; src.data.len()];

    if inverted == 0 {
        dst.data.copy_from_slice(&src.data);
    } else {
        for (dp, sp) in dst.data.iter_mut().zip(src.data.iter()) {
            *dp = !(*sp);
        }
    }
    Some(dst)
}
