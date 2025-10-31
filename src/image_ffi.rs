//! ZBar Barcode Scanning Library (Rust Port)
//!
//! This crate provides barcode scanning functionality, originally based on the C ZBar library.
//! The conversion to Rust is being done incrementally.

use std::mem::swap;

use libc::{c_int, c_uint};

use crate::img_scanner::zbar_symbol_set_t;

#[derive(Default)]
pub struct zbar_image_t {
    pub format: u32,
    pub width: c_uint,
    pub height: c_uint,
    pub data: Vec<u8>,
    pub refcnt: c_int,
    pub seq: c_uint,
    syms: Option<Box<zbar_symbol_set_t>>,
}

impl zbar_image_t {
    #[inline]
    pub(crate) fn set_syms(&mut self, syms: Box<zbar_symbol_set_t>) {
        self.syms = Some(syms);
    }

    #[inline]
    pub(crate) fn syms(&self) -> Option<&zbar_symbol_set_t> {
        self.syms.as_deref()
    }

    pub(crate) fn swap_symbols_with(&mut self, other: &mut Self) {
        swap(&mut self.syms, &mut other.syms);
    }
}

pub(crate) fn _zbar_image_copy(src: &zbar_image_t, inverted: c_int) -> Option<zbar_image_t> {
    const FOURCC_Y800: u32 = 0x30303859; // fourcc('Y', '8', '0', '0')
    const FOURCC_GREY: u32 = 0x59455247; // fourcc('G', 'R', 'E', 'Y')

    if inverted != 0 && src.format != FOURCC_Y800 && src.format != FOURCC_GREY {
        return None;
    }

    let mut dst = zbar_image_t {
        format: src.format,
        width: src.width,
        height: src.height,
        data: vec![0; src.data.len()],
        ..Default::default()
    };

    if inverted == 0 {
        dst.data.copy_from_slice(&src.data);
    } else {
        for (dp, sp) in dst.data.iter_mut().zip(src.data.iter()) {
            *dp = !(*sp);
        }
    }
    Some(dst)
}
