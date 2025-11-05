//! ZBar Barcode Scanning Library (Rust Port)
//!
//! This crate provides barcode scanning functionality, originally based on the C ZBar library.
//! The conversion to Rust is being done incrementally.

use std::mem::swap;


use crate::img_scanner::zbar_symbol_set_t;

#[derive(Default)]
pub struct zbar_image_t {
    pub width: u32,
    pub height: u32,
    pub data: Vec<u8>,
    pub refcnt: i32,
    pub seq: u32,
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

    pub(crate) fn copy(&self, inverted: bool) -> Option<Self> {
        let mut dst = Self {
            width: self.width,
            height: self.height,
            data: vec![0; self.data.len()],
            ..Default::default()
        };

        if !inverted {
            dst.data.copy_from_slice(&self.data);
        } else {
            for (dp, sp) in dst.data.iter_mut().zip(self.data.iter()) {
                *dp = !(*sp);
            }
        }
        Some(dst)
    }
}
