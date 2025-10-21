//! ZBar Barcode Scanning Library (Rust Port)
//!
//! This crate provides barcode scanning functionality, originally based on the C ZBar library.
//! The conversion to Rust is being done incrementally.

use std::ptr::{null, null_mut, NonNull};

use libc::{c_int, c_uint};

use crate::{
    ffi::refcnt,
    img_scanner::zbar_symbol_set_t,
    symbol::{zbar_symbol_set_ref, zbar_symbol_t},
};

#[inline]
unsafe fn image_alloc_zeroed() -> *mut zbar_image_t {
    let image = Box::new(zbar_image_t::default());
    Box::into_raw(image)
}

#[inline]
unsafe fn image_free_struct(img: *mut zbar_image_t) {
    drop(Box::from_raw(img))
}

#[derive(Default)]
#[allow(non_camel_case_types)]
pub struct zbar_image_t {
    pub format: u32,
    pub width: c_uint,
    pub height: c_uint,
    pub data: Vec<u8>,
    pub refcnt: c_int,
    pub seq: c_uint,
    pub syms: Option<NonNull<zbar_symbol_set_t>>,
}

impl zbar_image_t {
    #[inline]
    pub fn syms_ptr(&self) -> *mut zbar_symbol_set_t {
        self.syms.map_or(null_mut(), NonNull::as_ptr)
    }

    #[inline]
    pub fn set_syms_ptr(&mut self, ptr: *mut zbar_symbol_set_t) {
        self.syms = NonNull::new(ptr);
    }

    #[inline]
    pub fn take_syms_ptr(&mut self) -> *mut zbar_symbol_set_t {
        self.syms.take().map_or(null_mut(), NonNull::as_ptr)
    }
}

pub unsafe fn zbar_image_create() -> *mut zbar_image_t {
    let img = image_alloc_zeroed();
    if img.is_null() {
        return null_mut();
    }
    let img_ref = &mut *img;
    refcnt(&mut img_ref.refcnt, 1);
    img
}

pub unsafe fn _zbar_image_free(img: *mut zbar_image_t) {
    let image = &mut *img;
    if let Some(syms) = image.syms.take() {
        zbar_symbol_set_ref(syms.as_ptr(), -1);
    }
    image_free_struct(img);
}

pub unsafe fn zbar_image_destroy(img: *mut zbar_image_t) {
    let img_ref = &mut *img;
    if refcnt(&mut img_ref.refcnt, -1) == 0 {
        _zbar_image_free(img);
    }
}

pub unsafe fn zbar_image_first_symbol(img: *const zbar_image_t) -> *const zbar_symbol_t {
    let img = &*img;
    img.syms
        .map_or(null(), |syms| unsafe { (*syms.as_ptr()).head })
}

pub unsafe fn _zbar_image_swap_symbols(a: *mut zbar_image_t, b: *mut zbar_image_t) {
    let a = &mut *a;
    let b = &mut *b;
    std::mem::swap(&mut a.syms, &mut b.syms);
}

pub unsafe fn _zbar_image_copy(src: *const zbar_image_t, inverted: c_int) -> *mut zbar_image_t {
    const FOURCC_Y800: u32 = 0x30303859; // fourcc('Y', '8', '0', '0')
    const FOURCC_GREY: u32 = 0x59455247; // fourcc('G', 'R', 'E', 'Y')

    let src = &*src;

    if inverted != 0 && src.format != FOURCC_Y800 && src.format != FOURCC_GREY {
        return null_mut();
    }

    let dst = zbar_image_create();
    let dst_ref = &mut *dst;
    dst_ref.format = src.format;
    dst_ref.width = src.width;
    dst_ref.height = src.height;
    dst_ref.data = vec![0; src.data.len()];

    if inverted == 0 {
        dst_ref.data.copy_from_slice(&src.data);
    } else {
        for (dp, sp) in dst_ref.data.iter_mut().zip(src.data.iter()) {
            *dp = !(*sp);
        }
    }
    dst
}
