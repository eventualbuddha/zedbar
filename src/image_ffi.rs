//! ZBar Barcode Scanning Library (Rust Port)
//!
//! This crate provides barcode scanning functionality, originally based on the C ZBar library.
//! The conversion to Rust is being done incrementally.

use std::{
    ffi::c_void,
    mem::{size_of, transmute},
    ptr::{null, null_mut},
};

use libc::{c_char, c_int, c_uint, c_ulong};

use crate::{ffi::refcnt, img_scanner::zbar_symbol_set_t, symbol::zbar_symbol_set_ref};

#[inline]
unsafe fn image_alloc_zeroed() -> *mut zbar_image_t {
    libc::calloc(1, size_of::<zbar_image_t>()) as *mut zbar_image_t
}

#[inline]
unsafe fn image_free_struct(img: *mut zbar_image_t) {
    libc::free(img as *mut c_void);
}

#[inline]
unsafe fn image_alloc_data(size: usize) -> *mut c_void {
    libc::malloc(size)
}

#[inline]
unsafe fn image_free_data_ptr(ptr: *mut c_void) {
    libc::free(ptr);
}

#[allow(non_camel_case_types)]
pub struct zbar_image_t {
    pub format: u32,
    pub width: c_uint,
    pub height: c_uint,
    pub data: *mut c_void,
    pub datalen: c_ulong,
    pub userdata: *mut c_void,
    pub cleanup: *mut c_void,
    pub refcnt: c_int,
    pub srcidx: c_int,
    pub next: *mut zbar_image_t,
    pub seq: c_uint,
    pub syms: *mut zbar_symbol_set_t,
}

#[allow(non_camel_case_types)]
pub struct zbar_symbol_t {
    pub symbol_type: c_int,
    pub configs: c_uint,
    pub modifiers: c_uint,
    pub data_alloc: c_uint,
    pub datalen: c_uint,
    pub data: *mut c_char,
    pub pts_alloc: c_uint,
    pub npts: c_uint,
    pub pts: *mut c_void,
    pub orient: c_int,
    pub refcnt: c_int,
    pub next: *mut zbar_symbol_t,
    pub syms: *mut zbar_symbol_set_t,
    pub time: c_ulong,
    pub cache_count: c_int,
    pub quality: c_int,
}

pub unsafe fn zbar_image_create() -> *mut zbar_image_t {
    let img = image_alloc_zeroed();
    if img.is_null() {
        return null_mut();
    }
    refcnt(&mut (*img).refcnt, 1);
    (*img).srcidx = -1;
    img
}

pub unsafe fn _zbar_image_free(img: *mut zbar_image_t) {
    if !(*img).syms.is_null() {
        zbar_symbol_set_ref((*img).syms, -1);
        (*img).syms = null_mut();
    }
    image_free_struct(img);
}

pub unsafe fn zbar_image_destroy(img: *mut zbar_image_t) {
    if refcnt(&mut (*img).refcnt, -1) != 0 {
        if !(*img).cleanup.is_null() {
            let cleanup = transmute::<*mut c_void, fn(*mut zbar_image_t)>((*img).cleanup);
            cleanup(img);
        }
        _zbar_image_free(img);
    }
}

pub unsafe fn zbar_image_free_data(img: *mut zbar_image_t) {
    if img.is_null() {
        return;
    }
    if !(*img).cleanup.is_null() && !(*img).data.is_null() {
        if (*img).cleanup != zbar_image_free_data as *mut c_void {
            /* using function address to detect this case is a bad idea;
             * windows link libraries add an extra layer of indirection...
             * this works around that problem (bug #2796277)
             */
            let cleanup = transmute::<*mut c_void, fn(*mut zbar_image_t)>((*img).cleanup);
            (*img).cleanup = zbar_image_free_data as *mut c_void;
            cleanup(img);
        } else {
            image_free_data_ptr((*img).data);
        }
    }
    (*img).data = null_mut();
}

pub unsafe fn zbar_image_first_symbol(img: *const zbar_image_t) -> *const zbar_symbol_t {
    if (*img).syms.is_null() {
        null()
    } else {
        (*(*img).syms).head
    }
}

pub unsafe fn _zbar_image_swap_symbols(a: *mut zbar_image_t, b: *mut zbar_image_t) {
    std::mem::swap(&mut (*a).syms, &mut (*b).syms);
}

pub unsafe fn _zbar_image_copy(src: *const zbar_image_t, inverted: c_int) -> *mut zbar_image_t {
    const FOURCC_Y800: u32 = 0x30303859; // fourcc('Y', '8', '0', '0')
    const FOURCC_GREY: u32 = 0x59455247; // fourcc('G', 'R', 'E', 'Y')

    if inverted != 0 && (*src).format != FOURCC_Y800 && (*src).format != FOURCC_GREY {
        return null_mut();
    }

    let dst = zbar_image_create();
    (*dst).format = (*src).format;
    (*dst).width = (*src).width;
    (*dst).height = (*src).height;
    (*dst).datalen = (*src).datalen;
    (*dst).data = image_alloc_data((*src).datalen as usize);
    assert!(!(*dst).data.is_null());

    if inverted == 0 {
        libc::memcpy((*dst).data, (*src).data, (*src).datalen as usize);
    } else {
        let len = (*src).datalen as usize;
        let mut sp = (*src).data as *const usize;
        let mut dp = (*dst).data as *mut usize;

        // Do it word per word to speed up
        let word_len = len / std::mem::size_of::<usize>();
        for _ in 0..word_len {
            *dp = !(*sp);
            sp = sp.add(1);
            dp = dp.add(1);
        }

        // Deal with non-aligned remains
        let mut spc = sp as *const u8;
        let mut dpc = dp as *mut u8;
        let remain = len % std::mem::size_of::<usize>();
        for _ in 0..remain {
            *dpc = !(*spc);
            spc = spc.add(1);
            dpc = dpc.add(1);
        }
    }
    (*dst).cleanup = zbar_image_free_data as *mut c_void;
    dst
}
