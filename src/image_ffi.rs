//! ZBar Barcode Scanning Library (Rust Port)
//!
//! This crate provides barcode scanning functionality, originally based on the C ZBar library.
//! The conversion to Rust is being done incrementally.

use std::{
    ffi::c_void,
    mem::transmute,
    ptr::{null, null_mut},
};

use libc::c_int;

use crate::{refcnt, symbol::zbar_symbol_set_ref, zbar_image_t, zbar_symbol_t};

pub unsafe fn zbar_image_create() -> *mut zbar_image_t {
    let img = libc::calloc(1, std::mem::size_of::<zbar_image_t>()) as *mut zbar_image_t;
    refcnt(&mut (*img).refcnt, 1);
    (*img).srcidx = -1;
    img
}

pub unsafe fn _zbar_image_free(img: *mut zbar_image_t) {
    if !(*img).syms.is_null() {
        zbar_symbol_set_ref((*img).syms, -1);
        (*img).syms = null_mut();
    }
    libc::free(img as *mut c_void);
}

#[no_mangle]
pub unsafe extern "C" fn zbar_image_destroy(img: *mut zbar_image_t) {
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
            libc::free((*img).data);
        }
    }
    (*img).data = null_mut();
}

#[no_mangle]
pub unsafe extern "C" fn zbar_image_first_symbol(img: *const zbar_image_t) -> *const zbar_symbol_t {
    if (*img).syms.is_null() {
        null()
    } else {
        (*(*img).syms).head
    }
}

#[no_mangle]
pub unsafe extern "C" fn _zbar_image_swap_symbols(a: *mut zbar_image_t, b: *mut zbar_image_t) {
    std::mem::swap(&mut (*a).syms, &mut (*b).syms);
}

#[no_mangle]
pub unsafe extern "C" fn _zbar_image_copy_size(dst: *mut zbar_image_t, src: *const zbar_image_t) {
    (*dst).width = (*src).width;
    (*dst).height = (*src).height;
}

#[no_mangle]
pub unsafe extern "C" fn _zbar_image_copy(
    src: *const zbar_image_t,
    inverted: c_int,
) -> *mut zbar_image_t {
    const FOURCC_Y800: u32 = 0x30303859; // fourcc('Y', '8', '0', '0')
    const FOURCC_GREY: u32 = 0x59455247; // fourcc('G', 'R', 'E', 'Y')

    if inverted != 0 && (*src).format != FOURCC_Y800 && (*src).format != FOURCC_GREY {
        return null_mut();
    }

    let dst = zbar_image_create();
    (*dst).format = (*src).format;
    _zbar_image_copy_size(dst, src);
    (*dst).datalen = (*src).datalen;
    (*dst).data = libc::malloc((*src).datalen as usize);
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
