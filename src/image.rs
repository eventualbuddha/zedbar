//! Image handling and format support

use crate::symbol::{Symbol, SymbolSet};
use crate::{ffi, zbar_image_t};
use crate::{Error, Result};
use std::ptr;

/// Image formats supported by ZBar
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ImageFormat {
    /// 8-bit grayscale
    Y800,
    /// Same as Y800
    Gray,
}

impl ImageFormat {
    fn to_fourcc(self) -> u32 {
        match self {
            ImageFormat::Y800 | ImageFormat::Gray => fourcc(b'Y', b'8', b'0', b'0'),
        }
    }
}

/// Create a fourcc code from 4 bytes
fn fourcc(a: u8, b: u8, c: u8, d: u8) -> u32 {
    (a as u32) | ((b as u32) << 8) | ((c as u32) << 16) | ((d as u32) << 24)
}

/// An image containing barcode data
pub struct Image {
    ptr: *mut zbar_image_t,
}

impl Image {
    /// Create a new empty image
    pub fn new() -> Self {
        let ptr = unsafe { ffi::zbar_image_create() };
        Image { ptr }
    }

    /// Create an image from grayscale data
    pub fn from_gray(data: &[u8], width: u32, height: u32) -> Result<Self> {
        if (data.len() as u64) != (width as u64) * (height as u64) {
            return Err(Error::Invalid);
        }

        let image = Self::new();
        unsafe {
            ffi::zbar_image_set_format(image.ptr, ImageFormat::Gray.to_fourcc() as u64);
            ffi::zbar_image_set_size(image.ptr, width, height);
            ffi::zbar_image_set_data(
                image.ptr,
                data.as_ptr() as *const std::ffi::c_void,
                data.len() as u64,
                ptr::null(),
            );
        }
        Ok(image)
    }

    /// Get the image width
    pub fn width(&self) -> u32 {
        unsafe { ffi::zbar_image_get_width(self.ptr) }
    }

    /// Get the image height  
    pub fn height(&self) -> u32 {
        unsafe { ffi::zbar_image_get_height(self.ptr) }
    }

    /// Get the image format
    pub fn format(&self) -> u32 {
        unsafe { ffi::zbar_image_get_format(self.ptr) as u32 }
    }

    /// Get access to the raw image data
    pub fn data(&self) -> &[u8] {
        unsafe {
            let data_ptr = ffi::zbar_image_get_data(self.ptr);
            if data_ptr.is_null() {
                &[]
            } else {
                let len = (self.width() * self.height()) as usize;
                std::slice::from_raw_parts(data_ptr as *const u8, len)
            }
        }
    }

    /// Get the symbols found in this image (if it has been scanned)
    pub fn symbols(&self) -> SymbolSet {
        let first_symbol_ptr = unsafe { ffi::zbar_image_first_symbol(self.ptr) };
        let first_symbol = unsafe { Symbol::from_ptr(first_symbol_ptr) };
        SymbolSet::new(first_symbol)
    }

    /// Get a raw pointer to the underlying C image object
    pub(crate) fn as_ptr(&self) -> *mut zbar_image_t {
        self.ptr
    }
}

impl Drop for Image {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe {
                ffi::zbar_image_destroy(self.ptr);
            }
        }
    }
}

impl Default for Image {
    fn default() -> Self {
        Self::new()
    }
}
