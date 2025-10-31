//! Image handling and format support

use crate::image_ffi::zbar_image_t;
use crate::symbol::SymbolSet;
use crate::{Error, Result};

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
    image: zbar_image_t,
}

impl Image {
    /// Create a new empty image
    pub fn new() -> Self {
        Image {
            image: zbar_image_t::default(),
        }
    }

    pub(crate) fn as_mut_image(&mut self) -> &mut zbar_image_t {
        &mut self.image
    }

    /// Create an image from grayscale data
    pub fn from_gray(data: &[u8], width: u32, height: u32) -> Result<Self> {
        if (data.len() as u64) != (width as u64) * (height as u64) {
            return Err(Error::Invalid);
        }

        let mut image = Self::new();
        image.image.format = ImageFormat::Gray.to_fourcc();
        image.image.width = width;
        image.image.height = height;
        image.image.data.extend_from_slice(data);
        Ok(image)
    }

    /// Get the image width
    pub fn width(&self) -> u32 {
        self.image.width
    }

    /// Get the image height  
    pub fn height(&self) -> u32 {
        self.image.height
    }

    /// Get the image format
    pub fn format(&self) -> u32 {
        self.image.format
    }

    /// Get access to the raw image data
    pub fn data(&self) -> &[u8] {
        &self.image.data
    }

    /// Get the symbols found in this image (if it has been scanned)
    pub fn symbols(&self) -> SymbolSet<'_> {
        if let Some(syms) = self.image.syms() {
            SymbolSet::from_slice(&syms.symbols)
        } else {
            SymbolSet::from_slice(&[])
        }
    }
}

impl Default for Image {
    fn default() -> Self {
        Self::new()
    }
}
