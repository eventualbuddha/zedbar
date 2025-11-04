//! Image handling and format support

use crate::image_ffi::zbar_image_t;
use crate::symbol::SymbolSet;
use crate::{Error, Result};

/// An image containing barcode data
#[derive(Default)]
pub struct Image {
    image: zbar_image_t,
}

impl Image {
    pub(crate) fn as_mut_image(&mut self) -> &mut zbar_image_t {
        &mut self.image
    }

    /// Create an image from grayscale data
    pub fn from_gray(data: &[u8], width: u32, height: u32) -> Result<Self> {
        if (data.len() as u64) != (width as u64) * (height as u64) {
            return Err(Error::Invalid);
        }

        let mut image = Self::default();
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
