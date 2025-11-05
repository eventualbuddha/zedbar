//! Image handling and format support
//!
//! This module provides the [`Image`] type for holding barcode image data.
//! Images must be in grayscale format (8-bit luminance).
//!
//! # Example
//!
//! ```
//! use zbar::{Image, Scanner};
//!
//! // Create image from grayscale data
//! let width = 640;
//! let height = 480;
//! let data = vec![0u8; (width * height) as usize];
//! let mut image = Image::from_gray(&data, width, height).unwrap();
//!
//! // Scan the image
//! let mut scanner = Scanner::new();
//! let symbols = scanner.scan(&mut image);
//! ```

use crate::image_ffi::zbar_image_t;
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
}
