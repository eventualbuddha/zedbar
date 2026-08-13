//! Image handling and format support
//!
//! This module provides the [`Image`] type for holding barcode image data.
//! Images must be in grayscale format (8-bit luminance).
//!
//! # Example
//!
//! ```
//! use zedbar::{Image, Scanner};
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

use crate::image_data::ImageData;
use crate::{Error, Result};

/// An image containing barcode data
#[derive(Default)]
pub struct Image {
    image: ImageData,
}

impl Image {
    pub(crate) fn as_mut_image(&mut self) -> &mut ImageData {
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

    /// Crop a rectangular region from the image.
    ///
    /// Returns `None` if the region is out of bounds or empty.
    pub fn crop(&self, x: u32, y: u32, width: u32, height: u32) -> Option<Self> {
        self.image
            .crop(x, y, width, height)
            .map(|image| Image { image })
    }

    /// Upscale the image by an integer factor using bilinear interpolation.
    ///
    /// Returns `None` if `scale` < 2, the image is empty, or the new
    /// dimensions would overflow.
    pub fn upscale(&self, scale: u32) -> Option<Self> {
        self.image.upscale(scale).map(|image| Image { image })
    }
}

#[cfg(feature = "image")]
impl Image {
    /// Create an image from a decoded [`image::DynamicImage`], compositing
    /// any transparency over a white background.
    ///
    /// Prefer this over `img.to_luma8()`: converting an image with an alpha
    /// channel straight to grayscale discards the alpha, so a barcode drawn
    /// in black on a transparent background — the usual shape of one
    /// rendered from SVG — becomes black on black and cannot be found.
    /// Compositing over white first treats transparent areas as the quiet
    /// zone they are meant to represent.
    ///
    /// Images without an alpha channel take the direct conversion, so this
    /// is equivalent to `to_luma8()` for them.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use zedbar::{Image, Scanner};
    ///
    /// let img = image::open("barcode.png").unwrap();
    /// let mut img = Image::from_dynamic(&img).unwrap();
    /// let symbols = Scanner::new().scan(&mut img);
    /// ```
    pub fn from_dynamic(img: &::image::DynamicImage) -> Result<Self> {
        let gray = to_luma_over_white(img);
        Self::from_gray(gray.as_raw(), gray.width(), gray.height())
    }
}

/// Convert to grayscale, compositing any alpha channel over white.
#[cfg(feature = "image")]
fn to_luma_over_white(img: &::image::DynamicImage) -> ::image::GrayImage {
    use ::image::{DynamicImage, Rgb, RgbImage};

    if !img.color().has_alpha() {
        return img.to_luma8();
    }

    let rgba = img.to_rgba8();
    let mut rgb = RgbImage::new(rgba.width(), rgba.height());
    for (dst, src) in rgb.pixels_mut().zip(rgba.pixels()) {
        let [r, g, b, a] = src.0;
        // src-over white: c * a + 255 * (1 - a), in 8-bit fixed point.
        let over_white =
            |c: u8| ((c as u32 * a as u32 + 255 * (255 - a as u32) + 127) / 255).min(255) as u8;
        *dst = Rgb([over_white(r), over_white(g), over_white(b)]);
    }
    // Reuse the image crate's own luminance weights so opaque pixels convert
    // exactly as they would have.
    DynamicImage::ImageRgb8(rgb).to_luma8()
}
