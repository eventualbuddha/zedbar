//! Image Module
//!
//! This module provides image handling and barcode scanning functionality.

#[derive(Default)]
pub struct ImageData {
    pub width: u32,
    pub height: u32,
    pub data: Vec<u8>,
}

impl ImageData {
    pub(crate) fn copy(&self, inverted: bool) -> Option<Self> {
        let mut dst = Self {
            width: self.width,
            height: self.height,
            data: vec![0; self.data.len()],
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

    /// Upscales the image using bilinear interpolation.
    ///
    /// Returns None if the image is empty, scale factor is invalid, or
    /// dimensions would overflow.
    pub(crate) fn upscale(&self, scale: u32) -> Option<Self> {
        if scale < 2 || self.width == 0 || self.height == 0 {
            return None;
        }

        // Use checked arithmetic to prevent overflow
        let new_width = self.width.checked_mul(scale)?;
        let new_height = self.height.checked_mul(scale)?;
        let total_pixels = new_width.checked_mul(new_height)?;
        let mut data = vec![0u8; total_pixels as usize];

        let w = self.width as usize;
        let h = self.height as usize;
        let nw = new_width as usize;

        for ny in 0..new_height as usize {
            for nx in 0..nw {
                // Map back to source coordinates with sub-pixel precision
                // We use (ny + 0.5) / scale - 0.5 to center the mapping
                let sy_f = (ny as f32 + 0.5) / scale as f32 - 0.5;
                let sx_f = (nx as f32 + 0.5) / scale as f32 - 0.5;

                let sy0 = sy_f.floor().max(0.0) as usize;
                let sx0 = sx_f.floor().max(0.0) as usize;
                let sy1 = (sy0 + 1).min(h - 1);
                let sx1 = (sx0 + 1).min(w - 1);

                // Clamp interpolation weights to [0, 1] to avoid artifacts at borders
                let fy = (sy_f - sy0 as f32).clamp(0.0, 1.0);
                let fx = (sx_f - sx0 as f32).clamp(0.0, 1.0);

                // Bilinear interpolation
                let p00 = self.data[sy0 * w + sx0] as f32;
                let p10 = self.data[sy0 * w + sx1] as f32;
                let p01 = self.data[sy1 * w + sx0] as f32;
                let p11 = self.data[sy1 * w + sx1] as f32;

                let value = p00 * (1.0 - fx) * (1.0 - fy)
                    + p10 * fx * (1.0 - fy)
                    + p01 * (1.0 - fx) * fy
                    + p11 * fx * fy;

                data[ny * nw + nx] = value.round().clamp(0.0, 255.0) as u8;
            }
        }

        Some(Self {
            width: new_width,
            height: new_height,
            data,
        })
    }
}
