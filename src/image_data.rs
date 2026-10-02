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
    /// Copy the image, optionally inverting every pixel.
    pub(crate) fn copy(&self, inverted: bool) -> Self {
        let data = if inverted {
            self.data.iter().map(|p| !p).collect()
        } else {
            self.data.clone()
        };

        Self {
            width: self.width,
            height: self.height,
            data,
        }
    }

    /// Crops a rectangular region from the image.
    ///
    /// Returns None if the crop region is out of bounds or empty.
    pub(crate) fn crop(&self, x: u32, y: u32, w: u32, h: u32) -> Option<Self> {
        if w == 0 || h == 0 || x.checked_add(w)? > self.width || y.checked_add(h)? > self.height {
            return None;
        }

        let mut data = vec![0u8; (w as usize) * (h as usize)];
        for row in 0..h {
            let src_start = ((y + row) * self.width + x) as usize;
            let dst_start = (row * w) as usize;
            data[dst_start..dst_start + w as usize]
                .copy_from_slice(&self.data[src_start..src_start + w as usize]);
        }

        Some(Self {
            width: w,
            height: h,
            data,
        })
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

    /// Downscales the image by an integer factor with a box filter twice
    /// the size of the step, so each output pixel averages a `2f x 2f`
    /// window and adjacent windows overlap by half.
    ///
    /// Output pixel `(ox, oy)` is centered on source pixel
    /// `(ox * f + f, oy * f + f)`.
    ///
    /// Returns None if `factor` < 2 or either side is shorter than one
    /// window, `2 * factor`.
    pub(crate) fn downscale(&self, factor: u32) -> Option<Self> {
        if factor < 2 {
            return None;
        }
        let f = factor as usize;
        let w = self.width as usize;
        let h = self.height as usize;
        let new_width = (w / f).checked_sub(1)?;
        let new_height = (h / f).checked_sub(1)?;
        if new_width == 0 || new_height == 0 {
            return None;
        }

        let window = 2 * f;
        let area = (window * window) as u64;
        let half = area / 2;

        // Sum each row over a sliding horizontal window first, then sum
        // those over the vertical window: O(f) per output pixel instead of
        // O(f^2).
        let mut row_sums = vec![0u64; h * new_width];
        for y in 0..h {
            let row = &self.data[y * w..(y + 1) * w];
            let sums = &mut row_sums[y * new_width..(y + 1) * new_width];
            for (ox, sum) in sums.iter_mut().enumerate() {
                *sum = row[ox * f..ox * f + window].iter().map(|&p| p as u64).sum();
            }
        }

        let mut data = vec![0u8; new_width * new_height];
        let mut acc = vec![0u64; new_width];
        for oy in 0..new_height {
            acc.fill(0);
            for sy in oy * f..oy * f + window {
                let sums = &row_sums[sy * new_width..(sy + 1) * new_width];
                acc.iter_mut().zip(sums).for_each(|(a, &s)| *a += s);
            }
            let out = &mut data[oy * new_width..(oy + 1) * new_width];
            for (o, a) in out.iter_mut().zip(&acc) {
                *o = ((a + half) / area) as u8;
            }
        }

        Some(Self {
            width: new_width as u32,
            height: new_height as u32,
            data,
        })
    }
}
