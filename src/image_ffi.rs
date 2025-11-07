//! ZBar Image Module
//!
//! This module provides image handling and barcode scanning functionality.

#[derive(Default)]
pub struct zbar_image_t {
    pub width: u32,
    pub height: u32,
    pub data: Vec<u8>,
}

impl zbar_image_t {
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
}
