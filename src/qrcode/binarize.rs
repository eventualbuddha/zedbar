//! QR Code Image Binarization
//!
//! This module implements adaptive thresholding for converting grayscale images
//! to binary (black and white) for QR code detection.
//!
//! The algorithm compares each pixel value to the mean value of a large window
//! surrounding it. This simple approach works better than more complex methods
//! like Sauvola or Gatos for QR codes, as it doesn't over-shrink isolated black
//! dots inside the code.
//!
//! Copyright (C) 2008-2009 Timothy B. Terriberry (tterribe@xiph.org)
//! Licensed under LGPL 2.1 or later

use std::cmp::{max, min};

/// Binarizes a grayscale image using adaptive thresholding
///
/// This compares each pixel value to the mean value of a large window surrounding it.
/// The window size is chosen to be large enough that it doesn't fit completely inside
/// the center of a finder pattern of a version 1 QR code at full resolution.
///
/// # Arguments
/// * `img` - The input grayscale image data (row-major order)
/// * `width` - Width of the image in pixels
/// * `height` - Height of the image in pixels
///
/// # Returns
/// A binary mask where 0xFF represents foreground (black) and 0x00 represents background (white).
/// Returns an empty vector if width or height is <= 0.
pub fn binarize(img: &[u8], width: i32, height: i32) -> Vec<u8> {
    if width <= 0 || height <= 0 {
        return Vec::new();
    }

    let width = width as usize;
    let height = height as usize;
    let mut mask = vec![0u8; width * height];

    // Determine window size (power of 2, between 16 and 256)
    // The window should be roughly 1/8 of the image dimension
    let mut logwindw = 4;
    while logwindw < 8 && (1 << logwindw) < ((width + 7) >> 3) {
        logwindw += 1;
    }

    let mut logwindh = 4;
    while logwindh < 8 && (1 << logwindh) < ((height + 7) >> 3) {
        logwindh += 1;
    }

    let windw = 1 << logwindw;
    let windh = 1 << logwindh;

    // Column sums for efficient window sum calculation
    let mut col_sums = vec![0u32; width];

    // Initialize sums down each column
    for x in 0..width {
        let g = img[x] as u32;
        col_sums[x] = (g << (logwindh - 1)) + g;
    }

    // Add initial rows to column sums
    for y in 1..(windh >> 1) {
        let y1offs = min(y, height - 1) * width;
        for x in 0..width {
            let g = img[y1offs + x] as u32;
            col_sums[x] += g;
        }
    }

    // Process each row
    for y in 0..height {
        // Initialize the sum over the window for this row
        let mut m = (col_sums[0] << (logwindw - 1)) + col_sums[0];
        for x in 1..(windw >> 1) {
            let x1 = min(x, width - 1);
            m += col_sums[x1];
        }

        // Process each pixel in the row
        for x in 0..width {
            // Perform the test against the threshold T = (m/n) - D,
            // where n = windw * windh and D = 3
            let g = img[y * width + x] as u32;
            mask[y * width + x] = if ((g + 3) << (logwindw + logwindh)) < m {
                0xFF
            } else {
                0x00
            };

            // Update the window sum
            if x + 1 < width {
                let x0 = max(0, x as i32 - (windw as i32 >> 1)) as usize;
                let x1 = min(x + (windw >> 1), width - 1);
                m += col_sums[x1];
                m -= col_sums[x0];
            }
        }

        // Update the column sums for next row
        if y + 1 < height {
            let y0offs = max(0, y as i32 - (windh as i32 >> 1)) as usize * width;
            let y1offs = min(y + (windh >> 1), height - 1) * width;
            for x in 0..width {
                col_sums[x] -= img[y0offs + x] as u32;
                col_sums[x] += img[y1offs + x] as u32;
            }
        }
    }

    mask
}

// C FFI exports

/// C FFI wrapper for qr_binarize
///
/// # Safety
///
/// - `img` must point to a valid buffer of at least `width * height` bytes
/// - The returned pointer must be freed by the caller using `free()`
/// - Returns NULL if width or height is <= 0
#[no_mangle]
pub unsafe extern "C" fn qr_binarize(img: *const u8, width: i32, height: i32) -> *mut u8 {
    if img.is_null() || width <= 0 || height <= 0 {
        return std::ptr::null_mut();
    }

    let img_slice = std::slice::from_raw_parts(img, (width * height) as usize);
    let mask = binarize(img_slice, width, height);

    // Allocate C-compatible memory and copy the result
    let ptr = libc::malloc((width * height) as usize) as *mut u8;
    if !ptr.is_null() {
        std::ptr::copy_nonoverlapping(mask.as_ptr(), ptr, mask.len());
    }
    ptr
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_binarize_empty() {
        let img = vec![0u8; 0];
        let mask = binarize(&img, 0, 0);
        assert_eq!(mask.len(), 0);

        let mask = binarize(&img, -1, -1);
        assert_eq!(mask.len(), 0);
    }

    #[test]
    fn test_binarize_uniform() {
        // Uniform gray image - all pixels should be background (0x00)
        let img = vec![128u8; 100 * 100];
        let mask = binarize(&img, 100, 100);
        assert_eq!(mask.len(), 100 * 100);

        // All pixels should be the same (either all foreground or all background)
        let first = mask[0];
        assert!(mask.iter().all(|&p| p == first));
    }

    #[test]
    fn test_binarize_simple_pattern() {
        // Create a simple pattern: black square in white background
        let width = 50;
        let height = 50;
        let mut img = vec![255u8; width * height]; // White background

        // Add a black square in the center (20x20 to 30x30)
        for y in 20..30 {
            for x in 20..30 {
                img[y * width + x] = 0;
            }
        }

        let mask = binarize(&img, width as i32, height as i32);
        assert_eq!(mask.len(), width * height);

        // The black square should be detected as foreground (0xFF)
        // Check center of the black square
        assert_eq!(mask[25 * width + 25], 0xFF);
    }

    #[test]
    fn test_binarize_gradient() {
        // Test with a gradient image
        let width = 64;
        let height = 64;
        let mut img = vec![0u8; width * height];

        for y in 0..height {
            for x in 0..width {
                img[y * width + x] = ((x * 255) / (width - 1)) as u8;
            }
        }

        let mask = binarize(&img, width as i32, height as i32);
        assert_eq!(mask.len(), width * height);

        // Left side (dark) should be foreground, right side (light) should be background
        assert_eq!(mask[32 * width], 0xFF); // Left edge
        assert_eq!(mask[32 * width + 63], 0x00); // Right edge
    }

    #[test]
    fn test_binarize_small_image() {
        // Test with very small image
        let img = vec![0, 255, 255, 0]; // 2x2 checkerboard
        let mask = binarize(&img, 2, 2);
        assert_eq!(mask.len(), 4);
    }

    #[test]
    fn test_qr_binarize_c_ffi() {
        unsafe {
            // Test NULL pointer
            let result = qr_binarize(std::ptr::null(), 10, 10);
            assert!(result.is_null());

            // Test zero dimensions
            let img = [128u8; 100];
            let result = qr_binarize(img.as_ptr(), 0, 10);
            assert!(result.is_null());

            // Test valid input
            let img = [128u8; 100];
            let result = qr_binarize(img.as_ptr(), 10, 10);
            assert!(!result.is_null());

            // Verify the result and free it
            let mask_slice = std::slice::from_raw_parts(result, 100);
            assert_eq!(mask_slice.len(), 100);
            libc::free(result as *mut libc::c_void);
        }
    }
}
