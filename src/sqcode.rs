//! SQCode decoder module
//!
//! Copyright (C) 2018 Javier Serrano Polo <javier@jasp.net>
//! Rust port based on the C implementation

use crate::{image_ffi::zbar_image_t, img_scanner::zbar_image_scanner_t, SymbolType};
use libc::{c_int, size_t};
use std::io::Write;

#[derive(Debug, Copy, Clone, PartialEq)]
enum Shape {
    Dot,
    Corner,
    #[allow(dead_code)]
    Other,
    Void,
}

#[derive(Debug, Copy, Clone)]
struct Point {
    x: f32,
    y: f32,
}

#[derive(Debug, Copy, Clone)]
struct Dot {
    shape_type: Shape,
    x0: u32,
    y0: u32,
    width: u32,
    height: u32,
    center: Point,
}

pub struct SqReader {
    enabled: bool,
}
impl SqReader {
    pub(crate) fn new() -> Self {
        Self { enabled: true }
    }

    pub(crate) fn set_enabled(&mut self, enabled: bool) {
        self.enabled = enabled
    }

    pub(crate) fn reset(&mut self) {
        self.enabled = true
    }
}

const BASE64_TABLE: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

/// Base64 encode a buffer
fn base64_encode_buffer(s: &[u8]) -> Option<Vec<u8>> {
    let size = s.len();
    let encoded_size = size.div_ceil(3) * 4 + 1;
    let mut encoded = Vec::with_capacity(encoded_size);

    let mut i = 0;
    while i < size {
        let c = (s[i] >> 2) & 0x3f;
        encoded.push(BASE64_TABLE[c as usize]);

        let mut c = (s[i] << 4) & 0x30;
        i += 1;
        if i >= size {
            encoded.push(BASE64_TABLE[c as usize]);
            encoded.push(b'=');
            encoded.push(b'=');
            break;
        }

        c |= (s[i] >> 4) & 0x0f;
        encoded.push(BASE64_TABLE[c as usize]);
        c = (s[i] << 2) & 0x3c;
        i += 1;
        if i >= size {
            encoded.push(BASE64_TABLE[c as usize]);
            encoded.push(b'=');
            break;
        }

        c |= (s[i] >> 6) & 0x03;
        encoded.push(BASE64_TABLE[c as usize]);
        c = s[i] & 0x3f;
        encoded.push(BASE64_TABLE[c as usize]);
        i += 1;
    }

    encoded.push(0); // Null terminator
    Some(encoded)
}

/// Extract text from buffer and add to scanner results
fn sq_extract_text(iscn: &mut zbar_image_scanner_t, buf: &[u8], len: size_t) -> bool {
    let mut sym = iscn.alloc_sym(SymbolType::SqCode);

    let encoded = match base64_encode_buffer(&buf[..len]) {
        Some(e) => e,
        None => return true,
    };

    sym.data = encoded;
    iscn.add_symbol(sym);
    false
}

#[inline]
fn is_black_color(c: u8) -> bool {
    c <= 0x7f
}

unsafe fn is_black(img: &zbar_image_t, x: i32, y: i32) -> bool {
    if x < 0 || x >= img.width as i32 || y < 0 || y >= img.height as i32 {
        return false;
    }
    let data = img.data.as_ptr();
    let idx = y as usize * img.width as usize + x as usize;
    is_black_color(*data.add(idx))
}

fn set_dot_center(dot: &mut Dot, x: f32, y: f32) {
    dot.center.x = x;
    dot.center.y = y;
}

unsafe fn sq_scan_shape(img: &zbar_image_t, dot: &mut Dot, start_x: i32, start_y: i32) {
    if !is_black(img, start_x, start_y) {
        dot.shape_type = Shape::Void;
        dot.x0 = start_x as u32;
        dot.y0 = start_y as u32;
        dot.width = 0;
        dot.height = 0;
        set_dot_center(dot, start_x as f32, start_y as f32);
        return;
    }

    let mut x0 = start_x as u32;
    let mut y0 = start_y as u32;
    let mut width = 1u32;
    let mut height = 1u32;

    // Expand the shape to find all connected black pixels
    loop {
        let mut found_new = false;

        for x in (x0 as i32 - 1)..=(x0 as i32 + width as i32) {
            if is_black(img, x, y0 as i32 - 1) {
                y0 -= 1;
                height += 1;
                found_new = true;
                break;
            }
            if is_black(img, x, y0 as i32 + height as i32) {
                height += 1;
                found_new = true;
                break;
            }
        }

        if found_new {
            continue;
        }

        // Use checked arithmetic to avoid overflow
        let y_end = y0.saturating_add(height);
        for y in y0..y_end {
            if is_black(img, x0 as i32 - 1, y as i32) {
                x0 -= 1;
                width += 1;
                found_new = true;
                break;
            }
            if is_black(img, x0 as i32 + width as i32, y as i32) {
                width += 1;
                found_new = true;
                break;
            }
        }

        if !found_new {
            break;
        }
    }

    dot.x0 = x0;
    dot.y0 = y0;
    dot.width = width;
    dot.height = height;

    // Check if it's a corner shape
    if is_black(
        img,
        (x0 as f32 + 0.25 * width as f32) as i32,
        (y0 as f32 + 0.25 * height as f32) as i32,
    ) && !is_black(
        img,
        (x0 as f32 + 0.75 * width as f32) as i32,
        (y0 as f32 + 0.25 * height as f32) as i32,
    ) && !is_black(
        img,
        (x0 as f32 + 0.25 * width as f32) as i32,
        (y0 as f32 + 0.75 * height as f32) as i32,
    ) && is_black(
        img,
        (x0 as f32 + 0.75 * width as f32) as i32,
        (y0 as f32 + 0.75 * height as f32) as i32,
    ) {
        dot.shape_type = Shape::Corner;
        set_dot_center(
            dot,
            x0 as f32 + 0.5 * width as f32,
            y0 as f32 + 0.5 * height as f32,
        );
        return;
    }

    // Calculate weighted center for dot
    let data = img.data.as_ptr();
    let mut x_sum = 0u64;
    let mut y_sum = 0u64;
    let mut total_weight = 0u64;

    // Use checked arithmetic to avoid overflow
    let y_end = y0.saturating_add(height);
    let x_end = x0.saturating_add(width);

    for y in y0..y_end {
        for x in x0..x_end {
            if !is_black(img, x as i32, y as i32) {
                continue;
            }
            let idx = y as usize * img.width as usize + x as usize;
            let weight = (0xff - *data.add(idx)) as u64;
            x_sum += weight * x as u64;
            y_sum += weight * y as u64;
            total_weight += weight;
        }
    }

    dot.shape_type = Shape::Dot;
    set_dot_center(
        dot,
        x_sum as f32 / total_weight as f32 + 0.5,
        y_sum as f32 / total_weight as f32 + 0.5,
    );
}

unsafe fn find_left_dot(
    img: &zbar_image_t,
    dot: &Dot,
    found_x: &mut u32,
    found_y: &mut u32,
) -> bool {
    let y_end = dot.y0.saturating_add(dot.height);
    for y in dot.y0..y_end {
        for x in ((dot.x0 as i32 - 2 * dot.width as i32)..=(dot.x0 as i32 - 1)).rev() {
            if is_black(img, x, y as i32) {
                *found_x = x as u32;
                *found_y = y;
                return true;
            }
        }
    }
    false
}

unsafe fn find_right_dot(
    img: &zbar_image_t,
    dot: &Dot,
    found_x: &mut u32,
    found_y: &mut u32,
) -> bool {
    let y_end = dot.y0.saturating_add(dot.height);
    let x_start = dot.x0.saturating_add(dot.width);
    let x_end = dot.x0.saturating_add(3 * dot.width);

    for y in dot.y0..y_end {
        for x in x_start..x_end {
            if is_black(img, x as i32, y as i32) {
                *found_x = x;
                *found_y = y;
                return true;
            }
        }
    }
    false
}

unsafe fn find_bottom_dot(
    img: &zbar_image_t,
    dot: &Dot,
    found_x: &mut u32,
    found_y: &mut u32,
) -> bool {
    let x_end = dot
        .x0
        .checked_add(dot.width)
        .and_then(|v| v.checked_sub(1))
        .unwrap_or(0);
    let y_start = dot.y0.saturating_add(dot.height);
    let y_end = dot.y0.saturating_add(3 * dot.height);

    for x in (dot.x0..=x_end).rev() {
        for y in y_start..y_end {
            if is_black(img, x as i32, y as i32) {
                *found_x = x;
                *found_y = y;
                return true;
            }
        }
    }
    false
}

/// Main SQCode decoding function
///
/// # Safety
///
/// All pointers must be valid and properly initialized. The caller must ensure:
/// - `reader` points to a valid `SqReader` instance
/// - `iscn` points to a valid `zbar_image_scanner_t` instance
/// - `img` points to a valid `zbar_image_t` with properly initialized image data
pub(crate) unsafe fn sq_decode(
    reader: &mut SqReader,
    iscn: &mut zbar_image_scanner_t,
    img: &mut zbar_image_t,
) -> c_int {
    if !reader.enabled {
        return 0;
    }

    // Check image format (Y800 = fourcc('Y','8','0','0'))
    if img.format != 0x30303859 {
        let _ = writeln!(std::io::stderr(), "Unexpected image format");
        return 1;
    }

    // Find starting pixel
    let mut scan_x = 0u32;
    let mut scan_y = 0u32;
    let mut found_start = false;

    'outer: for y in 0..img.height {
        for x in 0..img.width {
            if is_black(img, x as i32, y as i32) {
                scan_x = x;
                scan_y = y;
                found_start = true;
                break 'outer;
            }
        }
    }

    if !found_start {
        return 1;
    }

    // Scan starting dot
    let mut start_dot = Dot {
        shape_type: Shape::Void,
        x0: 0,
        y0: 0,
        width: 0,
        height: 0,
        center: Point { x: 0.0, y: 0.0 },
    };
    sq_scan_shape(img, &mut start_dot, scan_x as i32, scan_y as i32);

    let start_corner = start_dot.shape_type == Shape::Corner;

    let mut top_border: Vec<Point>;
    let mut left_border: Vec<Point>;
    let mut right_border: Vec<Point>;
    let mut bottom_border: Vec<Point>;

    let mut border_len;
    if start_corner {
        border_len = 0;
        top_border = Vec::new();
    } else {
        border_len = 1;
        top_border = vec![start_dot.center];
    }

    // Scan left from starting dot
    let mut top_left_dot = start_dot;
    while find_left_dot(img, &top_left_dot, &mut scan_x, &mut scan_y) {
        sq_scan_shape(img, &mut top_left_dot, scan_x as i32, scan_y as i32);
        if top_left_dot.shape_type != Shape::Dot {
            return 1;
        }
        if border_len > 0 {
            border_len += 2;
            top_border.reserve(2);
            top_border.insert(0, top_left_dot.center);
            top_border.insert(1, Point { x: 0.0, y: 0.0 });
            let middle = Point {
                x: (top_border[0].x + top_border[2].x) / 2.0,
                y: (top_border[0].y + top_border[2].y) / 2.0,
            };
            top_border[1] = middle;
        } else {
            border_len = 1;
            top_border = vec![top_left_dot.center];
        }
    }
    if top_left_dot.shape_type != Shape::Dot {
        return 1;
    }

    // Scan right from starting dot
    let mut top_right_dot = start_dot;
    if !start_corner {
        while find_right_dot(img, &top_right_dot, &mut scan_x, &mut scan_y) {
            sq_scan_shape(img, &mut top_right_dot, scan_x as i32, scan_y as i32);
            if top_right_dot.shape_type == Shape::Corner {
                break;
            }
            if top_right_dot.shape_type != Shape::Dot {
                return 1;
            }
            border_len += 2;
            top_border.push(top_right_dot.center);
            let len = top_border.len();
            top_border.insert(len - 1, Point { x: 0.0, y: 0.0 });
            let len = top_border.len();
            let middle = Point {
                x: (top_border[len - 3].x + top_border[len - 1].x) / 2.0,
                y: (top_border[len - 3].y + top_border[len - 1].y) / 2.0,
            };
            top_border[len - 2] = middle;
        }
    }

    if border_len < 3 {
        return 1;
    }

    let inc_x = top_border[border_len - 1].x - top_border[border_len - 3].x;
    let inc_y = top_border[border_len - 1].y - top_border[border_len - 3].y;
    border_len += 3;
    top_border.resize(border_len, Point { x: 0.0, y: 0.0 });
    top_border[border_len - 3].x = top_border[border_len - 4].x + 0.5 * inc_x;
    top_border[border_len - 3].y = top_border[border_len - 4].y + 0.5 * inc_y;
    top_border[border_len - 2].x = top_border[border_len - 4].x + inc_x;
    top_border[border_len - 2].y = top_border[border_len - 4].y + inc_y;
    top_border[border_len - 1].x = top_border[border_len - 4].x + 1.5 * inc_x;
    top_border[border_len - 1].y = top_border[border_len - 4].y + 1.5 * inc_y;

    // Scan left border
    left_border = vec![Point { x: 0.0, y: 0.0 }; border_len];
    left_border[0] = top_border[0];

    let mut bottom_left_dot = top_left_dot;
    let mut cur_len = 1;
    while find_bottom_dot(img, &bottom_left_dot, &mut scan_x, &mut scan_y) {
        sq_scan_shape(img, &mut bottom_left_dot, scan_x as i32, scan_y as i32);
        if bottom_left_dot.shape_type == Shape::Corner {
            break;
        }
        if bottom_left_dot.shape_type != Shape::Dot {
            return 1;
        }
        cur_len += 2;
        if cur_len > border_len {
            return 1;
        }
        left_border[cur_len - 1] = bottom_left_dot.center;
        let middle = Point {
            x: (left_border[cur_len - 3].x + left_border[cur_len - 1].x) / 2.0,
            y: (left_border[cur_len - 3].y + left_border[cur_len - 1].y) / 2.0,
        };
        left_border[cur_len - 2] = middle;
    }

    if cur_len != border_len - 3 || bottom_left_dot.shape_type != Shape::Corner {
        return 1;
    }

    let inc_x = left_border[cur_len - 1].x - left_border[cur_len - 3].x;
    let inc_y = left_border[cur_len - 1].y - left_border[cur_len - 3].y;
    left_border[border_len - 3].x = left_border[border_len - 4].x + 0.5 * inc_x;
    left_border[border_len - 3].y = left_border[border_len - 4].y + 0.5 * inc_y;
    left_border[border_len - 2].x = left_border[border_len - 4].x + inc_x;
    left_border[border_len - 2].y = left_border[border_len - 4].y + inc_y;
    left_border[border_len - 1].x = left_border[border_len - 4].x + 1.5 * inc_x;
    left_border[border_len - 1].y = left_border[border_len - 4].y + 1.5 * inc_y;

    // Scan right border
    right_border = vec![Point { x: 0.0, y: 0.0 }; border_len];

    let mut bottom_right_dot = top_right_dot;
    let mut cur_len = 3;
    while find_bottom_dot(img, &bottom_right_dot, &mut scan_x, &mut scan_y) {
        sq_scan_shape(img, &mut bottom_right_dot, scan_x as i32, scan_y as i32);
        if bottom_right_dot.shape_type != Shape::Dot {
            return 1;
        }
        if cur_len == 3 {
            cur_len += 1;
            if cur_len > border_len {
                return 1;
            }
            right_border[cur_len - 1] = bottom_right_dot.center;
        } else {
            cur_len += 2;
            if cur_len > border_len {
                return 1;
            }
            right_border[cur_len - 1] = bottom_right_dot.center;
            let middle = Point {
                x: (right_border[cur_len - 3].x + right_border[cur_len - 1].x) / 2.0,
                y: (right_border[cur_len - 3].y + right_border[cur_len - 1].y) / 2.0,
            };
            right_border[cur_len - 2] = middle;
        }
    }

    if cur_len != border_len || border_len < 6 {
        return 1;
    }

    let inc_x = right_border[5].x - right_border[3].x;
    let inc_y = right_border[5].y - right_border[3].y;
    right_border[2].x = right_border[3].x - 0.5 * inc_x;
    right_border[2].y = right_border[3].y - 0.5 * inc_y;
    right_border[1].x = right_border[3].x - inc_x;
    right_border[1].y = right_border[3].y - inc_y;
    right_border[0].x = right_border[3].x - 1.5 * inc_x;
    right_border[0].y = right_border[3].y - 1.5 * inc_y;

    // Scan bottom border
    bottom_border = vec![Point { x: 0.0, y: 0.0 }; border_len];
    bottom_border[border_len - 1] = right_border[border_len - 1];

    let mut bottom_left2_dot = bottom_right_dot;
    let mut offset = border_len - 1;
    while find_left_dot(img, &bottom_left2_dot, &mut scan_x, &mut scan_y) {
        sq_scan_shape(img, &mut bottom_left2_dot, scan_x as i32, scan_y as i32);
        if bottom_left2_dot.shape_type == Shape::Corner {
            break;
        }
        if bottom_left2_dot.shape_type != Shape::Dot {
            return 1;
        }
        if offset < 2 {
            return 1;
        }
        offset -= 2;
        bottom_border[offset] = bottom_left2_dot.center;
        let middle = Point {
            x: (bottom_border[offset].x + bottom_border[offset + 2].x) / 2.0,
            y: (bottom_border[offset].y + bottom_border[offset + 2].y) / 2.0,
        };
        bottom_border[offset + 1] = middle;
    }

    if offset != 3 || bottom_left2_dot.shape_type != Shape::Corner {
        return 1;
    }

    let inc_x = bottom_border[5].x - bottom_border[3].x;
    let inc_y = bottom_border[5].y - bottom_border[3].y;
    bottom_border[2].x = bottom_border[3].x - 0.5 * inc_x;
    bottom_border[2].y = bottom_border[3].y - 0.5 * inc_y;
    bottom_border[1].x = bottom_border[3].x - inc_x;
    bottom_border[1].y = bottom_border[3].y - inc_y;
    bottom_border[0].x = bottom_border[3].x - 1.5 * inc_x;
    bottom_border[0].y = bottom_border[3].y - 1.5 * inc_y;

    // Size check
    if !(8 + 2 * (1 + 2)..=65535).contains(&border_len) {
        return 1;
    }

    let bit_side_len = border_len - 2 * (1 + 2);
    let bit_len = bit_side_len * bit_side_len;
    if bit_len % 8 != 0 {
        return 1;
    }

    let byte_len = bit_len / 8;
    let mut buf = vec![0u8; byte_len];

    let mut idx = 0;
    for y in 3..=border_len - 4 {
        for x in 3..=border_len - 4 {
            let bottom_weight = y as f32 / (border_len - 1) as f32;
            let top_weight = 1.0 - bottom_weight;
            let right_weight = x as f32 / (border_len - 1) as f32;
            let left_weight = 1.0 - right_weight;

            let top_left_source_x = top_border[x].x + left_border[y].x - left_border[0].x;
            let top_left_source_y = top_border[x].y + left_border[y].y - left_border[0].y;

            let bottom_right_source_x =
                bottom_border[x].x + right_border[y].x - right_border[border_len - 1].x;
            let bottom_right_source_y =
                bottom_border[x].y + right_border[y].y - right_border[border_len - 1].y;

            let data = img.data.as_ptr();
            let sample_x = top_left_source_x as usize;
            let sample_y = top_left_source_y as usize;
            let top_left_color = *data.add(sample_y * img.width as usize + sample_x);

            let sample_x = bottom_right_source_x as usize;
            let sample_y = bottom_right_source_y as usize;
            let bottom_right_color = *data.add(sample_y * img.width as usize + sample_x);

            let mixed_color = ((top_weight + left_weight) * top_left_color as f32
                + (bottom_weight + right_weight) * bottom_right_color as f32)
                / 2.0;

            if is_black_color(mixed_color as u8) {
                buf[idx / 8] |= 1 << (7 - (idx % 8));
            }
            idx += 1;
        }
    }

    let error = sq_extract_text(iscn, &buf, byte_len);
    if error {
        1
    } else {
        0
    }
}
