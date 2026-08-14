//! Regenerate the SQ Code fixtures.
//!
//! Run via: `cargo run --example generate_sqcode_fixture`
//!
//! Every other symbology here has a fixture produced by an off-the-shelf
//! encoder, but no SQ Code encoder appears to exist: zint, ZXing and zbar all
//! only read the format, and zbar ships no sample. So this derives the layout
//! from zbar's decoder in `zbar/sqcode.c` and renders it directly.
//!
//! Because the generator is written from the decoder it is verified against,
//! it would be worthless on its own. What makes the fixtures trustworthy is
//! that the *reference* `zbarimg` reads them and reports the payload we asked
//! for — an independent implementation agreeing on both the geometry and the
//! data. `tests/decode_examples.rs` re-checks that on every run.
//!
//! # Layout
//!
//! Positions are given in units of `H`, half the pitch between border dots,
//! with the top-left dot at the origin. A code whose data area is `S x S`
//! modules has `border_len = S + 6` grid indices, numbered `0..N-1`:
//!
//! ```text
//!   top row       dots at x = 0,2,..,N-4     y = 0
//!   left column   dots at y = 0,2,..,N-4     x = 0
//!   right column  dots at y = 3,5,..,N-1     x = (N-1)H
//!   bottom row    dots at x = 3,5,..,N-1     y = (N-1)H
//!   corners       top-right (N-1, 0), bottom-left (0, N-1)
//!   data          cells (x,y) for x,y in 3..=N-4, one module per H
//! ```
//!
//! A corner is two half-size squares filling the top-left and bottom-right
//! quadrants of its box, which is the shape `sq_scan_shape` tests for. Data
//! bits run row-major, most significant first, black = 1, and the decoded
//! payload is reported base64-encoded.

use image::{GrayImage, Luma};

/// Half the pitch between adjacent border dots, in pixels. Also the size of
/// one data module, since data cells sit on the half-pitch grid.
const H: u32 = 8;

/// Border dot / corner box size. Must exceed `H` so that the walk from the
/// last border dot to the corner (a gap of `3H`) stays inside the
/// `3 * width` search window `find_right_dot` and `find_bottom_dot` use, and
/// stay under `2H` so adjacent dots do not touch and merge.
const DOT: u32 = 12;

const QUIET: u32 = 24;

const FIXTURES: &[(&str, &[u8])] = &[
    ("examples/test-sqcode.png", b"zedbar!!"),
    ("examples/test-sqcode-large.png", b"zedbar SQ Code fix"),
];

struct Canvas {
    img: GrayImage,
    origin: u32,
}

impl Canvas {
    fn new(n: u32) -> Self {
        let size = (n - 1) * H + 2 * QUIET + DOT;
        Self {
            img: GrayImage::from_pixel(size, size, Luma([255])),
            origin: QUIET + DOT / 2,
        }
    }

    /// Pixel coordinate of grid index `i`.
    fn at(&self, i: u32) -> u32 {
        self.origin + i * H
    }

    /// Fill a `w x h` box centered on (`cx`, `cy`).
    fn fill(&mut self, cx: u32, cy: u32, w: u32, h: u32) {
        for y in cy.saturating_sub(h / 2)..cy.saturating_sub(h / 2) + h {
            for x in cx.saturating_sub(w / 2)..cx.saturating_sub(w / 2) + w {
                if x < self.img.width() && y < self.img.height() {
                    self.img.put_pixel(x, y, Luma([0]));
                }
            }
        }
    }

    fn dot(&mut self, ix: u32, iy: u32) {
        self.fill(self.at(ix), self.at(iy), DOT, DOT);
    }

    /// Two half squares meeting diagonally: black top-left and bottom-right,
    /// white top-right and bottom-left.
    fn corner(&mut self, ix: u32, iy: u32) {
        let half = DOT / 2;
        let (x0, y0) = (self.at(ix) - DOT / 2, self.at(iy) - DOT / 2);
        for y in y0..y0 + half {
            for x in x0..x0 + half {
                self.img.put_pixel(x, y, Luma([0]));
            }
        }
        for y in y0 + half..y0 + 2 * half {
            for x in x0 + half..x0 + 2 * half {
                self.img.put_pixel(x, y, Luma([0]));
            }
        }
    }
}

fn render(payload: &[u8]) -> GrayImage {
    let bits = payload.len() as u32 * 8;
    let side = (bits as f64).sqrt() as u32;
    assert_eq!(
        side * side,
        bits,
        "payload of {} bytes is {bits} bits, which is not a square data area",
        payload.len()
    );
    let n = side + 6;

    let mut c = Canvas::new(n);

    for i in (0..=n - 4).step_by(2) {
        c.dot(i, 0);
        c.dot(0, i);
    }
    for i in (3..n).step_by(2) {
        c.dot(n - 1, i);
        c.dot(i, n - 1);
    }
    c.corner(n - 1, 0);
    c.corner(0, n - 1);

    let mut idx = 0usize;
    for y in 3..=n - 4 {
        for x in 3..=n - 4 {
            if payload[idx / 8] & (1 << (7 - idx % 8)) != 0 {
                let (cx, cy) = (c.at(x), c.at(y));
                c.fill(cx, cy, H, H);
            }
            idx += 1;
        }
    }

    c.img
}

fn main() {
    for (path, payload) in FIXTURES {
        let img = render(payload);
        img.save(path)
            .unwrap_or_else(|e| panic!("save {path}: {e}"));
        println!(
            "{path}: {}x{} payload={:?}",
            img.width(),
            img.height(),
            String::from_utf8_lossy(payload)
        );
    }
    println!("\nVerify with: zbarimg -q examples/test-sqcode*.png");
}
