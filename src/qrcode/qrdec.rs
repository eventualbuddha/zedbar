//! QR Code decoder utilities
//!
//! This module provides low-level QR code decoding functions including
//! point geometry operations and error correction.

use std::ptr::null;

use libc::{c_int, c_uint};

use crate::{decoder_types::qr_finder_line, img_scanner::qr_reader, qrcode::util::qr_ilog};

use super::{
    isaac_init,
    qrdectxt::qr_point,
    rs::{rs_gf256_init, QR_PPOLY},
};

/// A point in QR code coordinate space: [x, y]
pub type QrPoint = [c_int; 2];

/// A line in QR code coordinate space: [A, B, C] for equation Ax + By + C = 0
pub type QrLine = [c_int; 3];

/// Number of bits in an int (typically 32)
const QR_INT_BITS: c_int = (std::mem::size_of::<c_int>() * 8) as c_int;

/// Affine transformation structure for QR code coordinate mapping
///
/// This structure represents an affine transformation between coordinate spaces,
/// including both forward and inverse transformations.
#[repr(C)]
pub struct QrAff {
    /// Forward transformation matrix [2][2]
    pub fwd: [[c_int; 2]; 2],
    /// Inverse transformation matrix [2][2]
    pub inv: [[c_int; 2]; 2],
    /// X offset
    pub x0: c_int,
    /// Y offset
    pub y0: c_int,
    /// Resolution bits
    pub res: c_int,
    /// Inverse resolution bits
    pub ires: c_int,
}

/// Helper function: maximum of two integers (branchless)
#[inline]
fn qr_maxi(a: c_int, b: c_int) -> c_int {
    a - ((a - b) & -((b > a) as c_int))
}

/// Helper function: flip sign of a if b is negative
#[inline]
fn qr_flipsigni(a: c_int, b: c_int) -> c_int {
    let mask = -((b < 0) as c_int);
    (a + mask) ^ mask
}

/// Helper function: divide with exact rounding
#[inline]
fn qr_divround(x: c_int, y: c_int) -> c_int {
    (x + qr_flipsigni(y >> 1, x)) / y
}

/// collection of finder lines
#[repr(C)]
pub struct qr_finder_lines {
    lines: *mut qr_finder_line,
    nlines: c_int,
    clines: c_int,
}

/// Initializes a client reader handle.
#[no_mangle]
pub unsafe extern "C" fn qr_reader_init(reader: *mut qr_reader) {
    isaac_init((&mut (*reader).isaac) as *mut _, null(), 0);
    rs_gf256_init((&mut (*reader).gf) as *mut _, QR_PPOLY);
}

/// Allocates a client reader handle.
#[no_mangle]
pub unsafe extern "C" fn _zbar_qr_create() -> *mut qr_reader {
    let reader = libc::calloc(1, size_of::<qr_reader>()) as *mut _;
    qr_reader_init(reader);
    reader
}

/// Frees a client reader handle.
#[no_mangle]
pub unsafe extern "C" fn _zbar_qr_destroy(reader: *mut qr_reader) {
    if !(*reader).finder_lines[0].lines.is_null() {
        libc::free((*reader).finder_lines[0].lines as *mut _);
    }

    if !(*reader).finder_lines[1].lines.is_null() {
        libc::free((*reader).finder_lines[1].lines as *mut _);
    }

    libc::free(reader as *mut _);
}

/// reset finder state between scans
#[no_mangle]
pub unsafe extern "C" fn _zbar_qr_reset(reader: *mut qr_reader) {
    (*reader).finder_lines[0].nlines = 0;
    (*reader).finder_lines[1].nlines = 0;
}

/// A cluster of lines crossing a finder pattern (all in the same direction).
pub struct qr_finder_cluster {
    /// Pointers to the lines crossing the pattern.
    lines: *mut *mut qr_finder_line,

    /// The number of lines in the cluster.
    nlines: c_int,
}

/// A point on the edge of a finder pattern. These are obtained from the
/// endpoints of the lines crossing this particular pattern.
struct qr_finder_edge_pt {
    /// The location of the edge point.
    pos: qr_point,

    /// A label classifying which edge this belongs to:
    /// 0: negative u edge (left)
    /// 1: positive u edge (right)
    /// 2: negative v edge (top)
    /// 3: positive v edge (bottom)*/
    edge: c_int,

    /// The (signed) perpendicular distance of the edge point from a line parallel
    /// to the edge passing through the finder center, in (u,v) coordinates.
    /// This is also re-used by RANSAC to store inlier flags.*/
    extent: c_int,
}

/*Determine if a horizontal line crosses a vertical line.
_hline: The horizontal line.
_vline: The vertical line.
Return: A non-zero value if the lines cross, or zero if they do not.*/
#[no_mangle]
pub unsafe extern "C" fn qr_finder_lines_are_crossing(
    _hline: *const qr_finder_line,
    _vline: *const qr_finder_line,
) -> c_int {
    c_int::from(
        (*_hline).pos[0] <= (*_vline).pos[0]
            && (*_vline).pos[0] < (*_hline).pos[0] + (*_hline).len
            && (*_vline).pos[1] <= (*_hline).pos[1]
            && (*_hline).pos[1] < (*_vline).pos[1] + (*_vline).len,
    )
}

/// Translate a point by the given offsets
///
/// Adds dx to the x coordinate and dy to the y coordinate.
#[no_mangle]
pub unsafe extern "C" fn qr_point_translate(point: *mut c_int, dx: c_int, dy: c_int) {
    if !point.is_null() {
        *point.offset(0) += dx;
        *point.offset(1) += dy;
    }
}

/// Calculate the squared distance between two points
///
/// Returns the squared Euclidean distance, which avoids the need for
/// expensive square root calculations when only relative distances matter.
#[no_mangle]
pub unsafe extern "C" fn qr_point_distance2(p1: *const c_int, p2: *const c_int) -> c_uint {
    let dx = *p1.offset(0) - *p2.offset(0);
    let dy = *p1.offset(1) - *p2.offset(1);
    (dx * dx + dy * dy) as c_uint
}

/// Check if three points are in counter-clockwise order
///
/// Returns the cross product of the vectors (p1-p0) and (p2-p0).
/// - Positive: points are in CCW order (in right-handed coordinate system)
/// - Zero: points are collinear
/// - Negative: points are in CW order
#[no_mangle]
pub unsafe extern "C" fn qr_point_ccw(
    p0: *const c_int,
    p1: *const c_int,
    p2: *const c_int,
) -> c_int {
    let p0x = *p0.offset(0);
    let p0y = *p0.offset(1);
    let p1x = *p1.offset(0);
    let p1y = *p1.offset(1);
    let p2x = *p2.offset(0);
    let p2y = *p2.offset(1);

    (p1x - p0x) * (p2y - p0y) - (p1y - p0y) * (p2x - p0x)
}

/// Evaluate a line equation at a point
///
/// Given a line defined by the equation A*x + B*y + C = 0,
/// this returns the value A*x + B*y + C for the given coordinates.
#[no_mangle]
pub unsafe extern "C" fn qr_line_eval(line: *const c_int, x: c_int, y: c_int) -> c_int {
    *line.offset(0) * x + *line.offset(1) * y + *line.offset(2)
}

/// Calculate step delta for moving along a line in affine space
///
/// This function computes how much the coordinate `dv` changes when
/// stepping by `du` along a line in an affine coordinate system.
///
/// # Arguments
/// - `aff`: The affine transformation
/// - `line`: The line equation coefficients [A, B, C]
/// - `v`: The coordinate index (0 for x, 1 for y)
/// - `du`: The step amount in the u direction
/// - `dv`: Output pointer for the computed step in v direction
///
/// # Returns
/// - 0 on success
/// - -1 if the line is too steep (>45 degrees from horizontal/vertical)
#[no_mangle]
pub unsafe extern "C" fn qr_aff_line_step(
    aff: *const QrAff,
    line: *const c_int,
    v: c_int,
    du: c_int,
    dv: *mut c_int,
) -> c_int {
    let aff = &*aff;
    let l0 = *line.offset(0);
    let l1 = *line.offset(1);

    // Compute n and d for the step calculation
    let mut n = aff.fwd[0][v as usize] * l0 + aff.fwd[1][v as usize] * l1;
    let mut d = aff.fwd[0][(1 - v) as usize] * l0 + aff.fwd[1][(1 - v) as usize] * l1;

    // Ensure d is positive
    if d < 0 {
        n = -n;
        d = -d;
    }

    // Calculate shift to prevent overflow
    let shift = qr_maxi(
        0,
        qr_ilog(du as u32) + qr_ilog(n.unsigned_abs()) + 3 - QR_INT_BITS,
    );
    let round = (1 << shift) >> 1;

    n = (n + round) >> shift;
    d = (d + round) >> shift;

    // The line should not be outside 45 degrees of horizontal/vertical
    // This helps ensure loop termination and avoids division by zero
    if n.abs() >= d {
        return -1;
    }

    n *= -du;
    let dv_result = qr_divround(n, d);

    if dv_result.abs() >= du {
        return -1;
    }

    *dv = dv_result;
    0
}

/// Calculate Hamming distance between two values
///
/// Counts the number of bit positions where the values differ,
/// up to a maximum of maxdiff.
#[no_mangle]
pub extern "C" fn qr_hamming_dist(y1: c_uint, y2: c_uint, maxdiff: c_int) -> c_int {
    let mut y = y1 ^ y2;
    let mut ret = 0;

    while ret < maxdiff && y != 0 {
        y &= y - 1;
        ret += 1;
    }

    ret
}

/// BCH(18,6,3) code words for version information
///
/// These codes are used for QR code version information,
/// which must be between 7 and 40 (inclusive).
const BCH18_6_CODES: [c_uint; 34] = [
    0x07C94, 0x085BC, 0x09A99, 0x0A4D3, 0x0BBF6, 0x0C762, 0x0D847, 0x0E60D, 0x0F928, 0x10B78,
    0x1145D, 0x12A17, 0x13532, 0x149A6, 0x15683, 0x168C9, 0x177EC, 0x18EC4, 0x191E1, 0x1AFAB,
    0x1B08E, 0x1CC1A, 0x1D33F, 0x1ED75, 0x1F250, 0x209D5, 0x216F0, 0x228BA, 0x2379F, 0x24B0B,
    0x2542E, 0x26A64, 0x27541, 0x28C69,
];

/// Correct a BCH(18,6,3) code word
///
/// Takes a code word and attempts to correct errors using the BCH(18,6,3) code.
/// The corrected value is written back to the input pointer.
///
/// Returns:
/// - The number of errors corrected (0-3)
/// - A negative value if more than 3 errors were detected (no correction performed)
#[no_mangle]
pub unsafe extern "C" fn bch18_6_correct(y: *mut c_uint) -> c_int {
    let y_val = *y;

    // Check the easy case first: see if the data bits were uncorrupted
    let x = y_val >> 12;
    if (7..=40).contains(&x) {
        let nerrs = qr_hamming_dist(y_val, BCH18_6_CODES[(x - 7) as usize], 4);
        if nerrs < 4 {
            *y = BCH18_6_CODES[(x - 7) as usize];
            return nerrs;
        }
    }

    // Exhaustive search is faster than field operations in GF(19)
    for (i, &code) in BCH18_6_CODES.iter().enumerate() {
        if i + 7 != (y_val >> 12) as usize {
            let nerrs = qr_hamming_dist(y_val, code, 4);
            if nerrs < 4 {
                *y = code;
                return nerrs;
            }
        }
    }

    -1
}
