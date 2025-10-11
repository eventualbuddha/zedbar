//! QR Code decoder utilities
//!
//! This module provides low-level QR code decoding functions including
//! point geometry operations and error correction.

use std::{
    cmp::Ordering,
    ffi::c_void,
    ptr::{null, null_mut},
};

use libc::{c_int, c_uchar, c_uint, calloc, free, malloc, memcpy, qsort, realloc, size_t};

use crate::{
    decoder_types::qr_finder_line,
    img_scanner::qr_reader,
    qrcode::util::{qr_ihypot, qr_ilog},
};

use super::{
    isaac_init,
    rs::{rs_gf256_init, QR_PPOLY},
};

/// A line in QR code coordinate space: [A, B, C] for equation Ax + By + C = 0
pub type qr_line = [c_int; 3];

/// A point in QR code coordinate space: [x, y]
pub type qr_point = [c_int; 2];

/// Number of bits in an int (typically 32)
const QR_INT_BITS: c_int = c_int::BITS as c_int;

/// Log2 of QR_INT_BITS (typically 5 for 32-bit ints)
const QR_INT_LOGBITS: c_int = 5; // qr_ilog(32) = 5

/// Number of bits of sub-module precision for alignment pattern search
const QR_ALIGN_SUBPREC: c_int = 2;

/// Number of bits of sub-module precision for finder pattern coordinates
const QR_FINDER_SUBPREC: c_int = 2;

/// Helper function: divide with exact rounding
///
/// Rounds towards positive infinity when x > 0, towards negative infinity when x < 0.
/// For x/y where the fractional part is exactly 0.5, rounds away from zero.
#[inline]
fn qr_divround(x: c_int, y: c_int) -> c_int {
    x.wrapping_add(x.signum().wrapping_mul(y >> 1)) / y
}

/// Fixed-point multiply with rounding and shift
///
/// Multiplies 32-bit numbers a and b, adds r, and takes bits [s, s+31] of the result.
/// This is used for fixed-point arithmetic to avoid overflow.
#[inline]
fn qr_fixmul(a: c_int, b: c_int, r: i64, s: c_int) -> c_int {
    ((a as i64 * b as i64 + r) >> s) as c_int
}

/// Returns the minimum of two integers
///
/// Uses bitwise arithmetic to avoid branches (matching C macro QR_MINI).
#[inline]
fn qr_mini(a: c_int, b: c_int) -> c_int {
    a + (((b - a) & -((b < a) as c_int)))
}

/// Returns the maximum of two integers
///
/// Uses bitwise arithmetic to avoid branches (matching C macro QR_MAXI).
#[inline]
fn qr_maxi(a: c_int, b: c_int) -> c_int {
    a - (((a - b) & -((b > a) as c_int)))
}

/// Extended multiply: multiplies 32-bit numbers a and b, adds r, and returns 64-bit result
///
/// This matches C macro QR_EXTMUL.
#[inline]
fn qr_extmul(a: c_int, b: c_int, r: i64) -> i64 {
    a as i64 * b as i64 + r
}

/// A cell in the sampling grid for homographic projection
///
/// Represents a mapping from a unit square to a quadrilateral in the image,
/// used for extracting QR code modules with perspective correction.
#[repr(C)]
pub struct qr_hom_cell {
    /// Forward transformation matrix [3][3]
    pub fwd: [[c_int; 3]; 3],
    /// X offset in image space
    pub x0: c_int,
    /// Y offset in image space
    pub y0: c_int,
    /// U offset in code space
    pub u0: c_int,
    /// V offset in code space
    pub v0: c_int,
}

/// Sampling grid for QR code module extraction
#[repr(C)]
pub struct qr_sampling_grid {
    /// Array of homography cells for mapping between code and image space
    pub cells: [*mut qr_hom_cell; 6],
    /// Mask indicating which modules are part of function patterns
    pub fpmask: *mut c_uint,
    /// Limits for each cell region
    pub cell_limits: [c_int; 6],
    /// Number of cells in use
    pub ncells: c_int,
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
#[repr(C)]
pub struct qr_finder_cluster {
    /// Pointers to the lines crossing the pattern.
    lines: *mut *mut qr_finder_line,

    /// The number of lines in the cluster.
    nlines: c_int,
}

/// A point on the edge of a finder pattern. These are obtained from the
/// endpoints of the lines crossing this particular pattern.
#[repr(C)]
pub struct qr_finder_edge_pt {
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

/// The center of a finder pattern obtained from the crossing of one or more
/// clusters of horizontal finder lines with one or more clusters of vertical
/// finder lines.
#[repr(C)]
pub struct qr_finder_center {
    /// The estimated location of the finder center.
    pos: qr_point,

    /// The list of edge points from the crossing lines.
    edge_pts: *mut qr_finder_edge_pt,

    /// The number of edge points from the crossing lines.
    nedge_pts: c_int,
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
pub unsafe extern "C" fn qr_line_eval(line: *const qr_line, x: c_int, y: c_int) -> c_int {
    (*line)[0] * x + (*line)[1] * y + (*line)[2]
}

#[no_mangle]
pub unsafe extern "C" fn qr_line_orient(_l: *mut qr_line, _x: c_int, _y: c_int) {
    if qr_line_eval(_l, _x, _y) < 0 {
        (*_l)[0] = -(*_l)[0];
        (*_l)[1] = -(*_l)[1];
        (*_l)[2] = -(*_l)[2];
    }
}

#[no_mangle]
pub unsafe extern "C" fn qr_line_isect(
    _p: *mut qr_point,
    _l0: *const qr_line,
    _l1: *const qr_line,
) -> c_int {
    let mut d = (*_l0)[0]
        .wrapping_mul((*_l1)[1])
        .wrapping_sub((*_l0)[1].wrapping_mul((*_l1)[0]));
    if d == 0 {
        return -1;
    }
    let mut x = (*_l0)[1]
        .wrapping_mul((*_l1)[2])
        .wrapping_sub((*_l1)[1].wrapping_mul((*_l0)[2]));
    let mut y = (*_l1)[0]
        .wrapping_mul((*_l0)[2])
        .wrapping_sub((*_l0)[0].wrapping_mul((*_l1)[2]));
    if d < 0 {
        x = -x;
        y = -y;
        d = -d;
    }
    (*_p)[0] = qr_divround(x, d);
    (*_p)[1] = qr_divround(y, d);
    0
}

/// Fit a line to covariance data using least-squares
///
/// # Parameters
/// - `_l`: Output line (Ax + By + C = 0)
/// - `_x0`, `_y0`: Centroid coordinates
/// - `_sxx`, `_sxy`, `_syy`: Covariance matrix values
/// - `_res`: Resolution bits for scaling
#[no_mangle]
pub unsafe extern "C" fn qr_line_fit(
    _l: *mut qr_line,
    _x0: c_int,
    _y0: c_int,
    _sxx: c_int,
    _sxy: c_int,
    _syy: c_int,
    _res: c_int,
) {
    let u = (_sxx - _syy).abs();
    let v = -_sxy << 1;
    let w = qr_ihypot(u, v) as c_int;

    // Compute shift factor to scale down into manageable range
    // Ensure product of any two of _l[0] and _l[1] fits within _res bits
    let dshift = 0.max(
        qr_ilog(u.max(v.abs()) as u32) + 1 - ((_res + 1) >> 1)
    );
    let dround = (1 << dshift) >> 1;

    if _sxx > _syy {
        (*_l)[0] = (v + dround) >> dshift;
        (*_l)[1] = (u + w + dround) >> dshift;
    } else {
        (*_l)[0] = (u + w + dround) >> dshift;
        (*_l)[1] = (v + dround) >> dshift;
    }
    (*_l)[2] = -(_x0.wrapping_mul((*_l)[0]).wrapping_add(_y0.wrapping_mul((*_l)[1])));
}

/// Perform a least-squares line fit to a list of points
///
/// At least two points are required.
///
/// # Parameters
/// - `_l`: Output line (Ax + By + C = 0)
/// - `_p`: Array of points
/// - `_np`: Number of points
/// - `_res`: Resolution bits for scaling
#[no_mangle]
pub unsafe extern "C" fn qr_line_fit_points(
    _l: *mut qr_line,
    _p: *mut qr_point,
    _np: c_int,
    _res: c_int,
) {
    let mut sx: c_int = 0;
    let mut sy: c_int = 0;
    let mut xmin = c_int::MAX;
    let mut xmax = c_int::MIN;
    let mut ymin = c_int::MAX;
    let mut ymax = c_int::MIN;

    // Compute centroid and bounds
    for i in 0.._np {
        let px = (*_p.offset(i as isize))[0];
        let py = (*_p.offset(i as isize))[1];
        sx = sx.wrapping_add(px);
        xmin = xmin.min(px);
        xmax = xmax.max(px);
        sy = sy.wrapping_add(py);
        ymin = ymin.min(py);
        ymax = ymax.max(py);
    }

    let xbar = (sx + (_np >> 1)) / _np;
    let ybar = (sy + (_np >> 1)) / _np;

    // Compute shift to prevent overflow in covariance calculation
    let sshift = 0.max(
        qr_ilog(
            (_np as u32).wrapping_mul(
                (xmax - xbar)
                    .max(xbar - xmin)
                    .max(ymax - ybar)
                    .max(ybar - ymin) as u32,
            ),
        ) - ((QR_INT_BITS - 1) >> 1),
    );
    let sround = (1 << sshift) >> 1;

    // Compute covariance matrix
    let mut sxx: c_int = 0;
    let mut sxy: c_int = 0;
    let mut syy: c_int = 0;
    for i in 0.._np {
        let dx = ((*_p.offset(i as isize))[0] - xbar + sround) >> sshift;
        let dy = ((*_p.offset(i as isize))[1] - ybar + sround) >> sshift;
        sxx = sxx.wrapping_add(dx.wrapping_mul(dx));
        sxy = sxy.wrapping_add(dx.wrapping_mul(dy));
        syy = syy.wrapping_add(dy.wrapping_mul(dy));
    }

    qr_line_fit(_l, xbar, ybar, sxx, sxy, syy, _res);
}

/// An affine homography.
/// This maps from the image (at subpel resolution) to a square domain with
/// power-of-two sides (of res bits) and back.
#[repr(C)]
pub struct qr_aff {
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

#[no_mangle]
pub unsafe extern "C" fn qr_aff_init(
    _aff: *mut qr_aff,
    _p0: *const qr_point,
    _p1: *const qr_point,
    _p2: *const qr_point,
    _res: c_int,
) {
    // det is ensured to be positive by our caller.
    let dx1 = (*_p1)[0] - (*_p0)[0];
    let dx2 = (*_p2)[0] - (*_p0)[0];
    let dy1 = (*_p1)[1] - (*_p0)[1];
    let dy2 = (*_p2)[1] - (*_p0)[1];
    let det = dx1 * dy2 - dy1 * dx2;
    let ires = c_int::max(((qr_ilog(det.unsigned_abs()) as u32 >> 1) - 2) as i32, 0);
    (*_aff).fwd[0][0] = dx1;
    (*_aff).fwd[0][1] = dx2;
    (*_aff).fwd[1][0] = dy1;
    (*_aff).fwd[1][1] = dy2;
    (*_aff).inv[0][0] = qr_divround(dy2 << _res, det >> ires);
    (*_aff).inv[0][1] = qr_divround(-dx2 << _res, det >> ires);
    (*_aff).inv[1][0] = qr_divround(-dy1 << _res, det >> ires);
    (*_aff).inv[1][1] = qr_divround(dx1 << _res, det >> ires);
    (*_aff).x0 = (*_p0)[0];
    (*_aff).y0 = (*_p0)[1];
    (*_aff).res = _res;
    (*_aff).ires = ires;
}

/// Map from the image (at subpel resolution) into the square domain.
#[no_mangle]
pub unsafe extern "C" fn qr_aff_unproject(
    _q: *mut qr_point,
    _aff: *const qr_aff,
    _x: c_int,
    _y: c_int,
) {
    (*_q)[0] = ((*_aff).inv[0][0]
        .wrapping_mul(_x.wrapping_sub((*_aff).x0))
        .wrapping_add((*_aff).inv[0][1].wrapping_mul(_y.wrapping_sub((*_aff).y0)))
        .wrapping_add((1 << (*_aff).ires) >> 1))
        >> (*_aff).ires;
    (*_q)[1] = ((*_aff).inv[1][0]
        .wrapping_mul(_x.wrapping_sub((*_aff).x0))
        .wrapping_add((*_aff).inv[1][1].wrapping_mul(_y.wrapping_sub((*_aff).y0)))
        .wrapping_add((1 << (*_aff).ires) >> 1))
        >> (*_aff).ires;
}

/// Map from the square domain into the image (at subpel resolution).
#[no_mangle]
pub unsafe extern "C" fn qr_aff_project(
    _p: *mut qr_point,
    _aff: *const qr_aff,
    _u: c_int,
    _v: c_int,
) {
    (*_p)[0] = (((*_aff).fwd[0][0]
        .wrapping_mul(_u)
        .wrapping_add((*_aff).fwd[0][1].wrapping_mul(_v))
        .wrapping_add(1 << ((*_aff).res - 1)))
        >> (*_aff).res)
        .wrapping_add((*_aff).x0);
    (*_p)[1] = (((*_aff).fwd[1][0]
        .wrapping_mul(_u)
        .wrapping_add((*_aff).fwd[1][1].wrapping_mul(_v))
        .wrapping_add(1 << ((*_aff).res - 1)))
        >> (*_aff).res)
        .wrapping_add((*_aff).y0);
}

/// A full homography.
/// Like the affine homography, this maps from the image (at subpel resolution)
/// to a square domain with power-of-two sides (of res bits) and back.
#[repr(C)]
pub struct qr_hom {
    fwd: [[c_int; 2]; 3],
    inv: [[c_int; 2]; 3],
    fwd22: c_int,
    inv22: c_int,
    x0: c_int,
    y0: c_int,
    res: c_int,
}

/// Initialize a homography mapping from four corner points
///
/// Computes both the forward and inverse homography transformations
/// from the given corner points to a square domain.
///
/// # Safety
/// This function is unsafe because it dereferences the raw _hom pointer.
#[no_mangle]
pub unsafe extern "C" fn qr_hom_init(
    _hom: *mut qr_hom,
    _x0: c_int,
    _y0: c_int,
    _x1: c_int,
    _y1: c_int,
    _x2: c_int,
    _y2: c_int,
    _x3: c_int,
    _y3: c_int,
    _res: c_int,
) {
    let dx10 = _x1 - _x0;
    let dx20 = _x2 - _x0;
    let dx30 = _x3 - _x0;
    let dx31 = _x3 - _x1;
    let dx32 = _x3 - _x2;
    let dy10 = _y1 - _y0;
    let dy20 = _y2 - _y0;
    let dy30 = _y3 - _y0;
    let dy31 = _y3 - _y1;
    let dy32 = _y3 - _y2;
    let a20 = dx32 * dy10 - dx10 * dy32;
    let a21 = dx20 * dy31 - dx31 * dy20;
    let a22 = dx32 * dy31 - dx31 * dy32;

    // Figure out if we need to downscale anything
    let b0 = qr_ilog(qr_maxi(dx10.abs(), dy10.abs()) as u32)
        + qr_ilog((a20 + a22).abs() as u32);
    let b1 = qr_ilog(qr_maxi(dx20.abs(), dy20.abs()) as u32)
        + qr_ilog((a21 + a22).abs() as u32);
    let b2 = qr_ilog(qr_maxi(qr_maxi(a20.abs(), a21.abs()), a22.abs()) as u32);
    let s1 = qr_maxi(0, _res + qr_maxi(qr_maxi(b0, b1), b2) - (QR_INT_BITS - 2));
    let r1 = (1i64 << s1) >> 1;

    // Compute the final coefficients of the forward transform
    // The 32x32->64 bit multiplies are really needed for accuracy with large versions
    (*_hom).fwd[0][0] = qr_fixmul(dx10, a20 + a22, r1, s1);
    (*_hom).fwd[0][1] = qr_fixmul(dx20, a21 + a22, r1, s1);
    (*_hom).x0 = _x0;
    (*_hom).fwd[1][0] = qr_fixmul(dy10, a20 + a22, r1, s1);
    (*_hom).fwd[1][1] = qr_fixmul(dy20, a21 + a22, r1, s1);
    (*_hom).y0 = _y0;
    (*_hom).fwd[2][0] = (a20 + r1 as c_int) >> s1;
    (*_hom).fwd[2][1] = (a21 + r1 as c_int) >> s1;
    (*_hom).fwd22 = if s1 > _res {
        (a22 + ((r1 >> _res) as c_int)) >> (s1 - _res)
    } else {
        a22 << (_res - s1)
    };

    // Now compute the inverse transform
    let b0 = qr_ilog(qr_maxi(qr_maxi(dx10.abs(), dx20.abs()), dx30.abs()) as u32)
        + qr_ilog(qr_maxi((*_hom).fwd[0][0].abs(), (*_hom).fwd[1][0].abs()) as u32);
    let b1 = qr_ilog(qr_maxi(qr_maxi(dy10.abs(), dy20.abs()), dy30.abs()) as u32)
        + qr_ilog(qr_maxi((*_hom).fwd[0][1].abs(), (*_hom).fwd[1][1].abs()) as u32);
    let b2 = qr_ilog(a22.abs() as u32) - s1;
    let s2 = qr_maxi(0, qr_maxi(b0, b1) + b2 - (QR_INT_BITS - 3));
    let r2 = (1i64 << s2) >> 1;
    let s1 = s1 + s2;
    let r1 = r1 << s2;

    // The 32x32->64 bit multiplies are really needed for accuracy with large versions
    (*_hom).inv[0][0] = qr_fixmul((*_hom).fwd[1][1], a22, r1, s1);
    (*_hom).inv[0][1] = qr_fixmul(-(*_hom).fwd[0][1], a22, r1, s1);
    (*_hom).inv[1][0] = qr_fixmul(-(*_hom).fwd[1][0], a22, r1, s1);
    (*_hom).inv[1][1] = qr_fixmul((*_hom).fwd[0][0], a22, r1, s1);
    (*_hom).inv[2][0] = qr_fixmul(
        (*_hom).fwd[1][0],
        (*_hom).fwd[2][1],
        -qr_extmul((*_hom).fwd[1][1], (*_hom).fwd[2][0], r2),
        s2,
    );
    (*_hom).inv[2][1] = qr_fixmul(
        (*_hom).fwd[0][1],
        (*_hom).fwd[2][0],
        -qr_extmul((*_hom).fwd[0][0], (*_hom).fwd[2][1], r2),
        s2,
    );
    (*_hom).inv22 = qr_fixmul(
        (*_hom).fwd[0][0],
        (*_hom).fwd[1][1],
        -qr_extmul((*_hom).fwd[0][1], (*_hom).fwd[1][0], r2),
        s2,
    );
    (*_hom).res = _res;
}

/// Finish a partial projection, converting from homogeneous coordinates to the
/// normal 2-D representation.
/// In loops, we can avoid many multiplies by computing the homogeneous _x, _y,
/// and _w incrementally, but we cannot avoid the divisions, done here.
#[no_mangle]
pub unsafe extern "C" fn qr_hom_fproject(
    _p: *mut qr_point,
    _hom: *const qr_hom,
    mut _x: c_int,
    mut _y: c_int,
    mut _w: c_int,
) {
    if _w == 0 {
        (*_p)[0] = if _x < 0 { c_int::MIN } else { c_int::MAX };
        (*_p)[1] = if _y < 0 { c_int::MIN } else { c_int::MAX };
    } else {
        if _w < 0 {
            _x = -_x;
            _y = -_y;
            _w = -_w;
        }
        (*_p)[0] = qr_divround(_x, _w) + (*_hom).x0;
        (*_p)[1] = qr_divround(_y, _w) + (*_hom).y0;
    }
}

/// All the information we've collected about a finder pattern in the current
/// configuration.
#[repr(C)]
pub struct qr_finder {
    /// The module size along each axis (in the square domain).
    size: [c_int; 2],

    /// The version estimated from the module size along each axis.
    eversion: [c_int; 2],

    /// The list of classified edge points for each edge.
    edge_pts: [*mut qr_finder_edge_pt; 4],

    /// The number of edge points classified as belonging to each edge.
    nedge_pts: [c_int; 4],

    /// The number of inliers found after running RANSAC on each edge.
    ninliers: [c_int; 4],

    /// The center of the finder pattern (in the square domain).
    o: qr_point,

    /// The finder center information from the original image.
    c: *mut qr_finder_center,
}

/// Map from the image (at subpel resolution) into the square domain.
/// Returns a negative value if the point went to infinity.
#[no_mangle]
pub unsafe extern "C" fn qr_hom_unproject(
    _q: *mut qr_point,
    _hom: *const qr_hom,
    mut _x: c_int,
    mut _y: c_int,
) -> c_int {
    _x = _x.wrapping_sub((*_hom).x0);
    _y = _y.wrapping_sub((*_hom).y0);
    let mut x = (*_hom).inv[0][0]
        .wrapping_mul(_x)
        .wrapping_add((*_hom).inv[0][1].wrapping_mul(_y));
    let mut y = (*_hom).inv[1][0]
        .wrapping_mul(_x)
        .wrapping_add((*_hom).inv[1][1].wrapping_mul(_y));
    let mut w = ((*_hom).inv[2][0]
        .wrapping_mul(_x)
        .wrapping_add((*_hom).inv[2][1].wrapping_mul(_y))
        .wrapping_add((*_hom).inv22)
        .wrapping_add(1 << ((*_hom).res - 1)))
        >> (*_hom).res;
    if w == 0 {
        (*_q)[0] = if x < 0 { c_int::MIN } else { c_int::MAX };
        (*_q)[1] = if y < 0 { c_int::MIN } else { c_int::MAX };
        return -1;
    } else {
        if w < 0 {
            x = -x;
            y = -y;
            w = -w;
        }
        (*_q)[0] = qr_divround(x, w);
        (*_q)[1] = qr_divround(y, w);
    }
    0
}

/// Bit reading buffer for QR code data
///
/// Bit reading code adapted from libogg/libtheora
/// Portions (C) Xiph.Org Foundation 1994-2008, BSD-style license.
#[repr(C)]
pub struct qr_pack_buf {
    buf: *const c_uchar,
    endbyte: c_int,
    endbit: c_int,
    storage: c_int,
}

/// Initialize a pack buffer for reading bits from data
#[no_mangle]
pub unsafe extern "C" fn qr_pack_buf_init(
    _b: *mut qr_pack_buf,
    _data: *const c_uchar,
    _ndata: c_int,
) {
    (*_b).buf = _data;
    (*_b).storage = _ndata;
    (*_b).endbyte = 0;
    (*_b).endbit = 0;
}

/// Read bits from the pack buffer
///
/// Assumes 0 <= _bits <= 16
/// Returns the read value, or -1 if there aren't enough bits available
#[no_mangle]
pub unsafe extern "C" fn qr_pack_buf_read(_b: *mut qr_pack_buf, _bits: c_int) -> c_int {
    let m = 16 - _bits;
    let bits = _bits + (*_b).endbit;
    let d = (*_b).storage - (*_b).endbyte;

    if d <= 2 {
        // Not the main path
        if d * 8 < bits {
            (*_b).endbyte += bits >> 3;
            (*_b).endbit = bits & 7;
            return -1;
        }
        // Special case to avoid reading p[0] below, which might be past the end of
        // the buffer; also skips some useless accounting
        else if bits == 0 {
            return 0;
        }
    }

    let p = (*_b).buf.add((*_b).endbyte as usize);
    let mut ret = (c_uint::from(*p) << (8 + (*_b).endbit)) as c_uint;
    if bits > 8 {
        ret |= c_uint::from(*p.add(1)) << (*_b).endbit;
        if bits > 16 {
            ret |= c_uint::from(*p.add(2)) >> (8 - (*_b).endbit);
        }
    }
    (*_b).endbyte += bits >> 3;
    (*_b).endbit = bits & 7;
    ((ret & 0xFFFF) >> m) as c_int
}

/// Get the number of bits available to read from the pack buffer
#[no_mangle]
pub unsafe extern "C" fn qr_pack_buf_avail(_b: *const qr_pack_buf) -> c_int {
    (((*_b).storage - (*_b).endbyte) << 3) - (*_b).endbit
}

/// Calculate the number of codewords in a QR code of a given version
///
/// This is a compact calculation that avoids a lookup table.
/// Returns the total number of data and error correction codewords.
#[no_mangle]
pub unsafe extern "C" fn qr_code_ncodewords(_version: c_uint) -> c_int {
    if _version == 1 {
        return 26;
    }
    let nalign = (_version / 7) + 2;
    (((_version << 4) * (_version + 8) - (5 * nalign) * (5 * nalign - 2)
        + 36 * c_uint::from(_version < 7)
        + 83)
        >> 3) as c_int
}

/// Comparison function for sorting vertical finder lines
///
/// Used by qsort to sort vertical lines in ascending order by X coordinate,
/// with ties broken by Y coordinate.
#[no_mangle]
pub unsafe extern "C" fn qr_finder_vline_cmp(_a: *const c_void, _b: *const c_void) -> c_int {
    let a = _a as *const qr_finder_line;
    let b = _b as *const qr_finder_line;
    ((c_int::from((*a).pos[0] > (*b).pos[0]) - c_int::from((*a).pos[0] < (*b).pos[0])) << 1)
        + c_int::from((*a).pos[1] > (*b).pos[1])
        - c_int::from((*a).pos[1] < (*b).pos[1])
}

/// Comparison function for sorting finder centers
///
/// Sorts primarily by number of edge points (descending), then by Y coordinate
/// (ascending), then by X coordinate (ascending).
#[no_mangle]
pub unsafe extern "C" fn qr_finder_center_cmp(_a: *const c_void, _b: *const c_void) -> c_int {
    let a = _a as *const qr_finder_center;
    let b = _b as *const qr_finder_center;
    ((c_int::from((*b).nedge_pts > (*a).nedge_pts)
        - c_int::from((*b).nedge_pts < (*a).nedge_pts))
        << 2)
        + ((c_int::from((*a).pos[1] > (*b).pos[1]) - c_int::from((*a).pos[1] < (*b).pos[1]))
            << 1)
        + c_int::from((*a).pos[0] > (*b).pos[0])
        - c_int::from((*a).pos[0] < (*b).pos[0])
}

/// Clusters adjacent lines into groups that are large enough to be crossing a
/// finder pattern (relative to their length).
///
/// # Parameters
/// - `_clusters`: The buffer in which to store the clusters found.
/// - `_neighbors`: The buffer used to store the lists of lines in each cluster.
/// - `_lines`: The list of lines to cluster.
///   Horizontal lines must be sorted in ascending order by Y coordinate, with ties broken by X coordinate.
///   Vertical lines must be sorted in ascending order by X coordinate, with ties broken by Y coordinate.
/// - `_nlines`: The number of lines in the set of lines to cluster.
/// - `_v`: 0 for horizontal lines, or 1 for vertical lines.
///
/// # Returns
/// The number of clusters found.
#[no_mangle]
pub unsafe extern "C" fn qr_finder_cluster_lines(
    _clusters: *mut qr_finder_cluster,
    _neighbors: *mut *mut qr_finder_line,
    _lines: *mut qr_finder_line,
    _nlines: c_int,
    _v: c_int,
) -> c_int {
    // Allocate mark array to track which lines have been clustered
    let mark = calloc(_nlines as usize, 1) as *mut c_uchar;
    let mut neighbors = _neighbors;
    let mut nclusters = 0;

    for i in 0..(_nlines - 1) {
        if *mark.offset(i as isize) != 0 {
            continue;
        }

        let mut nneighbors = 1;
        *neighbors.offset(0) = _lines.offset(i as isize);
        let mut len = (*_lines.offset(i as isize)).len;

        for j in (i + 1).._nlines {
            if *mark.offset(j as isize) != 0 {
                continue;
            }

            let a = *neighbors.offset((nneighbors - 1) as isize);
            let b = _lines.offset(j as isize);

            // The clustering threshold is proportional to the size of the lines,
            // since minor noise in large areas can interrupt patterns more easily
            // at high resolutions.
            let thresh = ((*a).len + 7) >> 2;

            // Check if lines are too far apart perpendicular to their direction
            if ((*a).pos[(1 - _v) as usize] - (*b).pos[(1 - _v) as usize]).abs() > thresh {
                break;
            }

            // Check if lines are too far apart along their direction
            if ((*a).pos[_v as usize] - (*b).pos[_v as usize]).abs() > thresh {
                continue;
            }

            // Check if line ends are too far apart
            if ((*a).pos[_v as usize] + (*a).len - (*b).pos[_v as usize] - (*b).len).abs()
                > thresh
            {
                continue;
            }

            // Check beginning offset alignment
            if (*a).boffs > 0
                && (*b).boffs > 0
                && ((*a).pos[_v as usize] - (*a).boffs - (*b).pos[_v as usize] + (*b).boffs).abs()
                    > thresh
            {
                continue;
            }

            // Check ending offset alignment
            if (*a).eoffs > 0
                && (*b).eoffs > 0
                && ((*a).pos[_v as usize] + (*a).len + (*a).eoffs - (*b).pos[_v as usize]
                    - (*b).len
                    - (*b).eoffs)
                    .abs()
                    > thresh
            {
                continue;
            }

            *neighbors.offset(nneighbors as isize) = _lines.offset(j as isize);
            nneighbors += 1;
            len += (*b).len;
        }

        // We require at least three lines to form a cluster, which eliminates a
        // large number of false positives, saving considerable decoding time.
        // This should still be sufficient for 1-pixel codes with no noise.
        if nneighbors < 3 {
            continue;
        }

        // The expected number of lines crossing a finder pattern is equal to their
        // average length.
        // We accept the cluster if size is at least 1/3 their average length (this
        // is a very small threshold, but was needed for some test images).
        len = ((len << 1) + nneighbors) / (nneighbors << 1);
        if nneighbors * (5 << QR_FINDER_SUBPREC) >= len {
            (*_clusters.offset(nclusters as isize)).lines = neighbors;
            (*_clusters.offset(nclusters as isize)).nlines = nneighbors;

            for j in 0..nneighbors {
                let line_offset = (*neighbors.offset(j as isize)).offset_from(_lines);
                *mark.offset(line_offset) = 1;
            }

            neighbors = neighbors.offset(nneighbors as isize);
            nclusters += 1;
        }
    }

    free(mark as *mut c_void);
    nclusters
}

/// Adds the coordinates of the edge points from the lines contained in the
/// given list of clusters to the list of edge points for a finder center.
///
/// Only the edge point position is initialized.
/// The edge label and extent are set by qr_finder_edge_pts_aff_classify()
/// or qr_finder_edge_pts_hom_classify().
///
/// # Parameters
/// - `_edge_pts`: The buffer in which to store the edge points.
/// - `_nedge_pts`: The current number of edge points in the buffer.
/// - `_neighbors`: The list of clusters.
/// - `_nneighbors`: The number of clusters in the list.
/// - `_v`: 0 for horizontal lines and 1 for vertical lines.
///
/// # Returns
/// The new total number of edge points.
#[no_mangle]
pub unsafe extern "C" fn qr_finder_edge_pts_fill(
    _edge_pts: *mut qr_finder_edge_pt,
    mut _nedge_pts: c_int,
    _neighbors: *mut *mut qr_finder_cluster,
    _nneighbors: c_int,
    _v: c_int,
) -> c_int {
    for i in 0.._nneighbors {
        let c = *_neighbors.offset(i as isize);
        for j in 0..(*c).nlines {
            let l = *(*c).lines.offset(j as isize);

            // Add beginning offset edge point if present
            if (*l).boffs > 0 {
                (*_edge_pts.offset(_nedge_pts as isize)).pos[0] = (*l).pos[0];
                (*_edge_pts.offset(_nedge_pts as isize)).pos[1] = (*l).pos[1];
                (*_edge_pts.offset(_nedge_pts as isize)).pos[_v as usize] -= (*l).boffs;
                _nedge_pts += 1;
            }

            // Add ending offset edge point if present
            if (*l).eoffs > 0 {
                (*_edge_pts.offset(_nedge_pts as isize)).pos[0] = (*l).pos[0];
                (*_edge_pts.offset(_nedge_pts as isize)).pos[1] = (*l).pos[1];
                (*_edge_pts.offset(_nedge_pts as isize)).pos[_v as usize] += (*l).len + (*l).eoffs;
                _nedge_pts += 1;
            }
        }
    }
    _nedge_pts
}

/// Finds horizontal clusters that cross corresponding vertical clusters,
/// presumably corresponding to a finder center.
///
/// # Parameters
/// - `_centers`: The buffer in which to store putative finder centers.
/// - `_edge_pts`: The buffer to use for the edge point lists for each finder center.
/// - `_hclusters`: The clusters of horizontal lines crossing finder patterns.
/// - `_nhclusters`: The number of horizontal line clusters.
/// - `_vclusters`: The clusters of vertical lines crossing finder patterns.
/// - `_nvclusters`: The number of vertical line clusters.
///
/// # Returns
/// The number of putative finder centers.
#[no_mangle]
pub unsafe extern "C" fn qr_finder_find_crossings(
    _centers: *mut qr_finder_center,
    mut _edge_pts: *mut qr_finder_edge_pt,
    _hclusters: *mut qr_finder_cluster,
    _nhclusters: c_int,
    _vclusters: *mut qr_finder_cluster,
    _nvclusters: c_int,
) -> c_int {
    let hneighbors = malloc((_nhclusters as usize) * size_of::<*mut qr_finder_cluster>())
        as *mut *mut qr_finder_cluster;
    let vneighbors = malloc((_nvclusters as usize) * size_of::<*mut qr_finder_cluster>())
        as *mut *mut qr_finder_cluster;
    let hmark = calloc(_nhclusters as usize, 1) as *mut c_uchar;
    let vmark = calloc(_nvclusters as usize, 1) as *mut c_uchar;
    let mut ncenters = 0;

    // TODO: This may need some re-working.
    // We should be finding groups of clusters such that _all_ horizontal lines in
    // _all_ horizontal clusters in the group cross _all_ vertical lines in _all_
    // vertical clusters in the group.
    // This is equivalent to finding the maximum bipartite clique in the
    // connectivity graph, which requires linear programming to solve efficiently.
    // In principle, that is easy to do, but a realistic implementation without
    // floating point is a lot of work (and computationally expensive).
    // Right now we are relying on a sufficient border around the finder patterns
    // to prevent false positives.

    for i in 0.._nhclusters {
        if *hmark.offset(i as isize) != 0 {
            continue;
        }

        let a = *(*_hclusters.offset(i as isize))
            .lines
            .offset(((*_hclusters.offset(i as isize)).nlines >> 1) as isize);
        let mut y = 0;
        let mut nvneighbors = 0;

        // Find vertical clusters that cross this horizontal cluster
        for j in 0.._nvclusters {
            if *vmark.offset(j as isize) != 0 {
                continue;
            }

            let b = *(*_vclusters.offset(j as isize))
                .lines
                .offset(((*_vclusters.offset(j as isize)).nlines >> 1) as isize);

            if qr_finder_lines_are_crossing(a, b) != 0 {
                *vmark.offset(j as isize) = 1;
                y += ((*b).pos[1] << 1) + (*b).len;
                if (*b).boffs > 0 && (*b).eoffs > 0 {
                    y += (*b).eoffs - (*b).boffs;
                }
                *vneighbors.offset(nvneighbors as isize) = _vclusters.offset(j as isize);
                nvneighbors += 1;
            }
        }

        if nvneighbors > 0 {
            let mut x = ((*a).pos[0] << 1) + (*a).len;
            if (*a).boffs > 0 && (*a).eoffs > 0 {
                x += (*a).eoffs - (*a).boffs;
            }
            *hneighbors.offset(0) = _hclusters.offset(i as isize);
            let mut nhneighbors = 1;

            let j_mid = nvneighbors >> 1;
            let b = *(*(*vneighbors.offset(j_mid as isize)))
                .lines
                .offset(((*(*vneighbors.offset(j_mid as isize))).nlines >> 1) as isize);

            // Find additional horizontal clusters that cross the vertical clusters
            for j in (i + 1).._nhclusters {
                if *hmark.offset(j as isize) != 0 {
                    continue;
                }

                let a = *(*_hclusters.offset(j as isize))
                    .lines
                    .offset(((*_hclusters.offset(j as isize)).nlines >> 1) as isize);

                if qr_finder_lines_are_crossing(a, b) != 0 {
                    *hmark.offset(j as isize) = 1;
                    x += ((*a).pos[0] << 1) + (*a).len;
                    if (*a).boffs > 0 && (*a).eoffs > 0 {
                        x += (*a).eoffs - (*a).boffs;
                    }
                    *hneighbors.offset(nhneighbors as isize) = _hclusters.offset(j as isize);
                    nhneighbors += 1;
                }
            }

            let c = _centers.offset(ncenters as isize);
            ncenters += 1;

            (*c).pos[0] = (x + nhneighbors) / (nhneighbors << 1);
            (*c).pos[1] = (y + nvneighbors) / (nvneighbors << 1);
            (*c).edge_pts = _edge_pts;

            let mut nedge_pts =
                qr_finder_edge_pts_fill(_edge_pts, 0, hneighbors, nhneighbors, 0);
            nedge_pts = qr_finder_edge_pts_fill(_edge_pts, nedge_pts, vneighbors, nvneighbors, 1);

            (*c).nedge_pts = nedge_pts;
            _edge_pts = _edge_pts.offset(nedge_pts as isize);
        }
    }

    free(vmark as *mut c_void);
    free(hmark as *mut c_void);
    free(vneighbors as *mut c_void);
    free(hneighbors as *mut c_void);

    // Sort the centers by decreasing numbers of edge points
    qsort(
        _centers as *mut c_void,
        ncenters as size_t,
        size_of::<qr_finder_center>(),
        Some(qr_finder_center_cmp),
    );

    ncenters
}

/// Locate QR finder pattern centers from scanned lines
///
/// Clusters horizontal and vertical lines that cross finder patterns,
/// then locates the centers where horizontal and vertical clusters intersect.
///
/// # Safety
/// This function is unsafe because it:
/// - Dereferences raw pointers
/// - Allocates and manages memory using C allocation functions
/// - Assumes the qr_reader structure has valid finder_lines arrays
#[no_mangle]
pub unsafe extern "C" fn qr_finder_centers_locate(
    _centers: *mut *mut qr_finder_center,
    _edge_pts: *mut *mut qr_finder_edge_pt,
    reader: *mut qr_reader,
    _width: c_int,
    _height: c_int,
) -> c_int {
    let hlines = (*reader).finder_lines[0].lines;
    let nhlines = (*reader).finder_lines[0].nlines;
    let vlines = (*reader).finder_lines[1].lines;
    let nvlines = (*reader).finder_lines[1].nlines;

    // Cluster the detected lines
    let hneighbors =
        malloc((nhlines as usize) * size_of::<*mut qr_finder_line>()) as *mut *mut qr_finder_line;

    // We require more than one line per cluster, so there are at most nhlines/2
    let hclusters = malloc(((nhlines >> 1) as usize) * size_of::<qr_finder_cluster>())
        as *mut qr_finder_cluster;

    let nhclusters = qr_finder_cluster_lines(hclusters, hneighbors, hlines, nhlines, 0);

    // We need vertical lines to be sorted by X coordinate, with ties broken by Y
    // coordinate, for clustering purposes.
    // We scan the image in the opposite order for cache efficiency, so sort the
    // lines we found here.
    qsort(
        vlines as *mut c_void,
        nvlines as size_t,
        size_of::<qr_finder_line>(),
        Some(qr_finder_vline_cmp),
    );

    let vneighbors =
        malloc((nvlines as usize) * size_of::<*mut qr_finder_line>()) as *mut *mut qr_finder_line;

    // We require more than one line per cluster, so there are at most nvlines/2
    let vclusters = malloc(((nvlines >> 1) as usize) * size_of::<qr_finder_cluster>())
        as *mut qr_finder_cluster;

    let nvclusters = qr_finder_cluster_lines(vclusters, vneighbors, vlines, nvlines, 1);

    // Find line crossings among the clusters
    let ncenters = if nhclusters >= 3 && nvclusters >= 3 {
        let mut nedge_pts = 0;
        for i in 0..nhclusters {
            nedge_pts += (*hclusters.offset(i as isize)).nlines;
        }
        for i in 0..nvclusters {
            nedge_pts += (*vclusters.offset(i as isize)).nlines;
        }
        nedge_pts <<= 1;

        let edge_pts =
            malloc((nedge_pts as usize) * size_of::<qr_finder_edge_pt>()) as *mut qr_finder_edge_pt;
        let centers = malloc((qr_mini(nhclusters, nvclusters) as usize) * size_of::<qr_finder_center>())
            as *mut qr_finder_center;

        let ncenters = qr_finder_find_crossings(
            centers, edge_pts, hclusters, nhclusters, vclusters, nvclusters,
        );

        *_centers = centers;
        *_edge_pts = edge_pts;
        ncenters
    } else {
        0
    };

    free(vclusters as *mut c_void);
    free(vneighbors as *mut c_void);
    free(hclusters as *mut c_void);
    free(hneighbors as *mut c_void);

    ncenters
}

/// Mark a given rectangular region as belonging to the function pattern
///
/// Function patterns include finder patterns, timing patterns, and alignment patterns.
/// We store bits column-wise, since that's how they're read out of the grid.
#[no_mangle]
pub unsafe extern "C" fn qr_sampling_grid_fp_mask_rect(
    _grid: *mut qr_sampling_grid,
    _dim: c_int,
    _u: c_int,
    _v: c_int,
    _w: c_int,
    _h: c_int,
) {
    let stride = (_dim + QR_INT_BITS - 1) >> QR_INT_LOGBITS;
    for j in _u..(_u + _w) {
        for i in _v..(_v + _h) {
            *(*_grid).fpmask.offset((j * stride + (i >> QR_INT_LOGBITS)) as isize) |=
                1 << (i & (QR_INT_BITS - 1));
        }
    }
}

/// Determine if a given grid location is inside a function pattern
///
/// Returns 1 if the location (_u, _v) is part of a function pattern, 0 otherwise.
#[no_mangle]
pub unsafe extern "C" fn qr_sampling_grid_is_in_fp(
    _grid: *const qr_sampling_grid,
    _dim: c_int,
    _u: c_int,
    _v: c_int,
) -> c_int {
    ((*(*_grid).fpmask.offset(
        (_u * ((_dim + QR_INT_BITS - 1) >> QR_INT_LOGBITS) + (_v >> QR_INT_LOGBITS)) as isize,
    )) >> (_v & (QR_INT_BITS - 1))
        & 1) as c_int
}

/// The spacing between alignment patterns after the second for versions >= 7
///
/// We could compact this more, but the code to access it would eliminate the gains.
#[no_mangle]
pub static QR_ALIGNMENT_SPACING: [c_uchar; 34] = [
    16, 18, 20, 22, 24, 26, 28, 20, 22, 24, 24, 26, 28, 28, 22, 24, 24, 26, 26, 28, 28, 24, 24,
    26, 26, 26, 28, 28, 24, 26, 26, 26, 28, 28,
];

/// Bulk data for the number of parity bytes per Reed-Solomon block
#[no_mangle]
pub static QR_RS_NPAR_VALS: [c_uchar; 71] = [
    // [ 0]
    7, 10, 13, 17, // [ 4]
    10, 16, 22, 28, 26, 26, 26, 22, 24, 22, 22, 26, 24, 18, 22, // [19]
    15, 26, 18, 22, 24, 30, 24, 20, 24, // [28]
    18, 16, 24, 28, 28, 28, 28, 30, 24, // [37]
    20, 18, 18, 26, 24, 28, 24, 30, 26, 28, 28, 26, 28, 30, 30, 22, 20, 24, // [55]
    20, 18, 26, 16, // [59]
    20, 30, 28, 24, 22, 26, 28, 26, 30, 28, 30, 30,
];

/// An offset into QR_RS_NPAR_VALS for each version that gives the number of
/// parity bytes per Reed-Solomon block for each error correction level
#[no_mangle]
pub static QR_RS_NPAR_OFFS: [c_uchar; 40] = [
    0, 4, 19, 55, 15, 28, 37, 12, 51, 39, 59, 62, 10, 24, 22, 41, 31, 44, 7, 65, 47, 33, 67, 67,
    48, 32, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67,
];

/// The number of Reed-Solomon blocks for each version and error correction level
#[no_mangle]
pub static QR_RS_NBLOCKS: [[c_uchar; 4]; 40] = [
    [1, 1, 1, 1],
    [1, 1, 1, 1],
    [1, 1, 2, 2],
    [1, 2, 2, 4],
    [1, 2, 4, 4],
    [2, 4, 4, 4],
    [2, 4, 6, 5],
    [2, 4, 6, 6],
    [2, 5, 8, 8],
    [4, 5, 8, 8],
    [4, 5, 8, 11],
    [4, 8, 10, 11],
    [4, 9, 12, 16],
    [4, 9, 16, 16],
    [6, 10, 12, 18],
    [6, 10, 17, 16],
    [6, 11, 16, 19],
    [6, 13, 18, 21],
    [7, 14, 21, 25],
    [8, 16, 20, 25],
    [8, 17, 23, 25],
    [9, 17, 23, 34],
    [9, 18, 25, 30],
    [10, 20, 27, 32],
    [12, 21, 29, 35],
    [12, 23, 34, 37],
    [12, 25, 34, 40],
    [13, 26, 35, 42],
    [14, 28, 38, 45],
    [15, 29, 40, 48],
    [16, 31, 43, 51],
    [17, 33, 45, 54],
    [18, 35, 48, 57],
    [19, 37, 51, 60],
    [19, 38, 53, 63],
    [20, 40, 56, 66],
    [21, 43, 59, 70],
    [22, 45, 62, 74],
    [24, 47, 65, 77],
    [25, 49, 68, 81],
];

pub fn qr_cmp_edge_pt(a: &qr_finder_edge_pt, b: &qr_finder_edge_pt) -> Ordering {
    match ((c_int::from(a.edge > b.edge) - c_int::from(a.edge < b.edge)) << 1)
        + c_int::from(a.extent > b.extent)
        - c_int::from(a.extent < b.extent)
    {
        ..0 => Ordering::Less,
        0 => Ordering::Equal,
        1.. => Ordering::Greater,
    }
}

/// Computes the index of the edge each edge point belongs to, and its (signed)
/// distance along the corresponding axis from the center of the finder pattern
/// (in the square domain).
/// The resulting list of edge points is sorted by edge index, with ties broken
/// by extent.
#[no_mangle]
pub unsafe extern "C" fn qr_finder_edge_pts_aff_classify(_f: *mut qr_finder, _aff: *const qr_aff) {
    let c = (*_f).c;
    for item in (*_f).nedge_pts.iter_mut() {
        *item = 0;
    }
    for i in 0..(*c).nedge_pts {
        let mut q = qr_point::default();
        let edge_pt = (*c).edge_pts.add(i as usize);
        qr_aff_unproject(
            (&mut q) as *mut _,
            _aff,
            (*edge_pt).pos[0],
            (*edge_pt).pos[1],
        );
        qr_point_translate(q.as_mut_ptr(), -(*_f).o[0], -(*_f).o[1]);
        let d = c_int::from(q[1].abs() > q[0].abs());
        let e = d << 1 | c_int::from(q[d as usize] >= 0);
        (*_f).nedge_pts[e as usize] += 1;
        (*edge_pt).edge = e;
        (*edge_pt).extent = q[d as usize];
    }

    let edge_pts = std::slice::from_raw_parts_mut((*c).edge_pts, (*c).nedge_pts as usize);
    edge_pts.sort_by(qr_cmp_edge_pt);
    (*_f).edge_pts[0] = (*c).edge_pts;
    for e in 1..(*_f).edge_pts.len() {
        (*_f).edge_pts[e] = (*_f).edge_pts[e - 1].add((*_f).nedge_pts[e - 1] as usize);
    }
}

/// Computes the index of the edge each edge point belongs to, and its (signed)
/// distance along the corresponding axis from the center of the finder pattern
/// (in the square domain), using homography projection.
///
/// This is similar to `qr_finder_edge_pts_aff_classify` but uses full homography
/// projection which can fail (return infinity). Failed projections are marked
/// with edge=4.
///
/// The resulting list of edge points is sorted by edge index, with ties broken
/// by extent.
#[no_mangle]
pub unsafe extern "C" fn qr_finder_edge_pts_hom_classify(_f: *mut qr_finder, _hom: *const qr_hom) {
    let c = (*_f).c;
    for item in (*_f).nedge_pts.iter_mut() {
        *item = 0;
    }
    for i in 0..(*c).nedge_pts {
        let mut q = qr_point::default();
        let edge_pt = (*c).edge_pts.add(i as usize);

        if qr_hom_unproject(
            (&mut q) as *mut _,
            _hom,
            (*edge_pt).pos[0],
            (*edge_pt).pos[1],
        ) >= 0
        {
            // Successful projection
            qr_point_translate(q.as_mut_ptr(), -(*_f).o[0], -(*_f).o[1]);
            let d = c_int::from(q[1].abs() > q[0].abs());
            let e = d << 1 | c_int::from(q[d as usize] >= 0);
            (*_f).nedge_pts[e as usize] += 1;
            (*edge_pt).edge = e;
            (*edge_pt).extent = q[d as usize];
        } else {
            // Projection failed (went to infinity)
            (*edge_pt).edge = 4;
            (*edge_pt).extent = q[0];
        }
    }

    let edge_pts = std::slice::from_raw_parts_mut((*c).edge_pts, (*c).nedge_pts as usize);
    edge_pts.sort_by(qr_cmp_edge_pt);
    (*_f).edge_pts[0] = (*c).edge_pts;
    for e in 1..(*_f).edge_pts.len() {
        (*_f).edge_pts[e] = (*_f).edge_pts[e - 1].add((*_f).nedge_pts[e - 1] as usize);
    }
}

#[no_mangle]
pub unsafe extern "C" fn qr_finder_quick_crossing_check(
    _img: *const c_uchar,
    _width: c_int,
    _height: c_int,
    _x0: c_int,
    _y0: c_int,
    _x1: c_int,
    _y1: c_int,
    _v: c_int,
) -> c_int {
    // The points must be inside the image, and have a !_v:_v:!_v pattern.
    // We don't scan the whole line initially, but quickly reject if the endpoints
    // aren't !_v, or the midpoint isn't _v.
    // If either end point is out of the image, or we don't encounter a _v pixel,
    // we return a negative value, indicating the region should be considered
    // empty.
    // Otherwise, we return a positive value to indicate it is non-empty.
    if _x0 < 0
        || _x0 >= _width
        || _y0 < 0
        || _y0 >= _height
        || _x1 < 0
        || _x1 >= _width
        || _y1 < 0
        || _y1 >= _height
    {
        return -1;
    }

    if (c_int::from(*_img.add((_y0 * _width + _x0) as usize) == 0)) != _v
        || (c_int::from(*_img.add((_y1 * _width + _x1) as usize) == 0)) != _v
    {
        return 1;
    }
    if (c_int::from(*_img.add((((_y0 + _y1) >> 1) * _width + ((_x0 + _x1) >> 1)) as usize) == 0))
        == _v
    {
        return -1;
    }
    0
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
    aff: *const qr_aff,
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
    let shift = c_int::max(
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

/// Initialize a homography cell for mapping between code and image space
///
/// This computes a homographic transformation from a quadrilateral in code space
/// (u0,v0)-(u1,v1)-(u2,v2)-(u3,v3) to a quadrilateral in image space
/// (x0,y0)-(x1,y1)-(x2,y2)-(x3,y3).
///
/// The transformation handles both affine and projective distortion, with careful
/// attention to numerical stability through dynamic range scaling.
#[no_mangle]
pub unsafe extern "C" fn qr_hom_cell_init(
    cell: *mut qr_hom_cell,
    u0: c_int,
    v0: c_int,
    u1: c_int,
    v1: c_int,
    u2: c_int,
    v2: c_int,
    u3: c_int,
    v3: c_int,
    x0: c_int,
    y0: c_int,
    x1: c_int,
    y1: c_int,
    x2: c_int,
    y2: c_int,
    x3: c_int,
    y3: c_int,
) {
    // Compute deltas for source points (code space)
    let du10 = u1 - u0;
    let du20 = u2 - u0;
    let du30 = u3 - u0;
    let du31 = u3 - u1;
    let du32 = u3 - u2;
    let dv10 = v1 - v0;
    let dv20 = v2 - v0;
    let dv30 = v3 - v0;
    let dv31 = v3 - v1;
    let dv32 = v3 - v2;

    // Compute coefficients of forward transform from unit square to source configuration
    let a20 = du32 * dv10 - du10 * dv32;
    let a21 = du20 * dv31 - du31 * dv20;
    let a22 = if a20 != 0 || a21 != 0 {
        // Non-affine arrangement
        du32 * dv31 - du31 * dv32
    } else {
        // Affine arrangement - no scaling needed for larger dynamic range
        1
    };
    let a00 = du10 * (a20 + a22);
    let a01 = du20 * (a21 + a22);
    let a10 = dv10 * (a20 + a22);
    let a11 = dv20 * (a21 + a22);

    // Compute inverse transform
    let i00_full = a11 * a22;
    let i01_full = -a01 * a22;
    let i10_full = -a10 * a22;
    let i11_full = a00 * a22;
    let i20_full = a10 * a21 - a11 * a20;
    let i21_full = a01 * a20 - a00 * a21;
    let i22 = a00 * a11 - a01 * a10;

    // Invert coefficients: divide i22 by all others
    // Zero means "infinity", used to signal when divisor is zero
    // QR_FLIPSIGNI flips sign of result if original value was negative
    let i00 = if i00_full != 0 {
        let result = qr_divround(i22, i00_full.abs());
        if i00_full < 0 {
            -result
        } else {
            result
        }
    } else {
        0
    };
    let i01 = if i01_full != 0 {
        let result = qr_divround(i22, i01_full.abs());
        if i01_full < 0 {
            -result
        } else {
            result
        }
    } else {
        0
    };
    let i10 = if i10_full != 0 {
        let result = qr_divround(i22, i10_full.abs());
        if i10_full < 0 {
            -result
        } else {
            result
        }
    } else {
        0
    };
    let i11 = if i11_full != 0 {
        let result = qr_divround(i22, i11_full.abs());
        if i11_full < 0 {
            -result
        } else {
            result
        }
    } else {
        0
    };
    let i20 = if i20_full != 0 {
        let result = qr_divround(i22, i20_full.abs());
        if i20_full < 0 {
            -result
        } else {
            result
        }
    } else {
        0
    };
    let i21 = if i21_full != 0 {
        let result = qr_divround(i22, i21_full.abs());
        if i21_full < 0 {
            -result
        } else {
            result
        }
    } else {
        0
    };

    // Compute map from unit square to image
    let dx10 = x1 - x0;
    let dx20 = x2 - x0;
    let dx30 = x3 - x0;
    let dx31 = x3 - x1;
    let dx32 = x3 - x2;
    let dy10 = y1 - y0;
    let dy20 = y2 - y0;
    let dy30 = y3 - y0;
    let dy31 = y3 - y1;
    let dy32 = y3 - y2;
    let a20 = dx32 * dy10 - dx10 * dy32;
    let a21 = dx20 * dy31 - dx31 * dy20;
    let a22 = dx32 * dy31 - dx31 * dy32;

    // Figure out if we need to downscale
    let b0 =
        qr_ilog(c_int::max(dx10.abs(), dy10.abs()) as u32) + qr_ilog((a20 + a22).unsigned_abs());
    let b1 =
        qr_ilog(c_int::max(dx20.abs(), dy20.abs()) as u32) + qr_ilog((a21 + a22).unsigned_abs());
    let b2 = qr_ilog(c_int::max(c_int::max(a20.abs(), a21.abs()), a22.abs()) as u32);
    let shift = c_int::max(
        0,
        c_int::max(c_int::max(b0, b1), b2) - (QR_INT_BITS - 3 - QR_ALIGN_SUBPREC),
    );
    let round = (1i64 << shift) >> 1;

    // Compute final coefficients of forward transform
    let a00 = qr_fixmul(dx10, a20 + a22, round, shift);
    let a01 = qr_fixmul(dx20, a21 + a22, round, shift);
    let a10 = qr_fixmul(dy10, a20 + a22, round, shift);
    let a11 = qr_fixmul(dy20, a21 + a22, round, shift);

    // Compose the two transforms (divide by inverted coefficients)
    (*cell).fwd[0][0] = (if i00 != 0 { qr_divround(a00, i00) } else { 0 })
        + (if i10 != 0 { qr_divround(a01, i10) } else { 0 });
    (*cell).fwd[0][1] = (if i01 != 0 { qr_divround(a00, i01) } else { 0 })
        + (if i11 != 0 { qr_divround(a01, i11) } else { 0 });
    (*cell).fwd[1][0] = (if i00 != 0 { qr_divround(a10, i00) } else { 0 })
        + (if i10 != 0 { qr_divround(a11, i10) } else { 0 });
    (*cell).fwd[1][1] = (if i01 != 0 { qr_divround(a10, i01) } else { 0 })
        + (if i11 != 0 { qr_divround(a11, i11) } else { 0 });
    (*cell).fwd[2][0] = ((if i00 != 0 { qr_divround(a20, i00) } else { 0 })
        + (if i10 != 0 { qr_divround(a21, i10) } else { 0 })
        + (if i20 != 0 { qr_divround(a22, i20) } else { 0 })
        + round as c_int) as c_int
        >> shift;
    (*cell).fwd[2][1] = ((if i01 != 0 { qr_divround(a20, i01) } else { 0 })
        + (if i11 != 0 { qr_divround(a21, i11) } else { 0 })
        + (if i21 != 0 { qr_divround(a22, i21) } else { 0 })
        + round as c_int) as c_int
        >> shift;
    (*cell).fwd[2][2] = ((a22 + round as c_int) >> shift) as c_int;

    // Compute offsets to distribute rounding error over whole range
    // (instead of concentrating it in the (u3,v3) corner)
    let mut x = (*cell).fwd[0][0] * du10 + (*cell).fwd[0][1] * dv10;
    let mut y = (*cell).fwd[1][0] * du10 + (*cell).fwd[1][1] * dv10;
    let mut w = (*cell).fwd[2][0] * du10 + (*cell).fwd[2][1] * dv10 + (*cell).fwd[2][2];
    let mut a02 = dx10 * w - x;
    let mut a12 = dy10 * w - y;

    x = (*cell).fwd[0][0] * du20 + (*cell).fwd[0][1] * dv20;
    y = (*cell).fwd[1][0] * du20 + (*cell).fwd[1][1] * dv20;
    w = (*cell).fwd[2][0] * du20 + (*cell).fwd[2][1] * dv20 + (*cell).fwd[2][2];
    a02 += dx20 * w - x;
    a12 += dy20 * w - y;

    x = (*cell).fwd[0][0] * du30 + (*cell).fwd[0][1] * dv30;
    y = (*cell).fwd[1][0] * du30 + (*cell).fwd[1][1] * dv30;
    w = (*cell).fwd[2][0] * du30 + (*cell).fwd[2][1] * dv30 + (*cell).fwd[2][2];
    a02 += dx30 * w - x;
    a12 += dy30 * w - y;

    (*cell).fwd[0][2] = (a02 + 2) >> 2;
    (*cell).fwd[1][2] = (a12 + 2) >> 2;
    (*cell).x0 = x0;
    (*cell).y0 = y0;
    (*cell).u0 = u0;
    (*cell).v0 = v0;
}

/// Get a bit from a binarized image
///
/// Samples a pixel from the binarized image, with coordinates in QR_FINDER_SUBPREC
/// subpixel units. Clamps coordinates to valid image bounds.
#[no_mangle]
pub unsafe extern "C" fn qr_img_get_bit(
    img: *const u8,
    width: c_int,
    height: c_int,
    mut x: c_int,
    mut y: c_int,
) -> c_int {
    x >>= QR_FINDER_SUBPREC;
    y >>= QR_FINDER_SUBPREC;
    let y_clamped = y.clamp(0, height - 1);
    let x_clamped = x.clamp(0, width - 1);
    let idx = y_clamped * width + x_clamped;
    if *img.offset(idx as isize) != 0 {
        1
    } else {
        0
    }
}

/// Finish a partial projection, converting from homogeneous coordinates to the
/// normal 2-D representation.
/// In loops, we can avoid many multiplies by computing the homogeneous _x, _y,
/// and _w incrementally, but we cannot avoid the divisions, done here.*/
#[no_mangle]
pub unsafe extern "C" fn qr_hom_cell_fproject(
    _p: *mut qr_point,
    _cell: *const qr_hom_cell,
    mut _x: c_int,
    mut _y: c_int,
    mut _w: c_int,
) {
    if _w == 0 {
        (*_p)[0] = if _x < 0 { c_int::MIN } else { c_int::MAX };
        (*_p)[1] = if _y < 0 { c_int::MIN } else { c_int::MAX };
    } else {
        if _w < 0 {
            _x = -_x;
            _y = -_y;
            _w = -_w;
        }
        (*_p)[0] = qr_divround(_x, _w) + (*_cell).x0;
        (*_p)[1] = qr_divround(_y, _w) + (*_cell).y0;
    }
}

#[no_mangle]
pub unsafe extern "C" fn qr_hom_cell_project(
    _p: *mut qr_point,
    _cell: *const qr_hom_cell,
    mut _u: c_int,
    mut _v: c_int,
    _res: c_int,
) {
    _u -= (*_cell).u0 << _res;
    _v -= (*_cell).v0 << _res;
    qr_hom_cell_fproject(
        _p,
        _cell,
        (*_cell).fwd[0][0] * _u + (*_cell).fwd[0][1] * _v + ((*_cell).fwd[0][2] << _res),
        (*_cell).fwd[1][0] * _u + (*_cell).fwd[1][1] * _v + ((*_cell).fwd[1][2] << _res),
        (*_cell).fwd[2][0] * _u + (*_cell).fwd[2][1] * _v + ((*_cell).fwd[2][2] << _res),
    );
}

/// Locate the crossing of a finder or alignment pattern along a line
///
/// Uses Bresenham's algorithm to trace along the line and find the exact
/// transitions from !_v to _v and back. Returns the midpoint of the segment.
///
/// Returns 0 on success, -1 if no crossing found.
#[no_mangle]
pub unsafe extern "C" fn qr_finder_locate_crossing(
    img: *const u8,
    width: c_int,
    _height: c_int,
    x0: c_int,
    y0: c_int,
    x1: c_int,
    y1: c_int,
    v: c_int,
    p: *mut qr_point,
) -> c_int {
    let mut x0_pos = [x0, y0];
    let mut x1_pos = [x1, y1];
    let dx = [(x1 - x0).abs(), (y1 - y0).abs()];
    let steep = if dx[1] > dx[0] { 1 } else { 0 };
    let step = [if x0 < x1 { 1 } else { -1 }, if y0 < y1 { 1 } else { -1 }];
    let derr = dx[1 - steep];

    // Find the first crossing from !v to v
    let mut err = 0;
    loop {
        // If we make it all the way to the other side, there's no crossing
        if x0_pos[steep] == x1_pos[steep] {
            return -1;
        }
        x0_pos[steep] += step[steep];
        err += derr;
        if (err << 1) > dx[steep] {
            x0_pos[1 - steep] += step[1 - steep];
            err -= dx[steep];
        }
        let pixel = *img.offset((x0_pos[1] * width + x0_pos[0]) as isize);
        if ((pixel == 0) as c_int) != v {
            break;
        }
    }

    // Find the last crossing from v to !v
    err = 0;
    loop {
        if x0_pos[steep] == x1_pos[steep] {
            break;
        }
        x1_pos[steep] -= step[steep];
        err += derr;
        if (err << 1) > dx[steep] {
            x1_pos[1 - steep] -= step[1 - steep];
            err -= dx[steep];
        }
        let pixel = *img.offset((x1_pos[1] * width + x1_pos[0]) as isize);
        if ((pixel == 0) as c_int) != v {
            break;
        }
    }

    // Return the midpoint of the v segment
    (*p)[0] = ((x0_pos[0] + x1_pos[0] + 1) << QR_FINDER_SUBPREC) >> 1;
    (*p)[1] = ((x0_pos[1] + x1_pos[1] + 1) << QR_FINDER_SUBPREC) >> 1;
    0
}

/// Fetch a 5x5 alignment pattern and return it as a 25-bit value
///
/// Samples a 5x5 grid of pixels offset by (x0, y0) from the template positions
/// in p, returning the result as a packed 25-bit unsigned integer.
unsafe fn qr_alignment_pattern_fetch(
    p: &[[qr_point; 5]; 5],
    x0: c_int,
    y0: c_int,
    img: *const u8,
    width: c_int,
    height: c_int,
) -> c_uint {
    let dx = x0 - p[2][2][0];
    let dy = y0 - p[2][2][1];
    let mut v = 0u32;
    let mut k = 0;
    for pi in p {
        for pij in pi {
            v |= (qr_img_get_bit(img, width, height, pij[0] + dx, pij[1] + dy) as c_uint) << k;
            k += 1;
        }
    }
    v
}

/// Search for an alignment pattern near the given location
///
/// Searches for an alignment pattern around the specified grid coordinates (u, v)
/// within a radius of r modules. Returns 0 on success, -1 if the best match is too poor.
#[no_mangle]
pub unsafe extern "C" fn qr_alignment_pattern_search(
    p: *mut qr_point,
    cell: *const qr_hom_cell,
    _u: c_int,
    _v: c_int,
    r: c_int,
    img: *const u8,
    width: c_int,
    height: c_int,
) -> c_int {
    let mut pattern: [[qr_point; 5]; 5] = [[Default::default(); 5]; 5];
    let mut nc = [0 as c_int; 4];
    let mut c: [qr_point; 4] = [[0; 2]; 4];
    let mut pc: qr_point = [0; 2];

    // Build up a basic template using cell to control shape and scale
    let u = (_u - 2) - (*cell).u0;
    let v = (_v - 2) - (*cell).v0;
    let mut x0 = (*cell).fwd[0][0] * u + (*cell).fwd[0][1] * v + (*cell).fwd[0][2];
    let mut y0 = (*cell).fwd[1][0] * u + (*cell).fwd[1][1] * v + (*cell).fwd[1][2];
    let mut w0 = (*cell).fwd[2][0] * u + (*cell).fwd[2][1] * v + (*cell).fwd[2][2];
    let dxdu = (*cell).fwd[0][0];
    let dydu = (*cell).fwd[1][0];
    let dwdu = (*cell).fwd[2][0];
    let dxdv = (*cell).fwd[0][1];
    let dydv = (*cell).fwd[1][1];
    let dwdv = (*cell).fwd[2][1];

    for item in pattern.iter_mut() {
        let mut x = x0;
        let mut y = y0;
        let mut w = w0;
        for subitem in item.iter_mut() {
            qr_hom_cell_fproject(subitem, cell, x, y, w);
            x += dxdu;
            y += dydu;
            w += dwdu;
        }
        x0 += dxdv;
        y0 += dydv;
        w0 += dwdv;
    }

    let mut bestx = pattern[2][2][0];
    let mut besty = pattern[2][2][1];
    let mut best_match = qr_alignment_pattern_fetch(&pattern, bestx, besty, img, width, height);
    let mut best_dist = qr_hamming_dist(best_match, 0x1F8D63F, 25);

    if best_dist > 0 {
        let u = _u - (*cell).u0;
        let v = _v - (*cell).v0;
        let mut x =
            ((*cell).fwd[0][0] * u + (*cell).fwd[0][1] * v + (*cell).fwd[0][2]) << QR_ALIGN_SUBPREC;
        let mut y =
            ((*cell).fwd[1][0] * u + (*cell).fwd[1][1] * v + (*cell).fwd[1][2]) << QR_ALIGN_SUBPREC;
        let mut w =
            ((*cell).fwd[2][0] * u + (*cell).fwd[2][1] * v + (*cell).fwd[2][2]) << QR_ALIGN_SUBPREC;

        // Search an area at most r modules around the target location, in concentric squares
        for i in 1..(r << QR_ALIGN_SUBPREC) {
            let side_len = (i << 1) - 1;
            x -= dxdu + dxdv;
            y -= dydu + dydv;
            w -= dwdu + dwdv;

            for j in 0..(4 * side_len) {
                qr_hom_cell_fproject(&mut pc, cell, x, y, w);
                let match_val =
                    qr_alignment_pattern_fetch(&pattern, pc[0], pc[1], img, width, height);
                let dist = qr_hamming_dist(match_val, 0x1F8D63F, best_dist + 1);
                if dist < best_dist {
                    best_match = match_val;
                    best_dist = dist;
                    bestx = pc[0];
                    besty = pc[1];
                }

                let dir = if j < 2 * side_len {
                    if j >= side_len {
                        1
                    } else {
                        0
                    }
                } else if j >= 3 * side_len {
                    1
                } else {
                    0
                };

                if j < 2 * side_len {
                    x += (*cell).fwd[0][dir];
                    y += (*cell).fwd[1][dir];
                    w += (*cell).fwd[2][dir];
                } else {
                    x -= (*cell).fwd[0][dir];
                    y -= (*cell).fwd[1][dir];
                    w -= (*cell).fwd[2][dir];
                }

                if best_dist == 0 {
                    break;
                }
            }
            if best_dist == 0 {
                break;
            }
        }
    }

    // If the best result we got was sufficiently bad, reject the match
    if best_dist > 6 {
        (*p)[0] = pattern[2][2][0];
        (*p)[1] = pattern[2][2][1];
        return -1;
    }

    // Now try to get a more accurate location of the pattern center
    let dx = bestx - pattern[2][2][0];
    let dy = besty - pattern[2][2][1];

    // We consider 8 lines across the finder pattern in turn
    const MASK_TESTS: [[c_uint; 2]; 8] = [
        [0x1040041, 0x1000001],
        [0x0041040, 0x0001000],
        [0x0110110, 0x0100010],
        [0x0011100, 0x0001000],
        [0x0420084, 0x0400004],
        [0x0021080, 0x0001000],
        [0x0006C00, 0x0004400],
        [0x0003800, 0x0001000],
    ];
    const MASK_COORDS: [[usize; 2]; 8] = [
        [0, 0],
        [1, 1],
        [4, 0],
        [3, 1],
        [2, 0],
        [2, 1],
        [0, 2],
        [1, 2],
    ];

    for i in 0..8 {
        if (best_match & MASK_TESTS[i][0]) == MASK_TESTS[i][1] {
            let x0 = (pattern[MASK_COORDS[i][1]][MASK_COORDS[i][0]][0] + dx) >> QR_FINDER_SUBPREC;
            if x0 < 0 || x0 >= width {
                continue;
            }
            let y0 = (pattern[MASK_COORDS[i][1]][MASK_COORDS[i][0]][1] + dy) >> QR_FINDER_SUBPREC;
            if y0 < 0 || y0 >= height {
                continue;
            }
            let x1 = (pattern[4 - MASK_COORDS[i][1]][4 - MASK_COORDS[i][0]][0] + dx)
                >> QR_FINDER_SUBPREC;
            if x1 < 0 || x1 >= width {
                continue;
            }
            let y1 = (pattern[4 - MASK_COORDS[i][1]][4 - MASK_COORDS[i][0]][1] + dy)
                >> QR_FINDER_SUBPREC;
            if y1 < 0 || y1 >= height {
                continue;
            }
            if qr_finder_locate_crossing(
                img,
                width,
                height,
                x0,
                y0,
                x1,
                y1,
                (i & 1) as c_int,
                &mut pc,
            ) == 0
            {
                let cx = pc[0] - bestx;
                let cy = pc[1] - besty;
                if (i & 1) != 0 {
                    // Weight crossings around the center dot more highly
                    nc[i >> 1] += 3;
                    c[i >> 1][0] += cx + (cx << 1);
                    c[i >> 1][1] += cy + (cy << 1);
                } else {
                    nc[i >> 1] += 1;
                    c[i >> 1][0] += cx;
                    c[i >> 1][1] += cy;
                }
            }
        }
    }

    // Sum offsets from lines in orthogonal directions
    for i in 0..2 {
        let a = nc[i << 1];
        let b = nc[(i << 1) | 1];
        if a != 0 && b != 0 {
            let w = c_int::max(a, b);
            c[i << 1][0] = qr_divround(w * (b * c[i << 1][0] + a * c[(i << 1) | 1][0]), a * b);
            c[i << 1][1] = qr_divround(w * (b * c[i << 1][1] + a * c[(i << 1) | 1][1]), a * b);
            nc[i << 1] = w << 1;
        } else {
            c[i << 1][0] += c[(i << 1) | 1][0];
            c[i << 1][1] += c[(i << 1) | 1][1];
            nc[i << 1] += b;
        }
    }

    // Average offsets from pairs of orthogonal lines
    c[0][0] += c[2][0];
    c[0][1] += c[2][1];
    nc[0] += nc[2];

    // If we actually found any such lines, apply the adjustment
    if nc[0] != 0 {
        let dx = qr_divround(c[0][0], nc[0]);
        let dy = qr_divround(c[0][1], nc[0]);
        // But only if it doesn't make things too much worse
        let match_val =
            qr_alignment_pattern_fetch(&pattern, bestx + dx, besty + dy, img, width, height);
        let dist = qr_hamming_dist(match_val, 0x1F8D63F, best_dist + 1);
        if dist <= best_dist + 1 {
            bestx += dx;
            besty += dy;
        }
    }

    (*p)[0] = bestx;
    (*p)[1] = besty;
    0
}

// QR Code data structures and functions

/// QR code data mode
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(non_camel_case_types)]
pub enum qr_mode {
    /// Numeric digits ('0'...'9')
    QR_MODE_NUM = 1,
    /// Alphanumeric characters ('0'...'9', 'A'...'Z', plus punctuation ' ', '$', '%', '*', '+', '-', '.', '/', ':')
    QR_MODE_ALNUM = 2,
    /// Structured-append header
    QR_MODE_STRUCT = 3,
    /// Raw 8-bit bytes
    QR_MODE_BYTE = 4,
    /// FNC1 marker in first position (GS1 formatting)
    QR_MODE_FNC1_1ST = 5,
    /// Extended Channel Interpretation code
    QR_MODE_ECI = 7,
    /// SJIS kanji characters
    QR_MODE_KANJI = 8,
    /// FNC1 marker in second position (industry application)
    QR_MODE_FNC1_2ND = 9,
}

/// Check if a mode has a data buffer associated with it
#[inline]
fn qr_mode_has_data(mode: qr_mode) -> bool {
    let mode_val = mode as c_uint;
    (mode_val & (mode_val - 1)) == 0
}

/// Data payload for a QR code data entry
#[repr(C)]
#[derive(Copy, Clone)]
pub union qr_code_data_payload {
    /// Data buffer for modes that have one
    pub data: qr_code_data_buffer,
    /// Decoded "Extended Channel Interpretation" data
    pub eci: c_uint,
    /// Decoded "Application Indicator" for FNC1 in 2nd position
    pub ai: c_int,
    /// Structured-append header data
    pub sa: qr_code_data_sa,
}

/// Data buffer
#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct qr_code_data_buffer {
    pub buf: *mut c_uchar,
    pub len: c_int,
}

/// Structured-append data
#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct qr_code_data_sa {
    pub sa_index: c_uchar,
    pub sa_size: c_uchar,
    pub sa_parity: c_uchar,
}

/// A single QR code data entry
#[repr(C)]
pub struct qr_code_data_entry {
    /// The mode of this data block
    pub mode: qr_mode,
    /// The payload (union)
    pub payload: qr_code_data_payload,
}

/// Low-level QR code data
#[repr(C)]
pub struct qr_code_data {
    /// The decoded data entries
    pub entries: *mut qr_code_data_entry,
    pub nentries: c_int,
    /// The code version (1...40)
    pub version: c_uchar,
    /// The ECC level (0...3, corresponding to 'L', 'M', 'Q', and 'H')
    pub ecc_level: c_uchar,
    /// The index of this code in the structured-append group
    pub sa_index: c_uchar,
    /// The size of the structured-append group, or 0 if there was no S-A header
    pub sa_size: c_uchar,
    /// The parity of the entire structured-append group
    pub sa_parity: c_uchar,
    /// The parity of this code
    pub self_parity: c_uchar,
    /// An approximate bounding box for the code
    pub bbox: [qr_point; 4],
}

/// List of QR code data
#[repr(C)]
pub struct qr_code_data_list {
    pub qrdata: *mut qr_code_data,
    pub nqrdata: c_int,
    pub cqrdata: c_int,
}

/// Clear a QR code data structure, freeing all allocated memory
#[no_mangle]
pub unsafe extern "C" fn qr_code_data_clear(qrdata: *mut qr_code_data) {
    for i in 0..(*qrdata).nentries {
        let entry = (*qrdata).entries.offset(i as isize);
        if qr_mode_has_data((*entry).mode) {
            free((*entry).payload.data.buf as *mut c_void);
        }
    }
    free((*qrdata).entries as *mut c_void);
}

/// Initialize a QR code data list
#[no_mangle]
pub unsafe extern "C" fn qr_code_data_list_init(qrlist: *mut qr_code_data_list) {
    (*qrlist).qrdata = null_mut();
    (*qrlist).nqrdata = 0;
    (*qrlist).cqrdata = 0;
}

/// Clear a QR code data list, freeing all allocated memory
#[no_mangle]
pub unsafe extern "C" fn qr_code_data_list_clear(qrlist: *mut qr_code_data_list) {
    for i in 0..(*qrlist).nqrdata {
        qr_code_data_clear((*qrlist).qrdata.offset(i as isize));
    }
    free((*qrlist).qrdata as *mut c_void);
    qr_code_data_list_init(qrlist);
}

/// Add a QR code data to the list
#[no_mangle]
pub unsafe extern "C" fn qr_code_data_list_add(
    qrlist: *mut qr_code_data_list,
    qrdata: *const qr_code_data,
) {
    if (*qrlist).nqrdata >= (*qrlist).cqrdata {
        (*qrlist).cqrdata = ((*qrlist).cqrdata << 1) | 1;
        (*qrlist).qrdata = realloc(
            (*qrlist).qrdata as *mut c_void,
            ((*qrlist).cqrdata as usize) * std::mem::size_of::<qr_code_data>(),
        ) as *mut qr_code_data;
    }
    memcpy(
        (*qrlist).qrdata.offset((*qrlist).nqrdata as isize) as *mut c_void,
        qrdata as *const c_void,
        std::mem::size_of::<qr_code_data>(),
    );
    (*qrlist).nqrdata += 1;
}
