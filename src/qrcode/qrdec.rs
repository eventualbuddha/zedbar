//! QR Code decoder utilities
//!
//! This module provides low-level QR code decoding functions including
//! point geometry operations and error correction.

use std::{cmp::Ordering, ptr::null};

use libc::{c_int, c_uchar, c_uint};

use crate::{decoder_types::qr_finder_line, img_scanner::qr_reader, qrcode::util::qr_ilog};

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

/// Number of bits of sub-module precision for alignment pattern search
const QR_ALIGN_SUBPREC: c_int = 2;

/// Helper function: divide with exact rounding
///
/// Rounds towards positive infinity when x > 0, towards negative infinity when x < 0.
/// For x/y where the fractional part is exactly 0.5, rounds away from zero.
#[inline]
fn qr_divround(x: c_int, y: c_int) -> c_int {
    (x + x.signum() * (y >> 1)) / y
}

/// Fixed-point multiply with rounding and shift
///
/// Multiplies 32-bit numbers a and b, adds r, and takes bits [s, s+31] of the result.
/// This is used for fixed-point arithmetic to avoid overflow.
#[inline]
fn qr_fixmul(a: c_int, b: c_int, r: i64, s: c_int) -> c_int {
    ((a as i64 * b as i64 + r) >> s) as c_int
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
    let mut d = (*_l0)[0] * (*_l1)[1] - (*_l0)[1] * (*_l1)[0];
    if d == 0 {
        return -1;
    }
    let mut x = (*_l0)[1] * (*_l1)[2] - (*_l1)[1] * (*_l0)[2];
    let mut y = (*_l1)[0] * (*_l0)[2] - (*_l0)[0] * (*_l1)[2];
    if d < 0 {
        x = -x;
        y = -y;
        d = -d;
    }
    (*_p)[0] = qr_divround(x, d);
    (*_p)[1] = qr_divround(y, d);
    0
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
    (*_q)[0] = ((*_aff).inv[0][0] * (_x - (*_aff).x0)
        + (*_aff).inv[0][1] * (_y - (*_aff).y0)
        + ((1 << (*_aff).ires) >> 1))
        >> (*_aff).ires;
    (*_q)[1] = ((*_aff).inv[1][0] * (_x - (*_aff).x0)
        + (*_aff).inv[1][1] * (_y - (*_aff).y0)
        + ((1 << (*_aff).ires) >> 1))
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
    (*_p)[0] = (((*_aff).fwd[0][0] * _u + (*_aff).fwd[0][1] * _v + (1 << ((*_aff).res - 1)))
        >> (*_aff).res)
        + (*_aff).x0;
    (*_p)[1] = (((*_aff).fwd[1][0] * _u + (*_aff).fwd[1][1] * _v + (1 << ((*_aff).res - 1)))
        >> (*_aff).res)
        + (*_aff).y0;
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
    _x -= (*_hom).x0;
    _y -= (*_hom).y0;
    let mut x = (*_hom).inv[0][0] * _x + (*_hom).inv[0][1] * _y;
    let mut y = (*_hom).inv[1][0] * _x + (*_hom).inv[1][1] * _y;
    let mut w = ((*_hom).inv[2][0] * _x
        + (*_hom).inv[2][1] * _y
        + (*_hom).inv22
        + (1 << ((*_hom).res - 1)))
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
        if i00_full < 0 { -result } else { result }
    } else {
        0
    };
    let i01 = if i01_full != 0 {
        let result = qr_divround(i22, i01_full.abs());
        if i01_full < 0 { -result } else { result }
    } else {
        0
    };
    let i10 = if i10_full != 0 {
        let result = qr_divround(i22, i10_full.abs());
        if i10_full < 0 { -result } else { result }
    } else {
        0
    };
    let i11 = if i11_full != 0 {
        let result = qr_divround(i22, i11_full.abs());
        if i11_full < 0 { -result } else { result }
    } else {
        0
    };
    let i20 = if i20_full != 0 {
        let result = qr_divround(i22, i20_full.abs());
        if i20_full < 0 { -result } else { result }
    } else {
        0
    };
    let i21 = if i21_full != 0 {
        let result = qr_divround(i22, i21_full.abs());
        if i21_full < 0 { -result } else { result }
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
    let b0 = qr_ilog(c_int::max(dx10.abs(), dy10.abs()) as u32) + qr_ilog((a20 + a22).abs() as u32);
    let b1 = qr_ilog(c_int::max(dx20.abs(), dy20.abs()) as u32) + qr_ilog((a21 + a22).abs() as u32);
    let b2 = qr_ilog(c_int::max(c_int::max(a20.abs(), a21.abs()), a22.abs()) as u32);
    let shift = c_int::max(0, c_int::max(c_int::max(b0, b1), b2) - (QR_INT_BITS - 3 - QR_ALIGN_SUBPREC));
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
