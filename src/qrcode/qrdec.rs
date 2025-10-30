//! QR Code decoder utilities
//!
//! This module provides low-level QR code decoding functions including
//! point geometry operations and error correction.

use std::{cmp::Ordering, mem::swap, slice::from_raw_parts};

use libc::{c_int, c_uchar, c_uint};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use reed_solomon::Decoder as RSDecoder;

use crate::{
    decoder::qr_finder_line,
    image_ffi::zbar_image_t,
    img_scanner::{qr_reader, zbar_image_scanner_t},
    qrcode::{
        binarize::binarize,
        qrdectxt::qr_code_data_list_extract_text,
        util::{qr_ihypot, qr_ilog, qr_isqrt},
    },
};

use super::bch15_5::bch15_5_correct;

/// A line in QR code coordinate space: [A, B, C] for equation Ax + By + C = 0
pub(crate) type qr_line = [c_int; 3];

/// A point in QR code coordinate space: [x, y]
pub(crate) type qr_point = [c_int; 2];

/// Number of bits in an int (typically 32)
const QR_INT_BITS: c_int = c_int::BITS as c_int;

/// Log2 of QR_INT_BITS (typically 5 for 32-bit ints)
const QR_INT_LOGBITS: c_int = 5; // qr_ilog(32) = 5

/// Number of bits of sub-module precision for alignment pattern search
const QR_ALIGN_SUBPREC: c_int = 2;

/// Number of bits of sub-module precision for finder pattern coordinates
const QR_FINDER_SUBPREC: c_int = 2;

/// A 14 bit resolution for a homography ensures that the ideal module size for a
/// version 40 code differs from that of a version 39 code by at least 2.
const QR_HOM_BITS: c_int = 14;

/// The amount that the estimated version numbers are allowed to differ from the
/// real version number and still be considered valid.
#[allow(dead_code)] // Used in functions not yet ported to Rust
const QR_SMALL_VERSION_SLACK: c_int = 1;

/// Since cell phone cameras can have severe radial distortion, the estimated
/// version for larger versions can be off by larger amounts.
const QR_LARGE_VERSION_SLACK: c_int = 3;

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

/// Extended multiply: multiplies 32-bit numbers a and b, adds r, and returns 64-bit result
///
/// This matches C macro QR_EXTMUL.
#[inline]
fn qr_extmul(a: c_int, b: c_int, r: i64) -> i64 {
    a as i64 * b as i64 + r
}

/// Sign mask: returns -1 if x < 0, else 0
///
/// This matches C macro QR_SIGNMASK.
#[inline]
fn qr_signmask(x: i64) -> i64 {
    -((x < 0) as i64)
}

/// Project alignment pattern center to corner of QR code
///
/// Given three corners and an alignment pattern center, compute the fourth corner
/// by geometric projection. Returns `Ok((brx, bry))` on success, `Err(-1)` if projection fails.
fn qr_hom_project_alignment_to_corner(
    p: &[qr_point],
    p3: &qr_point,
    dim: c_int,
) -> Result<(c_int, c_int), c_int> {
    let p0 = &p[0];
    let p1 = &p[1];
    let p2 = &p[2];

    let c21 = p2[0] * p1[1] - p2[1] * p1[0];
    let dx21 = p2[0] - p1[0];
    let dy21 = p2[1] - p1[1];

    let mut w = qr_extmul(
        dim - 7,
        c21,
        qr_extmul(
            dim - 13,
            p0[0] * dy21 - p0[1] * dx21,
            qr_extmul(6, p3[0] * dy21 - p3[1] * dx21, 0),
        ),
    );

    // The projection failed: invalid geometry
    if w == 0 {
        return Err(-1);
    }

    let mask = qr_signmask(w);
    w = (w + mask) ^ mask;

    // Inline division with rounding for i64 precision
    let brx_num = (qr_extmul(
        (dim - 7) * p0[0],
        p3[0] * dy21,
        qr_extmul(
            (dim - 13) * p3[0],
            c21 - p0[1] * dx21,
            qr_extmul(6 * p0[0], c21 - p3[1] * dx21, 0),
        ),
    ) + mask)
        ^ mask;
    let brx = ((brx_num + brx_num.signum() * (w >> 1)) / w) as c_int;

    let bry_num = (qr_extmul(
        (dim - 7) * p0[1],
        -p3[1] * dx21,
        qr_extmul(
            (dim - 13) * p3[1],
            c21 + p0[0] * dy21,
            qr_extmul(6 * p0[1], c21 + p3[0] * dy21, 0),
        ),
    ) + mask)
        ^ mask;
    let bry = ((bry_num + bry_num.signum() * (w >> 1)) / w) as c_int;

    Ok((brx, bry))
}

/// Fit a line to collected edge points, or use axis-aligned fallback if insufficient points
///
/// This is used for edges that have only one finder pattern, where we walk along the edge
/// collecting sample points. If we don't get enough points (> 1), we fall back to an
/// axis-aligned line in the affine coordinate system.
unsafe fn qr_hom_fit_edge_line(
    line: &mut qr_line,
    pts: &mut [qr_point],
    npts: usize,
    finder: &qr_finder,
    aff: &qr_aff,
    edge_axis: c_int,
) {
    if npts > 1 {
        qr_line_fit_points(line, pts, npts, aff.res);
    } else {
        // Project reference point from the finder pattern
        let p = if edge_axis == 1 {
            // Right edge: project from UR finder, extending 3 modules to the right
            qr_aff_project(aff, finder.o[0] + 3 * finder.size[0], finder.o[1])
        } else {
            // Bottom edge (axis 3): project from DL finder, extending 3 modules down
            qr_aff_project(aff, finder.o[0], finder.o[1] + 3 * finder.size[1])
        };

        // Calculate normalization shift (always uses column 1 of affine matrix)
        let shift = 0.max(
            qr_ilog((aff.fwd[0][1].abs()).max(aff.fwd[1][1].abs()) as c_uint)
                - ((aff.res + 1) >> 1),
        );
        let round = (1 << shift) >> 1;

        // Compute line coefficients using appropriate matrix column
        if edge_axis == 1 {
            // Right edge uses column 1 (vertical direction in affine space)
            line[0] = (aff.fwd[1][1] + round) >> shift;
            line[1] = (-aff.fwd[0][1] + round) >> shift;
        } else {
            // Bottom edge uses column 0 (horizontal direction in affine space)
            line[0] = (aff.fwd[1][0] + round) >> shift;
            line[1] = (-aff.fwd[0][0] + round) >> shift;
        }

        // Compute line constant term
        line[2] = -(line[0] * p[0] + line[1] * p[1]);
    }
}

/// A cell in the sampling grid for homographic projection
///
/// Represents a mapping from a unit square to a quadrilateral in the image,
/// used for extracting QR code modules with perspective correction.
#[derive(Copy, Clone, Default)]
struct qr_hom_cell {
    /// Forward transformation matrix [3][3]
    pub(crate) fwd: [[c_int; 3]; 3],
    /// X offset in image space
    pub(crate) x0: c_int,
    /// Y offset in image space
    pub(crate) y0: c_int,
    /// U offset in code space
    pub(crate) u0: c_int,
    /// V offset in code space
    pub(crate) v0: c_int,
}

/// Sampling grid for QR code module extraction
struct qr_sampling_grid {
    /// 2D array of homography cells for mapping between code and image space
    /// cells[i][j] represents the cell at row i, column j
    pub(crate) cells: Vec<Vec<qr_hom_cell>>,
    /// Mask indicating which modules are part of function patterns
    pub(crate) fpmask: Vec<c_uint>,
    /// Limits for each cell region
    pub(crate) cell_limits: [c_int; 6],
}

/// collection of finder lines
#[derive(Default)]
pub(crate) struct qr_finder_lines {
    lines: Vec<qr_finder_line>,
}

/// Allocates a client reader handle.
pub(crate) unsafe fn _zbar_qr_create() -> qr_reader {
    qr_reader {
        rng: ChaCha8Rng::from_seed([0u8; 32]),
        finder_lines: [qr_finder_lines::default(), qr_finder_lines::default()],
    }
}

/// reset finder state between scans
pub(crate) unsafe fn _zbar_qr_reset(reader: &mut qr_reader) {
    reader.finder_lines[0].lines.clear();
    reader.finder_lines[1].lines.clear();
}

/// A cluster of finder lines that have been grouped together.
/// Contains indices into the parent ClusteredLines structure.
#[derive(Clone, Default)]
pub(crate) struct Cluster {
    /// Indices into the lines array of the parent ClusteredLines structure.
    line_indices: Vec<usize>,
}

/// Encapsulates finder lines and their clusters.
/// This structure owns both the lines and the clusters, with clusters
/// referencing lines via indices rather than raw pointers.
#[derive(Clone, Default)]
pub(crate) struct ClusteredLines {
    lines: Vec<qr_finder_line>,
    clusters: Vec<Cluster>,
}

impl ClusteredLines {
    /// Create a new ClusteredLines from a vector of lines
    pub fn new(lines: Vec<qr_finder_line>) -> Self {
        Self {
            lines,
            clusters: Vec::new(),
        }
    }

    /// Access a line by index
    pub fn line(&self, idx: usize) -> &qr_finder_line {
        &self.lines[idx]
    }

    /// Get the number of clusters
    pub fn cluster_count(&self) -> usize {
        self.clusters.len()
    }

    /// Get a cluster by index
    pub fn cluster(&self, idx: usize) -> &Cluster {
        &self.clusters[idx]
    }

    /// Access a specific line within a cluster
    pub fn cluster_line(&self, cluster_idx: usize, line_idx: usize) -> &qr_finder_line {
        let line_index = self.clusters[cluster_idx].line_indices[line_idx];
        &self.lines[line_index]
    }

    /// Iterate over lines in a cluster
    pub fn cluster_lines(&self, cluster_idx: usize) -> impl Iterator<Item = &qr_finder_line> + '_ {
        self.clusters[cluster_idx]
            .line_indices
            .iter()
            .map(|&idx| &self.lines[idx])
    }
}

/// A point on the edge of a finder pattern. These are obtained from the
/// endpoints of the lines crossing this particular pattern.
#[derive(Copy, Clone, Default)]
pub(crate) struct qr_finder_edge_pt {
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
#[derive(Clone, Default)]
pub(crate) struct qr_finder_center {
    /// The estimated location of the finder center.
    pos: qr_point,

    /// The list of edge points from the crossing lines.
    edge_pts: Vec<qr_finder_edge_pt>,
}

/// Determine if a horizontal line crosses a vertical line.
///
/// # Returns
/// True if the lines cross, false otherwise.
pub(crate) fn qr_finder_lines_are_crossing(hline: &qr_finder_line, vline: &qr_finder_line) -> bool {
    hline.pos[0] <= vline.pos[0]
        && vline.pos[0] < hline.pos[0] + hline.len
        && vline.pos[1] <= hline.pos[1]
        && hline.pos[1] < vline.pos[1] + vline.len
}

/// Translate a point by the given offsets
///
/// Adds dx to the x coordinate and dy to the y coordinate.
pub(crate) fn qr_point_translate(point: &mut qr_point, dx: c_int, dy: c_int) {
    point[0] += dx;
    point[1] += dy;
}

/// Calculate the squared distance between two points
///
/// Returns the squared Euclidean distance, which avoids the need for
/// expensive square root calculations when only relative distances matter.
pub(crate) fn qr_point_distance2(p1: &qr_point, p2: &qr_point) -> c_uint {
    let dx = p1[0] - p2[0];
    let dy = p1[1] - p2[1];
    (dx * dx + dy * dy) as c_uint
}

/// Check if three points are in counter-clockwise order
///
/// Returns the cross product of the vectors (p1-p0) and (p2-p0).
/// - Positive: points are in CCW order (in right-handed coordinate system)
/// - Zero: points are collinear
/// - Negative: points are in CW order
pub(crate) fn qr_point_ccw(p0: &qr_point, p1: &qr_point, p2: &qr_point) -> c_int {
    let p0x = p0[0];
    let p0y = p0[1];
    let p1x = p1[0];
    let p1y = p1[1];
    let p2x = p2[0];
    let p2y = p2[1];

    (p1x - p0x) * (p2y - p0y) - (p1y - p0y) * (p2x - p0x)
}

/// Evaluate a line equation at a point
///
/// Given a line defined by the equation A*x + B*y + C = 0,
/// this returns the value A*x + B*y + C for the given coordinates.
pub(crate) fn qr_line_eval(line: &qr_line, x: c_int, y: c_int) -> c_int {
    line[0] * x + line[1] * y + line[2]
}

pub(crate) fn qr_line_orient(_l: &mut qr_line, _x: c_int, _y: c_int) {
    if qr_line_eval(_l, _x, _y) < 0 {
        _l[0] = -_l[0];
        _l[1] = -_l[1];
        _l[2] = -_l[2];
    }
}

pub(crate) unsafe fn qr_line_isect(_l0: &qr_line, _l1: &qr_line) -> Option<qr_point> {
    let mut d = (*_l0)[0]
        .wrapping_mul((*_l1)[1])
        .wrapping_sub((*_l0)[1].wrapping_mul((*_l1)[0]));
    if d == 0 {
        return None;
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
    Some([qr_divround(x, d), qr_divround(y, d)])
}

/// Fit a line to covariance data using least-squares
///
/// # Parameters
/// - `_x0`, `_y0`: Centroid coordinates
/// - `_sxx`, `_sxy`, `_syy`: Covariance matrix values
/// - `_res`: Resolution bits for scaling
///
/// # Returns
/// Line equation (Ax + By + C = 0)
pub(crate) fn qr_line_fit(
    _x0: c_int,
    _y0: c_int,
    _sxx: c_int,
    _sxy: c_int,
    _syy: c_int,
    _res: c_int,
) -> qr_line {
    let u = (_sxx - _syy).abs();
    let v = -_sxy << 1;
    let w = qr_ihypot(u, v) as c_int;

    // Compute shift factor to scale down into manageable range
    // Ensure product of any two of _l[0] and _l[1] fits within _res bits
    let dshift = 0.max(qr_ilog(u.max(v.abs()) as u32) + 1 - ((_res + 1) >> 1));
    let dround = (1 << dshift) >> 1;

    let mut line = qr_line::default();
    if _sxx > _syy {
        line[0] = (v + dround) >> dshift;
        line[1] = (u + w + dround) >> dshift;
    } else {
        line[0] = (u + w + dround) >> dshift;
        line[1] = (v + dround) >> dshift;
    }
    line[2] = -(_x0
        .wrapping_mul(line[0])
        .wrapping_add(_y0.wrapping_mul(line[1])));
    line
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
pub(crate) unsafe fn qr_line_fit_points(
    _l: &mut qr_line,
    _p: &mut [qr_point],
    _np: usize,
    _res: c_int,
) {
    let mut sx: c_int = 0;
    let mut sy: c_int = 0;
    let mut xmin = c_int::MAX;
    let mut xmax = c_int::MIN;
    let mut ymin = c_int::MAX;
    let mut ymax = c_int::MIN;

    // Compute centroid and bounds
    for point in _p.iter() {
        let px = point[0];
        let py = point[1];
        sx = sx.wrapping_add(px);
        xmin = xmin.min(px);
        xmax = xmax.max(px);
        sy = sy.wrapping_add(py);
        ymin = ymin.min(py);
        ymax = ymax.max(py);
    }

    let xbar = (sx + (_np >> 1) as c_int) / _np as c_int;
    let ybar = (sy + (_np >> 1) as c_int) / _np as c_int;

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
    for point in _p.iter().take(_np) {
        let dx = (point[0] - xbar + sround) >> sshift;
        let dy = (point[1] - ybar + sround) >> sshift;
        sxx = sxx.wrapping_add(dx.wrapping_mul(dx));
        sxy = sxy.wrapping_add(dx.wrapping_mul(dy));
        syy = syy.wrapping_add(dy.wrapping_mul(dy));
    }

    *_l = qr_line_fit(xbar, ybar, sxx, sxy, syy, _res);
}

/// An affine homography.
/// This maps from the image (at subpel resolution) to a square domain with
/// power-of-two sides (of res bits) and back.
#[derive(Copy, Clone, Default)]
pub(crate) struct qr_aff {
    /// Forward transformation matrix [2][2]
    pub(crate) fwd: [[c_int; 2]; 2],
    /// Inverse transformation matrix [2][2]
    pub(crate) inv: [[c_int; 2]; 2],
    /// X offset
    pub(crate) x0: c_int,
    /// Y offset
    pub(crate) y0: c_int,
    /// Resolution bits
    pub(crate) res: c_int,
    /// Inverse resolution bits
    pub(crate) ires: c_int,
}

pub(crate) fn qr_aff_init(
    _aff: &mut qr_aff,
    _p0: &qr_point,
    _p1: &qr_point,
    _p2: &qr_point,
    _res: c_int,
) {
    // det is ensured to be positive by our caller.
    let dx1 = _p1[0] - _p0[0];
    let dx2 = _p2[0] - _p0[0];
    let dy1 = _p1[1] - _p0[1];
    let dy2 = _p2[1] - _p0[1];
    let det = dx1 * dy2 - dy1 * dx2;
    let ires = c_int::max(((qr_ilog(det.unsigned_abs()) as u32 >> 1) - 2) as i32, 0);
    _aff.fwd[0][0] = dx1;
    _aff.fwd[0][1] = dx2;
    _aff.fwd[1][0] = dy1;
    _aff.fwd[1][1] = dy2;
    _aff.inv[0][0] = qr_divround(dy2 << _res, det >> ires);
    _aff.inv[0][1] = qr_divround(-dx2 << _res, det >> ires);
    _aff.inv[1][0] = qr_divround(-dy1 << _res, det >> ires);
    _aff.inv[1][1] = qr_divround(dx1 << _res, det >> ires);
    _aff.x0 = _p0[0];
    _aff.y0 = _p0[1];
    _aff.res = _res;
    _aff.ires = ires;
}

/// Map from the image (at subpel resolution) into the square domain.
pub(crate) fn qr_aff_unproject(_aff: &qr_aff, _x: c_int, _y: c_int) -> qr_point {
    let x = (_aff.inv[0][0]
        .wrapping_mul(_x.wrapping_sub(_aff.x0))
        .wrapping_add(_aff.inv[0][1].wrapping_mul(_y.wrapping_sub(_aff.y0)))
        .wrapping_add((1 << _aff.ires) >> 1))
        >> _aff.ires;
    let y = (_aff.inv[1][0]
        .wrapping_mul(_x.wrapping_sub(_aff.x0))
        .wrapping_add(_aff.inv[1][1].wrapping_mul(_y.wrapping_sub(_aff.y0)))
        .wrapping_add((1 << _aff.ires) >> 1))
        >> _aff.ires;
    [x, y]
}

/// Map from the square domain into the image (at subpel resolution).
pub(crate) fn qr_aff_project(_aff: &qr_aff, _u: c_int, _v: c_int) -> qr_point {
    let x = ((_aff.fwd[0][0]
        .wrapping_mul(_u)
        .wrapping_add(_aff.fwd[0][1].wrapping_mul(_v))
        .wrapping_add(1 << (_aff.res - 1)))
        >> _aff.res)
        .wrapping_add(_aff.x0);
    let y = ((_aff.fwd[1][0]
        .wrapping_mul(_u)
        .wrapping_add(_aff.fwd[1][1].wrapping_mul(_v))
        .wrapping_add(1 << (_aff.res - 1)))
        >> _aff.res)
        .wrapping_add(_aff.y0);
    [x, y]
}

/// A full homography.
/// Like the affine homography, this maps from the image (at subpel resolution)
/// to a square domain with power-of-two sides (of res bits) and back.
#[derive(Copy, Clone, Default)]
pub(crate) struct qr_hom {
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
#[allow(clippy::too_many_arguments)]
pub(crate) unsafe fn qr_hom_init(
    _hom: &mut qr_hom,
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
    let b0 =
        qr_ilog(c_int::max(dx10.abs(), dy10.abs()) as u32) + qr_ilog((a20 + a22).unsigned_abs());
    let b1 =
        qr_ilog(c_int::max(dx20.abs(), dy20.abs()) as u32) + qr_ilog((a21 + a22).unsigned_abs());
    let b2 = qr_ilog(c_int::max(c_int::max(a20.abs(), a21.abs()), a22.abs()) as u32);
    let s1 = c_int::max(
        0,
        _res + c_int::max(c_int::max(b0, b1), b2) - (QR_INT_BITS - 2),
    );
    let r1 = (1i64 << s1) >> 1;

    // Compute the final coefficients of the forward transform
    // The 32x32->64 bit multiplies are really needed for accuracy with large versions
    _hom.fwd[0][0] = qr_fixmul(dx10, a20 + a22, r1, s1);
    _hom.fwd[0][1] = qr_fixmul(dx20, a21 + a22, r1, s1);
    _hom.x0 = _x0;
    _hom.fwd[1][0] = qr_fixmul(dy10, a20 + a22, r1, s1);
    _hom.fwd[1][1] = qr_fixmul(dy20, a21 + a22, r1, s1);
    _hom.y0 = _y0;
    _hom.fwd[2][0] = (a20 + r1 as c_int) >> s1;
    _hom.fwd[2][1] = (a21 + r1 as c_int) >> s1;
    _hom.fwd22 = if s1 > _res {
        (a22 + ((r1 >> _res) as c_int)) >> (s1 - _res)
    } else {
        a22 << (_res - s1)
    };

    // Now compute the inverse transform
    let b0 = qr_ilog(c_int::max(c_int::max(dx10.abs(), dx20.abs()), dx30.abs()) as u32)
        + qr_ilog(c_int::max(_hom.fwd[0][0].abs(), _hom.fwd[1][0].abs()) as u32);
    let b1 = qr_ilog(c_int::max(c_int::max(dy10.abs(), dy20.abs()), dy30.abs()) as u32)
        + qr_ilog(c_int::max(_hom.fwd[0][1].abs(), _hom.fwd[1][1].abs()) as u32);
    let b2 = qr_ilog(a22.unsigned_abs()) - s1;
    let s2 = c_int::max(0, c_int::max(b0, b1) + b2 - (QR_INT_BITS - 3));
    let r2 = (1i64 << s2) >> 1;
    let s1 = s1 + s2;
    let r1 = r1 << s2;

    // The 32x32->64 bit multiplies are really needed for accuracy with large versions
    _hom.inv[0][0] = qr_fixmul(_hom.fwd[1][1], a22, r1, s1);
    _hom.inv[0][1] = qr_fixmul(-_hom.fwd[0][1], a22, r1, s1);
    _hom.inv[1][0] = qr_fixmul(-_hom.fwd[1][0], a22, r1, s1);
    _hom.inv[1][1] = qr_fixmul(_hom.fwd[0][0], a22, r1, s1);
    _hom.inv[2][0] = qr_fixmul(
        _hom.fwd[1][0],
        _hom.fwd[2][1],
        -qr_extmul(_hom.fwd[1][1], _hom.fwd[2][0], r2),
        s2,
    );
    _hom.inv[2][1] = qr_fixmul(
        _hom.fwd[0][1],
        _hom.fwd[2][0],
        -qr_extmul(_hom.fwd[0][0], _hom.fwd[2][1], r2),
        s2,
    );
    _hom.inv22 = qr_fixmul(
        _hom.fwd[0][0],
        _hom.fwd[1][1],
        -qr_extmul(_hom.fwd[0][1], _hom.fwd[1][0], r2),
        s2,
    );
    _hom.res = _res;
}

/// Fit a homography to correct large-scale perspective distortion
///
/// This function attempts to correct perspective distortion by fitting lines
/// to the edges of the QR code area and then building a homography transformation.
#[allow(clippy::too_many_arguments)]
pub(crate) unsafe fn qr_hom_fit(
    _hom: &mut qr_hom,
    _ul: &mut qr_finder,
    _ur: &mut qr_finder,
    _dl: &mut qr_finder,
    _p: &mut [qr_point],
    _aff: &qr_aff,
    rng: &mut ChaCha8Rng,
    img: &[u8],
    _width: c_int,
    _height: c_int,
) -> c_int {
    let mut l: [qr_line; 4] = [[0; 3]; 4];

    // We attempt to correct large-scale perspective distortion by fitting lines
    // to the edge of the code area.

    // Fitting lines is easy for the edges on which we have two finder patterns.
    // After the fit, UL is guaranteed to be on the proper side, but if either of
    // the other two finder patterns aren't, something is wrong.
    qr_finder_ransac(_ul, _aff, rng, 0);
    qr_finder_ransac(_dl, _aff, rng, 0);
    qr_line_fit_finder_pair(&mut l[0], _aff, _ul, _dl, 0);
    if qr_line_eval(
        &l[0],
        _dl.c.as_ref().unwrap().pos[0],
        _dl.c.as_ref().unwrap().pos[1],
    ) < 0
        || qr_line_eval(
            &l[0],
            _ur.c.as_ref().unwrap().pos[0],
            _ur.c.as_ref().unwrap().pos[1],
        ) < 0
    {
        return -1;
    }
    qr_finder_ransac(_ul, _aff, rng, 2);
    qr_finder_ransac(_ur, _aff, rng, 2);
    qr_line_fit_finder_pair(&mut l[2], _aff, _ul, _ur, 2);
    if qr_line_eval(
        &l[2],
        _dl.c.as_ref().unwrap().pos[0],
        _dl.c.as_ref().unwrap().pos[1],
    ) < 0
        || qr_line_eval(
            &l[2],
            _ur.c.as_ref().unwrap().pos[0],
            _ur.c.as_ref().unwrap().pos[1],
        ) < 0
    {
        return -1;
    }

    // The edges which only have one finder pattern are more difficult.
    let drv = _ur.size[1] >> 1;
    qr_finder_ransac(_ur, _aff, rng, 1);
    let mut dru = 0;
    if qr_line_fit_finder_edge(&mut l[1], _ur, 1, _aff.res) >= 0 {
        if qr_line_eval(
            &l[1],
            _ul.c.as_ref().unwrap().pos[0],
            _ul.c.as_ref().unwrap().pos[1],
        ) < 0
            || qr_line_eval(
                &l[1],
                _dl.c.as_ref().unwrap().pos[0],
                _dl.c.as_ref().unwrap().pos[1],
            ) < 0
        {
            return -1;
        }
        // Figure out the change in ru for a given change in rv when stepping along the fitted line
        dru = match qr_aff_line_step(_aff, &l[1], 1, drv) {
            Ok(v) => v,
            Err(_) => return -1,
        };
    }
    let mut ru = _ur.o[0] + 3 * _ur.size[0] - 2 * dru;
    let mut rv = _ur.o[1] - 2 * drv;

    let dbu = _dl.size[0] >> 1;
    qr_finder_ransac(_dl, _aff, rng, 3);
    let mut dbv = 0;
    if qr_line_fit_finder_edge(&mut l[3], _dl, 3, _aff.res) >= 0 {
        if qr_line_eval(
            &l[3],
            _ul.c.as_ref().unwrap().pos[0],
            _ul.c.as_ref().unwrap().pos[1],
        ) < 0
            || qr_line_eval(
                &l[3],
                _ur.c.as_ref().unwrap().pos[0],
                _ur.c.as_ref().unwrap().pos[1],
            ) < 0
        {
            return -1;
        }
        // Figure out the change in bv for a given change in bu when stepping along the fitted line
        dbv = match qr_aff_line_step(_aff, &l[3], 0, dbu) {
            Ok(v) => v,
            Err(_) => return -1,
        };
    }
    let mut bu = _dl.o[0] - 2 * dbu;
    let mut bv = _dl.o[1] + 3 * _dl.size[1] - 2 * dbv;

    // Set up the initial point lists
    let mut nr = _ur.ninliers[1];
    let mut rlastfit = nr;
    let mut cr = nr + (_dl.o[1] - rv + drv - 1) / drv;
    let mut r: Vec<qr_point> = Vec::with_capacity(cr as usize);
    for i in 0.._ur.ninliers[1] {
        r.push((*_ur.edge_pts[1].add(i as usize)).pos);
    }

    let mut nb = _dl.ninliers[3];
    let mut blastfit = nb;
    let mut cb = nb + (_ur.o[0] - bu + dbu - 1) / dbu;
    let mut b: Vec<qr_point> = Vec::with_capacity(cb as usize);
    for i in 0.._dl.ninliers[3] {
        b.push((*_dl.edge_pts[3].add(i as usize)).pos);
    }

    // Set up the step parameters for the affine projection
    let ox = (_aff.x0 << _aff.res) + (1 << (_aff.res - 1));
    let oy = (_aff.y0 << _aff.res) + (1 << (_aff.res - 1));
    let mut rx = _aff.fwd[0][0] * ru + _aff.fwd[0][1] * rv + ox;
    let mut ry = _aff.fwd[1][0] * ru + _aff.fwd[1][1] * rv + oy;
    let mut drxi = _aff.fwd[0][0] * dru + _aff.fwd[0][1] * drv;
    let mut dryi = _aff.fwd[1][0] * dru + _aff.fwd[1][1] * drv;
    let drxj = _aff.fwd[0][0] * _ur.size[0];
    let dryj = _aff.fwd[1][0] * _ur.size[0];
    let mut bx = _aff.fwd[0][0] * bu + _aff.fwd[0][1] * bv + ox;
    let mut by = _aff.fwd[1][0] * bu + _aff.fwd[1][1] * bv + oy;
    let mut dbxi = _aff.fwd[0][0] * dbu + _aff.fwd[0][1] * dbv;
    let mut dbyi = _aff.fwd[1][0] * dbu + _aff.fwd[1][1] * dbv;
    let dbxj = _aff.fwd[0][1] * _dl.size[1];
    let dbyj = _aff.fwd[1][1] * _dl.size[1];

    // Now step along the lines, looking for new sample points
    let mut nrempty = 0;
    let mut nbempty = 0;
    loop {
        let rdone = rv >= bv.min((_dl.o[1] + bv) >> 1) || nrempty > 14;
        let bdone = bu >= ru.min((_ur.o[0] + ru) >> 1) || nbempty > 14;

        if !rdone && (bdone || rv < bu) {
            let x0 = (rx + drxj) >> (_aff.res + QR_FINDER_SUBPREC);
            let y0 = (ry + dryj) >> (_aff.res + QR_FINDER_SUBPREC);
            let x1 = (rx - drxj) >> (_aff.res + QR_FINDER_SUBPREC);
            let y1 = (ry - dryj) >> (_aff.res + QR_FINDER_SUBPREC);

            if nr >= cr {
                cr = (cr << 1) | 1;
                r.reserve((cr - nr) as usize);
            }

            let mut ret = qr_finder_quick_crossing_check(img, _width, _height, x0, y0, x1, y1, 1);
            if ret == 0 {
                r.push([0; 2]);
                ret = qr_finder_locate_crossing(
                    img,
                    _width,
                    _height,
                    x0,
                    y0,
                    x1,
                    y1,
                    1,
                    &mut r[nr as usize],
                );
            }

            if ret >= 0 {
                if ret == 0 {
                    let q = qr_aff_unproject(_aff, r[nr as usize][0], r[nr as usize][1]);
                    // Move the current point halfway towards the crossing
                    ru = (ru + q[0]) >> 1;
                    // But ensure that rv monotonically increases
                    if q[1] + drv > rv {
                        rv = (rv + q[1]) >> 1;
                    }
                    rx = _aff.fwd[0][0] * ru + _aff.fwd[0][1] * rv + ox;
                    ry = _aff.fwd[1][0] * ru + _aff.fwd[1][1] * rv + oy;
                    nr += 1;
                    // Re-fit the line to update the step direction periodically
                    if nr > 1.max(rlastfit + (rlastfit >> 2)) {
                        qr_line_fit_points(&mut l[1], &mut r, nr as usize, _aff.res);
                        if let Ok(new_dru) = qr_aff_line_step(_aff, &l[1], 1, drv) {
                            dru = new_dru;
                            drxi = _aff.fwd[0][0] * dru + _aff.fwd[0][1] * drv;
                            dryi = _aff.fwd[1][0] * dru + _aff.fwd[1][1] * drv;
                        }
                        rlastfit = nr;
                    }
                }
                nrempty = 0;
            } else {
                nrempty += 1;
            }
            ru += dru;
            // Our final defense: if we overflow, stop
            if rv + drv > rv {
                rv += drv;
            } else {
                nrempty = c_int::MAX;
            }
            rx += drxi;
            ry += dryi;
        } else if !bdone {
            let x0 = (bx + dbxj) >> (_aff.res + QR_FINDER_SUBPREC);
            let y0 = (by + dbyj) >> (_aff.res + QR_FINDER_SUBPREC);
            let x1 = (bx - dbxj) >> (_aff.res + QR_FINDER_SUBPREC);
            let y1 = (by - dbyj) >> (_aff.res + QR_FINDER_SUBPREC);

            if nb >= cb {
                cb = (cb << 1) | 1;
                b.reserve((cb - nb) as usize);
            }

            let mut ret = qr_finder_quick_crossing_check(img, _width, _height, x0, y0, x1, y1, 1);
            if ret == 0 {
                b.push([0; 2]);
                ret = qr_finder_locate_crossing(
                    img,
                    _width,
                    _height,
                    x0,
                    y0,
                    x1,
                    y1,
                    1,
                    &mut b[nb as usize],
                );
            }

            if ret >= 0 {
                if ret == 0 {
                    let q = qr_aff_unproject(_aff, b[nb as usize][0], b[nb as usize][1]);
                    // Move the current point halfway towards the crossing
                    // But ensure that bu monotonically increases
                    if q[0] + dbu > bu {
                        bu = (bu + q[0]) >> 1;
                    }
                    bv = (bv + q[1]) >> 1;
                    bx = _aff.fwd[0][0] * bu + _aff.fwd[0][1] * bv + ox;
                    by = _aff.fwd[1][0] * bu + _aff.fwd[1][1] * bv + oy;
                    nb += 1;
                    // Re-fit the line to update the step direction periodically
                    if nb > 1.max(blastfit + (blastfit >> 2)) {
                        qr_line_fit_points(&mut l[3], &mut b, nb as usize, _aff.res);
                        if let Ok(new_dbv) = qr_aff_line_step(_aff, &l[3], 0, dbu) {
                            dbv = new_dbv;
                            dbxi = _aff.fwd[0][0] * dbu + _aff.fwd[0][1] * dbv;
                            dbyi = _aff.fwd[1][0] * dbu + _aff.fwd[1][1] * dbv;
                        }
                        blastfit = nb;
                    }
                }
                nbempty = 0;
            } else {
                nbempty += 1;
            }
            // Our final defense: if we overflow, stop
            if bu + dbu > bu {
                bu += dbu;
            } else {
                nbempty = c_int::MAX;
            }
            bv += dbv;
            bx += dbxi;
            by += dbyi;
        } else {
            break;
        }
    }

    // Fit the new lines
    qr_hom_fit_edge_line(&mut l[1], &mut r, nr as usize, _ur, _aff, 1);
    qr_hom_fit_edge_line(&mut l[3], &mut b, nb as usize, _dl, _aff, 3);

    // Compute line intersections
    for i in 0..4 {
        _p[i] = match qr_line_isect(&l[i & 1], &l[2 + (i >> 1)]) {
            Some(p) => p,
            None => return -1,
        };
        // It's plausible for points to be somewhat outside the image, but too far
        // and too much of the pattern will be gone for it to be decodable
        let p_i = &_p[i];
        if p_i[0] < (-_width << QR_FINDER_SUBPREC)
            || p_i[0] >= ((_width << QR_FINDER_SUBPREC) + 1)
            || p_i[1] < (-_height << QR_FINDER_SUBPREC)
            || p_i[1] >= ((_height << QR_FINDER_SUBPREC) + 1)
        {
            return -1;
        }
    }

    // By default, use the edge intersection point for the bottom-right corner
    let mut brx = _p[3][0];
    let mut bry = _p[3][1];

    // However, if our average version estimate is greater than 1, try to search for an alignment pattern
    let version4 = _ul.eversion[0] + _ul.eversion[1] + _ur.eversion[0] + _dl.eversion[1];
    if version4 > 4 {
        let mut cell = qr_hom_cell::default();
        let mut p3 = qr_point::default();
        let dim = 17 + version4;
        qr_hom_cell_init(
            &mut cell,
            0,
            0,
            dim - 1,
            0,
            0,
            dim - 1,
            dim - 1,
            dim - 1,
            _p[0][0],
            _p[0][1],
            _p[1][0],
            _p[1][1],
            _p[2][0],
            _p[2][1],
            _p[3][0],
            _p[3][1],
        );
        if qr_alignment_pattern_search(&mut p3, &cell, dim - 7, dim - 7, 4, img, _width, _height)
            >= 0
        {
            // We do need four points in a square to initialize our homography,
            // so project the point from the alignment center to the corner of the code area
            match qr_hom_project_alignment_to_corner(_p, &p3, dim) {
                Ok((x, y)) => {
                    brx = x;
                    bry = y;
                }
                Err(_) => return -1,
            }
        }
    }

    // Now we have four points that map to a square: initialize the projection
    qr_hom_init(
        _hom,
        _p[0][0],
        _p[0][1],
        _p[1][0],
        _p[1][1],
        _p[2][0],
        _p[2][1],
        brx,
        bry,
        QR_HOM_BITS,
    );

    0
}

/// Finish a partial projection, converting from homogeneous coordinates to the
/// normal 2-D representation.
/// In loops, we can avoid many multiplies by computing the homogeneous _x, _y,
/// and _w incrementally, but we cannot avoid the divisions, done here.
pub(crate) fn qr_hom_fproject(
    _hom: &qr_hom,
    mut _x: c_int,
    mut _y: c_int,
    mut _w: c_int,
) -> qr_point {
    if _w == 0 {
        [
            if _x < 0 { c_int::MIN } else { c_int::MAX },
            if _y < 0 { c_int::MIN } else { c_int::MAX },
        ]
    } else {
        if _w < 0 {
            _x = -_x;
            _y = -_y;
            _w = -_w;
        }
        [qr_divround(_x, _w) + _hom.x0, qr_divround(_y, _w) + _hom.y0]
    }
}

/// All the information we've collected about a finder pattern in the current
/// configuration.
#[derive(Default)]
pub(crate) struct qr_finder {
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

/// Estimates the size of a module after classifying the edge points.
///
/// _width:  The distance between UL and UR in the square domain.
/// _height: The distance between UL and DL in the square domain.
///
/// Returns 0 on success, or -1 if the module size or version could not be estimated.
///
/// # Safety
/// This function is unsafe because it dereferences the raw _f pointer.
pub(crate) unsafe fn qr_finder_estimate_module_size_and_version(
    _f: &mut qr_finder,
    _width: c_int,
    _height: c_int,
) -> c_int {
    let mut offs = qr_point::default();
    let mut sums: [c_int; 4] = [0; 4];
    let mut nsums: [c_int; 4] = [0; 4];

    for e in 0..4 {
        if _f.nedge_pts[e] > 0 {
            let edge_pts = _f.edge_pts[e];
            let mut n = _f.nedge_pts[e];

            // Average the samples for this edge, dropping the top and bottom 25%
            let mut sum: c_int = 0;
            for i in (n >> 2)..(n - (n >> 2)) {
                sum = sum.wrapping_add((*edge_pts.offset(i as isize)).extent);
            }
            n -= (n >> 2) << 1;
            let mean = qr_divround(sum, n);
            offs[e >> 1] = offs[e >> 1].wrapping_add(mean);
            sums[e] = sum;
            nsums[e] = n;
        } else {
            nsums[e] = 0;
            sums[e] = 0;
        }
    }

    // If we have samples on both sides of an axis, refine our idea of where the
    // unprojected finder center is located.
    if _f.nedge_pts[0] > 0 && _f.nedge_pts[1] > 0 {
        _f.o[0] = _f.o[0].wrapping_sub(offs[0] >> 1);
        sums[0] = sums[0].wrapping_sub((offs[0].wrapping_mul(nsums[0])) >> 1);
        sums[1] = sums[1].wrapping_sub((offs[0].wrapping_mul(nsums[1])) >> 1);
    }
    if _f.nedge_pts[2] > 0 && _f.nedge_pts[3] > 0 {
        _f.o[1] = _f.o[1].wrapping_sub(offs[1] >> 1);
        sums[2] = sums[2].wrapping_sub((offs[1].wrapping_mul(nsums[2])) >> 1);
        sums[3] = sums[3].wrapping_sub((offs[1].wrapping_mul(nsums[3])) >> 1);
    }

    // We must have _some_ samples along each axis... if we don't, our transform
    // must be pretty severely distorting the original square (e.g., with
    // coordinates so large as to cause overflow).
    let mut nusize = nsums[0].wrapping_add(nsums[1]);
    if nusize <= 0 {
        return -1;
    }

    // The module size is 1/3 the average edge extent.
    nusize = nusize.wrapping_mul(3);
    let mut usize = sums[1].wrapping_sub(sums[0]);
    usize = ((usize.wrapping_shl(1)).wrapping_add(nusize)) / (nusize.wrapping_shl(1));
    if usize <= 0 {
        return -1;
    }

    // Now estimate the version directly from the module size and the distance
    // between the finder patterns.
    // This is done independently using the extents along each axis.
    // If either falls significantly outside the valid range (1 to 40), reject the
    // configuration.
    let uversion = (_width.wrapping_sub(usize.wrapping_mul(8))) / (usize.wrapping_shl(2));
    if !(1..=40 + QR_LARGE_VERSION_SLACK).contains(&uversion) {
        return -1;
    }

    // Now do the same for the other axis.
    let mut nvsize = nsums[2].wrapping_add(nsums[3]);
    if nvsize <= 0 {
        return -1;
    }
    nvsize = nvsize.wrapping_mul(3);
    let mut vsize = sums[3].wrapping_sub(sums[2]);
    vsize = ((vsize.wrapping_shl(1)).wrapping_add(nvsize)) / (nvsize.wrapping_shl(1));
    if vsize <= 0 {
        return -1;
    }

    let vversion = (_height.wrapping_sub(vsize.wrapping_mul(8))) / (vsize.wrapping_shl(2));
    if !(1..=40 + QR_LARGE_VERSION_SLACK).contains(&vversion) {
        return -1;
    }

    // If the estimated version using extents along one axis is significantly
    // different than the estimated version along the other axis, then the axes
    // have significantly different scalings (relative to the grid).
    // This can happen, e.g., when we have multiple adjacent QR codes, and we've
    // picked two finder patterns from one and the third finder pattern from
    // another
    if uversion.wrapping_sub(vversion).abs() > QR_LARGE_VERSION_SLACK {
        return -1;
    }

    _f.size[0] = usize;
    _f.size[1] = vsize;

    // We intentionally do not compute an average version from the sizes along
    // both axes.
    // In the presence of projective distortion, one of them will be much more
    // accurate than the other.
    _f.eversion[0] = uversion;
    _f.eversion[1] = vversion;

    0
}

/// Eliminate outliers from the classified edge points with RANSAC.
///
/// Uses the RANSAC (RANdom SAmple Consensus) algorithm to identify inliers
/// among the edge points and eliminate outliers.
///
/// # Safety
/// This function is unsafe because it dereferences raw pointers.
pub(crate) unsafe fn qr_finder_ransac(
    _f: &mut qr_finder,
    _hom: &qr_aff,
    rng: &mut ChaCha8Rng,
    _e: c_int,
) {
    let edge_pts = _f.edge_pts[_e as usize];
    let n = _f.nedge_pts[_e as usize];
    let mut best_ninliers = 0;

    if n > 1 {
        // 17 iterations is enough to guarantee an outlier-free sample with more
        // than 99% probability given as many as 50% outliers.
        let mut max_iters = 17;
        let mut iter_count = 0;

        while iter_count < max_iters {
            iter_count += 1;

            // Pick two random points on this edge
            let p0i = rng.random_range(0..n) as usize;
            let mut p1i = rng.random_range(0..n - 1) as usize;
            if p1i >= p0i {
                p1i += 1;
            }

            let p0 = &(*edge_pts.add(p0i)).pos;
            let p1 = &(*edge_pts.add(p1i)).pos;

            // If the corresponding line is not within 45 degrees of the proper
            // orientation in the square domain, reject it outright.
            // This can happen, e.g., when highly skewed orientations cause points to
            // be misclassified into the wrong edge.
            let mut q0 = qr_aff_unproject(_hom, p0[0], p0[1]);
            let mut q1 = qr_aff_unproject(_hom, p1[0], p1[1]);
            qr_point_translate(&mut q0, -_f.o[0], -_f.o[1]);
            qr_point_translate(&mut q1, -_f.o[0], -_f.o[1]);

            if (q0[(_e >> 1) as usize] - q1[(_e >> 1) as usize]).abs()
                > (q0[(1 - (_e >> 1)) as usize] - q1[(1 - (_e >> 1)) as usize]).abs()
            {
                continue;
            }

            // Identify the other edge points which are inliers.
            // The squared distance should be distributed as a Chi^2 distribution
            // with one degree of freedom, which means for a 95% confidence the
            // point should lie within a factor 3.8414588 ~= 4 times the expected
            // variance of the point locations.
            // We grossly approximate the standard deviation as 1 pixel in one
            // direction, and 0.5 pixels in the other (because we average two
            // coordinates).
            let thresh =
                qr_isqrt(qr_point_distance2(p0, p1) << (2 * QR_FINDER_SUBPREC + 1)) as c_int;
            let mut ninliers = 0;

            for j in 0..n {
                if qr_point_ccw(p0, p1, &(*edge_pts.offset(j as isize)).pos).abs() <= thresh {
                    (*edge_pts.offset(j as isize)).extent |= 1;
                    ninliers += 1;
                } else {
                    (*edge_pts.offset(j as isize)).extent &= !1;
                }
            }

            if ninliers > best_ninliers {
                for j in 0..n {
                    (*edge_pts.offset(j as isize)).extent <<= 1;
                }
                best_ninliers = ninliers;

                // The actual number of iterations required is
                //   log(1-alpha)/log(1-r*r),
                // where alpha is the required probability of taking a sample with
                // no outliers (e.g., 0.99) and r is the estimated ratio of inliers
                // (e.g. ninliers/n).
                // This is just a rough (but conservative) approximation, but it
                // should be good enough to stop the iteration early when we find
                // a good set of inliers.
                if ninliers > n >> 1 {
                    max_iters = (67 * n - 63 * ninliers - 1) / (n << 1);
                }
            }
        }

        // Now collect all the inliers at the beginning of the list
        let mut i = 0;
        let mut j = 0;
        while j < best_ninliers {
            if (*edge_pts.offset(i as isize)).extent & 2 != 0 {
                if j < i {
                    let tmp = *edge_pts.offset(i as isize);
                    *edge_pts.offset(j as isize) = *edge_pts.offset(i as isize);
                    *edge_pts.offset(i as isize) = tmp;
                }
                j += 1;
            }
            i += 1;
        }
    }

    _f.ninliers[_e as usize] = best_ninliers;
}

/// Perform a least-squares line fit to an edge of a finder pattern using the inliers found by RANSAC.
///
/// # Safety
/// This function is unsafe because it dereferences raw pointers.
pub(crate) unsafe fn qr_line_fit_finder_edge(
    _l: &mut qr_line,
    _f: &qr_finder,
    _e: c_int,
    _res: c_int,
) -> c_int {
    let npts = _f.ninliers[_e as usize] as usize;
    if npts < 2 {
        return -1;
    }

    let mut pts = vec![qr_point::default(); npts];
    let edge_pts = _f.edge_pts[_e as usize];
    let edge_pts = from_raw_parts(edge_pts, npts);

    for (point, edge_point) in pts.iter_mut().zip(edge_pts) {
        point[0] = edge_point.pos[0];
        point[1] = edge_point.pos[1];
    }

    qr_line_fit_points(_l, &mut pts, npts, _res);

    // Make sure the center of the finder pattern lies in the positive halfspace of the line.
    qr_line_orient(_l, (*_f.c).pos[0], (*_f.c).pos[1]);

    0
}

/// Perform a least-squares line fit to a pair of common finder edges using the inliers found by RANSAC.
///
/// Unlike a normal edge fit, we guarantee that this one succeeds by creating at
/// least one point on each edge using the estimated module size if it has no inliers.
///
/// # Safety
/// This function is unsafe because it dereferences raw pointers.
pub(crate) unsafe fn qr_line_fit_finder_pair(
    _l: &mut qr_line,
    _aff: &qr_aff,
    _f0: &qr_finder,
    _f1: &qr_finder,
    _e: c_int,
) {
    let mut n0 = _f0.ninliers[_e as usize] as usize;
    let n1 = _f1.ninliers[_e as usize] as usize;

    // We could write a custom version of qr_line_fit_points that accesses
    // edge_pts directly, but this saves on code size and doesn't measurably slow things down.
    let npts = usize::max(n0, 1) + usize::max(n1, 1);
    let mut pts = vec![qr_point::default(); npts];

    if n0 > 0 {
        let edge_pts = _f0.edge_pts[_e as usize];
        let edge_pts = from_raw_parts(edge_pts, n0);
        for (point, edge_point) in pts.iter_mut().zip(edge_pts).take(n0) {
            point[0] = edge_point.pos[0];
            point[1] = edge_point.pos[1];
        }
    } else {
        let mut q: qr_point = [_f0.o[0], _f0.o[1]];
        q[(_e >> 1) as usize] += _f0.size[(_e >> 1) as usize] * (2 * (_e & 1) - 1);
        pts[0] = qr_aff_project(_aff, q[0], q[1]);
        n0 += 1;
    }

    if n1 > 0 {
        let edge_pts = _f1.edge_pts[_e as usize];
        let edge_pts = from_raw_parts(edge_pts, n1);
        for (point, edge_point) in pts[n0..].iter_mut().zip(edge_pts).take(n1) {
            point[0] = edge_point.pos[0];
            point[1] = edge_point.pos[1];
        }
    } else {
        let mut q: qr_point = [_f1.o[0], _f1.o[1]];
        q[(_e >> 1) as usize] += _f1.size[(_e >> 1) as usize] * (2 * (_e & 1) - 1);
        pts[n0] = qr_aff_project(_aff, q[0], q[1]);
    }

    qr_line_fit_points(_l, &mut pts, npts, _aff.res);

    // Make sure at least one finder center lies in the positive halfspace.
    qr_line_orient(_l, (*_f0.c).pos[0], (*_f0.c).pos[1]);
}

/// Map from the image (at subpel resolution) into the square domain.
/// Returns `Err(point)` with infinity values if the point went to infinity,
/// otherwise returns `Ok(point)` with the projected coordinates.
pub(crate) fn qr_hom_unproject(
    _hom: &qr_hom,
    mut _x: c_int,
    mut _y: c_int,
) -> Result<qr_point, qr_point> {
    _x = _x.wrapping_sub(_hom.x0);
    _y = _y.wrapping_sub(_hom.y0);
    let mut x = _hom.inv[0][0]
        .wrapping_mul(_x)
        .wrapping_add(_hom.inv[0][1].wrapping_mul(_y));
    let mut y = _hom.inv[1][0]
        .wrapping_mul(_x)
        .wrapping_add(_hom.inv[1][1].wrapping_mul(_y));
    let mut w = (_hom.inv[2][0]
        .wrapping_mul(_x)
        .wrapping_add(_hom.inv[2][1].wrapping_mul(_y))
        .wrapping_add(_hom.inv22)
        .wrapping_add(1 << (_hom.res - 1)))
        >> _hom.res;
    if w == 0 {
        Err([
            if x < 0 { c_int::MIN } else { c_int::MAX },
            if y < 0 { c_int::MIN } else { c_int::MAX },
        ])
    } else {
        if w < 0 {
            x = -x;
            y = -y;
            w = -w;
        }
        Ok([qr_divround(x, w), qr_divround(y, w)])
    }
}

/// Bit reading buffer for QR code data
///
/// Bit reading code adapted from libogg/libtheora
/// Portions (C) Xiph.Org Foundation 1994-2008, BSD-style license.
pub(crate) struct qr_pack_buf<'a> {
    buf: &'a [c_uchar],
    endbyte: c_int,
    endbit: c_int,
}

/// Read bits from the pack buffer
///
/// Assumes 0 <= _bits <= 16
/// Returns the read value, or `None` if there aren't enough bits available
pub(crate) fn qr_pack_buf_read(_b: &mut qr_pack_buf, _bits: c_int) -> Option<c_int> {
    let m = 16 - _bits;
    let bits = _bits + _b.endbit;
    let storage = _b.buf.len() as c_int;
    let d = storage - _b.endbyte;

    if d <= 2 {
        // Not the main path
        if d * 8 < bits {
            _b.endbyte += bits >> 3;
            _b.endbit = bits & 7;
            return None;
        }
        // Special case to avoid reading p[0] below, which might be past the end of
        // the buffer; also skips some useless accounting
        else if bits == 0 {
            return Some(0);
        }
    }

    let idx = _b.endbyte as usize;
    let mut ret = (c_uint::from(_b.buf[idx]) << (8 + _b.endbit)) as c_uint;
    if bits > 8 {
        ret |= c_uint::from(_b.buf[idx + 1]) << _b.endbit;
        if bits > 16 {
            ret |= c_uint::from(_b.buf[idx + 2]) >> (8 - _b.endbit);
        }
    }
    _b.endbyte += bits >> 3;
    _b.endbit = bits & 7;
    Some(((ret & 0xFFFF) >> m) as c_int)
}

/// Get the number of bits available to read from the pack buffer
pub(crate) fn qr_pack_buf_avail(_b: &qr_pack_buf) -> c_int {
    let storage = _b.buf.len() as c_int;
    ((storage - _b.endbyte) << 3) - _b.endbit
}

/// Calculate the number of codewords in a QR code of a given version
///
/// This is a compact calculation that avoids a lookup table.
/// Returns the total number of data and error correction codewords.
pub(crate) fn qr_code_ncodewords(_version: c_uint) -> usize {
    if _version == 1 {
        return 26;
    }
    let nalign = (_version / 7) + 2;
    (((_version << 4) * (_version + 8) - (5 * nalign) * (5 * nalign - 2)
        + 36 * c_uint::from(_version < 7)
        + 83)
        >> 3) as usize
}

/// Clusters adjacent lines into groups that are large enough to be crossing a
/// finder pattern (relative to their length).
///
/// # Parameters
/// - `lines`: The list of lines to cluster.
///   Horizontal lines must be sorted in ascending order by Y coordinate, with ties broken by X coordinate.
///   Vertical lines must be sorted in ascending order by X coordinate, with ties broken by Y coordinate.
/// - `_v`: 0 for horizontal lines, or 1 for vertical lines.
///
/// # Returns
/// A ClusteredLines structure containing the lines and their clusters.
fn qr_finder_cluster_lines(lines: Vec<qr_finder_line>, _v: c_int) -> ClusteredLines {
    let nlines = lines.len();
    let mut clustered = ClusteredLines::new(lines);

    if nlines < 2 {
        return clustered;
    }

    // Allocate mark array to track which lines have been clustered
    let mut mark = vec![0u8; nlines];
    // Buffer to store line indices for each cluster as we build it
    let mut neighbor_indices = Vec::with_capacity(nlines);

    for i in 0..(nlines - 1) {
        if mark[i] != 0 {
            continue;
        }

        let cluster_start = neighbor_indices.len();
        neighbor_indices.push(i);
        let mut len = clustered.line(i).len as usize;

        for (j, mj) in mark.iter().enumerate().take(nlines).skip(i + 1) {
            if *mj != 0 {
                continue;
            }

            let last_idx = *neighbor_indices.last().unwrap();
            let a = clustered.line(last_idx);
            let b = clustered.line(j);

            // The clustering threshold is proportional to the size of the lines,
            // since minor noise in large areas can interrupt patterns more easily
            // at high resolutions.
            let thresh = (a.len + 7) >> 2;

            // Check if lines are too far apart perpendicular to their direction
            if (a.pos[(1 - _v) as usize] - b.pos[(1 - _v) as usize]).abs() > thresh {
                break;
            }

            // Check if lines are too far apart along their direction
            if (a.pos[_v as usize] - b.pos[_v as usize]).abs() > thresh {
                continue;
            }

            // Check if line ends are too far apart
            if (a.pos[_v as usize] + a.len - b.pos[_v as usize] - b.len).abs() > thresh {
                continue;
            }

            // Check beginning offset alignment
            if a.boffs > 0
                && b.boffs > 0
                && (a.pos[_v as usize] - a.boffs - b.pos[_v as usize] + b.boffs).abs() > thresh
            {
                continue;
            }

            // Check ending offset alignment
            if a.eoffs > 0
                && b.eoffs > 0
                && (a.pos[_v as usize] + a.len + a.eoffs - b.pos[_v as usize] - b.len - b.eoffs)
                    .abs()
                    > thresh
            {
                continue;
            }

            neighbor_indices.push(j);
            len += b.len as usize;
        }

        let nneighbors = neighbor_indices.len() - cluster_start;

        // We require at least three lines to form a cluster, which eliminates a
        // large number of false positives, saving considerable decoding time.
        // This should still be sufficient for 1-pixel codes with no noise.
        if nneighbors < 3 {
            // Remove the indices we just added since this isn't a valid cluster
            neighbor_indices.truncate(cluster_start);
            continue;
        }

        // The expected number of lines crossing a finder pattern is equal to their
        // average length.
        // We accept the cluster if size is at least 1/3 their average length (this
        // is a very small threshold, but was needed for some test images).
        len = ((len << 1) + nneighbors) / (nneighbors << 1);
        if nneighbors * (5 << QR_FINDER_SUBPREC) >= len {
            // Mark all lines in this cluster
            for &idx in &neighbor_indices[cluster_start..] {
                mark[idx] = 1;
            }

            // Create the cluster
            let cluster = Cluster {
                line_indices: neighbor_indices[cluster_start..].to_vec(),
            };
            clustered.clusters.push(cluster);
        } else {
            // Remove the indices since this cluster didn't meet the threshold
            neighbor_indices.truncate(cluster_start);
        }
    }

    clustered
}

enum Direction {
    Horizontal,
    Vertical,
}

/// Gets the coordinates of the edge points based on the lines contained in the
/// given list of cluster indices from a ClusteredLines structure.
///
/// Only the edge point position is initialized.
/// The edge label and extent are set by qr_finder_edge_pts_aff_classify()
/// or qr_finder_edge_pts_hom_classify().
fn qr_finder_get_edge_pts(
    clustered: &ClusteredLines,
    cluster_indices: &[usize],
    dir: Direction,
) -> Vec<qr_finder_edge_pt> {
    let mut edge_pts = vec![];
    let v_idx = match dir {
        Direction::Horizontal => 0,
        Direction::Vertical => 1,
    };

    for &cluster_idx in cluster_indices {
        for l in clustered.cluster_lines(cluster_idx) {
            // Add beginning offset edge point if present
            if l.boffs > 0 {
                let mut pt = qr_finder_edge_pt::default();
                pt.pos[0] = l.pos[0];
                pt.pos[1] = l.pos[1];
                pt.pos[v_idx] -= l.boffs;
                edge_pts.push(pt);
            }

            // Add ending offset edge point if present
            if l.eoffs > 0 {
                let mut pt = qr_finder_edge_pt::default();
                pt.pos[0] = l.pos[0];
                pt.pos[1] = l.pos[1];
                pt.pos[v_idx] += l.len + l.eoffs;
                edge_pts.push(pt);
            }
        }
    }
    edge_pts
}

/// Finds horizontal clusters that cross corresponding vertical clusters,
/// presumably corresponding to a finder center.
///
/// # Parameters
/// - `hclustered`: The horizontal lines and their clusters.
/// - `vclustered`: The vertical lines and their clusters.
///
/// # Returns
/// Vec of putative finder centers, each owning its edge points.
pub(crate) fn qr_finder_find_crossings(
    hclustered: &ClusteredLines,
    vclustered: &ClusteredLines,
) -> Vec<qr_finder_center> {
    let nhclusters = hclustered.cluster_count();
    let nvclusters = vclustered.cluster_count();
    let mut hmark = vec![0u8; nhclusters];
    let mut vmark = vec![0u8; nvclusters];

    let mut centers = vec![];

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

    for i in 0..nhclusters {
        if hmark[i] != 0 {
            continue;
        }

        // Get the middle line from this horizontal cluster
        let h_mid_idx = hclustered.cluster(i).line_indices.len() >> 1;
        let a = hclustered.cluster_line(i, h_mid_idx);
        let mut y = 0;
        let mut vneighbor_indices = vec![];

        // Find vertical clusters that cross this horizontal cluster
        for (j, vj) in vmark.iter_mut().enumerate() {
            if *vj != 0 {
                continue;
            }

            // Get the middle line from this vertical cluster
            let v_mid_idx = vclustered.cluster(j).line_indices.len() >> 1;
            let b = vclustered.cluster_line(j, v_mid_idx);

            if qr_finder_lines_are_crossing(a, b) {
                *vj = 1;
                y += (b.pos[1] << 1) + b.len;
                if b.boffs > 0 && b.eoffs > 0 {
                    y += b.eoffs - b.boffs;
                }
                vneighbor_indices.push(j);
            }
        }

        if !vneighbor_indices.is_empty() {
            let mut x = (a.pos[0] << 1) + a.len;
            if a.boffs > 0 && a.eoffs > 0 {
                x += a.eoffs - a.boffs;
            }
            let mut hneighbor_indices = vec![i];

            let j_mid = vneighbor_indices.len() >> 1;
            let v_mid_line_idx = vclustered
                .cluster(vneighbor_indices[j_mid])
                .line_indices
                .len()
                >> 1;
            let b = vclustered.cluster_line(vneighbor_indices[j_mid], v_mid_line_idx);

            // Find additional horizontal clusters that cross the vertical clusters
            for (j, hj) in hmark.iter_mut().enumerate().take(nhclusters).skip(i + 1) {
                if *hj != 0 {
                    continue;
                }

                let h_mid_idx = hclustered.cluster(j).line_indices.len() >> 1;
                let a = hclustered.cluster_line(j, h_mid_idx);

                if qr_finder_lines_are_crossing(a, b) {
                    *hj = 1;
                    x += (a.pos[0] << 1) + a.len;
                    if a.boffs > 0 && a.eoffs > 0 {
                        x += a.eoffs - a.boffs;
                    }
                    hneighbor_indices.push(j);
                }
            }

            centers.push(qr_finder_center {
                pos: [
                    ((x as usize + hneighbor_indices.len()) / (hneighbor_indices.len() << 1))
                        as i32,
                    ((y as usize + vneighbor_indices.len()) / (vneighbor_indices.len() << 1))
                        as i32,
                ],
                edge_pts: [
                    qr_finder_get_edge_pts(hclustered, &hneighbor_indices, Direction::Horizontal),
                    qr_finder_get_edge_pts(vclustered, &vneighbor_indices, Direction::Vertical),
                ]
                .concat(),
            })
        }
    }

    centers.sort_by(|a, b| {
        b.edge_pts
            .len()
            .cmp(&a.edge_pts.len())
            .then_with(|| a.pos[1].cmp(&b.pos[1]))
            .then_with(|| a.pos[0].cmp(&b.pos[0]))
    });

    centers
}

/// Locate QR finder pattern centers from scanned lines
///
/// Clusters horizontal and vertical lines that cross finder patterns,
/// then locates the centers where horizontal and vertical clusters intersect.
pub(crate) fn qr_finder_centers_locate(
    reader: &mut qr_reader,
    _width: c_int,
    _height: c_int,
) -> Vec<qr_finder_center> {
    // Take ownership of the horizontal lines and cluster them
    let hlines = std::mem::take(&mut reader.finder_lines[0].lines);
    let hclustered = qr_finder_cluster_lines(hlines, 0);

    // We need vertical lines to be sorted by X coordinate, with ties broken by Y
    // coordinate, for clustering purposes.
    // We scan the image in the opposite order, so sort the lines we found here.
    let mut vlines = std::mem::take(&mut reader.finder_lines[1].lines);
    if vlines.len() > 1 {
        vlines.sort_by(|a, b| {
            a.pos[0]
                .cmp(&b.pos[0])
                .then_with(|| a.pos[1].cmp(&b.pos[1]))
        });
    }
    let vclustered = qr_finder_cluster_lines(vlines, 1);

    // Find line crossings among the clusters
    if hclustered.cluster_count() >= 3 && vclustered.cluster_count() >= 3 {
        qr_finder_find_crossings(&hclustered, &vclustered)
    } else {
        Vec::new()
    }
}

/// Mark a given rectangular region as belonging to the function pattern
///
/// Function patterns include finder patterns, timing patterns, and alignment patterns.
/// We store bits column-wise, since that's how they're read out of the grid.
fn qr_sampling_grid_fp_mask_rect(
    _grid: &mut qr_sampling_grid,
    _dim: c_int,
    _u: c_int,
    _v: c_int,
    _w: c_int,
    _h: c_int,
) {
    let stride = (_dim + QR_INT_BITS - 1) >> QR_INT_LOGBITS;
    for j in _u..(_u + _w) {
        for i in _v..(_v + _h) {
            let idx = (j * stride + (i >> QR_INT_LOGBITS)) as usize;
            _grid.fpmask[idx] |= 1 << (i & (QR_INT_BITS - 1));
        }
    }
}

/// Determine if a given grid location is inside a function pattern
///
/// Returns whether the location (_u, _v) is part of a function pattern.
fn qr_sampling_grid_is_in_fp(_grid: &qr_sampling_grid, _dim: usize, _u: c_int, _v: c_int) -> bool {
    let idx = _u as usize * ((_dim + QR_INT_BITS as usize - 1) >> QR_INT_LOGBITS as usize)
        + (_v as usize >> QR_INT_LOGBITS as usize);
    ((_grid.fpmask[idx]) >> (_v & (QR_INT_BITS - 1)) & 1) == 1
}

/// The spacing between alignment patterns after the second for versions >= 7
///
/// We could compact this more, but the code to access it would eliminate the gains.
pub(crate) static QR_ALIGNMENT_SPACING: [c_uchar; 34] = [
    16, 18, 20, 22, 24, 26, 28, 20, 22, 24, 24, 26, 28, 28, 22, 24, 24, 26, 26, 28, 28, 24, 24, 26,
    26, 26, 28, 28, 24, 26, 26, 26, 28, 28,
];

/// Bulk data for the number of parity bytes per Reed-Solomon block
pub(crate) static QR_RS_NPAR_VALS: [c_uchar; 71] = [
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
pub(crate) static QR_RS_NPAR_OFFS: [c_uchar; 40] = [
    0, 4, 19, 55, 15, 28, 37, 12, 51, 39, 59, 62, 10, 24, 22, 41, 31, 44, 7, 65, 47, 33, 67, 67,
    48, 32, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67,
];

/// The number of Reed-Solomon blocks for each version and error correction level
pub(crate) static QR_RS_NBLOCKS: [[c_uchar; 4]; 40] = [
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

pub(crate) fn qr_cmp_edge_pt(a: &qr_finder_edge_pt, b: &qr_finder_edge_pt) -> Ordering {
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
pub(crate) unsafe fn qr_finder_edge_pts_aff_classify(_f: &mut qr_finder, _aff: &qr_aff) {
    let c = _f.c;
    for item in _f.nedge_pts.iter_mut() {
        *item = 0;
    }
    for edge_pt in (*c).edge_pts.iter_mut() {
        let mut q = qr_aff_unproject(_aff, edge_pt.pos[0], edge_pt.pos[1]);
        qr_point_translate(&mut q, -_f.o[0], -_f.o[1]);
        let d = c_int::from(q[1].abs() > q[0].abs());
        let e = d << 1 | c_int::from(q[d as usize] >= 0);
        _f.nedge_pts[e as usize] += 1;
        edge_pt.edge = e;
        edge_pt.extent = q[d as usize];
    }

    (*c).edge_pts.sort_by(qr_cmp_edge_pt);
    _f.edge_pts[0] = (*c).edge_pts.as_mut_ptr();
    for e in 1.._f.edge_pts.len() {
        _f.edge_pts[e] = _f.edge_pts[e - 1].add(_f.nedge_pts[e - 1] as usize);
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
pub(crate) unsafe fn qr_finder_edge_pts_hom_classify(_f: &mut qr_finder, _hom: &qr_hom) {
    let c = _f.c;
    for item in _f.nedge_pts.iter_mut() {
        *item = 0;
    }
    for edge_pt in (*c).edge_pts.iter_mut() {
        match qr_hom_unproject(_hom, edge_pt.pos[0], edge_pt.pos[1]) {
            Ok(mut q) => {
                // Successful projection
                qr_point_translate(&mut q, -_f.o[0], -_f.o[1]);
                let d = c_int::from(q[1].abs() > q[0].abs());
                let e = d << 1 | c_int::from(q[d as usize] >= 0);
                _f.nedge_pts[e as usize] += 1;
                edge_pt.edge = e;
                edge_pt.extent = q[d as usize];
            }
            Err(q) => {
                // Projection failed (went to infinity)
                edge_pt.edge = 4;
                edge_pt.extent = q[0];
            }
        }
    }

    (*c).edge_pts.sort_by(qr_cmp_edge_pt);
    _f.edge_pts[0] = (*c).edge_pts.as_mut_ptr();
    for e in 1.._f.edge_pts.len() {
        _f.edge_pts[e] = _f.edge_pts[e - 1].add(_f.nedge_pts[e - 1] as usize);
    }
}

#[allow(clippy::too_many_arguments)]
pub(crate) unsafe fn qr_finder_quick_crossing_check(
    _img: &[u8],
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

    if (c_int::from(_img[(_y0 * _width + _x0) as usize] == 0)) != _v
        || (c_int::from(_img[(_y1 * _width + _x1) as usize] == 0)) != _v
    {
        return 1;
    }
    if (c_int::from(_img[(((_y0 + _y1) >> 1) * _width + ((_x0 + _x1) >> 1)) as usize] == 0)) == _v {
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
///
/// # Returns
/// - `Ok(dv)` with the computed step in v direction on success
/// - `Err(-1)` if the line is too steep (>45 degrees from horizontal/vertical)
pub(crate) fn qr_aff_line_step(
    aff: &qr_aff,
    line: &qr_line,
    v: c_int,
    du: c_int,
) -> Result<c_int, c_int> {
    let l0 = line[0];
    let l1 = line[1];

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
        return Err(-1);
    }

    n *= -du;
    let dv_result = qr_divround(n, d);

    if dv_result.abs() >= du {
        return Err(-1);
    }

    Ok(dv_result)
}

/// Calculate Hamming distance between two values
///
/// Counts the number of bit positions where the values differ,
/// up to a maximum of maxdiff.
pub(crate) fn qr_hamming_dist(y1: c_uint, y2: c_uint, maxdiff: c_int) -> c_int {
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

/// Decode the version number from a QR code's version information
///
/// Reads the 18-bit version information pattern (6 data bits + 12 parity bits)
/// from the specified direction and uses BCH error correction to recover the version.
///
/// # Safety
/// This function is unsafe because it dereferences raw pointers and accesses image data.
pub(crate) unsafe fn qr_finder_version_decode(
    _f: &qr_finder,
    _hom: &qr_hom,
    _img: &[u8],
    _width: c_int,
    _height: c_int,
    _dir: c_int,
) -> c_int {
    let mut q = qr_point::default();
    let mut v: c_uint = 0;

    // Calculate starting point for version info block
    q[_dir as usize] = _f.o[_dir as usize] - 7 * _f.size[_dir as usize];
    q[(1 - _dir) as usize] = _f.o[(1 - _dir) as usize] - 3 * _f.size[(1 - _dir) as usize];

    // Project the starting point through the homography
    let mut x0 = _hom.fwd[0][0] * q[0] + _hom.fwd[0][1] * q[1];
    let mut y0 = _hom.fwd[1][0] * q[0] + _hom.fwd[1][1] * q[1];
    let mut w0 = _hom.fwd[2][0] * q[0] + _hom.fwd[2][1] * q[1] + _hom.fwd22;

    // Calculate increments for stepping through version info grid
    let dxi = _hom.fwd[0][(1 - _dir) as usize] * _f.size[(1 - _dir) as usize];
    let dyi = _hom.fwd[1][(1 - _dir) as usize] * _f.size[(1 - _dir) as usize];
    let dwi = _hom.fwd[2][(1 - _dir) as usize] * _f.size[(1 - _dir) as usize];
    let dxj = _hom.fwd[0][_dir as usize] * _f.size[_dir as usize];
    let dyj = _hom.fwd[1][_dir as usize] * _f.size[_dir as usize];
    let dwj = _hom.fwd[2][_dir as usize] * _f.size[_dir as usize];

    // Read the 6x3 = 18 bits of version information
    let mut k = 0;
    for _i in 0..6 {
        let mut x = x0;
        let mut y = y0;
        let mut w = w0;
        for _j in 0..3 {
            let p = qr_hom_fproject(_hom, x, y, w);
            v |= (qr_img_get_bit(_img, _width, _height, p[0], p[1]) as c_uint) << k;
            k += 1;
            x += dxj;
            y += dyj;
            w += dwj;
        }
        x0 += dxi;
        y0 += dyi;
        w0 += dwi;
    }

    // Use BCH error correction to decode the version
    match bch18_6_correct(v) {
        Ok((corrected, _nerrs)) => (corrected >> 12) as c_int,
        Err(Bch18_6CorrectError::Unrecoverable) => -1,
    }
}

/// Decode the format information from a QR code
///
/// Reads the 15-bit format information pattern (5 data bits + 10 parity bits)
/// from around the three finder patterns and uses BCH(15,5) error correction
/// to decode it. The function tries all combinations of duplicate samples and
/// picks the most popular valid code.
///
/// # Safety
/// This function is unsafe because it dereferences raw pointers and accesses image data.
pub(crate) unsafe fn qr_finder_fmt_info_decode(
    _ul: &qr_finder,
    _ur: &qr_finder,
    _dl: &qr_finder,
    _hom: &qr_hom,
    _img: &[u8],
    _width: c_int,
    _height: c_int,
) -> c_int {
    let mut lo: [c_uint; 2] = [0; 2];
    let mut hi: [c_uint; 2] = [0; 2];

    // Read the bits around the UL corner
    lo[0] = 0;
    let mut u = _ul.o[0] + 5 * _ul.size[0];
    let mut v = _ul.o[1] - 3 * _ul.size[1];
    let mut x = _hom.fwd[0][0] * u + _hom.fwd[0][1] * v;
    let mut y = _hom.fwd[1][0] * u + _hom.fwd[1][1] * v;
    let mut w = _hom.fwd[2][0] * u + _hom.fwd[2][1] * v + _hom.fwd22;
    let mut dx = _hom.fwd[0][1] * _ul.size[1];
    let mut dy = _hom.fwd[1][1] * _ul.size[1];
    let mut dw = _hom.fwd[2][1] * _ul.size[1];

    let mut k = 0;
    let mut i = 0;
    loop {
        // Skip the timing pattern row
        if i != 6 {
            let p = qr_hom_fproject(_hom, x, y, w);
            lo[0] |= (qr_img_get_bit(_img, _width, _height, p[0], p[1]) as c_uint) << k;
            k += 1;
            // Don't advance q in the last iteration... we'll start the next loop from
            // the current position.
            if i >= 8 {
                break;
            }
        }
        x += dx;
        y += dy;
        w += dw;
        i += 1;
    }

    hi[0] = 0;
    dx = -_hom.fwd[0][0] * _ul.size[0];
    dy = -_hom.fwd[1][0] * _ul.size[0];
    dw = -_hom.fwd[2][0] * _ul.size[0];
    while i > 0 {
        i -= 1;
        x += dx;
        y += dy;
        w += dw;
        // Skip the timing pattern column
        if i != 6 {
            let p = qr_hom_fproject(_hom, x, y, w);
            hi[0] |= (qr_img_get_bit(_img, _width, _height, p[0], p[1]) as c_uint) << k;
            k += 1;
        }
    }

    // Read the bits next to the UR corner
    lo[1] = 0;
    u = _ur.o[0] + 3 * _ur.size[0];
    v = _ur.o[1] + 5 * _ur.size[1];
    x = _hom.fwd[0][0] * u + _hom.fwd[0][1] * v;
    y = _hom.fwd[1][0] * u + _hom.fwd[1][1] * v;
    w = _hom.fwd[2][0] * u + _hom.fwd[2][1] * v + _hom.fwd22;
    dx = -_hom.fwd[0][0] * _ur.size[0];
    dy = -_hom.fwd[1][0] * _ur.size[0];
    dw = -_hom.fwd[2][0] * _ur.size[0];
    for k in 0..8 {
        let p = qr_hom_fproject(_hom, x, y, w);
        lo[1] |= (qr_img_get_bit(_img, _width, _height, p[0], p[1]) as c_uint) << k;
        x += dx;
        y += dy;
        w += dw;
    }

    // Read the bits next to the DL corner
    hi[1] = 0;
    u = _dl.o[0] + 5 * _dl.size[0];
    v = _dl.o[1] - 3 * _dl.size[1];
    x = _hom.fwd[0][0] * u + _hom.fwd[0][1] * v;
    y = _hom.fwd[1][0] * u + _hom.fwd[1][1] * v;
    w = _hom.fwd[2][0] * u + _hom.fwd[2][1] * v + _hom.fwd22;
    dx = _hom.fwd[0][1] * _dl.size[1];
    dy = _hom.fwd[1][1] * _dl.size[1];
    dw = _hom.fwd[2][1] * _dl.size[1];
    for k in 8..15 {
        let p = qr_hom_fproject(_hom, x, y, w);
        hi[1] |= (qr_img_get_bit(_img, _width, _height, p[0], p[1]) as c_uint) << k;
        x += dx;
        y += dy;
        w += dw;
    }

    // For each group of bits we have two samples... try them in all combinations
    // and pick the most popular valid code, breaking ties using the number of
    // bit errors.
    let imax = 2 << (hi[0] != hi[1]) as c_int;
    let di = 1 + (lo[0] == lo[1]) as c_int;
    let mut fmt_info: [c_int; 4] = [0; 4];
    let mut count: [c_int; 4] = [0; 4];
    let mut nerrs: [c_int; 4] = [0; 4];
    let mut nfmt_info = 0;

    let mut i = 0;
    while i < imax {
        let mut v = (lo[(i & 1) as usize] | hi[(i >> 1) as usize]) ^ 0x5412;
        let mut ret = bch15_5_correct(&mut v);
        v >>= 10;
        if ret < 0 {
            ret = 4;
        }

        let mut j = 0;
        loop {
            if j >= nfmt_info {
                fmt_info[j] = v as c_int;
                count[j] = 1;
                nerrs[j] = ret;
                nfmt_info += 1;
                break;
            }
            if fmt_info[j] == v as c_int {
                count[j] += 1;
                if ret < nerrs[j] {
                    nerrs[j] = ret;
                }
                break;
            }
            j += 1;
        }
        i += di;
    }

    let mut besti = 0;
    for i in 1..nfmt_info {
        if (nerrs[besti] > 3 && nerrs[i] <= 3)
            || count[i] > count[besti]
            || (count[i] == count[besti] && nerrs[i] < nerrs[besti])
        {
            besti = i;
        }
    }

    if nerrs[besti] < 4 {
        fmt_info[besti]
    } else {
        -1
    }
}

/// Make a data mask for QR code decoding
///
/// Fills a buffer with the specified data mask pattern. The mask is stored
/// column-wise since that's how bits are read out of the QR code grid.
pub(crate) fn qr_make_data_mask(_dim: usize, _pattern: c_int) -> Vec<c_uint> {
    let data_word_count = _dim * ((_dim + QR_INT_BITS as usize - 1) >> QR_INT_LOGBITS as usize);
    let mut _mask = vec![0; data_word_count];
    let stride = (_dim + QR_INT_BITS as usize - 1) >> QR_INT_LOGBITS as usize;

    // Note that we store bits column-wise, since that's how they're read out of the grid.
    match _pattern {
        // Pattern 0: 10101010 (i+j+1&1)
        //            01010101
        //            10101010
        //            01010101
        0 => {
            let mut m: c_uint = 0x55;
            for j in 0.._dim {
                // Replicate byte value across all bytes of c_uint (like memset does)
                let replicated = m * 0x01010101;
                _mask[j * stride..][..stride].fill(replicated);
                m ^= 0xFF;
            }
        }
        // Pattern 1: 11111111 (i+1&1)
        //            00000000
        //            11111111
        //            00000000
        1 => {
            // Replicate byte value across all bytes of c_uint (like memset does)
            _mask[.._dim * stride].fill(0x55555555);
        }
        // Pattern 2: 10010010 ((j+1)%3&1)
        //            10010010
        //            10010010
        //            10010010
        2 => {
            let mut m: c_uint = 0xFF;
            for j in 0.._dim {
                // Replicate byte value across all bytes of c_uint (like memset does)
                let replicated = (m & 0xFF) * 0x01010101;
                _mask[j * stride..][..stride].fill(replicated);
                m = (m << 8) | (m >> 16);
            }
        }
        // Pattern 3: 10010010 ((i+j+1)%3&1)
        //            00100100
        //            01001001
        //            10010010
        3 => {
            let mut mj: c_uint = 0;
            for i in 0..((QR_INT_BITS + 2) / 3) {
                mj |= 1 << (3 * i);
            }
            for j in 0.._dim {
                let mut mi = mj;
                for i in 0..stride {
                    _mask[j * stride + i] = mi;
                    mi = (mi >> (QR_INT_BITS % 3)) | (mi << (3 - (QR_INT_BITS % 3)));
                }
                mj = (mj >> 1) | (mj << 2);
            }
        }
        // Pattern 4: 11100011 ((i>>1)+(j/3)+1&1)
        //            11100011
        //            00011100
        //            00011100
        4 => {
            let mut m: c_uint = 7;
            for j in 0.._dim {
                // Replicate byte value across all bytes of c_uint (like memset does)
                let replicated = ((0xCC ^ (m & 1).wrapping_neg()) & 0xFF) * 0x01010101;
                _mask[j * stride..][..stride].fill(replicated);
                m = (m >> 1) | (m << 5);
            }
        }
        // Pattern 5: 11111111 !((i*j)%6)
        //            10000010
        //            10010010
        //            10101010
        5 => {
            for j in 0.._dim {
                let mut m: c_uint = 0;
                for i in 0..6 {
                    m |= ((((i * j) % 6) == 0) as c_uint) << i;
                }
                let mut i = 6;
                while i < QR_INT_BITS {
                    m |= m << i;
                    i <<= 1;
                }
                for i in 0..stride {
                    _mask[j * stride + i] = m;
                    m = (m >> (QR_INT_BITS % 6)) | (m << (6 - (QR_INT_BITS % 6)));
                }
            }
        }
        // Pattern 6: 11111111 ((i*j)%3+i*j+1&1)
        //            11100011
        //            11011011
        //            10101010
        6 => {
            for j in 0.._dim {
                let mut m: c_uint = 0;
                for i in 0..6 {
                    m |= ((((i * j) % 3 + i * j + 1) & 1) as c_uint) << i;
                }
                let mut i = 6;
                while i < QR_INT_BITS {
                    m |= m << i;
                    i <<= 1;
                }
                for i in 0..stride {
                    _mask[j * stride + i] = m;
                    m = (m >> (QR_INT_BITS % 6)) | (m << (6 - (QR_INT_BITS % 6)));
                }
            }
        }
        // Pattern 7 (default): 10101010 ((i*j)%3+i+j+1&1)
        //                       00011100
        //                       10001110
        //                       01010101
        _ => {
            for j in 0.._dim {
                let mut m: c_uint = 0;
                for i in 0..6 {
                    m |= ((((i * j) % 3 + i + j + 1) & 1) as c_uint) << i;
                }
                let mut i = 6;
                while i < QR_INT_BITS {
                    m |= m << i;
                    i <<= 1;
                }
                for i in 0..stride {
                    _mask[j * stride + i] = m;
                    m = (m >> (QR_INT_BITS % 6)) | (m << (6 - (QR_INT_BITS % 6)));
                }
            }
        }
    }

    _mask
}

/// Initialize a QR sampling grid
///
/// Sets up the sampling grid for reading QR code data bits from an image.
/// This includes creating homography cells, masking function patterns,
/// and locating alignment patterns.
///
/// # Safety
/// This function is unsafe because it allocates memory and dereferences raw pointers.
#[allow(clippy::too_many_arguments)]
unsafe fn qr_sampling_grid_init(
    _grid: &mut qr_sampling_grid,
    _version: c_int,
    _ul_pos: &qr_point,
    _ur_pos: &qr_point,
    _dl_pos: &qr_point,
    _p: &mut [qr_point],
    _img: &[u8],
    _width: c_int,
    _height: c_int,
) {
    let mut base_cell = qr_hom_cell {
        fwd: [[0; 3]; 3],
        x0: 0,
        y0: 0,
        u0: 0,
        v0: 0,
    };
    let mut align_pos: [c_int; 7] = [0; 7];

    let dim = 17 + (_version << 2);
    let nalign = (_version / 7) + 2;

    // Create a base cell to bootstrap the alignment pattern search
    qr_hom_cell_init(
        &mut base_cell,
        0,
        0,
        dim - 1,
        0,
        0,
        dim - 1,
        dim - 1,
        dim - 1,
        _p[0][0],
        _p[0][1],
        _p[1][0],
        _p[1][1],
        _p[2][0],
        _p[2][1],
        _p[3][0],
        _p[3][1],
    );

    // Allocate the 2D array of cells
    let ncells = nalign as usize - 1;
    _grid.cells = vec![vec![qr_hom_cell::default(); ncells]; ncells];

    // Initialize the function pattern mask
    let stride = ((dim + QR_INT_BITS - 1) >> QR_INT_LOGBITS) as usize;
    _grid.fpmask = vec![0; dim as usize * stride];

    // Mask out the finder patterns (and separators and format info bits)
    qr_sampling_grid_fp_mask_rect(&mut *_grid, dim, 0, 0, 9, 9);
    qr_sampling_grid_fp_mask_rect(&mut *_grid, dim, 0, dim - 8, 9, 8);
    qr_sampling_grid_fp_mask_rect(&mut *_grid, dim, dim - 8, 0, 8, 9);

    // Mask out the version number bits
    if _version > 6 {
        qr_sampling_grid_fp_mask_rect(&mut *_grid, dim, 0, dim - 11, 6, 3);
        qr_sampling_grid_fp_mask_rect(&mut *_grid, dim, dim - 11, 0, 3, 6);
    }

    // Mask out the timing patterns
    qr_sampling_grid_fp_mask_rect(&mut *_grid, dim, 9, 6, dim - 17, 1);
    qr_sampling_grid_fp_mask_rect(&mut *_grid, dim, 6, 9, 1, dim - 17);

    // If we have no alignment patterns (e.g., this is a version 1 code), just use
    // the base cell and hope it's good enough
    if _version < 2 {
        _grid.cells[0][0] = base_cell;
    } else {
        let mut q = vec![qr_point::default(); (nalign * nalign) as usize];
        let mut p = vec![qr_point::default(); (nalign * nalign) as usize];

        // Initialize the alignment pattern position list
        align_pos[0] = 6;
        align_pos[(nalign - 1) as usize] = dim - 7;
        if _version > 6 {
            let d = QR_ALIGNMENT_SPACING[(_version - 7) as usize] as c_int;
            for i in (1..(nalign - 1)).rev() {
                align_pos[i as usize] = align_pos[(i + 1) as usize] - d;
            }
        }

        // Three of the corners use a finder pattern instead of a separate alignment pattern
        q[0][0] = 3;
        q[0][1] = 3;
        p[0][0] = _ul_pos[0];
        p[0][1] = _ul_pos[1];
        q[(nalign - 1) as usize][0] = dim - 4;
        q[(nalign - 1) as usize][1] = 3;
        p[(nalign - 1) as usize][0] = _ur_pos[0];
        p[(nalign - 1) as usize][1] = _ur_pos[1];
        q[((nalign - 1) * nalign) as usize][0] = 3;
        q[((nalign - 1) * nalign) as usize][1] = dim - 4;
        p[((nalign - 1) * nalign) as usize][0] = _dl_pos[0];
        p[((nalign - 1) * nalign) as usize][1] = _dl_pos[1];

        // Scan for alignment patterns using a diagonal sweep
        for k in 1..(2 * nalign - 1) {
            let jmax = c_int::min(k, nalign - 1) - (if k == nalign - 1 { 1 } else { 0 });
            let jmin = c_int::max(0, k - (nalign - 1)) + (if k == nalign - 1 { 1 } else { 0 });
            for j in jmin..=jmax {
                let i = jmax - (j - jmin);
                let k_idx = i * nalign + j;
                let u = align_pos[j as usize];
                let v = align_pos[i as usize];
                q[k_idx as usize][0] = u;
                q[k_idx as usize][1] = v;

                // Mask out the alignment pattern
                qr_sampling_grid_fp_mask_rect(&mut *_grid, dim, u - 2, v - 2, 5, 5);

                // Pick a cell to use to govern the alignment pattern search
                let cell = if i > 1 && j > 1 {
                    // Each predictor is basically a straight-line extrapolation from two
                    // neighboring alignment patterns
                    let mut p0 = qr_hom_cell_project(
                        &_grid.cells[(i - 2) as usize][(j - 1) as usize],
                        u,
                        v,
                        0,
                    );
                    let mut p1 = qr_hom_cell_project(
                        &_grid.cells[(i - 2) as usize][(j - 2) as usize],
                        u,
                        v,
                        0,
                    );
                    let mut p2 = qr_hom_cell_project(
                        &_grid.cells[(i - 1) as usize][(j - 2) as usize],
                        u,
                        v,
                        0,
                    );

                    // Take the median of the predictions as the search center
                    // QR_SORT2I implementation using swap
                    if p0[0] > p1[0] {
                        swap(&mut p0[0], &mut p1[0]);
                    }
                    if p0[1] > p1[1] {
                        swap(&mut p0[1], &mut p1[1]);
                    }
                    if p1[0] > p2[0] {
                        swap(&mut p1[0], &mut p2[0]);
                    }
                    if p1[1] > p2[1] {
                        swap(&mut p1[1], &mut p2[1]);
                    }
                    if p0[0] > p1[0] {
                        swap(&mut p0[0], &mut p1[0]);
                    }
                    if p0[1] > p1[1] {
                        swap(&mut p0[1], &mut p1[1]);
                    }

                    // We need a cell that has the target point at a known (u,v) location
                    qr_hom_cell_init(
                        &mut _grid.cells[(i - 1) as usize][(j - 1) as usize],
                        q[(k_idx - nalign - 1) as usize][0],
                        q[(k_idx - nalign - 1) as usize][1],
                        q[(k_idx - nalign) as usize][0],
                        q[(k_idx - nalign) as usize][1],
                        q[(k_idx - 1) as usize][0],
                        q[(k_idx - 1) as usize][1],
                        q[k_idx as usize][0],
                        q[k_idx as usize][1],
                        p[(k_idx - nalign - 1) as usize][0],
                        p[(k_idx - nalign - 1) as usize][1],
                        p[(k_idx - nalign) as usize][0],
                        p[(k_idx - nalign) as usize][1],
                        p[(k_idx - 1) as usize][0],
                        p[(k_idx - 1) as usize][1],
                        p1[0],
                        p1[1],
                    );
                    &_grid.cells[(i - 1) as usize][(j - 1) as usize]
                } else if i > 1 && j > 0 {
                    &_grid.cells[(i - 2) as usize][(j - 1) as usize]
                } else if i > 0 && j > 1 {
                    &_grid.cells[(i - 1) as usize][(j - 2) as usize]
                } else {
                    &base_cell
                };

                // Use a very small search radius
                qr_alignment_pattern_search(
                    &mut p[k_idx as usize],
                    cell,
                    u,
                    v,
                    2,
                    _img,
                    _width,
                    _height,
                );

                if i > 0 && j > 0 {
                    qr_hom_cell_init(
                        &mut _grid.cells[(i - 1) as usize][(j - 1) as usize],
                        q[(k_idx - nalign - 1) as usize][0],
                        q[(k_idx - nalign - 1) as usize][1],
                        q[(k_idx - nalign) as usize][0],
                        q[(k_idx - nalign) as usize][1],
                        q[(k_idx - 1) as usize][0],
                        q[(k_idx - 1) as usize][1],
                        q[k_idx as usize][0],
                        q[k_idx as usize][1],
                        p[(k_idx - nalign - 1) as usize][0],
                        p[(k_idx - nalign - 1) as usize][1],
                        p[(k_idx - nalign) as usize][0],
                        p[(k_idx - nalign) as usize][1],
                        p[(k_idx - 1) as usize][0],
                        p[(k_idx - 1) as usize][1],
                        p[k_idx as usize][0],
                        p[k_idx as usize][1],
                    );
                }
            }
        }
    }

    // Set the limits over which each cell is used
    let last_cell_idx = _grid.cells.len() - 1;
    _grid.cell_limits[..last_cell_idx].copy_from_slice(&align_pos[1..][..last_cell_idx]);
    _grid.cell_limits[last_cell_idx] = dim;

    // Produce a bounding square for the code
    _p[0] = qr_hom_cell_project(&_grid.cells[0][0], -1, -1, 1);
    _p[1] = qr_hom_cell_project(&_grid.cells[0][last_cell_idx], (dim << 1) - 1, -1, 1);
    _p[2] = qr_hom_cell_project(&_grid.cells[last_cell_idx][0], -1, (dim << 1) - 1, 1);
    _p[3] = qr_hom_cell_project(
        &_grid.cells[last_cell_idx][last_cell_idx],
        (dim << 1) - 1,
        (dim << 1) - 1,
        1,
    );

    // Clamp the points somewhere near the image
    for point in _p {
        point[0] = c_int::max(
            -(_width << QR_FINDER_SUBPREC),
            c_int::min(point[0], (_width << QR_FINDER_SUBPREC) + 1),
        );
        point[1] = c_int::max(
            -(_height << QR_FINDER_SUBPREC),
            c_int::min(point[1], (_height << QR_FINDER_SUBPREC) + 1),
        );
    }
}

/// Sample QR code data bits from an image using the sampling grid
///
/// Reads bits from the image using the homography cells in the sampling grid,
/// XORing them with the data mask pattern. Bits are stored column-wise.
fn qr_sampling_grid_sample(
    _grid: &qr_sampling_grid,
    _dim: usize,
    _fmt_info: c_int,
    _img: &[u8],
    _width: c_int,
    _height: c_int,
) -> Vec<c_uint> {
    // We initialize the buffer with the data mask and XOR bits into it as we read
    // them out of the image instead of unmasking in a separate step
    let mut _data_bits = qr_make_data_mask(_dim, _fmt_info & 7);
    let stride = (_dim + QR_INT_BITS as usize - 1) >> QR_INT_LOGBITS as usize;
    let mut u0 = 0;

    // We read data cell-by-cell to avoid having to constantly change which
    // projection we're using as we read each bit.
    // Note that bits are stored column-wise, since that's how we'll scan them.
    for j in 0.._grid.cells.len() {
        let u1 = _grid.cell_limits[j];
        let mut v0 = 0;
        for i in 0.._grid.cells.len() {
            let v1 = _grid.cell_limits[i];
            let cell = &_grid.cells[i][j];
            let du = u0 - cell.u0;
            let dv = v0 - cell.v0;
            let mut x0 = cell.fwd[0][0] * du + cell.fwd[0][1] * dv + cell.fwd[0][2];
            let mut y0 = cell.fwd[1][0] * du + cell.fwd[1][1] * dv + cell.fwd[1][2];
            let mut w0 = cell.fwd[2][0] * du + cell.fwd[2][1] * dv + cell.fwd[2][2];

            for u in u0..u1 {
                let mut x = x0;
                let mut y = y0;
                let mut w = w0;
                for v in v0..v1 {
                    // Skip doing all the divisions and bounds checks if the bit is in the
                    // function pattern
                    if !qr_sampling_grid_is_in_fp(_grid, _dim, u, v) {
                        let p = qr_hom_cell_fproject(cell, x, y, w);
                        _data_bits[(u as usize) * stride + ((v >> QR_INT_LOGBITS) as usize)] ^=
                            (qr_img_get_bit(_img, _width, _height, p[0], p[1]) as c_uint)
                                << (v & (QR_INT_BITS - 1));
                    }
                    x += cell.fwd[0][1];
                    y += cell.fwd[1][1];
                    w += cell.fwd[2][1];
                }
                x0 += cell.fwd[0][0];
                y0 += cell.fwd[1][0];
                w0 += cell.fwd[2][0];
            }
            v0 = v1;
        }
        u0 = u1;
    }
    _data_bits
}

/// Arrange sample bits into bytes and Reed-Solomon blocks
///
/// Takes the bit data read from the QR code and groups it into bytes,
/// distributing those bytes across the Reed-Solomon blocks.
/// The block pointers are advanced by this routine.
///
/// # Safety
/// This function is unsafe because it dereferences and modifies raw pointers.
pub(crate) unsafe fn qr_samples_unpack(
    block_data: &mut [u8],
    mut block_positions: Vec<usize>,
    _nshort_data: usize,
    mut _nshort_blocks: usize,
    _data_bits: &[c_uint],
    _fp_mask: &[c_uint],
    _dim: usize,
) {
    let _nblocks = block_positions.len();
    let stride = (_dim + QR_INT_BITS as usize - 1) >> QR_INT_LOGBITS as usize;

    // If _all_ the blocks are short, don't skip anything
    if _nshort_blocks >= _nblocks {
        _nshort_blocks = 0;
    }

    let mut bits: c_uint = 0;
    let mut biti = 0;
    let mut blocki = 0;
    let mut blockj = 0;
    let mut j = _dim - 1;

    // Scan columns in pairs from right to left
    while j > 0 {
        // Scan up a pair of columns
        let mut nbits = ((_dim - 1) & (QR_INT_BITS as usize - 1)) + 1;
        let l = j * stride;

        for i in (0..stride).rev() {
            let data1 = _data_bits[l + i];
            let fp_mask1 = _fp_mask[l + i];
            let data2 = _data_bits[l + i - stride];
            let fp_mask2 = _fp_mask[l + i - stride];

            let mut nbits_inner = nbits;
            while nbits_inner > 0 {
                nbits_inner -= 1;
                // Pull a bit from the right column
                if ((fp_mask1 >> nbits_inner) & 1) == 0 {
                    bits = (bits << 1) | ((data1 >> nbits_inner) & 1);
                    biti += 1;
                }
                // Pull a bit from the left column
                if ((fp_mask2 >> nbits_inner) & 1) == 0 {
                    bits = (bits << 1) | ((data2 >> nbits_inner) & 1);
                    biti += 1;
                }
                // If we finished a byte, drop it in a block
                if biti >= 8 {
                    biti -= 8;
                    let pos = block_positions[blocki];
                    block_data[pos] = (bits >> biti) as c_uchar;
                    block_positions[blocki] += 1;
                    blocki += 1;

                    if blocki >= _nblocks {
                        blockj += 1;
                        blocki = if blockj == _nshort_data {
                            _nshort_blocks
                        } else {
                            0
                        };
                    }
                }
            }
            nbits = QR_INT_BITS as usize;
        }

        if j < 2 {
            break;
        }
        j -= 2;
        // Skip the column with the vertical timing pattern
        if j == 6 {
            j -= 1;
        }

        // Scan down a pair of columns
        let l = j * stride;
        for i in 0..stride {
            let mut data1 = _data_bits[l + i];
            let mut fp_mask1 = _fp_mask[l + i];
            let mut data2 = _data_bits[l + i - stride];
            let mut fp_mask2 = _fp_mask[l + i - stride];
            let mut nbits = usize::min(_dim - (i << QR_INT_LOGBITS as usize), QR_INT_BITS as usize);

            while nbits > 0 {
                nbits -= 1;
                // Pull a bit from the right column
                if (fp_mask1 & 1) == 0 {
                    bits = (bits << 1) | (data1 & 1);
                    biti += 1;
                }
                data1 >>= 1;
                fp_mask1 >>= 1;

                // Pull a bit from the left column
                if (fp_mask2 & 1) == 0 {
                    bits = (bits << 1) | (data2 & 1);
                    biti += 1;
                }
                data2 >>= 1;
                fp_mask2 >>= 1;

                // If we finished a byte, drop it in a block
                if biti >= 8 {
                    biti -= 8;
                    let pos = block_positions[blocki];
                    block_data[pos] = (bits >> biti) as c_uchar;
                    block_positions[blocki] += 1;
                    blocki += 1;

                    if blocki >= _nblocks {
                        blockj += 1;
                        blocki = if blockj == _nshort_data {
                            _nshort_blocks
                        } else {
                            0
                        };
                    }
                }
            }
        }

        if j < 2 {
            break;
        }
        j -= 2;
    }
}

/// The characters available in QR_MODE_ALNUM
const QR_ALNUM_TABLE: &[u8; 45] = b"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ $%*+-./:";

/// Parse QR code data from corrected codewords
///
/// Decodes the various data modes (numeric, alphanumeric, byte, Kanji, etc.)
/// and populates the qr_code_data structure.
///
/// # Safety
/// This function is unsafe because it dereferences raw pointers and performs memory allocation.
pub(crate) unsafe fn qr_code_data_parse(
    _qrdata: &mut qr_code_data,
    _version: c_int,
    data: &[u8],
) -> c_int {
    // The number of bits used to encode the character count for each version range and data mode
    const LEN_BITS: [[c_int; 4]; 3] = [
        [10, 9, 8, 8],    // Versions 1-9
        [12, 11, 16, 10], // Versions 10-26
        [14, 13, 16, 12], // Versions 27-40
    ];

    let mut self_parity: c_uint = 0;

    // Entries are stored directly in the struct during parsing
    _qrdata.entries = vec![];
    _qrdata.sa_size = 0;

    // The versions are divided into 3 ranges that each use a different number of bits for length fields
    let len_bits_idx = (if _version > 9 { 1 } else { 0 }) + (if _version > 26 { 1 } else { 0 });

    let mut qpb = qr_pack_buf {
        buf: data,
        endbyte: 0,
        endbit: 0,
    };

    // While we have enough bits to read a mode...
    while qr_pack_buf_avail(&qpb) >= 4 {
        let Some(mode) = qr_pack_buf_read(&mut qpb, 4) else {
            return -1;
        };

        // Mode 0 is a terminator
        if mode == 0 {
            break;
        }

        let Ok(mode) = qr_mode::try_from(mode) else {
            // Unknown mode - we can't skip it, so fail
            return -1;
        };

        let payload = match mode {
            qr_mode::Num => {
                let Some(len) = qr_pack_buf_read(&mut qpb, LEN_BITS[len_bits_idx][0]) else {
                    return -1;
                };

                let count = len / 3;
                let rem = len % 3;
                if qr_pack_buf_avail(&qpb) < 10 * count + 7 * ((rem >> 1) & 1) + 4 * (rem & 1) {
                    return -1;
                }

                let mut data = vec![0; len as usize];
                let mut buf = data.as_mut_slice();

                // Read groups of 3 digits encoded in 10 bits
                for _ in 0..count {
                    let Some(bits) = qr_pack_buf_read(&mut qpb, 10) else {
                        return -1;
                    };
                    let bits = bits as c_uint;
                    if bits >= 1000 {
                        return -1;
                    }
                    let c = b'0' + (bits / 100) as u8;
                    self_parity ^= c as c_uint;
                    buf[0] = c;
                    buf = &mut buf[1..];

                    let bits = bits % 100;
                    let c = b'0' + (bits / 10) as u8;
                    self_parity ^= c as c_uint;
                    buf[0] = c;
                    buf = &mut buf[1..];

                    let c = b'0' + (bits % 10) as u8;
                    self_parity ^= c as c_uint;
                    buf[0] = c;
                    buf = &mut buf[1..];
                }

                // Read the last two digits encoded in 7 bits
                if rem > 1 {
                    let Some(bits) = qr_pack_buf_read(&mut qpb, 7) else {
                        return -1;
                    };
                    let bits = bits as c_uint;
                    if bits >= 100 {
                        return -1;
                    }
                    let c = b'0' + (bits / 10) as u8;
                    self_parity ^= c as c_uint;
                    buf[0] = c;
                    buf = &mut buf[1..];

                    let c = b'0' + (bits % 10) as u8;
                    self_parity ^= c as c_uint;
                    buf[0] = c;
                }
                // Or the last one digit encoded in 4 bits
                else if rem != 0 {
                    let Some(bits) = qr_pack_buf_read(&mut qpb, 4) else {
                        return -1;
                    };
                    let bits = bits as c_uint;
                    if bits >= 10 {
                        return -1;
                    }
                    let c = b'0' + bits as u8;
                    self_parity ^= c as c_uint;
                    buf[0] = c;
                }

                qr_code_data_payload::Numeric(data)
            }
            qr_mode::Alnum => {
                let Some(len) = qr_pack_buf_read(&mut qpb, LEN_BITS[len_bits_idx][1]) else {
                    return -1;
                };

                let count = len >> 1;
                let rem = len & 1;
                if qr_pack_buf_avail(&qpb) < 11 * count + 6 * rem {
                    return -1;
                }

                let mut data = vec![0; len as usize];
                let mut buf = data.as_mut_slice();

                // Read groups of two characters encoded in 11 bits
                for _ in 0..count {
                    let Some(bits) = qr_pack_buf_read(&mut qpb, 11) else {
                        return -1;
                    };
                    let bits = bits as c_uint;
                    if bits >= 2025 {
                        return -1;
                    }
                    let c = QR_ALNUM_TABLE[(bits / 45) as usize];
                    self_parity ^= c as c_uint;
                    buf[0] = c;
                    buf = &mut buf[1..];

                    let c = QR_ALNUM_TABLE[(bits % 45) as usize];
                    self_parity ^= c as c_uint;
                    buf[0] = c;
                    buf = &mut buf[1..];
                }

                // Read the last character encoded in 6 bits
                if rem != 0 {
                    let Some(bits) = qr_pack_buf_read(&mut qpb, 6) else {
                        return -1;
                    };
                    let bits = bits as c_uint;
                    if bits >= 45 {
                        return -1;
                    }
                    let c = QR_ALNUM_TABLE[bits as usize];
                    self_parity ^= c as c_uint;
                    buf[0] = c;
                }

                qr_code_data_payload::Alphanumeric(data)
            }
            qr_mode::Struct => {
                // Structured-append header
                let Some(bits) = qr_pack_buf_read(&mut qpb, 16) else {
                    return -1;
                };

                let mut sa = qr_code_data_sa::default();

                // If this is the first S-A header, save it
                if _qrdata.sa_size == 0 {
                    _qrdata.sa_index = ((bits >> 12) & 0xF) as c_uchar;
                    sa.sa_index = _qrdata.sa_index;
                    _qrdata.sa_size = (((bits >> 8) & 0xF) + 1) as c_uchar;
                    sa.sa_size = _qrdata.sa_size;
                    _qrdata.sa_parity = (bits & 0xFF) as c_uchar;
                    sa.sa_parity = _qrdata.sa_parity;
                }

                qr_code_data_payload::StructuredAppendedHeaderData(sa)
            }
            qr_mode::Byte => {
                let Some(len) = qr_pack_buf_read(&mut qpb, LEN_BITS[len_bits_idx][2]) else {
                    return -1;
                };
                if len < 0 {
                    return -1;
                }

                if qr_pack_buf_avail(&qpb) < len << 3 {
                    return -1;
                }

                let mut data = vec![0; len as usize];

                for b in data.iter_mut() {
                    let Some(c) = qr_pack_buf_read(&mut qpb, 8) else {
                        return -1;
                    };
                    let c = c as c_uchar;
                    self_parity ^= c as c_uint;
                    *b = c;
                }

                qr_code_data_payload::Bytes(data)
            }
            qr_mode::Fnc1_1st => {
                // FNC1 first position marker
                // No data to read
                qr_code_data_payload::Fnc1FirstPositionMarker
            }
            qr_mode::Eci => {
                // Extended Channel Interpretation
                let Some(bits) = qr_pack_buf_read(&mut qpb, 8) else {
                    return -1;
                };

                let val = if (bits & 0x80) == 0 {
                    // One byte
                    bits as c_uint
                } else if (bits & 0x40) == 0 {
                    // Two bytes
                    let mut val = ((bits & 0x3F) as c_uint) << 8;
                    let Some(bits) = qr_pack_buf_read(&mut qpb, 8) else {
                        return -1;
                    };
                    val |= bits as c_uint;
                    val
                } else if (bits & 0x20) == 0 {
                    // Three bytes
                    let mut val = ((bits & 0x1F) as c_uint) << 16;
                    let Some(bits) = qr_pack_buf_read(&mut qpb, 16) else {
                        return -1;
                    };
                    val |= bits as c_uint;
                    if val >= 1000000 {
                        return -1;
                    }
                    val
                } else {
                    // Invalid lead byte
                    return -1;
                };

                qr_code_data_payload::ExtendedChannelInterpretation(val)
            }
            qr_mode::Kanji => {
                let Some(len) = qr_pack_buf_read(&mut qpb, LEN_BITS[len_bits_idx][3]) else {
                    return -1;
                };
                if len < 0 {
                    return -1;
                }

                if qr_pack_buf_avail(&qpb) < 13 * len {
                    return -1;
                }

                let mut data = vec![0; (2 * len) as usize];

                // Decode 2-byte SJIS characters encoded in 13 bits
                for i in 0..(len as usize) {
                    let Some(bits) = qr_pack_buf_read(&mut qpb, 13) else {
                        return -1;
                    };
                    let mut bits = bits as c_uint;
                    bits = (((bits / 0xC0) << 8) | (bits % 0xC0)) + 0x8140;
                    if bits >= 0xA000 {
                        bits += 0x4000;
                    }
                    self_parity ^= bits;
                    data[i * 2] = (bits >> 8) as c_uchar;
                    data[i * 2 + 1] = (bits & 0xFF) as c_uchar;
                }

                qr_code_data_payload::Kanji(data)
            }
            qr_mode::Fnc1_2nd => {
                // FNC1 second position marker
                let Some(bits) = qr_pack_buf_read(&mut qpb, 8) else {
                    return -1;
                };
                if !((0..100).contains(&bits)
                    || (165..191).contains(&bits)
                    || (197..223).contains(&bits))
                {
                    return -1;
                }
                qr_code_data_payload::ApplicationIndicator(bits)
            }
        };

        _qrdata.entries.push(qr_code_data_entry { payload });
    }

    // Store the parity of the data from this code
    _qrdata.self_parity = (((self_parity >> 8) ^ self_parity) & 0xFF) as c_uchar;

    0
}

/// Decode a QR code from an image
///
/// This is the main decoding function that:
/// 1. Initializes the sampling grid
/// 2. Samples data bits from the image
/// 3. Groups bits into Reed-Solomon codewords
/// 4. Performs error correction on each block
/// 5. Parses the corrected data
///
/// # Safety
/// This function is unsafe because it dereferences raw pointers and performs memory allocation.
#[allow(clippy::too_many_arguments)]
unsafe fn qr_code_decode(
    _qrdata: &mut qr_code_data,
    _ul_pos: &qr_point,
    _ur_pos: &qr_point,
    _dl_pos: &qr_point,
    _version: c_int,
    _fmt_info: c_int,
    _img: &[u8],
    _width: c_int,
    _height: c_int,
) -> c_int {
    let mut grid: qr_sampling_grid = qr_sampling_grid {
        cells: Vec::new(),
        fpmask: Vec::new(),
        cell_limits: [0; 6],
    };

    // Read the bits out of the image
    qr_sampling_grid_init(
        &mut grid,
        _version,
        _ul_pos,
        _ur_pos,
        _dl_pos,
        &mut _qrdata.bbox,
        _img,
        _width,
        _height,
    );

    let dim = 17 + (_version << 2) as usize;
    let data_bits = qr_sampling_grid_sample(&grid, dim, _fmt_info, _img, _width, _height);

    // Group those bits into Reed-Solomon codewords
    let ecc_level = (_fmt_info >> 3) ^ 1;
    let nblocks = QR_RS_NBLOCKS[(_version - 1) as usize][ecc_level as usize] as usize;
    let npar = QR_RS_NPAR_VALS
        [(QR_RS_NPAR_OFFS[(_version - 1) as usize] as usize) + (ecc_level as usize)]
        as usize;
    let ncodewords = qr_code_ncodewords(_version as c_uint);
    let block_sz = ncodewords / nblocks;
    let nshort_blocks = nblocks - (ncodewords % nblocks);

    let mut block_data = vec![0u8; ncodewords];
    let mut block_positions = vec![0usize; nblocks];

    // Initialize block starting positions
    block_positions[0] = 0;
    for i in 1..nblocks {
        block_positions[i] = block_positions[i - 1] + block_sz + (i > nshort_blocks) as usize;
    }

    qr_samples_unpack(
        &mut block_data,
        block_positions,
        block_sz - npar,
        nshort_blocks,
        &data_bits,
        &grid.fpmask,
        dim,
    );

    // Perform the error correction using reed-solomon crate
    let mut ndata = 0;
    let mut ncodewords_processed = 0;
    let mut ret = 0;

    for i in 0..nblocks {
        let block_szi = block_sz + (i >= nshort_blocks) as usize;

        // Use reed-solomon crate for error correction
        // QR codes always use m0=0, so we can use the crate directly
        let block_slice = &mut block_data[ncodewords_processed..][..block_szi];

        // Save original for error counting
        let original = block_slice.to_vec();

        let decoder = RSDecoder::new(npar);
        ret = match decoder.correct(&original, None) {
            Ok(corrected) => {
                // Copy corrected data back
                let corrected_data = corrected.to_vec();
                block_slice.copy_from_slice(&corrected_data);
                // Count errors by comparing original with corrected
                let mut error_count = 0;
                for j in 0..block_szi {
                    if original[j] != corrected_data[j] {
                        error_count += 1;
                    }
                }
                if error_count > 0 {
                    error_count
                } else {
                    0
                }
            }
            Err(_) => -1, // Correction failed
        };

        // For version 1 symbols and version 2-L and 3-L symbols, we aren't allowed
        // to use all the parity bytes for correction.
        // They are instead used to improve detection.
        if ret < 0
            || (_version == 1 && ret > ((ecc_level + 1) << 1))
            || (_version == 2 && ecc_level == 0 && ret > 4)
        {
            ret = -1;
            break;
        }

        let ndatai = block_szi - npar;
        block_data.copy_within(ncodewords_processed..ncodewords_processed + ndatai, ndata);
        ncodewords_processed += block_szi;
        ndata += ndatai;
    }

    // Parse the corrected bitstream
    if ret >= 0 {
        ret = qr_code_data_parse(_qrdata, _version, &block_data[..ndata]);
        if ret < 0 {
            qr_code_data_clear(_qrdata);
        }
        _qrdata.version = _version as c_uchar;
        _qrdata.ecc_level = ecc_level as c_uchar;
    }

    ret
}

enum Bch18_6CorrectError {
    Unrecoverable,
}

/// Correct a BCH(18,6,3) code word
///
/// Takes a code word and attempts to correct errors using the BCH(18,6,3) code.
///
/// Returns:
/// - `Ok((corrected_value, num_errors))` where `num_errors` is 0-3
/// - `Err(Bch18_6CorrectError::Unrecoverable)` if more than 3 errors were detected
fn bch18_6_correct(y: c_uint) -> Result<(c_uint, c_int), Bch18_6CorrectError> {
    // Check the easy case first: see if the data bits were uncorrupted
    let x = y >> 12;
    if (7..=40).contains(&x) {
        let nerrs = qr_hamming_dist(y, BCH18_6_CODES[(x - 7) as usize], 4);
        if nerrs < 4 {
            return Ok((BCH18_6_CODES[(x - 7) as usize], nerrs));
        }
    }

    // Exhaustive search is faster than field operations in GF(19)
    for (i, &code) in BCH18_6_CODES.iter().enumerate() {
        if i + 7 != (y >> 12) as usize {
            let nerrs = qr_hamming_dist(y, code, 4);
            if nerrs < 4 {
                return Ok((code, nerrs));
            }
        }
    }

    Err(Bch18_6CorrectError::Unrecoverable)
}

/// Initialize a homography cell for mapping between code and image space
///
/// This computes a homographic transformation from a quadrilateral in code space
/// (u0,v0)-(u1,v1)-(u2,v2)-(u3,v3) to a quadrilateral in image space
/// (x0,y0)-(x1,y1)-(x2,y2)-(x3,y3).
///
/// The transformation handles both affine and projective distortion, with careful
/// attention to numerical stability through dynamic range scaling.
#[allow(clippy::too_many_arguments)]
fn qr_hom_cell_init(
    cell: &mut qr_hom_cell,
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
    cell.fwd[0][0] = (if i00 != 0 { qr_divround(a00, i00) } else { 0 })
        + (if i10 != 0 { qr_divround(a01, i10) } else { 0 });
    cell.fwd[0][1] = (if i01 != 0 { qr_divround(a00, i01) } else { 0 })
        + (if i11 != 0 { qr_divround(a01, i11) } else { 0 });
    cell.fwd[1][0] = (if i00 != 0 { qr_divround(a10, i00) } else { 0 })
        + (if i10 != 0 { qr_divround(a11, i10) } else { 0 });
    cell.fwd[1][1] = (if i01 != 0 { qr_divround(a10, i01) } else { 0 })
        + (if i11 != 0 { qr_divround(a11, i11) } else { 0 });
    cell.fwd[2][0] = ((if i00 != 0 { qr_divround(a20, i00) } else { 0 })
        + (if i10 != 0 { qr_divround(a21, i10) } else { 0 })
        + (if i20 != 0 { qr_divround(a22, i20) } else { 0 })
        + round as c_int) as c_int
        >> shift;
    cell.fwd[2][1] = ((if i01 != 0 { qr_divround(a20, i01) } else { 0 })
        + (if i11 != 0 { qr_divround(a21, i11) } else { 0 })
        + (if i21 != 0 { qr_divround(a22, i21) } else { 0 })
        + round as c_int) as c_int
        >> shift;
    cell.fwd[2][2] = ((a22 + round as c_int) >> shift) as c_int;

    // Compute offsets to distribute rounding error over whole range
    // (instead of concentrating it in the (u3,v3) corner)
    let mut x = cell.fwd[0][0] * du10 + cell.fwd[0][1] * dv10;
    let mut y = cell.fwd[1][0] * du10 + cell.fwd[1][1] * dv10;
    let mut w = cell.fwd[2][0] * du10 + cell.fwd[2][1] * dv10 + cell.fwd[2][2];
    let mut a02 = dx10 * w - x;
    let mut a12 = dy10 * w - y;

    x = cell.fwd[0][0] * du20 + cell.fwd[0][1] * dv20;
    y = cell.fwd[1][0] * du20 + cell.fwd[1][1] * dv20;
    w = cell.fwd[2][0] * du20 + cell.fwd[2][1] * dv20 + cell.fwd[2][2];
    a02 += dx20 * w - x;
    a12 += dy20 * w - y;

    x = cell.fwd[0][0] * du30 + cell.fwd[0][1] * dv30;
    y = cell.fwd[1][0] * du30 + cell.fwd[1][1] * dv30;
    w = cell.fwd[2][0] * du30 + cell.fwd[2][1] * dv30 + cell.fwd[2][2];
    a02 += dx30 * w - x;
    a12 += dy30 * w - y;

    cell.fwd[0][2] = (a02 + 2) >> 2;
    cell.fwd[1][2] = (a12 + 2) >> 2;
    cell.x0 = x0;
    cell.y0 = y0;
    cell.u0 = u0;
    cell.v0 = v0;
}

/// Get a bit from a binarized image
///
/// Samples a pixel from the binarized image, with coordinates in QR_FINDER_SUBPREC
/// subpixel units. Clamps coordinates to valid image bounds.
pub(crate) fn qr_img_get_bit(
    img: &[u8],
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
    if img[idx as usize] != 0 {
        1
    } else {
        0
    }
}

/// Finish a partial projection, converting from homogeneous coordinates to the
/// normal 2-D representation.
/// In loops, we can avoid many multiplies by computing the homogeneous _x, _y,
/// and _w incrementally, but we cannot avoid the divisions, done here.
fn qr_hom_cell_fproject(
    _cell: &qr_hom_cell,
    mut _x: c_int,
    mut _y: c_int,
    mut _w: c_int,
) -> qr_point {
    if _w == 0 {
        [
            if _x < 0 { c_int::MIN } else { c_int::MAX },
            if _y < 0 { c_int::MIN } else { c_int::MAX },
        ]
    } else {
        if _w < 0 {
            _x = -_x;
            _y = -_y;
            _w = -_w;
        }
        [
            qr_divround(_x, _w) + _cell.x0,
            qr_divround(_y, _w) + _cell.y0,
        ]
    }
}

fn qr_hom_cell_project(_cell: &qr_hom_cell, mut _u: c_int, mut _v: c_int, _res: c_int) -> qr_point {
    _u -= _cell.u0 << _res;
    _v -= _cell.v0 << _res;
    qr_hom_cell_fproject(
        _cell,
        _cell.fwd[0][0] * _u + _cell.fwd[0][1] * _v + (_cell.fwd[0][2] << _res),
        _cell.fwd[1][0] * _u + _cell.fwd[1][1] * _v + (_cell.fwd[1][2] << _res),
        _cell.fwd[2][0] * _u + _cell.fwd[2][1] * _v + (_cell.fwd[2][2] << _res),
    )
}

/// Locate the crossing of a finder or alignment pattern along a line
///
/// Uses Bresenham's algorithm to trace along the line and find the exact
/// transitions from !_v to _v and back. Returns the midpoint of the segment.
///
/// Returns 0 on success, -1 if no crossing found.
#[allow(clippy::too_many_arguments)]
pub(crate) unsafe fn qr_finder_locate_crossing(
    img: &[u8],
    width: c_int,
    _height: c_int,
    x0: c_int,
    y0: c_int,
    x1: c_int,
    y1: c_int,
    v: c_int,
    p: &mut qr_point,
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
        let pixel = img[(x0_pos[1] * width + x0_pos[0]) as usize];
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
        let pixel = img[(x1_pos[1] * width + x1_pos[0]) as usize];
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
fn qr_alignment_pattern_fetch(
    p: &[[qr_point; 5]; 5],
    x0: c_int,
    y0: c_int,
    img: &[u8],
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
#[allow(clippy::too_many_arguments)]
unsafe fn qr_alignment_pattern_search(
    p: *mut qr_point,
    cell: &qr_hom_cell,
    _u: c_int,
    _v: c_int,
    r: c_int,
    img: &[u8],
    width: c_int,
    height: c_int,
) -> c_int {
    let mut pattern: [[qr_point; 5]; 5] = [[Default::default(); 5]; 5];
    let mut nc = [0 as c_int; 4];
    let mut c: [qr_point; 4] = [[0; 2]; 4];
    let mut pc = qr_point::default();

    // Build up a basic template using cell to control shape and scale
    let u = (_u - 2) - cell.u0;
    let v = (_v - 2) - cell.v0;
    let mut x0 = cell.fwd[0][0] * u + cell.fwd[0][1] * v + cell.fwd[0][2];
    let mut y0 = cell.fwd[1][0] * u + cell.fwd[1][1] * v + cell.fwd[1][2];
    let mut w0 = cell.fwd[2][0] * u + cell.fwd[2][1] * v + cell.fwd[2][2];
    let dxdu = cell.fwd[0][0];
    let dydu = cell.fwd[1][0];
    let dwdu = cell.fwd[2][0];
    let dxdv = cell.fwd[0][1];
    let dydv = cell.fwd[1][1];
    let dwdv = cell.fwd[2][1];

    for item in pattern.iter_mut() {
        let mut x = x0;
        let mut y = y0;
        let mut w = w0;
        for subitem in item.iter_mut() {
            *subitem = qr_hom_cell_fproject(cell, x, y, w);
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
        let u = _u - cell.u0;
        let v = _v - cell.v0;
        let mut x = (cell.fwd[0][0] * u + cell.fwd[0][1] * v + cell.fwd[0][2]) << QR_ALIGN_SUBPREC;
        let mut y = (cell.fwd[1][0] * u + cell.fwd[1][1] * v + cell.fwd[1][2]) << QR_ALIGN_SUBPREC;
        let mut w = (cell.fwd[2][0] * u + cell.fwd[2][1] * v + cell.fwd[2][2]) << QR_ALIGN_SUBPREC;

        // Search an area at most r modules around the target location, in concentric squares
        for i in 1..(r << QR_ALIGN_SUBPREC) {
            let side_len = (i << 1) - 1;
            x -= dxdu + dxdv;
            y -= dydu + dydv;
            w -= dwdu + dwdv;

            for j in 0..(4 * side_len) {
                pc = qr_hom_cell_fproject(cell, x, y, w);
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
                    x += cell.fwd[0][dir];
                    y += cell.fwd[1][dir];
                    w += cell.fwd[2][dir];
                } else {
                    x -= cell.fwd[0][dir];
                    y -= cell.fwd[1][dir];
                    w -= cell.fwd[2][dir];
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
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum qr_mode {
    /// Numeric digits ('0'...'9')
    Num = 1,
    /// Alphanumeric characters ('0'...'9', 'A'...'Z', plus punctuation ' ', '$', '%', '*', '+', '-', '.', '/', ':')
    Alnum = 2,
    /// Structured-append header
    Struct = 3,
    /// Raw 8-bit bytes
    Byte = 4,
    /// FNC1 marker in first position (GS1 formatting)
    Fnc1_1st = 5,
    /// Extended Channel Interpretation code
    Eci = 7,
    /// SJIS kanji characters
    Kanji = 8,
    /// FNC1 marker in second position (industry application)
    Fnc1_2nd = 9,
}

impl TryFrom<c_int> for qr_mode {
    type Error = ();

    fn try_from(value: c_int) -> Result<Self, Self::Error> {
        Ok(match value {
            1 => Self::Num,
            2 => Self::Alnum,
            3 => Self::Struct,
            4 => Self::Byte,
            5 => Self::Fnc1_1st,
            7 => Self::Eci,
            8 => Self::Kanji,
            9 => Self::Fnc1_2nd,
            _ => return Err(()),
        })
    }
}

/// Data payload for a QR code data entry
#[derive(Clone)]
pub(crate) enum qr_code_data_payload {
    /// Numeric digits ('0'...'9')
    Numeric(Vec<u8>),
    /// Alphanumeric characters ('0'...'9', 'A'...'Z', plus punctuation ' ', '$', '%', '*', '+', '-', '.', '/', ':')
    Alphanumeric(Vec<u8>),
    /// Structured-append header
    #[allow(dead_code)]
    StructuredAppendedHeaderData(qr_code_data_sa),
    /// Raw 8-bit bytes
    Bytes(Vec<u8>),
    /// FNC1 marker in first position (GS1 formatting)
    Fnc1FirstPositionMarker,
    /// Extended Channel Interpretation code
    ExtendedChannelInterpretation(c_uint),
    /// SJIS kanji characters
    Kanji(Vec<u8>),
    /// FNC1 marker in second position (industry application)
    ApplicationIndicator(c_int),
}

/// Structured-append data
#[derive(Debug, Copy, Clone, Default)]
pub(crate) struct qr_code_data_sa {
    pub(crate) sa_index: c_uchar,
    pub(crate) sa_size: c_uchar,
    pub(crate) sa_parity: c_uchar,
}

/// A single QR code data entry
pub(crate) struct qr_code_data_entry {
    /// The payload
    pub(crate) payload: qr_code_data_payload,
}

/// Low-level QR code data
#[derive(Default)]
pub(crate) struct qr_code_data {
    /// The decoded data entries
    pub(crate) entries: Vec<qr_code_data_entry>,
    /// The code version (1...40)
    pub(crate) version: c_uchar,
    /// The ECC level (0...3, corresponding to 'L', 'M', 'Q', and 'H')
    pub(crate) ecc_level: c_uchar,
    /// The index of this code in the structured-append group
    pub(crate) sa_index: c_uchar,
    /// The size of the structured-append group, or 0 if there was no S-A header
    pub(crate) sa_size: c_uchar,
    /// The parity of the entire structured-append group
    pub(crate) sa_parity: c_uchar,
    /// The parity of this code
    pub(crate) self_parity: c_uchar,
    /// An approximate bounding box for the code
    pub(crate) bbox: [qr_point; 4],
}

/// List of QR code data
#[derive(Default)]
pub(crate) struct qr_code_data_list {
    pub(crate) qrdata: Vec<qr_code_data>,
}

/// Clear a QR code data structure.
pub(crate) unsafe fn qr_code_data_clear(qrdata: &mut qr_code_data) {
    for entry in qrdata.entries.iter_mut() {
        match &mut entry.payload {
            qr_code_data_payload::Numeric(items)
            | qr_code_data_payload::Alphanumeric(items)
            | qr_code_data_payload::Bytes(items)
            | qr_code_data_payload::Kanji(items) => items.clear(),
            qr_code_data_payload::StructuredAppendedHeaderData(_)
            | qr_code_data_payload::Fnc1FirstPositionMarker
            | qr_code_data_payload::ExtendedChannelInterpretation(_)
            | qr_code_data_payload::ApplicationIndicator(_) => {}
        }
    }
}

/// Try to decode a QR code with the given configuration of three finder patterns.
///
/// This function tries different orderings of the three finder pattern centers
/// to find a valid QR code configuration, estimates module size and version,
/// fits a homography transformation, and attempts to decode the QR code.
///
/// Returns the version number if successful, -1 otherwise.
pub(crate) unsafe fn qr_reader_try_configuration(
    reader: &mut qr_reader,
    qrdata: &mut qr_code_data,
    img: &[u8],
    _width: c_int,
    _height: c_int,
    centers: [*mut qr_finder_center; 3],
) -> c_int {
    let mut ci: [usize; 7] = [0; 7];
    let mut maxd: c_uint;

    let mut i0: usize;

    // Sort the points in counter-clockwise order
    let ccw: c_int = qr_point_ccw(&(*centers[0]).pos, &(*centers[1]).pos, &(*centers[2]).pos);

    // Colinear points can't be the corners of a quadrilateral
    if ccw == 0 {
        return -1;
    }

    // Include a few extra copies of the cyclical list to avoid mods
    ci[6] = 0;
    ci[3] = 0;
    ci[0] = 0;
    ci[4] = if ccw < 0 { 2 } else { 1 };
    ci[1] = ci[4];
    ci[5] = if ccw < 0 { 1 } else { 2 };
    ci[2] = ci[5];

    // Assume the points farthest from each other are the opposite corners,
    // and find the top-left point
    maxd = qr_point_distance2(&(*centers[1]).pos, &(*centers[2]).pos);
    i0 = 0;
    for i in 1..3 {
        let d = qr_point_distance2(&(*centers[ci[i + 1]]).pos, &(*centers[ci[i + 2]]).pos);
        if d > maxd {
            i0 = i;
            maxd = d;
        }
    }

    // However, try all three possible orderings, just to be sure (a severely
    // skewed projection could move opposite corners closer than adjacent)
    for i in i0..i0 + 3 {
        let mut aff = qr_aff::default();
        let mut hom = qr_hom::default();
        let mut ul = qr_finder::default();
        let mut ur = qr_finder::default();
        let mut dl = qr_finder::default();
        let mut bbox: [qr_point; 4] = Default::default();

        let ur_version: c_int;
        let mut fmt_info: c_int;

        ul.c = centers[ci[i]];
        ur.c = centers[ci[i + 1]];
        dl.c = centers[ci[i + 2]];

        // Estimate the module size and version number from the two opposite corners.
        // The module size is not constant in the image, so we compute an affine
        // projection from the three points we have to a square domain, and
        // estimate it there.
        // Although it should be the same along both axes, we keep separate
        // estimates to account for any remaining projective distortion.
        let res: c_int = QR_INT_BITS
            - 2
            - QR_FINDER_SUBPREC
            - qr_ilog((c_int::max(_width, _height) - 1) as c_uint);
        qr_aff_init(&mut aff, &(*ul.c).pos, &(*ur.c).pos, &(*dl.c).pos, res);
        ur.o = qr_aff_unproject(&aff, (*ur.c).pos[0], (*ur.c).pos[1]);
        qr_finder_edge_pts_aff_classify(&mut ur, &aff);
        if qr_finder_estimate_module_size_and_version(&mut ur, 1 << res, 1 << res) < 0 {
            continue;
        }
        dl.o = qr_aff_unproject(&aff, (*dl.c).pos[0], (*dl.c).pos[1]);
        qr_finder_edge_pts_aff_classify(&mut dl, &aff);
        if qr_finder_estimate_module_size_and_version(&mut dl, 1 << res, 1 << res) < 0 {
            continue;
        }

        // If the estimated versions are significantly different, reject the
        // configuration
        if (ur.eversion[1] - dl.eversion[0]).abs() > QR_LARGE_VERSION_SLACK {
            continue;
        }

        ul.o = qr_aff_unproject(&aff, (*ul.c).pos[0], (*ul.c).pos[1]);
        qr_finder_edge_pts_aff_classify(&mut ul, &aff);
        if qr_finder_estimate_module_size_and_version(&mut ul, 1 << res, 1 << res) < 0
            || (ul.eversion[1] - ur.eversion[1]).abs() > QR_LARGE_VERSION_SLACK
            || (ul.eversion[0] - dl.eversion[0]).abs() > QR_LARGE_VERSION_SLACK
        {
            continue;
        }

        // If we made it this far, upgrade the affine homography to a full
        // homography
        if qr_hom_fit(
            &mut hom,
            &mut ul,
            &mut ur,
            &mut dl,
            &mut bbox,
            &aff,
            &mut reader.rng,
            img,
            _width,
            _height,
        ) < 0
        {
            continue;
        }

        qrdata.bbox = bbox;

        ul.o = qr_hom_unproject(&hom, (*ul.c).pos[0], (*ul.c).pos[1]).unwrap_or_default();
        ur.o = qr_hom_unproject(&hom, (*ur.c).pos[0], (*ur.c).pos[1]).unwrap_or_default();
        dl.o = qr_hom_unproject(&hom, (*dl.c).pos[0], (*dl.c).pos[1]).unwrap_or_default();
        qr_finder_edge_pts_hom_classify(&mut ur, &hom);
        let width = ur.o[0] - ul.o[0];
        let height = ur.o[0] - ul.o[0];
        if qr_finder_estimate_module_size_and_version(&mut ur, width, height) < 0 {
            continue;
        }
        qr_finder_edge_pts_hom_classify(&mut dl, &hom);
        let width = dl.o[1] - ul.o[1];
        let height = dl.o[1] - ul.o[1];
        if qr_finder_estimate_module_size_and_version(&mut dl, width, height) < 0 {
            continue;
        }

        // If we have a small version (less than 7), there's no encoded version
        // information. If the estimated version on the two corners matches and is
        // sufficiently small, we assume this is the case.
        if ur.eversion[1] == dl.eversion[0] && ur.eversion[1] < 7 {
            ur_version = ur.eversion[1];
        } else {
            // If the estimated versions are significantly different, reject the
            // configuration
            if (ur.eversion[1] - dl.eversion[0]).abs() > QR_LARGE_VERSION_SLACK {
                continue;
            }

            // Otherwise we try to read the actual version data from the image.
            // If the real version is not sufficiently close to our estimated version,
            // then we assume there was an unrecoverable decoding error (so many bit
            // errors we were within 3 errors of another valid code), and throw that
            // value away.
            // If no decoded version could be sufficiently close, we don't even try.
            let ur_version_tmp = if ur.eversion[1] >= 7 - QR_LARGE_VERSION_SLACK {
                let ver = qr_finder_version_decode(&ur, &hom, img, _width, _height, 0);
                if (ver - ur.eversion[1]).abs() > QR_LARGE_VERSION_SLACK {
                    -1
                } else {
                    ver
                }
            } else {
                -1
            };

            let dl_version_tmp = if dl.eversion[0] >= 7 - QR_LARGE_VERSION_SLACK {
                let ver = qr_finder_version_decode(&dl, &hom, img, _width, _height, 1);
                if (ver - dl.eversion[0]).abs() > QR_LARGE_VERSION_SLACK {
                    -1
                } else {
                    ver
                }
            } else {
                -1
            };

            // If we got at least one valid version, or we got two and they match,
            // then we found a valid configuration
            if ur_version_tmp >= 0 {
                if dl_version_tmp >= 0 && dl_version_tmp != ur_version_tmp {
                    continue;
                }
                ur_version = ur_version_tmp;
            } else if dl_version_tmp < 0 {
                continue;
            } else {
                ur_version = dl_version_tmp;
            }
        }

        qr_finder_edge_pts_hom_classify(&mut ul, &hom);
        let width = ur.o[0] - dl.o[0];
        let height = dl.o[1] - ul.o[1];
        if qr_finder_estimate_module_size_and_version(&mut ul, width, height) < 0
            || (ul.eversion[1] - ur.eversion[1]).abs() > QR_SMALL_VERSION_SLACK
            || (ul.eversion[0] - dl.eversion[0]).abs() > QR_SMALL_VERSION_SLACK
        {
            continue;
        }

        fmt_info = qr_finder_fmt_info_decode(&ul, &ur, &dl, &hom, img, _width, _height);
        if fmt_info < 0
            || qr_code_decode(
                qrdata,
                &(*ul.c).pos,
                &(*ur.c).pos,
                &(*dl.c).pos,
                ur_version,
                fmt_info,
                img,
                _width,
                _height,
            ) < 0
        {
            // The code may be flipped.
            // Try again, swapping the UR and DL centers.
            // We should get a valid version either way, so it's relatively cheap to
            // check this, as we've already filtered out a lot of invalid
            // configurations.
            // Swap using temporary variables to avoid borrow checker issues
            let t = hom.inv[0][0];
            hom.inv[0][0] = hom.inv[1][0];
            hom.inv[1][0] = t;
            let t = hom.inv[0][1];
            hom.inv[0][1] = hom.inv[1][1];
            hom.inv[1][1] = t;
            hom.fwd[0].swap(0, 1);
            hom.fwd[1].swap(0, 1);
            hom.fwd[2].swap(0, 1);
            ul.o.swap(0, 1);
            ul.size.swap(0, 1);
            ur.o.swap(0, 1);
            ur.size.swap(0, 1);
            dl.o.swap(0, 1);
            dl.size.swap(0, 1);

            fmt_info = qr_finder_fmt_info_decode(&ul, &dl, &ur, &hom, img, _width, _height);
            if fmt_info < 0 {
                continue;
            }

            let t = bbox[1][0];
            bbox[1][0] = bbox[2][0];
            bbox[2][0] = t;
            let t = bbox[1][1];
            bbox[1][1] = bbox[2][1];
            bbox[2][1] = t;
            qrdata.bbox = bbox;

            if qr_code_decode(
                qrdata,
                &(*ul.c).pos,
                &(*dl.c).pos,
                &(*ur.c).pos,
                ur_version,
                fmt_info,
                img,
                _width,
                _height,
            ) < 0
            {
                continue;
            }
        }

        return ur_version;
    }

    -1
}

/// Add a found finder line to the reader's line list
pub(crate) fn _zbar_qr_found_line(
    reader: &mut qr_reader,
    dir: c_int,
    line: &qr_finder_line,
) -> c_int {
    let lines = &mut reader.finder_lines[dir as usize];
    lines.lines.push(*line);
    0
}

/// Match finder centers and decode QR codes
unsafe fn qr_reader_match_centers(
    reader: &mut qr_reader,
    _qrlist: &mut qr_code_data_list,
    _centers: &mut [qr_finder_center],
    _ncenters: usize,
    img: &[u8],
    _width: c_int,
    _height: c_int,
) {
    // The number of centers should be small, so an O(n^3) exhaustive search of
    // which ones go together should be reasonable.
    let mut mark = vec![0u8; _ncenters];
    let nfailures_max = c_int::max(8192, (_width * _height) >> 9);
    let mut nfailures = 0;

    for i in 0.._ncenters {
        // TODO: We might be able to accelerate this step significantly by
        // considering the remaining finder centers in a more intelligent order,
        // based on the first finder center we just chose.
        for j in i + 1.._ncenters {
            if mark[i] != 0 {
                break;
            }

            for k in j + 1.._ncenters {
                if mark[j] != 0 {
                    break;
                }

                if mark[k] == 0 {
                    let mut qrdata = qr_code_data::default();
                    let version = qr_reader_try_configuration(
                        reader,
                        &mut qrdata,
                        img,
                        _width,
                        _height,
                        [
                            &mut _centers[i] as *mut qr_finder_center,
                            &mut _centers[j] as *mut qr_finder_center,
                            &mut _centers[k] as *mut qr_finder_center,
                        ],
                    );

                    if version >= 0 {
                        let mut ninside: usize;

                        // Convert the bounding box we're returning to the user to normal
                        // image coordinates
                        for point in qrdata.bbox.as_mut_slice() {
                            point[0] >>= QR_FINDER_SUBPREC;
                            point[1] >>= QR_FINDER_SUBPREC;
                        }

                        // Mark these centers as used
                        mark[i] = 1;
                        mark[j] = 1;
                        mark[k] = 1;

                        // Find any other finder centers located inside this code
                        ninside = 0;
                        for l in 0.._ncenters {
                            if mark[l] == 0
                                && qr_point_ccw(&qrdata.bbox[0], &qrdata.bbox[1], &_centers[l].pos)
                                    >= 0
                                && qr_point_ccw(&qrdata.bbox[1], &qrdata.bbox[3], &_centers[l].pos)
                                    >= 0
                                && qr_point_ccw(&qrdata.bbox[3], &qrdata.bbox[2], &_centers[l].pos)
                                    >= 0
                                && qr_point_ccw(&qrdata.bbox[2], &qrdata.bbox[0], &_centers[l].pos)
                                    >= 0
                            {
                                mark[l] = 2;
                                ninside += 1;
                            }
                        }

                        if ninside >= 3 {
                            // We might have a "Double QR": a code inside a code.
                            // Copy the relevant centers to a new array and do a search confined
                            // to that subset.
                            let mut inside = vec![];
                            for l in 0.._ncenters {
                                if mark[l] == 2 {
                                    inside.push(_centers[l].clone());
                                }
                            }
                            ninside = inside.len();
                            qr_reader_match_centers(
                                reader,
                                _qrlist,
                                &mut inside,
                                ninside,
                                img,
                                _width,
                                _height,
                            );
                        }

                        // Mark _all_ such centers used: codes cannot partially overlap
                        for m in mark.iter_mut().take(_ncenters) {
                            if *m == 2 {
                                *m = 1;
                            }
                        }

                        nfailures = 0;

                        // Add the data to the list
                        _qrlist.qrdata.push(qrdata);
                    } else {
                        nfailures += 1;
                        if nfailures > nfailures_max {
                            // Give up.
                            // We're unlikely to find a valid code in all this clutter, and we
                            // could spend quite a lot of time trying.
                            // mark will be automatically freed by Drop
                            return;
                        }
                    }
                }
            }
        }
    }
}

/// Decode QR codes from an image
pub(crate) unsafe fn qr_decode(
    reader: &mut qr_reader,
    iscn: &mut zbar_image_scanner_t,
    img: &mut zbar_image_t,
) -> c_int {
    if reader.finder_lines[0].lines.len() < 9 || reader.finder_lines[1].lines.len() < 9 {
        return 0;
    }

    let mut centers = qr_finder_centers_locate(reader, 0, 0);

    if centers.len() >= 3 {
        let bin = binarize(&img.data, img.width as usize, img.height as usize);

        let mut qrlist = qr_code_data_list::default();

        let ncenters = centers.len();
        qr_reader_match_centers(
            reader,
            &mut qrlist,
            &mut centers,
            ncenters,
            &bin,
            img.width as c_int,
            img.height as c_int,
        );

        let nqrdata = if !qrlist.qrdata.is_empty() {
            qr_code_data_list_extract_text(&qrlist, iscn)
        } else {
            0
        };

        qrlist.qrdata.clear();
        nqrdata
    } else {
        0
    }
}
