/*Copyright (C) 2008-2009  Timothy B. Terriberry (tterribe@xiph.org)
  You can redistribute this library and/or modify it under the terms of the
   GNU Lesser General Public License as published by the Free Software
   Foundation; either version 2.1 of the License, or (at your option) any later
   version.*/

#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "bch15_5.h"
#include "binarize.h"
#include "image.h"
#include "isaac.h"
#include "qrcode.h"
#include "rs.h"
#include "util.h"

#include "qrdec.h"

typedef int qr_line[3];

typedef struct qr_finder_cluster qr_finder_cluster;
typedef struct qr_finder_edge_pt qr_finder_edge_pt;
typedef struct qr_finder_center qr_finder_center;

typedef struct qr_aff qr_aff;
typedef struct qr_hom qr_hom;

typedef struct qr_finder qr_finder;

typedef struct qr_hom_cell qr_hom_cell;
typedef struct qr_sampling_grid qr_sampling_grid;
typedef struct qr_pack_buf qr_pack_buf;

/*The number of bits in an int.
  Note the cast to (int): this prevents this value from "promoting" whole
   expressions to an (unsigned) size_t.*/
#define QR_INT_BITS    ((int)sizeof(int) * CHAR_BIT)
#define QR_INT_LOGBITS (QR_ILOG(QR_INT_BITS))

/*A 14 bit resolution for a homography ensures that the ideal module size for a
   version 40 code differs from that of a version 39 code by at least 2.*/
#define QR_HOM_BITS (14)

/*The number of bits of sub-module precision to use when searching for
   alignment patterns.
  Two bits allows an alignment pattern to be found even if the modules have
   been eroded by up to 50% (due to blurring, etc.).
  This must be at least one, since it affects the dynamic range of the
   transforms, and we sample at half-module resolution to compute a bounding
   quadrilateral for the code.*/
#define QR_ALIGN_SUBPREC (2)

/* collection of finder lines */
typedef struct qr_finder_lines {
    qr_finder_line *lines;
    int nlines, clines;
} qr_finder_lines;

struct qr_reader {
    /*The GF(256) representation used in Reed-Solomon decoding.*/
    rs_gf256 gf;
    /*The random number generator used by RANSAC.*/
    isaac_ctx isaac;
    /* current finder state, horizontal and vertical lines */
    qr_finder_lines finder_lines[2];
};

/*Initializes a client reader handle.*/
extern void qr_reader_init(qr_reader *reader);

/*Allocates a client reader handle.*/
extern qr_reader *_zbar_qr_create(void);

/*Frees a client reader handle.*/
extern void _zbar_qr_destroy(qr_reader *reader);

/* reset finder state between scans */
extern void _zbar_qr_reset(qr_reader *reader);

/*A cluster of lines crossing a finder pattern (all in the same direction).*/
struct qr_finder_cluster {
    /*Pointers to the lines crossing the pattern.*/
    qr_finder_line **lines;
    /*The number of lines in the cluster.*/
    int nlines;
};

/*A point on the edge of a finder pattern.
  These are obtained from the endpoints of the lines crossing this particular
   pattern.*/
struct qr_finder_edge_pt {
    /*The location of the edge point.*/
    qr_point pos;
    /*A label classifying which edge this belongs to:
  0: negative u edge (left)
  1: positive u edge (right)
  2: negative v edge (top)
  3: positive v edge (bottom)*/
    int edge;
    /*The (signed) perpendicular distance of the edge point from a line parallel
   to the edge passing through the finder center, in (u,v) coordinates.
  This is also re-used by RANSAC to store inlier flags.*/
    int extent;
};

/*The center of a finder pattern obtained from the crossing of one or more
   clusters of horizontal finder lines with one or more clusters of vertical
   finder lines.*/
struct qr_finder_center {
    /*The estimated location of the finder center.*/
    qr_point pos;
    /*The list of edge points from the crossing lines.*/
    qr_finder_edge_pt *edge_pts;
    /*The number of edge points from the crossing lines.*/
    int nedge_pts;
};

/* Implemented in Rust (src/qrcode/qrdec.rs) */
extern int qr_finder_vline_cmp(const void *_a, const void *_b);

/*Clusters adjacent lines into groups that are large enough to be crossing a
   finder pattern (relative to their length).
  _clusters:  The buffer in which to store the clusters found.
  _neighbors: The buffer used to store the lists of lines in each cluster.
  _lines:     The list of lines to cluster.
              Horizontal lines must be sorted in ascending order by Y
               coordinate, with ties broken by X coordinate.
              Vertical lines must be sorted in ascending order by X coordinate,
               with ties broken by Y coordinate.
  _nlines:    The number of lines in the set of lines to cluster.
  _v:         0 for horizontal lines, or 1 for vertical lines.
  Return: The number of clusters.
  Implemented in Rust (src/qrcode/qrdec.rs) */
extern int qr_finder_cluster_lines(qr_finder_cluster *_clusters,
				   qr_finder_line **_neighbors,
				   qr_finder_line *_lines, int _nlines, int _v);

/*Adds the coordinates of the edge points from the lines contained in the
   given list of clusters to the list of edge points for a finder center.
  Only the edge point position is initialized.
  The edge label and extent are set by qr_finder_edge_pts_aff_classify()
   or qr_finder_edge_pts_hom_classify().
  _edge_pts:   The buffer in which to store the edge points.
  _nedge_pts:  The current number of edge points in the buffer.
  _neighbors:  The list of lines in the cluster.
  _nneighbors: The number of lines in the list of lines in the cluster.
  _v:          0 for horizontal lines and 1 for vertical lines.
  Return: The new total number of edge points.
  Implemented in Rust (src/qrcode/qrdec.rs) */
extern int qr_finder_edge_pts_fill(qr_finder_edge_pt *_edge_pts, int _nedge_pts,
				   qr_finder_cluster **_neighbors,
				   int _nneighbors, int _v);

/* Implemented in Rust (src/qrcode/qrdec.rs) */
extern int qr_finder_center_cmp(const void *_a, const void *_b);

/*Determine if a horizontal line crosses a vertical line.
  _hline: The horizontal line.
  _vline: The vertical line.
  Return: A non-zero value if the lines cross, or zero if they do not.*/
extern int qr_finder_lines_are_crossing(const qr_finder_line *_hline,
					const qr_finder_line *_vline);

/*Finds horizontal clusters that cross corresponding vertical clusters,
   presumably corresponding to a finder center.
  _center:     The buffer in which to store putative finder centers.
  _edge_pts:   The buffer to use for the edge point lists for each finder
                center.
  _hclusters:  The clusters of horizontal lines crossing finder patterns.
  _nhclusters: The number of horizontal line clusters.
  _vclusters:  The clusters of vertical lines crossing finder patterns.
  _nvclusters: The number of vertical line clusters.
  Return: The number of putative finder centers.
  Implemented in Rust (src/qrcode/qrdec.rs) */
extern int qr_finder_find_crossings(qr_finder_center *_centers,
				    qr_finder_edge_pt *_edge_pts,
				    qr_finder_cluster *_hclusters,
				    int _nhclusters,
				    qr_finder_cluster *_vclusters,
				    int _nvclusters);

/*Locates a set of putative finder centers in the image.
  First we search for horizontal and vertical lines that have
   (dark:light:dark:light:dark) runs with size ratios of roughly (1:1:3:1:1).
  Then we cluster them into groups such that each subsequent pair of endpoints
   is close to the line before it in the cluster.
  This will locate many line clusters that don't cross a finder pattern, but
   qr_finder_find_crossings() will filter most of them out.
  Where horizontal and vertical clusters cross, a prospective finder center is
   returned.
  _centers:  Returns a pointer to a freshly-allocated list of finder centers.
             This must be freed by the caller.
  _edge_pts: Returns a pointer to a freshly-allocated list of edge points
              around those centers.
             This must be freed by the caller.
  _img:      The binary image to search.
  _width:    The width of the image.
  _height:   The height of the image.
  Return: The number of putative finder centers located.*/
/*Locate QR finder pattern centers from scanned lines.
  Clusters horizontal and vertical lines that cross finder patterns,
  then locates the centers where horizontal and vertical clusters intersect.
  Implemented in Rust (src/qrcode/qrdec.rs) */
extern int qr_finder_centers_locate(qr_finder_center **_centers,
				     qr_finder_edge_pt **_edge_pts,
				     qr_reader *reader, int _width, int _height);

/* Implemented in Rust (src/qrcode/qrdec.rs) */
extern void qr_point_translate(qr_point _point, int _dx, int _dy);
extern unsigned qr_point_distance2(const qr_point _p1, const qr_point _p2);

/*Returns the cross product of the three points, which is positive if they are
   in CCW order (in a right-handed coordinate system), and 0 if they're
   colinear.
   Implemented in Rust (src/qrcode/qrdec.rs) */
extern int qr_point_ccw(const qr_point _p0, const qr_point _p1,
			const qr_point _p2);

/*Evaluates a line equation at a point.
  _line: The line to evaluate.
  _x:    The X coordinate of the point.
  _y:    The y coordinate of the point.
  Return: The value of the line equation _line[0]*_x+_line[1]*_y+_line[2].
  Implemented in Rust (src/qrcode/qrdec.rs) */
extern int qr_line_eval(qr_line _line, int _x, int _y);

/*Computes a line passing through the given point using the specified second
   order statistics.
  Given a line defined by the equation
    A*x+B*y+C = 0 ,
   the least squares fit to n points (x_i,y_i) must satisfy the two equations
    A^2 + (Syy - Sxx)/Sxy*A*B - B^2 = 0 ,
    C = -(xbar*A+ybar*B) ,
   where
    xbar = sum(x_i)/n ,
    ybar = sum(y_i)/n ,
    Sxx = sum((x_i-xbar)**2) ,
    Sxy = sum((x_i-xbar)*(y_i-ybar)) ,
    Syy = sum((y_i-ybar)**2) .
  The quadratic can be solved for the ratio (A/B) or (B/A):
    A/B = (Syy + sqrt((Sxx-Syy)**2 + 4*Sxy**2) - Sxx)/(-2*Sxy) ,
    B/A = (Sxx + sqrt((Sxx-Syy)**2 + 4*Sxy**2) - Syy)/(-2*Sxy) .
  We pick the one that leads to the larger ratio to avoid destructive
   cancellation (and e.g., 0/0 for horizontal or vertical lines).
  The above solutions correspond to the actual minimum.
  The other solution of the quadratic corresponds to a saddle point of the
   least squares objective function.
  _l:   Returns the fitted line values A, B, and C.
  _x0:  The X coordinate of the point the line is supposed to pass through.
  _y0:  The Y coordinate of the point the line is supposed to pass through.
  _sxx: The sum Sxx.
  _sxy: The sum Sxy.
  _syy: The sum Syy.
  _res: The maximum number of bits occupied by the product of any two of
         _l[0] or _l[1].
        Smaller numbers give less angular resolution, but allow more overhead
         room for computations.*/
/* Implemented in Rust (src/qrcode/qrdec.rs) */
extern void qr_line_fit(qr_line *_l, int _x0, int _y0, int _sxx, int _sxy,
			int _syy, int _res);

/*Perform a least-squares line fit to a list of points.
  At least two points are required.*/
/* Implemented in Rust (src/qrcode/qrdec.rs) */
extern void qr_line_fit_points(qr_line *_l, qr_point *_p, int _np, int _res);

extern void qr_line_orient(qr_line *_l, int _x, int _y);

extern int qr_line_isect(qr_point *_p, const qr_line *_l0, const qr_line *_l1);

/*An affine homography.
  This maps from the image (at subpel resolution) to a square domain with
   power-of-two sides (of res bits) and back.*/
struct qr_aff {
    int fwd[2][2];
    int inv[2][2];
    int x0;
    int y0;
    int res;
    int ires;
};

extern void qr_aff_init(qr_aff *_aff, const qr_point *_p0, const qr_point *_p1,
			const qr_point *_p2, int _res);

/*Map from the image (at subpel resolution) into the square domain.*/
extern void qr_aff_unproject(qr_point *_q, const qr_aff *_aff, int _x, int _y);

/*Map from the square domain into the image (at subpel resolution).*/
extern void qr_aff_project(qr_point *_p, const qr_aff *_aff, int _u, int _v);

/*A full homography.
  Like the affine homography, this maps from the image (at subpel resolution)
   to a square domain with power-of-two sides (of res bits) and back.*/
struct qr_hom {
    int fwd[3][2];
    int inv[3][2];
    int fwd22;
    int inv22;
    int x0;
    int y0;
    int res;
};

/*Initialize a homography mapping from four corner points.
  Computes both the forward and inverse homography transformations
  from the given corner points to a square domain.
  Implemented in Rust (src/qrcode/qrdec.rs) */
extern void qr_hom_init(qr_hom *_hom, int _x0, int _y0, int _x1, int _y1,
			int _x2, int _y2, int _x3, int _y3, int _res);

/*Map from the image (at subpel resolution) into the square domain.
  Returns a negative value if the point went to infinity.*/
extern int qr_hom_unproject(qr_point *_q, const qr_hom *_hom, int _x, int _y);

/*Finish a partial projection, converting from homogeneous coordinates to the
   normal 2-D representation.
  In loops, we can avoid many multiplies by computing the homogeneous _x, _y,
   and _w incrementally, but we cannot avoid the divisions, done here.*/
extern void qr_hom_fproject(qr_point *_p, const qr_hom *_hom, int _x, int _y,
			    int _w);

/*All the information we've collected about a finder pattern in the current
   configuration.*/
struct qr_finder {
    /*The module size along each axis (in the square domain).*/
    int size[2];
    /*The version estimated from the module size along each axis.*/
    int eversion[2];
    /*The list of classified edge points for each edge.*/
    qr_finder_edge_pt *edge_pts[4];
    /*The number of edge points classified as belonging to each edge.*/
    int nedge_pts[4];
    /*The number of inliers found after running RANSAC on each edge.*/
    int ninliers[4];
    /*The center of the finder pattern (in the square domain).*/
    qr_point o;
    /*The finder center information from the original image.*/
    qr_finder_center *c;
};

/*Computes the index of the edge each edge point belongs to, and its (signed)
   distance along the corresponding axis from the center of the finder pattern
   (in the square domain).
  The resulting list of edge points is sorted by edge index, with ties broken
   by extent.
  Implemented in Rust (src/qrcode/qrdec.rs) */
extern void qr_finder_edge_pts_aff_classify(qr_finder *_f, const qr_aff *_aff);

/*Computes the index of the edge each edge point belongs to, and its (signed)
   distance along the corresponding axis from the center of the finder pattern
   (in the square domain), using homography projection.
  The resulting list of edge points is sorted by edge index, with ties broken
   by extent.
  Implemented in Rust (src/qrcode/qrdec.rs) */
extern void qr_finder_edge_pts_hom_classify(qr_finder *_f, const qr_hom *_hom);

/*TODO: Perhaps these thresholds should be on the module size instead?
  Unfortunately, I'd need real-world images of codes with larger versions to
   see if these thresholds are still effective, but such versions aren't used
   often.*/

/*The amount that the estimated version numbers are allowed to differ from the
   real version number and still be considered valid.*/
#define QR_SMALL_VERSION_SLACK (1)
/*Since cell phone cameras can have severe radial distortion, the estimated
   version for larger versions can be off by larger amounts.*/
#define QR_LARGE_VERSION_SLACK (3)

/*Estimates the size of a module after classifying the edge points.
  _width:  The distance between UL and UR in the square domain.
  _height: The distance between UL and DL in the square domain.
  Returns 0 on success, or -1 if the module size or version could not be estimated.
  Implemented in Rust (src/qrcode/qrdec.rs) */
extern int qr_finder_estimate_module_size_and_version(qr_finder *_f, int _width,
						       int _height);

/*Eliminate outliers from the classified edge points with RANSAC.
  Uses the RANSAC (RANdom SAmple Consensus) algorithm to identify inliers
  among the edge points and eliminate outliers.
  Implemented in Rust (src/qrcode/qrdec.rs) */
extern void qr_finder_ransac(qr_finder *_f, const qr_aff *_hom,
			      isaac_ctx *_isaac, int _e);

/*Perform a least-squares line fit to an edge of a finder pattern using the
  inliers found by RANSAC.
  Implemented in Rust (src/qrcode/qrdec.rs) */
extern int qr_line_fit_finder_edge(qr_line *_l, const qr_finder *_f, int _e,
				    int _res);

/*Perform a least-squares line fit to a pair of common finder edges using the
  inliers found by RANSAC.
  Unlike a normal edge fit, we guarantee that this one succeeds by creating at
  least one point on each edge using the estimated module size if it has no inliers.
  Implemented in Rust (src/qrcode/qrdec.rs) */
extern void qr_line_fit_finder_pair(qr_line *_l, const qr_aff *_aff,
				     const qr_finder *_f0, const qr_finder *_f1,
				     int _e);

extern int qr_finder_quick_crossing_check(const unsigned char *_img, int _width,
					  int _height, int _x0, int _y0,
					  int _x1, int _y1, int _v);

/*Locate the midpoint of a _v segment along a !_v:_v:!_v line from (_x0,_y0) to
   (_x1,_y1).
  All coordinates, which are NOT in subpel resolution, must lie inside the
   image, and the endpoints are already assumed to have the value !_v.
  The returned value is in subpel resolution.*/
/* Locate the crossing of a finder or alignment pattern along a line.
   Implemented in Rust (src/qrcode/qrdec.rs) */
extern int qr_finder_locate_crossing(const unsigned char *_img, int _width,
				     int _height, int _x0, int _y0, int _x1,
				     int _y1, int _v, qr_point _p);

/* Calculate step delta for moving along a line in affine space.
   Implemented in Rust (src/qrcode/qrdec.rs) */
extern int qr_aff_line_step(const qr_aff *_aff, qr_line _l, int _v, int _du,
			    int *_dv);

/*Computes the Hamming distance between two bit patterns (the number of bits
   that differ).
  May stop counting after _maxdiff differences.*/
/* Implemented in Rust (src/qrcode/qrdec.rs) */
extern int qr_hamming_dist(unsigned _y1, unsigned _y2, int _maxdiff);

/*Retrieve a bit (guaranteed to be 0 or 1) from the image, given coordinates in
   subpel resolution which have not been bounds checked.*/
extern int qr_img_get_bit(const unsigned char *_img, int _width, int _height,
			  int _x, int _y);

/*A homography from one region of the grid back to the image.
  Unlike a qr_hom, this does not include an inverse transform and maps directly
   from the grid points, not a square with power-of-two sides.*/
struct qr_hom_cell {
    int fwd[3][3];
    int x0;
    int y0;
    int u0;
    int v0;
};

/* Initialize a homography cell for mapping between code and image space.
   Implemented in Rust (src/qrcode/qrdec.rs) */
extern void qr_hom_cell_init(qr_hom_cell *_cell, int _u0, int _v0, int _u1,
			     int _v1, int _u2, int _v2, int _u3, int _v3,
			     int _x0, int _y0, int _x1, int _y1, int _x2,
			     int _y2, int _x3, int _y3);

/*Finish a partial projection, converting from homogeneous coordinates to the
   normal 2-D representation.
  In loops, we can avoid many multiplies by computing the homogeneous _x, _y,
   and _w incrementally, but we cannot avoid the divisions, done here.*/
extern void qr_hom_cell_fproject(qr_point *_p, const qr_hom_cell *_cell, int _x,
				 int _y, int _w);

extern void qr_hom_cell_project(qr_point *_p, const qr_hom_cell *_cell, int _u,
				int _v, int _res);

/*Retrieves the bits corresponding to the alignment pattern template centered
   at the given location in the original image (at subpel precision).*/
/* Search for an alignment pattern near the given location.
   Implemented in Rust (src/qrcode/qrdec.rs) */
extern int qr_alignment_pattern_search(qr_point _p, const qr_hom_cell *_cell,
				       int _u, int _v, int _r,
				       const unsigned char *_img, int _width,
				       int _height);

/*Project alignment pattern center to corner of QR code.
  Given three corners and an alignment pattern center, compute the fourth corner
  by geometric projection. Returns 0 on success, -1 if projection fails.
  Implemented in Rust (src/qrcode/qrdec.rs) */
extern int qr_hom_project_alignment_to_corner(int *brx, int *bry,
					       const qr_point *p,
					       const qr_point *p3, int dim);

/*Fit a line to collected edge points, or use axis-aligned fallback if insufficient points.
  This is used for edges that have only one finder pattern, where we walk along the edge
  collecting sample points. If we don't get enough points (> 1), we fall back to an
  axis-aligned line in the affine coordinate system.*/
static void qr_hom_fit_edge_line(qr_line *line, qr_point *pts, int npts,
				  const qr_finder *finder, const qr_aff *aff,
				  int edge_axis)
{
    if (npts > 1) {
	qr_line_fit_points(line, pts, npts, aff->res);
    } else {
	qr_point p;
	int shift, round;

	/* Project reference point from the finder pattern */
	if (edge_axis == 1) {
	    /* Right edge: project from UR finder, extending 3 modules to the right */
	    qr_aff_project(&p, aff, finder->o[0] + 3 * finder->size[0],
			   finder->o[1]);
	} else {
	    /* Bottom edge (axis 3): project from DL finder, extending 3 modules down */
	    qr_aff_project(&p, aff, finder->o[0],
			   finder->o[1] + 3 * finder->size[1]);
	}

	/* Calculate normalization shift (always uses column 1 of affine matrix) */
	shift = QR_MAXI(0, qr_ilog(QR_MAXI(abs(aff->fwd[0][1]),
					    abs(aff->fwd[1][1]))) -
			((aff->res + 1) >> 1));
	round = (1 << shift) >> 1;

	/* Compute line coefficients using appropriate matrix column */
	if (edge_axis == 1) {
	    /* Right edge uses column 1 (vertical direction in affine space) */
	    (*line)[0] = (aff->fwd[1][1] + round) >> shift;
	    (*line)[1] = (-aff->fwd[0][1] + round) >> shift;
	} else {
	    /* Bottom edge uses column 0 (horizontal direction in affine space) */
	    (*line)[0] = (aff->fwd[1][0] + round) >> shift;
	    (*line)[1] = (-aff->fwd[0][0] + round) >> shift;
	}
	/* Compute line constant term */
	(*line)[2] = -((*line)[0] * p[0] + (*line)[1] * p[1]);
    }
}

static int qr_hom_fit(qr_hom *_hom, qr_finder *_ul, qr_finder *_ur,
		      qr_finder *_dl, qr_point _p[4], const qr_aff *_aff,
		      isaac_ctx *_isaac, const unsigned char *_img, int _width,
		      int _height)
{
    qr_point *b;
    int nb;
    int cb;
    qr_point *r;
    int nr;
    int cr;
    qr_line l[4];
    qr_point q;
    int ox;
    int oy;
    int ru;
    int rv;
    int dru;
    int drv;
    int bu;
    int bv;
    int dbu;
    int dbv;
    int rx;
    int ry;
    int drxi;
    int dryi;
    int drxj;
    int dryj;
    int rdone;
    int nrempty;
    int rlastfit;
    int bx;
    int by;
    int dbxi;
    int dbyi;
    int dbxj;
    int dbyj;
    int bdone;
    int nbempty;
    int blastfit;
    int version4;
    int brx;
    int bry;
    int i;
    /*We attempt to correct large-scale perspective distortion by fitting lines
   to the edge of the code area.
  We could also look for an alignment pattern now, but that wouldn't work for
   version 1 codes, which have no alignment pattern.
  Even if the code is supposed to have one, there's go guarantee we'd find it
   intact.*/
    /*Fitting lines is easy for the edges on which we have two finder patterns.
  After the fit, UL is guaranteed to be on the proper side, but if either of
   the other two finder patterns aren't, something is wrong.*/
    qr_finder_ransac(_ul, _aff, _isaac, 0);
    qr_finder_ransac(_dl, _aff, _isaac, 0);
    qr_line_fit_finder_pair(&l[0], _aff, _ul, _dl, 0);
    if (qr_line_eval(l[0], _dl->c->pos[0], _dl->c->pos[1]) < 0 ||
	qr_line_eval(l[0], _ur->c->pos[0], _ur->c->pos[1]) < 0) {
	return -1;
    }
    qr_finder_ransac(_ul, _aff, _isaac, 2);
    qr_finder_ransac(_ur, _aff, _isaac, 2);
    qr_line_fit_finder_pair(&l[2], _aff, _ul, _ur, 2);
    if (qr_line_eval(l[2], _dl->c->pos[0], _dl->c->pos[1]) < 0 ||
	qr_line_eval(l[2], _ur->c->pos[0], _ur->c->pos[1]) < 0) {
	return -1;
    }
    /*The edges which only have one finder pattern are more difficult.
  We start by fitting a line to the edge of the one finder pattern we do
   have.
  This can fail due to an insufficient number of sample points, and even if
   it succeeds can be fairly inaccurate, because all of the points are
   clustered in one corner of the QR code.
  If it fails, we just use an axis-aligned line in the affine coordinate
   system.
  Then we walk along the edge of the entire code, looking for
   light:dark:light patterns perpendicular to the edge.
  Wherever we find one, we take the center of the dark portion as an
   additional sample point.
  At the end, we re-fit the line using all such sample points found.*/
    drv = _ur->size[1] >> 1;
    qr_finder_ransac(_ur, _aff, _isaac, 1);
    if (qr_line_fit_finder_edge(&l[1], _ur, 1, _aff->res) >= 0) {
	if (qr_line_eval(l[1], _ul->c->pos[0], _ul->c->pos[1]) < 0 ||
	    qr_line_eval(l[1], _dl->c->pos[0], _dl->c->pos[1]) < 0) {
	    return -1;
	}
	/*Figure out the change in ru for a given change in rv when stepping along
   the fitted line.*/
	if (qr_aff_line_step(_aff, l[1], 1, drv, &dru) < 0)
	    return -1;
    } else
	dru = 0;
    ru	= _ur->o[0] + 3 * _ur->size[0] - 2 * dru;
    rv	= _ur->o[1] - 2 * drv;
    dbu = _dl->size[0] >> 1;
    qr_finder_ransac(_dl, _aff, _isaac, 3);
    if (qr_line_fit_finder_edge(&l[3], _dl, 3, _aff->res) >= 0) {
	if (qr_line_eval(l[3], _ul->c->pos[0], _ul->c->pos[1]) < 0 ||
	    qr_line_eval(l[3], _ur->c->pos[0], _ur->c->pos[1]) < 0) {
	    return -1;
	}
	/*Figure out the change in bv for a given change in bu when stepping along
   the fitted line.*/
	if (qr_aff_line_step(_aff, l[3], 0, dbu, &dbv) < 0)
	    return -1;
    } else
	dbv = 0;
    bu = _dl->o[0] - 2 * dbu;
    bv = _dl->o[1] + 3 * _dl->size[1] - 2 * dbv;
    /*Set up the initial point lists.*/
    nr = rlastfit = _ur->ninliers[1];
    cr		  = nr + (_dl->o[1] - rv + drv - 1) / drv;
    r		  = (qr_point *)malloc(cr * sizeof(*r));
    for (i = 0; i < _ur->ninliers[1]; i++) {
	memcpy(r[i], _ur->edge_pts[1][i].pos, sizeof(r[i]));
    }
    nb = blastfit = _dl->ninliers[3];
    cb		  = nb + (_ur->o[0] - bu + dbu - 1) / dbu;
    b		  = (qr_point *)malloc(cb * sizeof(*b));
    for (i = 0; i < _dl->ninliers[3]; i++) {
	memcpy(b[i], _dl->edge_pts[3][i].pos, sizeof(b[i]));
    }
    /*Set up the step parameters for the affine projection.*/
    ox	 = (_aff->x0 << _aff->res) + (1 << (_aff->res - 1));
    oy	 = (_aff->y0 << _aff->res) + (1 << (_aff->res - 1));
    rx	 = _aff->fwd[0][0] * ru + _aff->fwd[0][1] * rv + ox;
    ry	 = _aff->fwd[1][0] * ru + _aff->fwd[1][1] * rv + oy;
    drxi = _aff->fwd[0][0] * dru + _aff->fwd[0][1] * drv;
    dryi = _aff->fwd[1][0] * dru + _aff->fwd[1][1] * drv;
    drxj = _aff->fwd[0][0] * _ur->size[0];
    dryj = _aff->fwd[1][0] * _ur->size[0];
    bx	 = _aff->fwd[0][0] * bu + _aff->fwd[0][1] * bv + ox;
    by	 = _aff->fwd[1][0] * bu + _aff->fwd[1][1] * bv + oy;
    dbxi = _aff->fwd[0][0] * dbu + _aff->fwd[0][1] * dbv;
    dbyi = _aff->fwd[1][0] * dbu + _aff->fwd[1][1] * dbv;
    dbxj = _aff->fwd[0][1] * _dl->size[1];
    dbyj = _aff->fwd[1][1] * _dl->size[1];
    /*Now step along the lines, looking for new sample points.*/
    nrempty = nbempty = 0;
    for (;;) {
	int ret;
	int x0;
	int y0;
	int x1;
	int y1;
	/*If we take too many steps without encountering a non-zero pixel, assume
   we have wandered off the edge and stop looking before we hit the other
   side of the quiet region.
  Otherwise, stop when the lines cross (if they do so inside the affine
   region) or come close to crossing (outside the affine region).
  TODO: We don't have any way of detecting when we've wandered into the
   code interior; we could stop if the outside sample ever shows up dark,
   but this could happen because of noise in the quiet region, too.*/
	rdone = rv >= QR_MINI(bv, (_dl->o[1] + bv) >> 1) || nrempty > 14;
	bdone = bu >= QR_MINI(ru, (_ur->o[0] + ru) >> 1) || nbempty > 14;
	if (!rdone && (bdone || rv < bu)) {
	    x0 = (rx + drxj) >> (_aff->res + QR_FINDER_SUBPREC);
	    y0 = (ry + dryj) >> (_aff->res + QR_FINDER_SUBPREC);
	    x1 = (rx - drxj) >> (_aff->res + QR_FINDER_SUBPREC);
	    y1 = (ry - dryj) >> (_aff->res + QR_FINDER_SUBPREC);
	    if (nr >= cr) {
		cr = cr << 1 | 1;
		r  = (qr_point *)realloc(r, cr * sizeof(*r));
	    }
	    ret = qr_finder_quick_crossing_check(_img, _width, _height, x0, y0,
						 x1, y1, 1);
	    if (!ret) {
		ret = qr_finder_locate_crossing(_img, _width, _height, x0, y0,
						x1, y1, 1, r[nr]);
	    }
	    if (ret >= 0) {
		if (!ret) {
		    qr_aff_unproject(&q, _aff, r[nr][0], r[nr][1]);
		    /*Move the current point halfway towards the crossing.
  We don't move the whole way to give us some robustness to noise.*/
		    ru = (ru + q[0]) >> 1;
		    /*But ensure that rv monotonically increases.*/
		    if (q[1] + drv > rv)
			rv = (rv + q[1]) >> 1;
		    rx = _aff->fwd[0][0] * ru + _aff->fwd[0][1] * rv + ox;
		    ry = _aff->fwd[1][0] * ru + _aff->fwd[1][1] * rv + oy;
		    nr++;
		    /*Re-fit the line to update the step direction periodically.*/
		    if (nr > QR_MAXI(1, rlastfit + (rlastfit >> 2))) {
			qr_line_fit_points(&l[1], r, nr, _aff->res);
			if (qr_aff_line_step(_aff, l[1], 1, drv, &dru) >= 0) {
			    drxi =
				_aff->fwd[0][0] * dru + _aff->fwd[0][1] * drv;
			    dryi =
				_aff->fwd[1][0] * dru + _aff->fwd[1][1] * drv;
			}
			rlastfit = nr;
		    }
		}
		nrempty = 0;
	    } else
		nrempty++;
	    ru += dru;
	    /*Our final defense: if we overflow, stop.*/
	    if (rv + drv > rv)
		rv += drv;
	    else
		nrempty = INT_MAX;
	    rx += drxi;
	    ry += dryi;
	} else if (!bdone) {
	    x0 = (bx + dbxj) >> (_aff->res + QR_FINDER_SUBPREC);
	    y0 = (by + dbyj) >> (_aff->res + QR_FINDER_SUBPREC);
	    x1 = (bx - dbxj) >> (_aff->res + QR_FINDER_SUBPREC);
	    y1 = (by - dbyj) >> (_aff->res + QR_FINDER_SUBPREC);
	    if (nb >= cb) {
		cb = cb << 1 | 1;
		b  = (qr_point *)realloc(b, cb * sizeof(*b));
	    }
	    ret = qr_finder_quick_crossing_check(_img, _width, _height, x0, y0,
						 x1, y1, 1);
	    if (!ret) {
		ret = qr_finder_locate_crossing(_img, _width, _height, x0, y0,
						x1, y1, 1, b[nb]);
	    }
	    if (ret >= 0) {
		if (!ret) {
		    qr_aff_unproject(&q, _aff, b[nb][0], b[nb][1]);
		    /*Move the current point halfway towards the crossing.
  We don't move the whole way to give us some robustness to noise.*/
		    /*But ensure that bu monotonically increases.*/
		    if (q[0] + dbu > bu)
			bu = (bu + q[0]) >> 1;
		    bv = (bv + q[1]) >> 1;
		    bx = _aff->fwd[0][0] * bu + _aff->fwd[0][1] * bv + ox;
		    by = _aff->fwd[1][0] * bu + _aff->fwd[1][1] * bv + oy;
		    nb++;
		    /*Re-fit the line to update the step direction periodically.*/
		    if (nb > QR_MAXI(1, blastfit + (blastfit >> 2))) {
			qr_line_fit_points(&l[3], b, nb, _aff->res);
			if (qr_aff_line_step(_aff, l[3], 0, dbu, &dbv) >= 0) {
			    dbxi =
				_aff->fwd[0][0] * dbu + _aff->fwd[0][1] * dbv;
			    dbyi =
				_aff->fwd[1][0] * dbu + _aff->fwd[1][1] * dbv;
			}
			blastfit = nb;
		    }
		}
		nbempty = 0;
	    } else
		nbempty++;
	    /*Our final defense: if we overflow, stop.*/
	    if (bu + dbu > bu)
		bu += dbu;
	    else
		nbempty = INT_MAX;
	    bv += dbv;
	    bx += dbxi;
	    by += dbyi;
	} else
	    break;
    }
    /*Fit the new lines.
  If we _still_ don't have enough sample points, then just use an
   axis-aligned line from the affine coordinate system (e.g., one parallel
   to the opposite edge in the image).*/
    qr_hom_fit_edge_line(&l[1], r, nr, _ur, _aff, 1);
    free(r);
    qr_hom_fit_edge_line(&l[3], b, nb, _dl, _aff, 3);
    free(b);
    for (i = 0; i < 4; i++) {
	if (qr_line_isect(&_p[i], &l[i & 1], &l[2 + (i >> 1)]) < 0)
	    return -1;
	/*It's plausible for points to be somewhat outside the image, but too far
   and too much of the pattern will be gone for it to be decodable.*/
	if (_p[i][0] < (-_width << QR_FINDER_SUBPREC) ||
	    _p[i][0] >= ((_width << QR_FINDER_SUBPREC) + 1) ||
	    _p[i][1] < (-_height << QR_FINDER_SUBPREC) ||
	    _p[i][1] >= ((_height << QR_FINDER_SUBPREC) + 1)) {
	    return -1;
	}
    }
    /*By default, use the edge intersection point for the bottom-right corner.*/
    brx = _p[3][0];
    bry = _p[3][1];
    /*However, if our average version estimate is greater than 1, NOW we try to
   search for an alignment pattern.
  We get a much better success rate by doing this after our initial attempt
   to promote the transform to a homography than before.
  You might also think it would be more reliable to use the interior finder
   pattern edges, since the outer ones may be obscured or damaged, and it
   would save us a reprojection below, since they would form a nice square
   with the location of the alignment pattern, but this turns out to be a bad
   idea.
  Non-linear distortion is usually maximal on the outside edge, and thus
   estimating the grid position from points on the interior means we might
   get mis-aligned by the time we reach the edge.*/
    version4 = _ul->eversion[0] + _ul->eversion[1] + _ur->eversion[0] +
	       _dl->eversion[1];
    if (version4 > 4) {
	qr_hom_cell cell;
	qr_point p3;
	int dim;
	dim = 17 + version4;
	qr_hom_cell_init(&cell, 0, 0, dim - 1, 0, 0, dim - 1, dim - 1, dim - 1,
			 _p[0][0], _p[0][1], _p[1][0], _p[1][1], _p[2][0],
			 _p[2][1], _p[3][0], _p[3][1]);
	if (qr_alignment_pattern_search(p3, &cell, dim - 7, dim - 7, 4, _img,
					_width, _height) >= 0) {
	    /*There's no real need to update the bounding box corner, and in fact we
   actively perform worse if we do.
  Clearly it was good enough for us to find this alignment pattern, so
   it should be good enough to use for grid initialization.
  The point of doing the search was to get more accurate version
   estimates and a better chance of decoding the version and format info.
  This is particularly important for small versions that have no encoded
   version info, since any mismatch in version renders the code
   undecodable.*/
	    /*We do, however, need four points in a square to initialize our
   homography, so project the point from the alignment center to the
   corner of the code area.*/
	    if (qr_hom_project_alignment_to_corner(&brx, &bry, _p, &p3, dim) < 0)
		return -1;
	}
    }
    /*Now we have four points that map to a square: initialize the projection.*/
    qr_hom_init(_hom, _p[0][0], _p[0][1], _p[1][0], _p[1][1], _p[2][0],
		_p[2][1], brx, bry, QR_HOM_BITS);
    return 0;
}

/*The BCH(18,6,3) codes and correction function are implemented in Rust.
  See src/qrcode/qrdec.rs */

/*Corrects a BCH(18,6,3) code word.
  _y: Contains the code word to be checked on input, and the corrected value on
       output.
  Return: The number of errors.
          If more than 3 errors are detected, returns a negative value and
           performs no correction.
  Implemented in Rust (src/qrcode/qrdec.rs) */
extern int bch18_6_correct(unsigned *_y);

/*Reads the version bits near a finder module and decodes the version number.
  Implemented in Rust (src/qrcode/qrdec.rs) */
extern int qr_finder_version_decode(qr_finder *_f, const qr_hom *_hom,
				    const unsigned char *_img, int _width,
				    int _height, int _dir);

/*Reads the format info bits near the finder modules and decodes them.
  Implemented in Rust (src/qrcode/qrdec.rs) */
extern int qr_finder_fmt_info_decode(qr_finder *_ul, qr_finder *_ur,
				     qr_finder *_dl, const qr_hom *_hom,
				     const unsigned char *_img, int _width,
				     int _height);

/*The grid used to sample the image bits.
  The grid is divided into separate cells bounded by finder patterns and/or
   alignment patterns, and a separate map back to the original image is
   constructed for each cell.
  All of these structural elements, as well as the timing patterns, version
   info, and format info, are marked in fpmask so they can easily be skipped
   during decode.*/
struct qr_sampling_grid {
    qr_hom_cell *cells[6];
    unsigned *fpmask;
    int cell_limits[6];
    int ncells;
};

/*Mark a given region as belonging to the function pattern.*/
/* Implemented in Rust (src/qrcode/qrdec.rs) */
extern void qr_sampling_grid_fp_mask_rect(qr_sampling_grid *_grid, int _dim,
					  int _u, int _v, int _w, int _h);

/*Determine if a given grid location is inside the function pattern.
  Implemented in Rust (src/qrcode/qrdec.rs) */
extern int qr_sampling_grid_is_in_fp(const qr_sampling_grid *_grid, int _dim,
				     int _u, int _v);

/*The spacing between alignment patterns after the second for versions >= 7.
  We could compact this more, but the code to access it would eliminate the
   gains.
  Implemented in Rust (src/qrcode/qrdec.rs) */
extern const unsigned char QR_ALIGNMENT_SPACING[34];

/*Initialize the sampling grid for each region of the code.
  _version:  The (decoded) version number.
  _ul_pos:   The location of the UL finder pattern.
  _ur_pos:   The location of the UR finder pattern.
  _dl_pos:   The location of the DL finder pattern.
  _p:        On input, contains estimated positions of the four corner modules.
             On output, contains a bounding quadrilateral for the code.
  _img:      The binary input image.
  _width:    The width of the input image.
  _height:   The height of the input image.
  Return: 0 on success, or a negative value on error.*/
/*Initialize a QR sampling grid.
  Implemented in Rust (src/qrcode/qrdec.rs) */
extern void qr_sampling_grid_init(qr_sampling_grid *_grid, int _version,
                                   const qr_point _ul_pos,
                                   const qr_point _ur_pos,
                                   const qr_point _dl_pos, qr_point _p[4],
                                   const unsigned char *_img, int _width,
                                   int _height);

/*Clear a QR sampling grid.
  Implemented in Rust (src/qrcode/qrdec.rs) */
extern void qr_sampling_grid_clear(qr_sampling_grid *_grid);

/*Generate the data mask corresponding to the given mask pattern.*/
/*Fill a data mask for QR code decoding.
  Implemented in Rust (src/qrcode/qrdec.rs) */
extern void qr_data_mask_fill(unsigned *_mask, int _dim, int _pattern);

/*Sample QR code data bits from the image using the sampling grid.
  Implemented in Rust (src/qrcode/qrdec.rs) */
extern void qr_sampling_grid_sample(const qr_sampling_grid *_grid,
				     unsigned *_data_bits, int _dim,
				     int _fmt_info, const unsigned char *_img,
				     int _width, int _height);

/*Arranges the sample bits read by qr_sampling_grid_sample() into bytes and
   groups those bytes into Reed-Solomon blocks.
  Implemented in Rust (src/qrcode/qrdec.rs) */
extern void qr_samples_unpack(unsigned char **_blocks, int _nblocks,
			       int _nshort_data, int _nshort_blocks,
			       const unsigned *_data_bits,
			       const unsigned *_fp_mask, int _dim);

/*Bit reading code blatantly stolen^W^Wadapted from libogg/libtheora (because
   I've already debugged it and I know it works).
  Portions (C) Xiph.Org Foundation 1994-2008, BSD-style license.
  Implemented in Rust (src/qrcode/qrdec.rs) */
struct qr_pack_buf {
    const unsigned char *buf;
    int endbyte;
    int endbit;
    int storage;
};

/* Implemented in Rust (src/qrcode/qrdec.rs) */
extern void qr_pack_buf_init(qr_pack_buf *_b, const unsigned char *_data,
			     int _ndata);

/* Implemented in Rust (src/qrcode/qrdec.rs) */
/*Assumes 0<=_bits<=16.*/
extern int qr_pack_buf_read(qr_pack_buf *_b, int _bits);

/* Implemented in Rust (src/qrcode/qrdec.rs) */
extern int qr_pack_buf_avail(const qr_pack_buf *_b);

/*Parse QR code data from corrected codewords.
  Implemented in Rust (src/qrcode/qrdec.rs) */
extern int qr_code_data_parse(qr_code_data *_qrdata, int _version,
			       const unsigned char *_data, int _ndata);

/* Clear a QR code data structure, freeing all allocated memory.
   Implemented in Rust (src/qrcode/qrdec.rs) */
extern void qr_code_data_clear(qr_code_data *_qrdata);

/* Initialize a QR code data list.
   Implemented in Rust (src/qrcode/qrdec.rs) */
extern void qr_code_data_list_init(qr_code_data_list *_qrlist);

/* Clear a QR code data list, freeing all allocated memory.
   Implemented in Rust (src/qrcode/qrdec.rs) */
extern void qr_code_data_list_clear(qr_code_data_list *_qrlist);

/* Add a QR code data to the list.
   Implemented in Rust (src/qrcode/qrdec.rs) */
extern void qr_code_data_list_add(qr_code_data_list *_qrlist,
				  qr_code_data *_qrdata);

/*The total number of codewords in a QR code.
  Implemented in Rust (src/qrcode/qrdec.rs) */
extern int qr_code_ncodewords(unsigned _version);

/*Bulk data for the number of parity bytes per Reed-Solomon block.
  Implemented in Rust (src/qrcode/qrdec.rs) */
extern const unsigned char QR_RS_NPAR_VALS[71];

/*An offset into QR_RS_NPAR_VALS for each version that gives the number of
   parity bytes per Reed-Solomon block for each error correction level.
  Implemented in Rust (src/qrcode/qrdec.rs) */
extern const unsigned char QR_RS_NPAR_OFFS[40];

/*The number of Reed-Solomon blocks for each version and error correction
   level.
  Implemented in Rust (src/qrcode/qrdec.rs) */
extern const unsigned char QR_RS_NBLOCKS[40][4];

/*Decode a QR code from an image.
  Implemented in Rust (src/qrcode/qrdec.rs) */
extern int qr_code_decode(qr_code_data *_qrdata, const rs_gf256 *_gf,
			   const qr_point *_ul_pos, const qr_point *_ur_pos,
			   const qr_point *_dl_pos, int _version, int _fmt_info,
			   const unsigned char *_img, int _width, int _height);

/*Searches for an arrangement of these three finder centers that yields a valid
   configuration.
  _c: On input, the three finder centers to consider in any order.
  Return: The detected version number, or a negative value on error.*/
static int qr_reader_try_configuration(qr_reader *_reader,
				       qr_code_data *_qrdata,
				       const unsigned char *_img, int _width,
				       int _height, qr_finder_center *_c[3])
{
    int ci[7];
    unsigned maxd;
    int ccw;
    int i0;
    int i;
    /*Sort the points in counter-clockwise order.*/
    ccw = qr_point_ccw(_c[0]->pos, _c[1]->pos, _c[2]->pos);
    /*Colinear points can't be the corners of a quadrilateral.*/
    if (!ccw)
	return -1;
    /*Include a few extra copies of the cyclical list to avoid mods.*/
    ci[6] = ci[3] = ci[0] = 0;
    ci[4] = ci[1] = 1 + (ccw < 0);
    ci[5] = ci[2] = 2 - (ccw < 0);
    /*Assume the points farthest from each other are the opposite corners, and
   find the top-left point.*/
    maxd = qr_point_distance2(_c[1]->pos, _c[2]->pos);
    i0	 = 0;
    for (i = 1; i < 3; i++) {
	unsigned d;
	d = qr_point_distance2(_c[ci[i + 1]]->pos, _c[ci[i + 2]]->pos);
	if (d > maxd) {
	    i0	 = i;
	    maxd = d;
	}
    }
    /*However, try all three possible orderings, just to be sure (a severely
   skewed projection could move opposite corners closer than adjacent).*/
    for (i = i0; i < i0 + 3; i++) {
	qr_aff aff;
	qr_hom hom;
	qr_finder ul;
	qr_finder ur;
	qr_finder dl;
	qr_point bbox[4];
	int res;
	int ur_version;
	int dl_version;
	int fmt_info;
	ul.c = _c[ci[i]];
	ur.c = _c[ci[i + 1]];
	dl.c = _c[ci[i + 2]];
	/*Estimate the module size and version number from the two opposite corners.
  The module size is not constant in the image, so we compute an affine
   projection from the three points we have to a square domain, and
   estimate it there.
  Although it should be the same along both axes, we keep separate
   estimates to account for any remaining projective distortion.*/
	res = QR_INT_BITS - 2 - QR_FINDER_SUBPREC -
	      qr_ilog(QR_MAXI(_width, _height) - 1);
	qr_aff_init(&aff, &ul.c->pos, &ur.c->pos, &dl.c->pos, res);
	qr_aff_unproject(&ur.o, &aff, ur.c->pos[0], ur.c->pos[1]);
	qr_finder_edge_pts_aff_classify(&ur, &aff);
	if (qr_finder_estimate_module_size_and_version(&ur, 1 << res,
						       1 << res) < 0)
	    continue;
	qr_aff_unproject(&dl.o, &aff, dl.c->pos[0], dl.c->pos[1]);
	qr_finder_edge_pts_aff_classify(&dl, &aff);
	if (qr_finder_estimate_module_size_and_version(&dl, 1 << res,
						       1 << res) < 0)
	    continue;
	/*If the estimated versions are significantly different, reject the
   configuration.*/
	if (abs(ur.eversion[1] - dl.eversion[0]) > QR_LARGE_VERSION_SLACK)
	    continue;
	qr_aff_unproject(&ul.o, &aff, ul.c->pos[0], ul.c->pos[1]);
	qr_finder_edge_pts_aff_classify(&ul, &aff);
	if (qr_finder_estimate_module_size_and_version(&ul, 1 << res,
						       1 << res) < 0 ||
	    abs(ul.eversion[1] - ur.eversion[1]) > QR_LARGE_VERSION_SLACK ||
	    abs(ul.eversion[0] - dl.eversion[0]) > QR_LARGE_VERSION_SLACK) {
	    continue;
	}
#if defined(QR_DEBUG)
	qr_finder_dump_aff_undistorted(&ul, &ur, &dl, &aff, _img, _width,
				       _height);
#endif
	/*If we made it this far, upgrade the affine homography to a full
   homography.*/
	if (qr_hom_fit(&hom, &ul, &ur, &dl, bbox, &aff, &_reader->isaac, _img,
		       _width, _height) < 0) {
	    continue;
	}
	memcpy(_qrdata->bbox, bbox, sizeof(bbox));
	qr_hom_unproject(&ul.o, &hom, ul.c->pos[0], ul.c->pos[1]);
	qr_hom_unproject(&ur.o, &hom, ur.c->pos[0], ur.c->pos[1]);
	qr_hom_unproject(&dl.o, &hom, dl.c->pos[0], dl.c->pos[1]);
	qr_finder_edge_pts_hom_classify(&ur, &hom);
	if (qr_finder_estimate_module_size_and_version(&ur, ur.o[0] - ul.o[0],
						       ur.o[0] - ul.o[0]) < 0) {
	    continue;
	}
	qr_finder_edge_pts_hom_classify(&dl, &hom);
	if (qr_finder_estimate_module_size_and_version(&dl, dl.o[1] - ul.o[1],
						       dl.o[1] - ul.o[1]) < 0) {
	    continue;
	}
#if defined(QR_DEBUG)
	qr_finder_dump_hom_undistorted(&ul, &ur, &dl, &hom, _img, _width,
				       _height);
#endif
	/*If we have a small version (less than 7), there's no encoded version
   information.
  If the estimated version on the two corners matches and is sufficiently
   small, we assume this is the case.*/
	if (ur.eversion[1] == dl.eversion[0] && ur.eversion[1] < 7) {
	    ur_version = ur.eversion[1];
	} else {
	    /*If the estimated versions are significantly different, reject the
   configuration.*/
	    if (abs(ur.eversion[1] - dl.eversion[0]) > QR_LARGE_VERSION_SLACK)
		continue;
	    /*Otherwise we try to read the actual version data from the image.
  If the real version is not sufficiently close to our estimated version,
   then we assume there was an unrecoverable decoding error (so many bit
   errors we were within 3 errors of another valid code), and throw that
   value away.
  If no decoded version could be sufficiently close, we don't even try.*/
	    if (ur.eversion[1] >= 7 - QR_LARGE_VERSION_SLACK) {
		ur_version = qr_finder_version_decode(&ur, &hom, _img, _width,
						      _height, 0);
		if (abs(ur_version - ur.eversion[1]) > QR_LARGE_VERSION_SLACK)
		    ur_version = -1;
	    } else
		ur_version = -1;
	    if (dl.eversion[0] >= 7 - QR_LARGE_VERSION_SLACK) {
		dl_version = qr_finder_version_decode(&dl, &hom, _img, _width,
						      _height, 1);
		if (abs(dl_version - dl.eversion[0]) > QR_LARGE_VERSION_SLACK)
		    dl_version = -1;
	    } else
		dl_version = -1;
	    /*If we got at least one valid version, or we got two and they match,
   then we found a valid configuration.*/
	    if (ur_version >= 0) {
		if (dl_version >= 0 && dl_version != ur_version)
		    continue;
	    } else if (dl_version < 0)
		continue;
	    else
		ur_version = dl_version;
	}
	qr_finder_edge_pts_hom_classify(&ul, &hom);
	if (qr_finder_estimate_module_size_and_version(&ul, ur.o[0] - dl.o[0],
						       dl.o[1] - ul.o[1]) < 0 ||
	    abs(ul.eversion[1] - ur.eversion[1]) > QR_SMALL_VERSION_SLACK ||
	    abs(ul.eversion[0] - dl.eversion[0]) > QR_SMALL_VERSION_SLACK) {
	    continue;
	}
	fmt_info = qr_finder_fmt_info_decode(&ul, &ur, &dl, &hom, _img, _width,
					     _height);
	if (fmt_info < 0 ||
	    qr_code_decode(_qrdata, &_reader->gf, &ul.c->pos, &ur.c->pos,
			   &dl.c->pos, ur_version, fmt_info, _img, _width,
			   _height) < 0) {
	    /*The code may be flipped.
  Try again, swapping the UR and DL centers.
  We should get a valid version either way, so it's relatively cheap to
   check this, as we've already filtered out a lot of invalid
   configurations.*/
	    QR_SWAP2I(hom.inv[0][0], hom.inv[1][0]);
	    QR_SWAP2I(hom.inv[0][1], hom.inv[1][1]);
	    QR_SWAP2I(hom.fwd[0][0], hom.fwd[0][1]);
	    QR_SWAP2I(hom.fwd[1][0], hom.fwd[1][1]);
	    QR_SWAP2I(hom.fwd[2][0], hom.fwd[2][1]);
	    QR_SWAP2I(ul.o[0], ul.o[1]);
	    QR_SWAP2I(ul.size[0], ul.size[1]);
	    QR_SWAP2I(ur.o[0], ur.o[1]);
	    QR_SWAP2I(ur.size[0], ur.size[1]);
	    QR_SWAP2I(dl.o[0], dl.o[1]);
	    QR_SWAP2I(dl.size[0], dl.size[1]);
#if defined(QR_DEBUG)
	    qr_finder_dump_hom_undistorted(&ul, &dl, &ur, &hom, _img, _width,
					   _height);
#endif
	    fmt_info = qr_finder_fmt_info_decode(&ul, &dl, &ur, &hom, _img,
						 _width, _height);
	    if (fmt_info < 0)
		continue;
	    QR_SWAP2I(bbox[1][0], bbox[2][0]);
	    QR_SWAP2I(bbox[1][1], bbox[2][1]);
	    memcpy(_qrdata->bbox, bbox, sizeof(bbox));
	    if (qr_code_decode(_qrdata, &_reader->gf, &ul.c->pos, &dl.c->pos,
			       &ur.c->pos, ur_version, fmt_info, _img, _width,
			       _height) < 0) {
		continue;
	    }
	}
	return ur_version;
    }
    return -1;
}

void qr_reader_match_centers(qr_reader *_reader, qr_code_data_list *_qrlist,
			     qr_finder_center *_centers, int _ncenters,
			     const unsigned char *_img, int _width, int _height)
{
    /*The number of centers should be small, so an O(n^3) exhaustive search of
   which ones go together should be reasonable.*/
    unsigned char *mark;
    int nfailures_max;
    int nfailures;
    int i;
    int j;
    int k;
    mark	  = (unsigned char *)calloc(_ncenters, sizeof(*mark));
    nfailures_max = QR_MAXI(8192, _width * _height >> 9);
    nfailures	  = 0;
    for (i = 0; i < _ncenters; i++) {
	/*TODO: We might be able to accelerate this step significantly by
   considering the remaining finder centers in a more intelligent order,
   based on the first finder center we just chose.*/
	for (j = i + 1; i < _ncenters && !mark[i] && j < _ncenters; j++) {
	    for (k = j + 1; j < _ncenters && !mark[j] && k < _ncenters; k++)
		if (!mark[k]) {
		    qr_finder_center *c[3];
		    qr_code_data qrdata;
		    int version;
		    c[0]    = _centers + i;
		    c[1]    = _centers + j;
		    c[2]    = _centers + k;
		    version = qr_reader_try_configuration(_reader, &qrdata,
							  _img, _width, _height,
							  c);
		    if (version >= 0) {
			int ninside;
			int l;
			/*Add the data to the list.*/
			qr_code_data_list_add(_qrlist, &qrdata);
			/*Convert the bounding box we're returning to the user to normal
 image coordinates.*/
			for (l = 0; l < 4; l++) {
			    _qrlist->qrdata[_qrlist->nqrdata - 1].bbox[l][0] >>=
				QR_FINDER_SUBPREC;
			    _qrlist->qrdata[_qrlist->nqrdata - 1].bbox[l][1] >>=
				QR_FINDER_SUBPREC;
			}
			/*Mark these centers as used.*/
			mark[i] = mark[j] = mark[k] = 1;
			/*Find any other finder centers located inside this code.*/
			for (l = ninside = 0; l < _ncenters; l++)
			    if (!mark[l]) {
				if (qr_point_ccw(qrdata.bbox[0], qrdata.bbox[1],
						 _centers[l].pos) >= 0 &&
				    qr_point_ccw(qrdata.bbox[1], qrdata.bbox[3],
						 _centers[l].pos) >= 0 &&
				    qr_point_ccw(qrdata.bbox[3], qrdata.bbox[2],
						 _centers[l].pos) >= 0 &&
				    qr_point_ccw(qrdata.bbox[2], qrdata.bbox[0],
						 _centers[l].pos) >= 0) {
				    mark[l] = 2;
				    ninside++;
				}
			    }
			if (ninside >= 3) {
			    /*We might have a "Double QR": a code inside a code.
Copy the relevant centers to a new array and do a search confined
 to that subset.*/
			    qr_finder_center *inside;
			    inside = (qr_finder_center *)malloc(
				ninside * sizeof(*inside));
			    for (l = ninside = 0; l < _ncenters; l++) {
				if (mark[l] == 2)
				    *&inside[ninside++] = *&_centers[l];
			    }
			    qr_reader_match_centers(_reader, _qrlist, inside,
						    ninside, _img, _width,
						    _height);
			    free(inside);
			}
			/*Mark _all_ such centers used: codes cannot partially overlap.*/
			for (l = 0; l < _ncenters; l++)
			    if (mark[l] == 2)
				mark[l] = 1;
			nfailures = 0;
		    } else if (++nfailures > nfailures_max) {
			/*Give up.
We're unlikely to find a valid code in all this clutter, and we
 could spent quite a lot of time trying.*/
			i = j = k = _ncenters;
		    }
		}
	}
    }
    free(mark);
}

int _zbar_qr_found_line(qr_reader *reader, int dir, const qr_finder_line *line)
{
    /* minimally intrusive brute force version */
    qr_finder_lines *lines = &reader->finder_lines[dir];

    if (lines->nlines >= lines->clines) {
	lines->clines *= 2;
	lines->lines =
	    realloc(lines->lines, ++lines->clines * sizeof(*lines->lines));
    }

    memcpy(lines->lines + lines->nlines++, line, sizeof(*line));

    return (0);
}

int _zbar_qr_decode(qr_reader *reader, zbar_image_scanner_t *iscn,
		    zbar_image_t *img)
{
    int nqrdata			= 0, ncenters;
    qr_finder_edge_pt *edge_pts = NULL;
    qr_finder_center *centers	= NULL;

    if (reader->finder_lines[0].nlines < 9 ||
	reader->finder_lines[1].nlines < 9)
	return (0);

    ncenters = qr_finder_centers_locate(&centers, &edge_pts, reader, 0, 0);

    if (ncenters >= 3) {
	void *bin = qr_binarize(img->data, img->width, img->height);

	qr_code_data_list qrlist;
	qr_code_data_list_init(&qrlist);

	qr_reader_match_centers(reader, &qrlist, centers, ncenters, bin,
				img->width, img->height);

	if (qrlist.nqrdata > 0)
	    nqrdata = qr_code_data_list_extract_text(&qrlist, iscn, img);

	qr_code_data_list_clear(&qrlist);
	free(bin);
    }

    if (centers)
	free(centers);
    if (edge_pts)
	free(edge_pts);
    return (nqrdata);
}
