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
  axis-aligned line in the affine coordinate system.
  Implemented in Rust (src/qrcode/qrdec.rs) */
extern void qr_hom_fit_edge_line(qr_line *line, qr_point *pts, int npts,
				  const qr_finder *finder, const qr_aff *aff,
				  int edge_axis);

/*Fit a homography to correct large-scale perspective distortion.
  Implemented in Rust (src/qrcode/qrdec.rs) */
extern int qr_hom_fit(qr_hom *_hom, qr_finder *_ul, qr_finder *_ur,
		      qr_finder *_dl, qr_point _p[4], const qr_aff *_aff,
		      isaac_ctx *_isaac, const unsigned char *_img, int _width,
		      int _height);

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
/*Try to decode a QR code with the given configuration of three finder patterns.
  Implemented in Rust (src/qrcode/qrdec.rs) */
extern int qr_reader_try_configuration(qr_reader *_reader,
				       qr_code_data *_qrdata,
				       const unsigned char *_img, int _width,
				       int _height, qr_finder_center *_c[3]);

/* Match finder centers and decode QR codes.
   Implemented in Rust (src/qrcode/qrdec.rs) */
extern void qr_reader_match_centers(qr_reader *_reader, qr_code_data_list *_qrlist,
				    qr_finder_center *_centers, int _ncenters,
				    const unsigned char *_img, int _width, int _height);

/* Add a found finder line to the reader's line list.
   Implemented in Rust (src/qrcode/qrdec.rs) */
extern int _zbar_qr_found_line(qr_reader *reader, int dir, const qr_finder_line *line);

/* Decode QR codes from an image.
   Implemented in Rust (src/qrcode/qrdec.rs) */
extern int _zbar_qr_decode(qr_reader *reader, zbar_image_scanner_t *iscn,
			   zbar_image_t *img);
