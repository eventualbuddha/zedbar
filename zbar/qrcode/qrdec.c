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
  Return: The number of clusters.*/
static int qr_finder_cluster_lines(qr_finder_cluster *_clusters,
				   qr_finder_line **_neighbors,
				   qr_finder_line *_lines, int _nlines, int _v)
{
    unsigned char *mark;
    qr_finder_line **neighbors;
    int nneighbors;
    int nclusters;
    int i;
    /*TODO: Kalman filters!*/
    mark      = (unsigned char *)calloc(_nlines, sizeof(*mark));
    neighbors = _neighbors;
    nclusters = 0;
    for (i = 0; i < _nlines - 1; i++)
	if (!mark[i]) {
	    int len;
	    int j;
	    nneighbors	 = 1;
	    neighbors[0] = _lines + i;
	    len		 = _lines[i].len;
	    for (j = i + 1; j < _nlines; j++)
		if (!mark[j]) {
		    const qr_finder_line *a;
		    const qr_finder_line *b;
		    int thresh;
		    a = neighbors[nneighbors - 1];
		    b = _lines + j;
		    /*The clustering threshold is proportional to the size of the lines,
since minor noise in large areas can interrupt patterns more easily
at high resolutions.*/
		    thresh = (a->len + 7) >> 2;
		    if (abs(a->pos[1 - _v] - b->pos[1 - _v]) > thresh)
			break;
		    if (abs(a->pos[_v] - b->pos[_v]) > thresh)
			continue;
		    if (abs(a->pos[_v] + a->len - b->pos[_v] - b->len) > thresh)
			continue;
		    if (a->boffs > 0 && b->boffs > 0 &&
			abs(a->pos[_v] - a->boffs - b->pos[_v] + b->boffs) >
			    thresh) {
			continue;
		    }
		    if (a->eoffs > 0 && b->eoffs > 0 &&
			abs(a->pos[_v] + a->len + a->eoffs - b->pos[_v] -
			    b->len - b->eoffs) > thresh) {
			continue;
		    }
		    neighbors[nneighbors++] = _lines + j;
		    len += b->len;
		}
	    /*We require at least three lines to form a cluster, which eliminates a
 large number of false positives, saving considerable decoding time.
This should still be sufficient for 1-pixel codes with no noise.*/
	    if (nneighbors < 3)
		continue;
	    /*The expected number of lines crossing a finder pattern is equal to their
 average length.
We accept the cluster if size is at least 1/3 their average length (this
 is a very small threshold, but was needed for some test images).*/
	    len = ((len << 1) + nneighbors) / (nneighbors << 1);
	    if (nneighbors * (5 << QR_FINDER_SUBPREC) >= len) {
		_clusters[nclusters].lines  = neighbors;
		_clusters[nclusters].nlines = nneighbors;
		for (j = 0; j < nneighbors; j++)
		    mark[neighbors[j] - _lines] = 1;
		neighbors += nneighbors;
		nclusters++;
	    }
	}
    free(mark);
    return nclusters;
}

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
  Return: The new total number of edge points.*/
static int qr_finder_edge_pts_fill(qr_finder_edge_pt *_edge_pts, int _nedge_pts,
				   qr_finder_cluster **_neighbors,
				   int _nneighbors, int _v)
{
    int i;
    for (i = 0; i < _nneighbors; i++) {
	qr_finder_cluster *c;
	int j;
	c = _neighbors[i];
	for (j = 0; j < c->nlines; j++) {
	    qr_finder_line *l;
	    l = c->lines[j];
	    if (l->boffs > 0) {
		_edge_pts[_nedge_pts].pos[0] = l->pos[0];
		_edge_pts[_nedge_pts].pos[1] = l->pos[1];
		_edge_pts[_nedge_pts].pos[_v] -= l->boffs;
		_nedge_pts++;
	    }
	    if (l->eoffs > 0) {
		_edge_pts[_nedge_pts].pos[0] = l->pos[0];
		_edge_pts[_nedge_pts].pos[1] = l->pos[1];
		_edge_pts[_nedge_pts].pos[_v] += l->len + l->eoffs;
		_nedge_pts++;
	    }
	}
    }
    return _nedge_pts;
}

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
  Return: The number of putative finder centers.*/
static int qr_finder_find_crossings(qr_finder_center *_centers,
				    qr_finder_edge_pt *_edge_pts,
				    qr_finder_cluster *_hclusters,
				    int _nhclusters,
				    qr_finder_cluster *_vclusters,
				    int _nvclusters)
{
    qr_finder_cluster **hneighbors;
    qr_finder_cluster **vneighbors;
    unsigned char *hmark;
    unsigned char *vmark;
    int ncenters;
    int i;
    int j;
    hneighbors =
	(qr_finder_cluster **)malloc(_nhclusters * sizeof(*hneighbors));
    vneighbors =
	(qr_finder_cluster **)malloc(_nvclusters * sizeof(*vneighbors));
    hmark    = (unsigned char *)calloc(_nhclusters, sizeof(*hmark));
    vmark    = (unsigned char *)calloc(_nvclusters, sizeof(*vmark));
    ncenters = 0;
    /*TODO: This may need some re-working.
  We should be finding groups of clusters such that _all_ horizontal lines in
   _all_ horizontal clusters in the group cross _all_ vertical lines in _all_
   vertical clusters in the group.
  This is equivalent to finding the maximum bipartite clique in the
   connectivity graph, which requires linear progamming to solve efficiently.
  In principle, that is easy to do, but a realistic implementation without
   floating point is a lot of work (and computationally expensive).
  Right now we are relying on a sufficient border around the finder patterns
   to prevent false positives.*/
    for (i = 0; i < _nhclusters; i++)
	if (!hmark[i]) {
	    qr_finder_line *a;
	    qr_finder_line *b;
	    int nvneighbors;
	    int nedge_pts;
	    int y;
	    a = _hclusters[i].lines[_hclusters[i].nlines >> 1];
	    y = nvneighbors = 0;
	    for (j = 0; j < _nvclusters; j++)
		if (!vmark[j]) {
		    b = _vclusters[j].lines[_vclusters[j].nlines >> 1];
		    if (qr_finder_lines_are_crossing(a, b)) {
			vmark[j] = 1;
			y += (b->pos[1] << 1) + b->len;
			if (b->boffs > 0 && b->eoffs > 0)
			    y += b->eoffs - b->boffs;
			vneighbors[nvneighbors++] = _vclusters + j;
		    }
		}
	    if (nvneighbors > 0) {
		qr_finder_center *c;
		int nhneighbors;
		int x;
		x = (a->pos[0] << 1) + a->len;
		if (a->boffs > 0 && a->eoffs > 0)
		    x += a->eoffs - a->boffs;
		hneighbors[0] = _hclusters + i;
		nhneighbors   = 1;
		j	      = nvneighbors >> 1;
		b = vneighbors[j]->lines[vneighbors[j]->nlines >> 1];
		for (j = i + 1; j < _nhclusters; j++)
		    if (!hmark[j]) {
			a = _hclusters[j].lines[_hclusters[j].nlines >> 1];
			if (qr_finder_lines_are_crossing(a, b)) {
			    hmark[j] = 1;
			    x += (a->pos[0] << 1) + a->len;
			    if (a->boffs > 0 && a->eoffs > 0)
				x += a->eoffs - a->boffs;
			    hneighbors[nhneighbors++] = _hclusters + j;
			}
		    }
		c	     = _centers + ncenters++;
		c->pos[0]    = (x + nhneighbors) / (nhneighbors << 1);
		c->pos[1]    = (y + nvneighbors) / (nvneighbors << 1);
		c->edge_pts  = _edge_pts;
		nedge_pts    = qr_finder_edge_pts_fill(_edge_pts, 0, hneighbors,
						       nhneighbors, 0);
		nedge_pts    = qr_finder_edge_pts_fill(_edge_pts, nedge_pts,
						       vneighbors, nvneighbors, 1);
		c->nedge_pts = nedge_pts;
		_edge_pts += nedge_pts;
	    }
	}
    free(vmark);
    free(hmark);
    free(vneighbors);
    free(hneighbors);
    /*Sort the centers by decreasing numbers of edge points.*/
    qsort(_centers, ncenters, sizeof(*_centers), qr_finder_center_cmp);
    return ncenters;
}

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
static int qr_finder_centers_locate(qr_finder_center **_centers,
				    qr_finder_edge_pt **_edge_pts,
				    qr_reader *reader,
				    int _width __attribute__((unused)),
				    int _height __attribute__((unused)))
{
    qr_finder_line *hlines = reader->finder_lines[0].lines;
    int nhlines		   = reader->finder_lines[0].nlines;
    qr_finder_line *vlines = reader->finder_lines[1].lines;
    int nvlines		   = reader->finder_lines[1].nlines;

    qr_finder_line **hneighbors;
    qr_finder_cluster *hclusters;
    int nhclusters;
    qr_finder_line **vneighbors;
    qr_finder_cluster *vclusters;
    int nvclusters;
    int ncenters;

    /*Cluster the detected lines.*/
    hneighbors = (qr_finder_line **)malloc(nhlines * sizeof(*hneighbors));
    /*We require more than one line per cluster, so there are at most nhlines/2.*/
    hclusters =
	(qr_finder_cluster *)malloc((nhlines >> 1) * sizeof(*hclusters));
    nhclusters =
	qr_finder_cluster_lines(hclusters, hneighbors, hlines, nhlines, 0);
    /*We need vertical lines to be sorted by X coordinate, with ties broken by Y
   coordinate, for clustering purposes.
  We scan the image in the opposite order for cache efficiency, so sort the
   lines we found here.*/
    qsort(vlines, nvlines, sizeof(*vlines), qr_finder_vline_cmp);
    vneighbors = (qr_finder_line **)malloc(nvlines * sizeof(*vneighbors));
    /*We require more than one line per cluster, so there are at most nvlines/2.*/
    vclusters =
	(qr_finder_cluster *)malloc((nvlines >> 1) * sizeof(*vclusters));
    nvclusters =
	qr_finder_cluster_lines(vclusters, vneighbors, vlines, nvlines, 1);
    /*Find line crossings among the clusters.*/
    if (nhclusters >= 3 && nvclusters >= 3) {
	qr_finder_edge_pt *edge_pts;
	qr_finder_center *centers;
	int nedge_pts;
	int i;
	nedge_pts = 0;
	for (i = 0; i < nhclusters; i++)
	    nedge_pts += hclusters[i].nlines;
	for (i = 0; i < nvclusters; i++)
	    nedge_pts += vclusters[i].nlines;
	nedge_pts <<= 1;
	edge_pts  = (qr_finder_edge_pt *)malloc(nedge_pts * sizeof(*edge_pts));
	centers	  = (qr_finder_center *)malloc(QR_MINI(nhclusters, nvclusters) *
					       sizeof(*centers));
	ncenters  = qr_finder_find_crossings(centers, edge_pts, hclusters,
					     nhclusters, vclusters, nvclusters);
	*_centers = centers;
	*_edge_pts = edge_pts;
    } else
	ncenters = 0;
    free(vclusters);
    free(vneighbors);
    free(hclusters);
    free(hneighbors);
    return ncenters;
}

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

static void qr_hom_init(qr_hom *_hom, int _x0, int _y0, int _x1, int _y1,
			int _x2, int _y2, int _x3, int _y3, int _res)
{
    int dx10;
    int dx20;
    int dx30;
    int dx31;
    int dx32;
    int dy10;
    int dy20;
    int dy30;
    int dy31;
    int dy32;
    int a20;
    int a21;
    int a22;
    int b0;
    int b1;
    int b2;
    int s1;
    int s2;
    int r1;
    int r2;
    dx10 = _x1 - _x0;
    dx20 = _x2 - _x0;
    dx30 = _x3 - _x0;
    dx31 = _x3 - _x1;
    dx32 = _x3 - _x2;
    dy10 = _y1 - _y0;
    dy20 = _y2 - _y0;
    dy30 = _y3 - _y0;
    dy31 = _y3 - _y1;
    dy32 = _y3 - _y2;
    a20	 = dx32 * dy10 - dx10 * dy32;
    a21	 = dx20 * dy31 - dx31 * dy20;
    a22	 = dx32 * dy31 - dx31 * dy32;
    /*Figure out if we need to downscale anything.*/
    b0 = qr_ilog(QR_MAXI(abs(dx10), abs(dy10))) + qr_ilog(abs(a20 + a22));
    b1 = qr_ilog(QR_MAXI(abs(dx20), abs(dy20))) + qr_ilog(abs(a21 + a22));
    b2 = qr_ilog(QR_MAXI(QR_MAXI(abs(a20), abs(a21)), abs(a22)));
    s1 = QR_MAXI(0, _res + QR_MAXI(QR_MAXI(b0, b1), b2) - (QR_INT_BITS - 2));
    r1 = (1 << s1) >> 1;
    /*Compute the final coefficients of the forward transform.
  The 32x32->64 bit multiplies are really needed for accuracy with large
   versions.*/
    _hom->fwd[0][0] = QR_FIXMUL(dx10, a20 + a22, r1, s1);
    _hom->fwd[0][1] = QR_FIXMUL(dx20, a21 + a22, r1, s1);
    _hom->x0	    = _x0;
    _hom->fwd[1][0] = QR_FIXMUL(dy10, a20 + a22, r1, s1);
    _hom->fwd[1][1] = QR_FIXMUL(dy20, a21 + a22, r1, s1);
    _hom->y0	    = _y0;
    _hom->fwd[2][0] = (a20 + r1) >> s1;
    _hom->fwd[2][1] = (a21 + r1) >> s1;
    _hom->fwd22	    = s1 > _res ? (a22 + (r1 >> _res)) >> (s1 - _res) :
				  a22 << (_res - s1);
    /*Now compute the inverse transform.*/
    b0 = qr_ilog(QR_MAXI(QR_MAXI(abs(dx10), abs(dx20)), abs(dx30))) +
	 qr_ilog(QR_MAXI(abs(_hom->fwd[0][0]), abs(_hom->fwd[1][0])));
    b1 = qr_ilog(QR_MAXI(QR_MAXI(abs(dy10), abs(dy20)), abs(dy30))) +
	 qr_ilog(QR_MAXI(abs(_hom->fwd[0][1]), abs(_hom->fwd[1][1])));
    b2 = qr_ilog(abs(a22)) - s1;
    s2 = QR_MAXI(0, QR_MAXI(b0, b1) + b2 - (QR_INT_BITS - 3));
    r2 = (1 << s2) >> 1;
    s1 += s2;
    r1 <<= s2;
    /*The 32x32->64 bit multiplies are really needed for accuracy with large
   versions.*/
    _hom->inv[0][0] = QR_FIXMUL(_hom->fwd[1][1], a22, r1, s1);
    _hom->inv[0][1] = QR_FIXMUL(-_hom->fwd[0][1], a22, r1, s1);
    _hom->inv[1][0] = QR_FIXMUL(-_hom->fwd[1][0], a22, r1, s1);
    _hom->inv[1][1] = QR_FIXMUL(_hom->fwd[0][0], a22, r1, s1);
    _hom->inv[2][0] =
	QR_FIXMUL(_hom->fwd[1][0], _hom->fwd[2][1],
		  -QR_EXTMUL(_hom->fwd[1][1], _hom->fwd[2][0], r2), s2);
    _hom->inv[2][1] =
	QR_FIXMUL(_hom->fwd[0][1], _hom->fwd[2][0],
		  -QR_EXTMUL(_hom->fwd[0][0], _hom->fwd[2][1], r2), s2);
    _hom->inv22 = QR_FIXMUL(_hom->fwd[0][0], _hom->fwd[1][1],
			    -QR_EXTMUL(_hom->fwd[0][1], _hom->fwd[1][0], r2),
			    s2);
    _hom->res	= _res;
}

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
  _height: The distance between UL and DL in the square domain.*/
static int qr_finder_estimate_module_size_and_version(qr_finder *_f, int _width,
						      int _height)
{
    qr_point offs;
    int sums[4];
    int nsums[4];
    int usize;
    int nusize;
    int vsize;
    int nvsize;
    int uversion;
    int vversion;
    int e;
    offs[0] = offs[1] = 0;
    for (e = 0; e < 4; e++)
	if (_f->nedge_pts[e] > 0) {
	    qr_finder_edge_pt *edge_pts;
	    int sum;
	    int mean;
	    int n;
	    int i;
	    /*Average the samples for this edge, dropping the top and bottom 25%.*/
	    edge_pts = _f->edge_pts[e];
	    n	     = _f->nedge_pts[e];
	    sum	     = 0;
	    for (i = (n >> 2); i < n - (n >> 2); i++)
		sum += edge_pts[i].extent;
	    n	 = n - ((n >> 2) << 1);
	    mean = QR_DIVROUND(sum, n);
	    offs[e >> 1] += mean;
	    sums[e]  = sum;
	    nsums[e] = n;
	} else
	    nsums[e] = sums[e] = 0;
    /*If we have samples on both sides of an axis, refine our idea of where the
   unprojected finder center is located.*/
    if (_f->nedge_pts[0] > 0 && _f->nedge_pts[1] > 0) {
	_f->o[0] -= offs[0] >> 1;
	sums[0] -= offs[0] * nsums[0] >> 1;
	sums[1] -= offs[0] * nsums[1] >> 1;
    }
    if (_f->nedge_pts[2] > 0 && _f->nedge_pts[3] > 0) {
	_f->o[1] -= offs[1] >> 1;
	sums[2] -= offs[1] * nsums[2] >> 1;
	sums[3] -= offs[1] * nsums[3] >> 1;
    }
    /*We must have _some_ samples along each axis... if we don't, our transform
   must be pretty severely distorting the original square (e.g., with
   coordinates so large as to cause overflow).*/
    nusize = nsums[0] + nsums[1];
    if (nusize <= 0)
	return -1;
    /*The module size is 1/3 the average edge extent.*/
    nusize *= 3;
    usize = sums[1] - sums[0];
    usize = ((usize << 1) + nusize) / (nusize << 1);
    if (usize <= 0)
	return -1;
    /*Now estimate the version directly from the module size and the distance
   between the finder patterns.
  This is done independently using the extents along each axis.
  If either falls significantly outside the valid range (1 to 40), reject the
   configuration.*/
    uversion = (_width - 8 * usize) / (usize << 2);
    if (uversion < 1 || uversion > 40 + QR_LARGE_VERSION_SLACK)
	return -1;
    /*Now do the same for the other axis.*/
    nvsize = nsums[2] + nsums[3];
    if (nvsize <= 0)
	return -1;
    nvsize *= 3;
    vsize = sums[3] - sums[2];
    vsize = ((vsize << 1) + nvsize) / (nvsize << 1);
    if (vsize <= 0)
	return -1;
    vversion = (_height - 8 * vsize) / (vsize << 2);
    if (vversion < 1 || vversion > 40 + QR_LARGE_VERSION_SLACK)
	return -1;
    /*If the estimated version using extents along one axis is significantly
   different than the estimated version along the other axis, then the axes
   have significantly different scalings (relative to the grid).
  This can happen, e.g., when we have multiple adjacent QR codes, and we've
   picked two finder patterns from one and the third finder pattern from
   another, e.g.:
    X---DL UL---X
    |....   |....
    X....  UR....
  Such a configuration might even pass any other geometric checks if we
   didn't reject it here.*/
    if (abs(uversion - vversion) > QR_LARGE_VERSION_SLACK)
	return -1;
    _f->size[0] = usize;
    _f->size[1] = vsize;
    /*We intentionally do not compute an average version from the sizes along
   both axes.
  In the presence of projective distortion, one of them will be much more
   accurate than the other.*/
    _f->eversion[0] = uversion;
    _f->eversion[1] = vversion;
    return 0;
}

/*Eliminate outliers from the classified edge points with RANSAC.*/
static void qr_finder_ransac(qr_finder *_f, const qr_aff *_hom,
			     isaac_ctx *_isaac, int _e)
{
    qr_finder_edge_pt *edge_pts;
    int best_ninliers;
    int n;
    edge_pts	  = _f->edge_pts[_e];
    n		  = _f->nedge_pts[_e];
    best_ninliers = 0;
    if (n > 1) {
	int max_iters;
	int i;
	int j;
	/*17 iterations is enough to guarantee an outlier-free sample with more
   than 99% probability given as many as 50% outliers.*/
	max_iters = 17;
	for (i = 0; i < max_iters; i++) {
	    qr_point q0;
	    qr_point q1;
	    int ninliers;
	    int thresh;
	    int p0i;
	    int p1i;
	    int *p0;
	    int *p1;
	    int j;
	    /*Pick two random points on this edge.*/
	    p0i = isaac_next_uint(_isaac, n);
	    p1i = isaac_next_uint(_isaac, n - 1);
	    if (p1i >= p0i)
		p1i++;
	    p0 = edge_pts[p0i].pos;
	    p1 = edge_pts[p1i].pos;
	    /*If the corresponding line is not within 45 degrees of the proper
   orientation in the square domain, reject it outright.
  This can happen, e.g., when highly skewed orientations cause points to
   be misclassified into the wrong edge.
  The irony is that using such points might produce a line which _does_
   pass the corresponding validity checks.*/
	    qr_aff_unproject(&q0, _hom, p0[0], p0[1]);
	    qr_aff_unproject(&q1, _hom, p1[0], p1[1]);
	    qr_point_translate(q0, -_f->o[0], -_f->o[1]);
	    qr_point_translate(q1, -_f->o[0], -_f->o[1]);
	    if (abs(q0[_e >> 1] - q1[_e >> 1]) >
		abs(q0[1 - (_e >> 1)] - q1[1 - (_e >> 1)]))
		continue;
	    /*Identify the other edge points which are inliers.
  The squared distance should be distributed as a \Chi^2 distribution
   with one degree of freedom, which means for a 95% confidence the
   point should lie within a factor 3.8414588 ~= 4 times the expected
   variance of the point locations.
  We grossly approximate the standard deviation as 1 pixel in one
   direction, and 0.5 pixels in the other (because we average two
   coordinates).*/
	    thresh   = qr_isqrt(qr_point_distance2(p0, p1)
				<< (2 * QR_FINDER_SUBPREC + 1));
	    ninliers = 0;
	    for (j = 0; j < n; j++) {
		if (abs(qr_point_ccw(p0, p1, edge_pts[j].pos)) <= thresh) {
		    edge_pts[j].extent |= 1;
		    ninliers++;
		} else
		    edge_pts[j].extent &= ~1;
	    }
	    if (ninliers > best_ninliers) {
		for (j = 0; j < n; j++)
		    edge_pts[j].extent <<= 1;
		best_ninliers = ninliers;
		/*The actual number of iterations required is
    log(1-\alpha)/log(1-r*r),
   where \alpha is the required probability of taking a sample with
    no outliers (e.g., 0.99) and r is the estimated ratio of inliers
    (e.g. ninliers/n).
  This is just a rough (but conservative) approximation, but it
   should be good enough to stop the iteration early when we find
   a good set of inliers.*/
		if (ninliers > n >> 1)
		    max_iters = (67 * n - 63 * ninliers - 1) / (n << 1);
	    }
	}
	/*Now collect all the inliers at the beginning of the list.*/
	for (i = j = 0; j < best_ninliers; i++)
	    if (edge_pts[i].extent & 2) {
		if (j < i) {
		    qr_finder_edge_pt tmp;
		    *&tmp	    = *(edge_pts + i);
		    *(edge_pts + j) = *(edge_pts + i);
		    *(edge_pts + i) = *&tmp;
		}
		j++;
	    }
    }
    _f->ninliers[_e] = best_ninliers;
}

/*Perform a least-squares line fit to an edge of a finder pattern using the
   inliers found by RANSAC.*/
static int qr_line_fit_finder_edge(qr_line *_l, const qr_finder *_f, int _e,
				   int _res)
{
    qr_finder_edge_pt *edge_pts;
    qr_point *pts;
    int npts;
    int i;
    npts = _f->ninliers[_e];
    if (npts < 2)
	return -1;
    /*We could write a custom version of qr_line_fit_points that accesses
   edge_pts directly, but this saves on code size and doesn't measurably slow
   things down.*/
    pts	     = (qr_point *)malloc(npts * sizeof(*pts));
    edge_pts = _f->edge_pts[_e];
    for (i = 0; i < npts; i++) {
	pts[i][0] = edge_pts[i].pos[0];
	pts[i][1] = edge_pts[i].pos[1];
    }
    qr_line_fit_points(_l, pts, npts, _res);
    /*Make sure the center of the finder pattern lies in the positive halfspace
   of the line.*/
    qr_line_orient(_l, _f->c->pos[0], _f->c->pos[1]);
    free(pts);
    return 0;
}

/*Perform a least-squares line fit to a pair of common finder edges using the
   inliers found by RANSAC.
  Unlike a normal edge fit, we guarantee that this one succeeds by creating at
   least one point on each edge using the estimated module size if it has no
   inliers.*/
static void qr_line_fit_finder_pair(qr_line *_l, const qr_aff *_aff,
				    const qr_finder *_f0, const qr_finder *_f1,
				    int _e)
{
    qr_point *pts;
    int npts;
    qr_finder_edge_pt *edge_pts;
    qr_point q;
    int n0;
    int n1;
    int i;
    n0 = _f0->ninliers[_e];
    n1 = _f1->ninliers[_e];
    /*We could write a custom version of qr_line_fit_points that accesses
   edge_pts directly, but this saves on code size and doesn't measurably slow
   things down.*/
    npts = QR_MAXI(n0, 1) + QR_MAXI(n1, 1);
    pts	 = (qr_point *)malloc(npts * sizeof(*pts));
    if (n0 > 0) {
	edge_pts = _f0->edge_pts[_e];
	for (i = 0; i < n0; i++) {
	    pts[i][0] = edge_pts[i].pos[0];
	    pts[i][1] = edge_pts[i].pos[1];
	}
    } else {
	q[0] = _f0->o[0];
	q[1] = _f0->o[1];
	q[_e >> 1] += _f0->size[_e >> 1] * (2 * (_e & 1) - 1);
	qr_aff_project(&pts[0], _aff, q[0], q[1]);
	n0++;
    }
    if (n1 > 0) {
	edge_pts = _f1->edge_pts[_e];
	for (i = 0; i < n1; i++) {
	    pts[n0 + i][0] = edge_pts[i].pos[0];
	    pts[n0 + i][1] = edge_pts[i].pos[1];
	}
    } else {
	q[0] = _f1->o[0];
	q[1] = _f1->o[1];
	q[_e >> 1] += _f1->size[_e >> 1] * (2 * (_e & 1) - 1);
	qr_aff_project(&pts[n0], _aff, q[0], q[1]);
	n1++;
    }
    qr_line_fit_points(_l, pts, npts, _aff->res);
    /*Make sure at least one finder center lies in the positive halfspace.*/
    qr_line_orient(_l, _f0->c->pos[0], _f0->c->pos[1]);
    free(pts);
}

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
    qr_point p;
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
    int shift;
    int round;
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
    if (nr > 1)
	qr_line_fit_points(&l[1], r, nr, _aff->res);
    else {
	qr_aff_project(&p, _aff, _ur->o[0] + 3 * _ur->size[0], _ur->o[1]);
	shift	= QR_MAXI(0, qr_ilog(QR_MAXI(abs(_aff->fwd[0][1]),
					     abs(_aff->fwd[1][1]))) -
				 ((_aff->res + 1) >> 1));
	round	= (1 << shift) >> 1;
	l[1][0] = (_aff->fwd[1][1] + round) >> shift;
	l[1][1] = (-_aff->fwd[0][1] + round) >> shift;
	l[1][2] = -(l[1][0] * p[0] + l[1][1] * p[1]);
    }
    free(r);
    if (nb > 1)
	qr_line_fit_points(&l[3], b, nb, _aff->res);
    else {
	qr_aff_project(&p, _aff, _dl->o[0], _dl->o[1] + 3 * _dl->size[1]);
	shift	= QR_MAXI(0, qr_ilog(QR_MAXI(abs(_aff->fwd[0][1]),
					     abs(_aff->fwd[1][1]))) -
				 ((_aff->res + 1) >> 1));
	round	= (1 << shift) >> 1;
	l[3][0] = (_aff->fwd[1][0] + round) >> shift;
	l[3][1] = (-_aff->fwd[0][0] + round) >> shift;
	l[3][2] = -(l[1][0] * p[0] + l[1][1] * p[1]);
    }
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
	    long long w;
	    long long mask;
	    int c21;
	    int dx21;
	    int dy21;
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
	    c21	 = _p[2][0] * _p[1][1] - _p[2][1] * _p[1][0];
	    dx21 = _p[2][0] - _p[1][0];
	    dy21 = _p[2][1] - _p[1][1];
	    w	 = QR_EXTMUL(dim - 7, c21,
			     QR_EXTMUL(dim - 13, _p[0][0] * dy21 - _p[0][1] * dx21,
				       QR_EXTMUL(6, p3[0] * dy21 - p3[1] * dx21,
						 0)));
	    /*The projection failed: invalid geometry.*/
	    if (w == 0)
		return -1;
	    mask = QR_SIGNMASK(w);
	    w	 = (w + mask) ^ mask;
	    brx	 = (int)QR_DIVROUND(
		 (QR_EXTMUL((dim - 7) * _p[0][0], p3[0] * dy21,
			    QR_EXTMUL((dim - 13) * p3[0], c21 - _p[0][1] * dx21,
				      QR_EXTMUL(6 * _p[0][0], c21 - p3[1] * dx21,
						0))) +
		  mask) ^
		     mask,
		 w);
	    bry = (int)QR_DIVROUND(
		(QR_EXTMUL((dim - 7) * _p[0][1], -p3[1] * dx21,
			   QR_EXTMUL((dim - 13) * p3[1], c21 + _p[0][0] * dy21,
				     QR_EXTMUL(6 * _p[0][1], c21 + p3[0] * dy21,
					       0))) +
		 mask) ^
		    mask,
		w);
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

/*Reads the version bits near a finder module and decodes the version number.*/
static int qr_finder_version_decode(qr_finder *_f, const qr_hom *_hom,
				    const unsigned char *_img, int _width,
				    int _height, int _dir)
{
    qr_point q;
    unsigned v;
    int x0;
    int y0;
    int w0;
    int dxi;
    int dyi;
    int dwi;
    int dxj;
    int dyj;
    int dwj;
    int ret;
    int i;
    int j;
    int k;
    v		= 0;
    q[_dir]	= _f->o[_dir] - 7 * _f->size[_dir];
    q[1 - _dir] = _f->o[1 - _dir] - 3 * _f->size[1 - _dir];
    x0		= _hom->fwd[0][0] * q[0] + _hom->fwd[0][1] * q[1];
    y0		= _hom->fwd[1][0] * q[0] + _hom->fwd[1][1] * q[1];
    w0		= _hom->fwd[2][0] * q[0] + _hom->fwd[2][1] * q[1] + _hom->fwd22;
    dxi		= _hom->fwd[0][1 - _dir] * _f->size[1 - _dir];
    dyi		= _hom->fwd[1][1 - _dir] * _f->size[1 - _dir];
    dwi		= _hom->fwd[2][1 - _dir] * _f->size[1 - _dir];
    dxj		= _hom->fwd[0][_dir] * _f->size[_dir];
    dyj		= _hom->fwd[1][_dir] * _f->size[_dir];
    dwj		= _hom->fwd[2][_dir] * _f->size[_dir];
    for (k = i = 0; i < 6; i++) {
	int x;
	int y;
	int w;
	x = x0;
	y = y0;
	w = w0;
	for (j = 0; j < 3; j++, k++) {
	    qr_point p;
	    qr_hom_fproject(&p, _hom, x, y, w);
	    v |= qr_img_get_bit(_img, _width, _height, p[0], p[1]) << k;
	    x += dxj;
	    y += dyj;
	    w += dwj;
	}
	x0 += dxi;
	y0 += dyi;
	w0 += dwi;
    }
    ret = bch18_6_correct(&v);
    /*TODO: I seem to have an image with the version bits in a different order
   (the transpose of the standard order).
  Even if I change the order here so I can parse the version on this image,
   I can't decode the rest of the code.
  If this is really needed, we should just re-order the bits.*/
    return ret >= 0 ? (int)(v >> 12) : ret;
}

/*Reads the format info bits near the finder modules and decodes them.*/
static int qr_finder_fmt_info_decode(qr_finder *_ul, qr_finder *_ur,
				     qr_finder *_dl, const qr_hom *_hom,
				     const unsigned char *_img, int _width,
				     int _height)
{
    qr_point p;
    unsigned lo[2];
    unsigned hi[2];
    int u;
    int v;
    int x;
    int y;
    int w;
    int dx;
    int dy;
    int dw;
    int fmt_info[4];
    int count[4];
    int nerrs[4];
    int nfmt_info;
    int besti;
    int imax;
    int di;
    int i;
    int k;
    /*Read the bits around the UL corner.*/
    lo[0] = 0;
    u	  = _ul->o[0] + 5 * _ul->size[0];
    v	  = _ul->o[1] - 3 * _ul->size[1];
    x	  = _hom->fwd[0][0] * u + _hom->fwd[0][1] * v;
    y	  = _hom->fwd[1][0] * u + _hom->fwd[1][1] * v;
    w	  = _hom->fwd[2][0] * u + _hom->fwd[2][1] * v + _hom->fwd22;
    dx	  = _hom->fwd[0][1] * _ul->size[1];
    dy	  = _hom->fwd[1][1] * _ul->size[1];
    dw	  = _hom->fwd[2][1] * _ul->size[1];
    for (k = i = 0;; i++) {
	/*Skip the timing pattern row.*/
	if (i != 6) {
	    qr_hom_fproject(&p, _hom, x, y, w);
	    lo[0] |= qr_img_get_bit(_img, _width, _height, p[0], p[1]) << k++;
	    /*Don't advance q in the last iteration... we'll start the next loop from
   the current position.*/
	    if (i >= 8)
		break;
	}
	x += dx;
	y += dy;
	w += dw;
    }
    hi[0] = 0;
    dx	  = -_hom->fwd[0][0] * _ul->size[0];
    dy	  = -_hom->fwd[1][0] * _ul->size[0];
    dw	  = -_hom->fwd[2][0] * _ul->size[0];
    while (i-- > 0) {
	x += dx;
	y += dy;
	w += dw;
	/*Skip the timing pattern column.*/
	if (i != 6) {
	    qr_hom_fproject(&p, _hom, x, y, w);
	    hi[0] |= qr_img_get_bit(_img, _width, _height, p[0], p[1]) << k++;
	}
    }
    /*Read the bits next to the UR corner.*/
    lo[1] = 0;
    u	  = _ur->o[0] + 3 * _ur->size[0];
    v	  = _ur->o[1] + 5 * _ur->size[1];
    x	  = _hom->fwd[0][0] * u + _hom->fwd[0][1] * v;
    y	  = _hom->fwd[1][0] * u + _hom->fwd[1][1] * v;
    w	  = _hom->fwd[2][0] * u + _hom->fwd[2][1] * v + _hom->fwd22;
    dx	  = -_hom->fwd[0][0] * _ur->size[0];
    dy	  = -_hom->fwd[1][0] * _ur->size[0];
    dw	  = -_hom->fwd[2][0] * _ur->size[0];
    for (k = 0; k < 8; k++) {
	qr_hom_fproject(&p, _hom, x, y, w);
	lo[1] |= qr_img_get_bit(_img, _width, _height, p[0], p[1]) << k;
	x += dx;
	y += dy;
	w += dw;
    }
    /*Read the bits next to the DL corner.*/
    hi[1] = 0;
    u	  = _dl->o[0] + 5 * _dl->size[0];
    v	  = _dl->o[1] - 3 * _dl->size[1];
    x	  = _hom->fwd[0][0] * u + _hom->fwd[0][1] * v;
    y	  = _hom->fwd[1][0] * u + _hom->fwd[1][1] * v;
    w	  = _hom->fwd[2][0] * u + _hom->fwd[2][1] * v + _hom->fwd22;
    dx	  = _hom->fwd[0][1] * _dl->size[1];
    dy	  = _hom->fwd[1][1] * _dl->size[1];
    dw	  = _hom->fwd[2][1] * _dl->size[1];
    for (k = 8; k < 15; k++) {
	qr_hom_fproject(&p, _hom, x, y, w);
	hi[1] |= qr_img_get_bit(_img, _width, _height, p[0], p[1]) << k;
	x += dx;
	y += dy;
	w += dw;
    }
    /*For each group of bits we have two samples... try them in all combinations
   and pick the most popular valid code, breaking ties using the number of
   bit errors.*/
    imax      = 2 << (hi[0] != hi[1]);
    di	      = 1 + (lo[0] == lo[1]);
    nfmt_info = 0;
    for (i = 0; i < imax; i += di) {
	unsigned v;
	int ret;
	int j;
	v   = (lo[i & 1] | hi[i >> 1]) ^ 0x5412;
	ret = bch15_5_correct(&v);
	v >>= 10;
	if (ret < 0)
	    ret = 4;
	for (j = 0;; j++) {
	    if (j >= nfmt_info) {
		fmt_info[j] = v;
		count[j]    = 1;
		nerrs[j]    = ret;
		nfmt_info++;
		break;
	    }
	    if (fmt_info[j] == (int)v) {
		count[j]++;
		if (ret < nerrs[j])
		    nerrs[j] = ret;
		break;
	    }
	}
    }
    besti = 0;
    for (i = 1; i < nfmt_info; i++) {
	if ((nerrs[besti] > 3 && nerrs[i] <= 3) || count[i] > count[besti] ||
	    (count[i] == count[besti] && nerrs[i] < nerrs[besti])) {
	    besti = i;
	}
    }
    return nerrs[besti] < 4 ? fmt_info[besti] : -1;
}

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
static void qr_sampling_grid_init(qr_sampling_grid *_grid, int _version,
				  const qr_point _ul_pos,
				  const qr_point _ur_pos,
				  const qr_point _dl_pos, qr_point _p[4],
				  const unsigned char *_img, int _width,
				  int _height)
{
    qr_hom_cell base_cell;
    int align_pos[7];
    int dim;
    int nalign;
    int i;
    dim	   = 17 + (_version << 2);
    nalign = (_version / 7) + 2;
    /*Create a base cell to bootstrap the alignment pattern search.*/
    qr_hom_cell_init(&base_cell, 0, 0, dim - 1, 0, 0, dim - 1, dim - 1, dim - 1,
		     _p[0][0], _p[0][1], _p[1][0], _p[1][1], _p[2][0], _p[2][1],
		     _p[3][0], _p[3][1]);
    /*Allocate the array of cells.*/
    _grid->ncells   = nalign - 1;
    _grid->cells[0] = (qr_hom_cell *)malloc((nalign - 1) * (nalign - 1) *
					    sizeof(*_grid->cells[0]));
    for (i = 1; i < _grid->ncells; i++)
	_grid->cells[i] = _grid->cells[i - 1] + _grid->ncells;
    /*Initialize the function pattern mask.*/
    _grid->fpmask =
	(unsigned *)calloc(dim, ((dim + QR_INT_BITS - 1) >> QR_INT_LOGBITS) *
				    sizeof(*_grid->fpmask));
    /*Mask out the finder patterns (and separators and format info bits).*/
    qr_sampling_grid_fp_mask_rect(_grid, dim, 0, 0, 9, 9);
    qr_sampling_grid_fp_mask_rect(_grid, dim, 0, dim - 8, 9, 8);
    qr_sampling_grid_fp_mask_rect(_grid, dim, dim - 8, 0, 8, 9);
    /*Mask out the version number bits.*/
    if (_version > 6) {
	qr_sampling_grid_fp_mask_rect(_grid, dim, 0, dim - 11, 6, 3);
	qr_sampling_grid_fp_mask_rect(_grid, dim, dim - 11, 0, 3, 6);
    }
    /*Mask out the timing patterns.*/
    qr_sampling_grid_fp_mask_rect(_grid, dim, 9, 6, dim - 17, 1);
    qr_sampling_grid_fp_mask_rect(_grid, dim, 6, 9, 1, dim - 17);
    /*If we have no alignment patterns (e.g., this is a version 1 code), just use
   the base cell and hope it's good enough.*/
    if (_version < 2)
	memcpy(_grid->cells[0], &base_cell, sizeof(base_cell));
    else {
	qr_point *q;
	qr_point *p;
	int j;
	int k;
	q = (qr_point *)malloc(nalign * nalign * sizeof(*q));
	p = (qr_point *)malloc(nalign * nalign * sizeof(*p));
	/*Initialize the alignment pattern position list.*/
	align_pos[0]	      = 6;
	align_pos[nalign - 1] = dim - 7;
	if (_version > 6) {
	    int d;
	    d = QR_ALIGNMENT_SPACING[_version - 7];
	    for (i = nalign - 1; i-- > 1;)
		align_pos[i] = align_pos[i + 1] - d;
	}
	/*Three of the corners use a finder pattern instead of a separate
   alignment pattern.*/
	q[0][0]			    = 3;
	q[0][1]			    = 3;
	p[0][0]			    = _ul_pos[0];
	p[0][1]			    = _ul_pos[1];
	q[nalign - 1][0]	    = dim - 4;
	q[nalign - 1][1]	    = 3;
	p[nalign - 1][0]	    = _ur_pos[0];
	p[nalign - 1][1]	    = _ur_pos[1];
	q[(nalign - 1) * nalign][0] = 3;
	q[(nalign - 1) * nalign][1] = dim - 4;
	p[(nalign - 1) * nalign][0] = _dl_pos[0];
	p[(nalign - 1) * nalign][1] = _dl_pos[1];
	/*Scan for alignment patterns using a diagonal sweep.*/
	for (k = 1; k < 2 * nalign - 1; k++) {
	    int jmin;
	    int jmax;
	    jmax = QR_MINI(k, nalign - 1) - (k == nalign - 1);
	    jmin = QR_MAXI(0, k - (nalign - 1)) + (k == nalign - 1);
	    for (j = jmin; j <= jmax; j++) {
		qr_hom_cell *cell;
		int u;
		int v;
		int k;
		i	= jmax - (j - jmin);
		k	= i * nalign + j;
		u	= align_pos[j];
		v	= align_pos[i];
		q[k][0] = u;
		q[k][1] = v;
		/*Mask out the alignment pattern.*/
		qr_sampling_grid_fp_mask_rect(_grid, dim, u - 2, v - 2, 5, 5);
		/*Pick a cell to use to govern the alignment pattern search.*/
		if (i > 1 && j > 1) {
		    qr_point p0;
		    qr_point p1;
		    qr_point p2;
		    /*Each predictor is basically a straight-line extrapolation from two
   neighboring alignment patterns (except possibly near the opposing
   finder patterns).*/
		    qr_hom_cell_project(&p0, _grid->cells[i - 2] + j - 1, u, v,
					0);
		    qr_hom_cell_project(&p1, _grid->cells[i - 2] + j - 2, u, v,
					0);
		    qr_hom_cell_project(&p2, _grid->cells[i - 1] + j - 2, u, v,
					0);
		    /*Take the median of the predictions as the search center.*/
		    QR_SORT2I(p0[0], p1[0]);
		    QR_SORT2I(p0[1], p1[1]);
		    QR_SORT2I(p1[0], p2[0]);
		    QR_SORT2I(p1[1], p2[1]);
		    QR_SORT2I(p0[0], p1[0]);
		    QR_SORT2I(p0[1], p1[1]);
		    /*We need a cell that has the target point at a known (u,v) location.
  Since our cells don't have inverses, just construct one from our
   neighboring points.*/
		    cell = _grid->cells[i - 1] + j - 1;
		    qr_hom_cell_init(cell, q[k - nalign - 1][0],
				     q[k - nalign - 1][1], q[k - nalign][0],
				     q[k - nalign][1], q[k - 1][0], q[k - 1][1],
				     q[k][0], q[k][1], p[k - nalign - 1][0],
				     p[k - nalign - 1][1], p[k - nalign][0],
				     p[k - nalign][1], p[k - 1][0], p[k - 1][1],
				     p1[0], p1[1]);
		} else if (i > 1 && j > 0)
		    cell = _grid->cells[i - 2] + j - 1;
		else if (i > 0 && j > 1)
		    cell = _grid->cells[i - 1] + j - 2;
		else
		    cell = &base_cell;
		/*Use a very small search radius.
  A large displacement here usually means a false positive (e.g., when
   the real alignment pattern is damaged or missing), which can
   severely distort the projection.*/
		qr_alignment_pattern_search(p[k], cell, u, v, 2, _img, _width,
					    _height);
		if (i > 0 && j > 0) {
		    qr_hom_cell_init(_grid->cells[i - 1] + j - 1,
				     q[k - nalign - 1][0], q[k - nalign - 1][1],
				     q[k - nalign][0], q[k - nalign][1],
				     q[k - 1][0], q[k - 1][1], q[k][0], q[k][1],
				     p[k - nalign - 1][0], p[k - nalign - 1][1],
				     p[k - nalign][0], p[k - nalign][1],
				     p[k - 1][0], p[k - 1][1], p[k][0],
				     p[k][1]);
		}
	    }
	}
	free(q);
	free(p);
    }
    /*Set the limits over which each cell is used.*/
    memcpy(_grid->cell_limits, align_pos + 1,
	   (_grid->ncells - 1) * sizeof(*_grid->cell_limits));
    _grid->cell_limits[_grid->ncells - 1] = dim;
    /*Produce a bounding square for the code (to mark finder centers with).
  Because of non-linear distortion, this might not actually bound the code,
   but it should be good enough.
  I don't think it's worth computing a convex hull or anything silly like
   that.*/
    qr_hom_cell_project(&_p[0], _grid->cells[0] + 0, -1, -1, 1);
    qr_hom_cell_project(&_p[1], _grid->cells[0] + _grid->ncells - 1,
			(dim << 1) - 1, -1, 1);
    qr_hom_cell_project(&_p[2], _grid->cells[_grid->ncells - 1] + 0, -1,
			(dim << 1) - 1, 1);
    qr_hom_cell_project(&_p[3],
			_grid->cells[_grid->ncells - 1] + _grid->ncells - 1,
			(dim << 1) - 1, (dim << 1) - 1, 1);
    /*Clamp the points somewhere near the image (this is really just in case a
   corner is near the plane at infinity).*/
    for (i = 0; i < 4; i++) {
	_p[i][0] = QR_CLAMPI((-_width << QR_FINDER_SUBPREC), _p[i][0],
			     ((_width << QR_FINDER_SUBPREC) + 1));
	_p[i][1] = QR_CLAMPI((-_height << QR_FINDER_SUBPREC), _p[i][1],
			     ((_height << QR_FINDER_SUBPREC) + 1));
    }
    /*TODO: Make fine adjustments using the timing patterns.
  Possible strategy: scan the timing pattern at QR_ALIGN_SUBPREC (or finer)
   resolution, use dynamic programming to match midpoints between
   transitions to the ideal grid locations.*/
}

static void qr_sampling_grid_clear(qr_sampling_grid *_grid)
{
    free(_grid->fpmask);
    free(_grid->cells[0]);
}

/*Generate the data mask corresponding to the given mask pattern.*/
static void qr_data_mask_fill(unsigned *_mask, int _dim, int _pattern)
{
    int stride;
    int i;
    int j;
    stride = (_dim + QR_INT_BITS - 1) >> QR_INT_LOGBITS;
    /*Note that we store bits column-wise, since that's how they're read out of
   the grid.*/
    switch (_pattern) {
    /*10101010 i+j+1&1
    01010101
    10101010
    01010101*/
    case 0: {
	int m;
	m = 0x55;
	for (j = 0; j < _dim; j++) {
	    memset(_mask + j * stride, m, stride * sizeof(*_mask));
	    m ^= 0xFF;
	}
    } break;
    /*11111111 i+1&1
    00000000
    11111111
    00000000*/
    case 1:
	memset(_mask, 0x55, _dim * stride * sizeof(*_mask));
	break;
    /*10010010 (j+1)%3&1
    10010010
    10010010
    10010010*/
    case 2: {
	unsigned m;
	m = 0xFF;
	for (j = 0; j < _dim; j++) {
	    memset(_mask + j * stride, m & 0xFF, stride * sizeof(*_mask));
	    m = m << 8 | m >> 16;
	}
    } break;
    /*10010010 (i+j+1)%3&1
    00100100
    01001001
    10010010*/
    case 3: {
	unsigned mi;
	unsigned mj;
	mj = 0;
	for (i = 0; i < (QR_INT_BITS + 2) / 3; i++)
	    mj |= 1 << 3 * i;
	for (j = 0; j < _dim; j++) {
	    mi = mj;
	    for (i = 0; i < stride; i++) {
		_mask[j * stride + i] = mi;
		mi		      = (mi >> (QR_INT_BITS % 3)) |
		     (mi << (3 - (QR_INT_BITS % 3)));
	    }
	    mj = mj >> 1 | mj << 2;
	}
    } break;
    /*11100011 (i>>1)+(j/3)+1&1
    11100011
    00011100
    00011100*/
    case 4: {
	unsigned m;
	m = 7;
	for (j = 0; j < _dim; j++) {
	    memset(_mask + j * stride, (0xCC ^ -(m & 1)) & 0xFF,
		   stride * sizeof(*_mask));
	    m = m >> 1 | m << 5;
	}
    } break;
    /*11111111 !((i*j)%6)
    10000010
    10010010
    10101010*/
    case 5: {
	for (j = 0; j < _dim; j++) {
	    unsigned m;
	    m = 0;
	    for (i = 0; i < 6; i++)
		m |= (!((i * j) % 6)) << i;
	    for (i = 6; i < QR_INT_BITS; i <<= 1)
		m |= m << i;
	    for (i = 0; i < stride; i++) {
		_mask[j * stride + i] = m;
		m = (m >> (QR_INT_BITS % 6)) | (m << (6 - (QR_INT_BITS % 6)));
	    }
	}
    } break;
    /*11111111 (i*j)%3+i*j+1&1
    11100011
    11011011
    10101010*/
    case 6: {
	for (j = 0; j < _dim; j++) {
	    unsigned m;
	    m = 0;
	    for (i = 0; i < 6; i++)
		m |= ((((i * j) % 3 + i * j + 1) & 1)) << i;
	    for (i = 6; i < QR_INT_BITS; i <<= 1)
		m |= m << i;
	    for (i = 0; i < stride; i++) {
		_mask[j * stride + i] = m;
		m = (m >> (QR_INT_BITS % 6)) | (m << (6 - (QR_INT_BITS % 6)));
	    }
	}
    } break;
    /*10101010 (i*j)%3+i+j+1&1
    00011100
    10001110
    01010101*/
    default: {
	for (j = 0; j < _dim; j++) {
	    unsigned m;
	    m = 0;
	    for (i = 0; i < 6; i++)
		m |= ((((i * j) % 3 + i + j + 1) & 1)) << i;
	    for (i = 6; i < QR_INT_BITS; i <<= 1)
		m |= m << i;
	    for (i = 0; i < stride; i++) {
		_mask[j * stride + i] = m;
		m = (m >> (QR_INT_BITS % 6)) | (m << (6 - (QR_INT_BITS % 6)));
	    }
	}
    } break;
    }
}

static void qr_sampling_grid_sample(const qr_sampling_grid *_grid,
				    unsigned *_data_bits, int _dim,
				    int _fmt_info, const unsigned char *_img,
				    int _width, int _height)
{
    int stride;
    int u0;
    int u1;
    int j;
    /*We initialize the buffer with the data mask and XOR bits into it as we read
   them out of the image instead of unmasking in a separate step.*/
    qr_data_mask_fill(_data_bits, _dim, _fmt_info & 7);
    stride = (_dim + QR_INT_BITS - 1) >> QR_INT_LOGBITS;
    u0	   = 0;
    /*We read data cell-by-cell to avoid having to constantly change which
   projection we're using as we read each bit.
  This (and the position-dependent data mask) is the reason we buffer the
   bits we read instead of converting them directly to codewords here.
  Note that bits are stored column-wise, since that's how we'll scan them.*/
    for (j = 0; j < _grid->ncells; j++) {
	int i;
	int v0;
	int v1;
	u1 = _grid->cell_limits[j];
	v0 = 0;
	for (i = 0; i < _grid->ncells; i++) {
	    qr_hom_cell *cell;
	    int x0;
	    int y0;
	    int w0;
	    int u;
	    int du;
	    int dv;
	    v1	 = _grid->cell_limits[i];
	    cell = _grid->cells[i] + j;
	    du	 = u0 - cell->u0;
	    dv	 = v0 - cell->v0;
	    x0 = cell->fwd[0][0] * du + cell->fwd[0][1] * dv + cell->fwd[0][2];
	    y0 = cell->fwd[1][0] * du + cell->fwd[1][1] * dv + cell->fwd[1][2];
	    w0 = cell->fwd[2][0] * du + cell->fwd[2][1] * dv + cell->fwd[2][2];
	    for (u = u0; u < u1; u++) {
		int x;
		int y;
		int w;
		int v;
		x = x0;
		y = y0;
		w = w0;
		for (v = v0; v < v1; v++) {
		    /*Skip doing all the divisions and bounds checks if the bit is in the
   function pattern.*/
		    if (!qr_sampling_grid_is_in_fp(_grid, _dim, u, v)) {
			qr_point p;
			qr_hom_cell_fproject(&p, cell, x, y, w);
			_data_bits[u * stride + (v >> QR_INT_LOGBITS)] ^=
			    qr_img_get_bit(_img, _width, _height, p[0], p[1])
			    << (v & (QR_INT_BITS - 1));
		    }
		    x += cell->fwd[0][1];
		    y += cell->fwd[1][1];
		    w += cell->fwd[2][1];
		}
		x0 += cell->fwd[0][0];
		y0 += cell->fwd[1][0];
		w0 += cell->fwd[2][0];
	    }
	    v0 = v1;
	}
	u0 = u1;
    }
}

/*Arranges the sample bits read by qr_sampling_grid_sample() into bytes and
   groups those bytes into Reed-Solomon blocks.
  The individual block pointers are destroyed by this routine.*/
static void qr_samples_unpack(unsigned char **_blocks, int _nblocks,
			      int _nshort_data, int _nshort_blocks,
			      const unsigned *_data_bits,
			      const unsigned *_fp_mask, int _dim)
{
    unsigned bits;
    int biti;
    int stride;
    int blocki;
    int blockj;
    int i;
    int j;
    stride = (_dim + QR_INT_BITS - 1) >> QR_INT_LOGBITS;
    /*If _all_ the blocks are short, don't skip anything (see below).*/
    if (_nshort_blocks >= _nblocks)
	_nshort_blocks = 0;
    /*Scan columns in pairs from right to left.*/
    bits = 0;
    for (blocki = blockj = biti = 0, j = _dim - 1; j > 0; j -= 2) {
	unsigned data1;
	unsigned data2;
	unsigned fp_mask1;
	unsigned fp_mask2;
	int nbits;
	int l;
	/*Scan up a pair of columns.*/
	nbits = ((_dim - 1) & (QR_INT_BITS - 1)) + 1;
	l     = j * stride;
	for (i = stride; i-- > 0;) {
	    data1    = _data_bits[l + i];
	    fp_mask1 = _fp_mask[l + i];
	    data2    = _data_bits[l + i - stride];
	    fp_mask2 = _fp_mask[l + i - stride];
	    while (nbits-- > 0) {
		/*Pull a bit from the right column.*/
		if (!((fp_mask1 >> nbits) & 1)) {
		    bits = (bits << 1) | ((data1 >> nbits) & 1);
		    biti++;
		}
		/*Pull a bit from the left column.*/
		if (!((fp_mask2 >> nbits) & 1)) {
		    bits = (bits << 1) | ((data2 >> nbits) & 1);
		    biti++;
		}
		/*If we finished a byte, drop it in a block.*/
		if (biti >= 8) {
		    biti -= 8;
		    *_blocks[blocki++]++ = (unsigned char)(bits >> biti);
		    /*For whatever reason, the long blocks are at the _end_ of the list,
   instead of the beginning.
  Even worse, the extra bytes they get come at the end of the data
   bytes, before the parity bytes.
  Hence the logic here: when we've filled up the data portion of the
   short blocks, skip directly to the long blocks for the next byte.
  It's also the reason we increment _blocks[blocki] on each store,
   instead of just indexing with blockj (after this iteration the
   number of bytes in each block differs).*/
		    if (blocki >= _nblocks)
			blocki = ++blockj == _nshort_data ? _nshort_blocks : 0;
		}
	    }
	    nbits = QR_INT_BITS;
	}
	j -= 2;
	/*Skip the column with the vertical timing pattern.*/
	if (j == 6)
	    j--;
	/*Scan down a pair of columns.*/
	l = j * stride;
	for (i = 0; i < stride; i++) {
	    data1    = _data_bits[l + i];
	    fp_mask1 = _fp_mask[l + i];
	    data2    = _data_bits[l + i - stride];
	    fp_mask2 = _fp_mask[l + i - stride];
	    nbits    = QR_MINI(_dim - (i << QR_INT_LOGBITS), QR_INT_BITS);
	    while (nbits-- > 0) {
		/*Pull a bit from the right column.*/
		if (!(fp_mask1 & 1)) {
		    bits = (bits << 1) | (data1 & 1);
		    biti++;
		}
		data1 >>= 1;
		fp_mask1 >>= 1;
		/*Pull a bit from the left column.*/
		if (!(fp_mask2 & 1)) {
		    bits = (bits << 1) | (data2 & 1);
		    biti++;
		}
		data2 >>= 1;
		fp_mask2 >>= 1;
		/*If we finished a byte, drop it in a block.*/
		if (biti >= 8) {
		    biti -= 8;
		    *_blocks[blocki++]++ = (unsigned char)(bits >> biti);
		    /*See comments on the "up" loop for the reason behind this mess.*/
		    if (blocki >= _nblocks)
			blocki = ++blockj == _nshort_data ? _nshort_blocks : 0;
		}
	    }
	}
    }
}

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

/*The characters available in QR_MODE_ALNUM.*/
static const unsigned char QR_ALNUM_TABLE[45] = {
    '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E',
    'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
    'U', 'V', 'W', 'X', 'Y', 'Z', ' ', '$', '%', '*', '+', '-', '.', '/', ':'
};

static int qr_code_data_parse(qr_code_data *_qrdata, int _version,
			      const unsigned char *_data, int _ndata)
{
    qr_pack_buf qpb;
    unsigned self_parity;
    int centries;
    int len_bits_idx;
    /*Entries are stored directly in the struct during parsing.
  Caller cleans up any allocated data on failure.*/
    _qrdata->entries  = NULL;
    _qrdata->nentries = 0;
    _qrdata->sa_size  = 0;
    self_parity	      = 0;
    centries	      = 0;
    /*The versions are divided into 3 ranges that each use a different number of
   bits for length fields.*/
    len_bits_idx = (_version > 9) + (_version > 26);
    qr_pack_buf_init(&qpb, _data, _ndata);
    /*While we have enough bits to read a mode...*/
    while (qr_pack_buf_avail(&qpb) >= 4) {
	qr_code_data_entry *entry;
	int mode;
	mode = qr_pack_buf_read(&qpb, 4);
	/*Mode 0 is a terminator.*/
	if (!mode)
	    break;
	if (_qrdata->nentries >= centries) {
	    centries	     = centries << 1 | 1;
	    _qrdata->entries = (qr_code_data_entry *)realloc(
		_qrdata->entries, centries * sizeof(*_qrdata->entries));
	}
	entry	    = _qrdata->entries + _qrdata->nentries++;
	entry->mode = mode;
	/*Set the buffer to NULL, because if parsing fails, we might try to free it
   on clean-up.*/
	entry->payload.data.buf = NULL;
	switch (mode) {
	    /*The number of bits used to encode the character count for each version
   range and each data mode.*/
	    static const unsigned char LEN_BITS[3][4] = { { 10, 9, 8, 8 },
							  { 12, 11, 16, 10 },
							  { 14, 13, 16, 12 } };
	case QR_MODE_NUM: {
	    unsigned char *buf;
	    unsigned bits;
	    unsigned c;
	    int len;
	    int count;
	    int rem;
	    len = qr_pack_buf_read(&qpb, LEN_BITS[len_bits_idx][0]);
	    if (len < 0)
		return -1;
	    /*Check to see if there are enough bits left now, so we don't have to
     in the decode loop.*/
	    count = len / 3;
	    rem	  = len % 3;
	    if (qr_pack_buf_avail(&qpb) <
		10 * count + 7 * (rem >> 1 & 1) + 4 * (rem & 1))
		return -1;
	    entry->payload.data.buf = buf =
		(unsigned char *)malloc(len * sizeof(*buf));
	    entry->payload.data.len = len;
	    /*Read groups of 3 digits encoded in 10 bits.*/
	    while (count-- > 0) {
		bits = qr_pack_buf_read(&qpb, 10);
		if (bits >= 1000)
		    return -1;
		c = '0' + bits / 100;
		self_parity ^= c;
		*buf++ = (unsigned char)c;
		bits %= 100;
		c = '0' + bits / 10;
		self_parity ^= c;
		*buf++ = (unsigned char)c;
		c      = '0' + bits % 10;
		self_parity ^= c;
		*buf++ = (unsigned char)c;
	    }
	    /*Read the last two digits encoded in 7 bits.*/
	    if (rem > 1) {
		bits = qr_pack_buf_read(&qpb, 7);
		if (bits >= 100)
		    return -1;
		c = '0' + bits / 10;
		self_parity ^= c;
		*buf++ = (unsigned char)c;
		c      = '0' + bits % 10;
		self_parity ^= c;
		*buf++ = (unsigned char)c;
	    }
	    /*Or the last one digit encoded in 4 bits.*/
	    else if (rem) {
		bits = qr_pack_buf_read(&qpb, 4);
		if (bits >= 10)
		    return -1;
		c = '0' + bits;
		self_parity ^= c;
		*buf++ = (unsigned char)c;
	    }
	} break;
	case QR_MODE_ALNUM: {
	    unsigned char *buf;
	    unsigned bits;
	    unsigned c;
	    int len;
	    int count;
	    int rem;
	    len = qr_pack_buf_read(&qpb, LEN_BITS[len_bits_idx][1]);
	    if (len < 0)
		return -1;
	    /*Check to see if there are enough bits left now, so we don't have to
     in the decode loop.*/
	    count = len >> 1;
	    rem	  = len & 1;
	    if (qr_pack_buf_avail(&qpb) < 11 * count + 6 * rem)
		return -1;
	    entry->payload.data.buf = buf =
		(unsigned char *)malloc(len * sizeof(*buf));
	    entry->payload.data.len = len;
	    /*Read groups of two characters encoded in 11 bits.*/
	    while (count-- > 0) {
		bits = qr_pack_buf_read(&qpb, 11);
		if (bits >= 2025)
		    return -1;
		c = QR_ALNUM_TABLE[bits / 45];
		self_parity ^= c;
		*buf++ = (unsigned char)c;
		c      = QR_ALNUM_TABLE[bits % 45];
		self_parity ^= c;
		*buf++ = (unsigned char)c;
		len -= 2;
	    }
	    /*Read the last character encoded in 6 bits.*/
	    if (rem) {
		bits = qr_pack_buf_read(&qpb, 6);
		if (bits >= 45)
		    return -1;
		c = QR_ALNUM_TABLE[bits];
		self_parity ^= c;
		*buf++ = (unsigned char)c;
	    }
	} break;
	/*Structured-append header.*/
	case QR_MODE_STRUCT: {
	    int bits;
	    bits = qr_pack_buf_read(&qpb, 16);
	    if (bits < 0)
		return -1;
	    /*We save a copy of the data in _qrdata for easy reference when
     grouping structured-append codes.
    If for some reason the code has multiple S-A headers, first one wins,
     since it is supposed to come before everything else (TODO: should we
     return an error instead?).*/
	    if (_qrdata->sa_size == 0) {
		_qrdata->sa_index = entry->payload.sa.sa_index =
		    (unsigned char)(bits >> 12 & 0xF);
		_qrdata->sa_size = entry->payload.sa.sa_size =
		    (unsigned char)((bits >> 8 & 0xF) + 1);
		_qrdata->sa_parity = entry->payload.sa.sa_parity =
		    (unsigned char)(bits & 0xFF);
	    }
	} break;
	case QR_MODE_BYTE: {
	    unsigned char *buf;
	    unsigned c;
	    int len;
	    len = qr_pack_buf_read(&qpb, LEN_BITS[len_bits_idx][2]);
	    if (len < 0)
		return -1;
	    /*Check to see if there are enough bits left now, so we don't have to
     in the decode loop.*/
	    if (qr_pack_buf_avail(&qpb) < len << 3)
		return -1;
	    entry->payload.data.buf = buf =
		(unsigned char *)malloc(len * sizeof(*buf));
	    entry->payload.data.len = len;
	    while (len-- > 0) {
		c = qr_pack_buf_read(&qpb, 8);
		self_parity ^= c;
		*buf++ = (unsigned char)c;
	    }
	} break;
	/*FNC1 first position marker.*/
	case QR_MODE_FNC1_1ST:
	    break;
	/*Extended Channel Interpretation data.*/
	case QR_MODE_ECI: {
	    unsigned val;
	    int bits;
	    /*ECI uses a variable-width encoding similar to UTF-8*/
	    bits = qr_pack_buf_read(&qpb, 8);
	    if (bits < 0)
		return -1;
	    /*One byte:*/
	    if (!(bits & 0x80))
		val = bits;
	    /*Two bytes:*/
	    else if (!(bits & 0x40)) {
		val  = bits & 0x3F << 8;
		bits = qr_pack_buf_read(&qpb, 8);
		if (bits < 0)
		    return -1;
		val |= bits;
	    }
	    /*Three bytes:*/
	    else if (!(bits & 0x20)) {
		val  = bits & 0x1F << 16;
		bits = qr_pack_buf_read(&qpb, 16);
		if (bits < 0)
		    return -1;
		val |= bits;
		/*Valid ECI values are 0...999999.*/
		if (val >= 1000000)
		    return -1;
	    }
	    /*Invalid lead byte.*/
	    else
		return -1;
	    entry->payload.eci = val;
	} break;
	case QR_MODE_KANJI: {
	    unsigned char *buf;
	    unsigned bits;
	    int len;
	    len = qr_pack_buf_read(&qpb, LEN_BITS[len_bits_idx][3]);
	    if (len < 0)
		return -1;
	    /*Check to see if there are enough bits left now, so we don't have to
     in the decode loop.*/
	    if (qr_pack_buf_avail(&qpb) < 13 * len)
		return -1;
	    entry->payload.data.buf = buf =
		(unsigned char *)malloc(2 * len * sizeof(*buf));
	    entry->payload.data.len = 2 * len;
	    /*Decode 2-byte SJIS characters encoded in 13 bits.*/
	    while (len-- > 0) {
		bits = qr_pack_buf_read(&qpb, 13);
		bits = (bits / 0xC0 << 8 | bits % 0xC0) + 0x8140;
		if (bits >= 0xA000)
		    bits += 0x4000;
		/*TODO: Are values 0xXX7F, 0xXXFD...0xXXFF always invalid?
    Should we reject them here?*/
		self_parity ^= bits;
		*buf++ = (unsigned char)(bits >> 8);
		*buf++ = (unsigned char)(bits & 0xFF);
	    }
	} break;
	/*FNC1 second position marker.*/
	case QR_MODE_FNC1_2ND: {
	    int bits;
	    /*FNC1 in the 2nd position encodes an Application Indicator in one
     byte, which is either a letter (A...Z or a...z) or a 2-digit number.
    The letters are encoded with their ASCII value plus 100, the numbers
     are encoded directly with their numeric value.
    Values 100...164, 191...196, and 223...255 are invalid, so we reject
     them here.*/
	    bits = qr_pack_buf_read(&qpb, 8);
	    if (!((bits >= 0 && bits < 100) || (bits >= 165 && bits < 191) ||
		  (bits >= 197 && bits < 223))) {
		return -1;
	    }
	    entry->payload.ai = bits;
	} break;
	/*Unknown mode number:*/
	default: {
	    /*Unfortunately, because we have to understand the format of a mode to
     know how many bits it occupies, we can't skip unknown modes.
    Therefore we have to fail.*/
	    return -1;
	} break;
	}
    }
    /*Store the parity of the data from this code, for S-A.
  The final parity is the 8-bit XOR of all the decoded bytes of literal data.
  We don't combine the 2-byte kanji codes into one byte in the loops above,
   because we can just do it here instead.*/
    _qrdata->self_parity = ((self_parity >> 8) ^ self_parity) & 0xFF;
    /*Success.*/
    _qrdata->entries = (qr_code_data_entry *)realloc(
	_qrdata->entries, _qrdata->nentries * sizeof(*_qrdata->entries));
    return 0;
}

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

/*Bulk data for the number of parity bytes per Reed-Solomon block.*/
static const unsigned char QR_RS_NPAR_VALS[71] = {
    /*[ 0]*/ 7,	 10, 13, 17,
    /*[ 4]*/ 10, 16, 22, 28, 26, 26, 26, 22, 24, 22, 22, 26,
    24,		 18, 22,
    /*[19]*/ 15, 26, 18, 22, 24, 30, 24, 20, 24,
    /*[28]*/ 18, 16, 24, 28, 28, 28, 28, 30, 24,
    /*[37]*/ 20, 18, 18, 26, 24, 28, 24, 30, 26, 28, 28, 26,
    28,		 30, 30, 22, 20, 24,
    /*[55]*/ 20, 18, 26, 16,
    /*[59]*/ 20, 30, 28, 24, 22, 26, 28, 26, 30, 28, 30, 30
};

/*An offset into QR_RS_NPAR_DATA for each version that gives the number of
   parity bytes per Reed-Solomon block for each error correction level.*/
static const unsigned char QR_RS_NPAR_OFFS[40] = {
    0,	4,  19, 55, 15, 28, 37, 12, 51, 39, 59, 62, 10, 24,
    22, 41, 31, 44, 7,	65, 47, 33, 67, 67, 48, 32, 67, 67,
    67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67
};

/*The number of Reed-Solomon blocks for each version and error correction
   level.*/
static const unsigned char QR_RS_NBLOCKS[40][4] = {
    { 1, 1, 1, 1 },	{ 1, 1, 1, 1 },	    { 1, 1, 2, 2 },
    { 1, 2, 2, 4 },	{ 1, 2, 4, 4 },	    { 2, 4, 4, 4 },
    { 2, 4, 6, 5 },	{ 2, 4, 6, 6 },	    { 2, 5, 8, 8 },
    { 4, 5, 8, 8 },	{ 4, 5, 8, 11 },    { 4, 8, 10, 11 },
    { 4, 9, 12, 16 },	{ 4, 9, 16, 16 },   { 6, 10, 12, 18 },
    { 6, 10, 17, 16 },	{ 6, 11, 16, 19 },  { 6, 13, 18, 21 },
    { 7, 14, 21, 25 },	{ 8, 16, 20, 25 },  { 8, 17, 23, 25 },
    { 9, 17, 23, 34 },	{ 9, 18, 25, 30 },  { 10, 20, 27, 32 },
    { 12, 21, 29, 35 }, { 12, 23, 34, 37 }, { 12, 25, 34, 40 },
    { 13, 26, 35, 42 }, { 14, 28, 38, 45 }, { 15, 29, 40, 48 },
    { 16, 31, 43, 51 }, { 17, 33, 45, 54 }, { 18, 35, 48, 57 },
    { 19, 37, 51, 60 }, { 19, 38, 53, 63 }, { 20, 40, 56, 66 },
    { 21, 43, 59, 70 }, { 22, 45, 62, 74 }, { 24, 47, 65, 77 },
    { 25, 49, 68, 81 }
};

/*Attempts to fully decode a QR code.
  _qrdata:   Returns the parsed code data.
  _gf:       Used for Reed-Solomon error correction.
  _ul_pos:   The location of the UL finder pattern.
  _ur_pos:   The location of the UR finder pattern.
  _dl_pos:   The location of the DL finder pattern.
  _version:  The (decoded) version number.
  _fmt_info: The decoded format info.
  _img:      The binary input image.
  _width:    The width of the input image.
  _height:   The height of the input image.
  Return: 0 on success, or a negative value on error.*/
static int qr_code_decode(qr_code_data *_qrdata, const rs_gf256 *_gf,
			  const qr_point _ul_pos, const qr_point _ur_pos,
			  const qr_point _dl_pos, int _version, int _fmt_info,
			  const unsigned char *_img, int _width, int _height)
{
    qr_sampling_grid grid;
    unsigned *data_bits;
    unsigned char **blocks;
    unsigned char *block_data;
    int nblocks;
    int nshort_blocks;
    int ncodewords;
    int block_sz;
    int ecc_level;
    int ndata;
    int npar;
    int dim;
    int ret;
    int i;
    /*Read the bits out of the image.*/
    qr_sampling_grid_init(&grid, _version, _ul_pos, _ur_pos, _dl_pos,
			  _qrdata->bbox, _img, _width, _height);
#if defined(QR_DEBUG)
    qr_sampling_grid_dump(&grid, _version, _img, _width, _height);
#endif
    dim	      = 17 + (_version << 2);
    data_bits = (unsigned *)malloc(
	dim * ((dim + QR_INT_BITS - 1) >> QR_INT_LOGBITS) * sizeof(*data_bits));
    qr_sampling_grid_sample(&grid, data_bits, dim, _fmt_info, _img, _width,
			    _height);
    /*Group those bits into Reed-Solomon codewords.*/
    ecc_level  = (_fmt_info >> 3) ^ 1;
    nblocks    = QR_RS_NBLOCKS[_version - 1][ecc_level];
    npar       = *(QR_RS_NPAR_VALS + QR_RS_NPAR_OFFS[_version - 1] + ecc_level);
    ncodewords = qr_code_ncodewords(_version);
    block_sz   = ncodewords / nblocks;
    nshort_blocks = nblocks - (ncodewords % nblocks);
    blocks	  = (unsigned char **)malloc(nblocks * sizeof(*blocks));
    block_data	  = (unsigned char *)malloc(ncodewords * sizeof(*block_data));
    blocks[0]	  = block_data;
    for (i = 1; i < nblocks; i++)
	blocks[i] = blocks[i - 1] + block_sz + (i > nshort_blocks);
    qr_samples_unpack(blocks, nblocks, block_sz - npar, nshort_blocks,
		      data_bits, grid.fpmask, dim);
    qr_sampling_grid_clear(&grid);
    free(blocks);
    free(data_bits);
    /*Perform the error correction.*/
    ndata      = 0;
    ncodewords = 0;
    ret	       = 0;
    for (i = 0; i < nblocks; i++) {
	int block_szi;
	int ndatai;
	block_szi = block_sz + (i >= nshort_blocks);
	ret = rs_correct(_gf, QR_M0, block_data + ncodewords, block_szi, npar,
			 NULL, 0);
	/*For version 1 symbols and version 2-L and 3-L symbols, we aren't allowed
   to use all the parity bytes for correction.
  They are instead used to improve detection.
  Version 1-L reserves 3 parity bytes for detection.
  Versions 1-M and 2-L reserve 2 parity bytes for detection.
  Versions 1-Q, 1-H, and 3-L reserve 1 parity byte for detection.
  We can ignore the version 3-L restriction because it has an odd number of
   parity bytes, and we don't support erasure detection.*/
	if (ret < 0 || (_version == 1 && ret > ((ecc_level + 1) << 1)) ||
	    (_version == 2 && ecc_level == 0 && ret > 4)) {
	    ret = -1;
	    break;
	}
	ndatai = block_szi - npar;
	memmove(block_data + ndata, block_data + ncodewords,
		ndatai * sizeof(*block_data));
	ncodewords += block_szi;
	ndata += ndatai;
    }
    /*Parse the corrected bitstream.*/
    if (ret >= 0) {
	ret = qr_code_data_parse(_qrdata, _version, block_data, ndata);
	/*We could return any partially decoded data, but then we'd have to have
   API support for that; a mode ignoring ECC errors might also be useful.*/
	if (ret < 0)
	    qr_code_data_clear(_qrdata);
	_qrdata->version   = _version;
	_qrdata->ecc_level = ecc_level;
    }
    free(block_data);
    return ret;
}

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
	    qr_code_decode(_qrdata, &_reader->gf, ul.c->pos, ur.c->pos,
			   dl.c->pos, ur_version, fmt_info, _img, _width,
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
	    if (qr_code_decode(_qrdata, &_reader->gf, ul.c->pos, dl.c->pos,
			       ur.c->pos, ur_version, fmt_info, _img, _width,
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
