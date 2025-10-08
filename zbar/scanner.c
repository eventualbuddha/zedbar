/*------------------------------------------------------------------------
 *  Copyright 2007-2010 (c) Jeff Brown <spadix@users.sourceforge.net>
 *
 *  This file is part of the ZBar Bar Code Reader.
 *
 *  The ZBar Bar Code Reader is free software; you can redistribute it
 *  and/or modify it under the terms of the GNU Lesser Public License as
 *  published by the Free Software Foundation; either version 2.1 of
 *  the License, or (at your option) any later version.
 *
 *  The ZBar Bar Code Reader is distributed in the hope that it will be
 *  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser Public License
 *  along with the ZBar Bar Code Reader; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 *  Boston, MA  02110-1301  USA
 *
 *  http://sourceforge.net/projects/zbar
 *------------------------------------------------------------------------*/

#include <stddef.h>
#include <stdlib.h> /* malloc, free, abs */
#include <string.h> /* memset */

#include "zbar.h"

#ifndef ZBAR_FIXED
#define ZBAR_FIXED 5
#endif
#define ROUND (1 << (ZBAR_FIXED - 1))

/* FIXME add runtime config API for these */
#ifndef ZBAR_SCANNER_THRESH_MIN
#define ZBAR_SCANNER_THRESH_MIN 4
#endif

#ifndef ZBAR_SCANNER_THRESH_INIT_WEIGHT
#define ZBAR_SCANNER_THRESH_INIT_WEIGHT .44
#endif
#define THRESH_INIT                                                          \
    ((unsigned)((ZBAR_SCANNER_THRESH_INIT_WEIGHT * (1 << (ZBAR_FIXED + 1)) + \
		 1) /                                                        \
		2))

#ifndef ZBAR_SCANNER_THRESH_FADE
#define ZBAR_SCANNER_THRESH_FADE 8
#endif

#ifndef ZBAR_SCANNER_EWMA_WEIGHT
#define ZBAR_SCANNER_EWMA_WEIGHT .78
#endif
#define EWMA_WEIGHT \
    ((unsigned)((ZBAR_SCANNER_EWMA_WEIGHT * (1 << (ZBAR_FIXED + 1)) + 1) / 2))

/* scanner state */
struct zbar_scanner_s {
    zbar_decoder_t *decoder; /* associated bar width decoder */
    unsigned y1_min_thresh;  /* minimum threshold */

    unsigned x; /* relative scan position of next sample */
    int y0[4];	/* short circular buffer of average intensities */

    int y1_sign;	/* slope at last crossing */
    unsigned y1_thresh; /* current slope threshold */

    unsigned cur_edge;	/* interpolated position of tracking edge */
    unsigned last_edge; /* interpolated position of last located edge */
    unsigned width;	/* last element width */
};

zbar_scanner_t *zbar_scanner_create(zbar_decoder_t *dcode)
{
    zbar_scanner_t *scn = malloc(sizeof(zbar_scanner_t));
    scn->decoder	= dcode;
    scn->y1_min_thresh	= ZBAR_SCANNER_THRESH_MIN;
    zbar_scanner_reset(scn);
    return (scn);
}

// zbar_scanner_destroy() implemented in Rust (src/line_scanner.rs)

zbar_symbol_type_t zbar_scanner_reset(zbar_scanner_t *scn)
{
    memset(&scn->x, 0, sizeof(zbar_scanner_t) - offsetof(zbar_scanner_t, x));
    scn->y1_thresh = scn->y1_min_thresh;
    if (scn->decoder)
	zbar_decoder_reset(scn->decoder);
    return (ZBAR_NONE);
}

// zbar_scanner_get_width() implemented in Rust (src/line_scanner.rs)
// zbar_scanner_get_edge() implemented in Rust (src/line_scanner.rs)
// zbar_scanner_get_color() implemented in Rust (src/line_scanner.rs)

// calc_thresh() implemented in Rust (src/line_scanner.rs)
// process_edge() implemented in Rust (src/line_scanner.rs)
// zbar_scanner_flush() implemented in Rust (src/line_scanner.rs)
// zbar_scanner_new_scan() implemented in Rust (src/line_scanner.rs)
// zbar_scan_y() implemented in Rust (src/line_scanner.rs)
