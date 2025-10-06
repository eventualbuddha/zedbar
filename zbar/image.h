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
#ifndef _IMAGE_H_
#define _IMAGE_H_

#include "config.h"
#include <assert.h>
#include <inttypes.h>
#include <stdlib.h>

#include <zbar.h>
#include "error.h"
#include "refcnt.h"
#include "symbol.h"

#define fourcc zbar_fourcc

/* unpack size/location of component */
#define RGB_SIZE(c)   ((c) >> 5)
#define RGB_OFFSET(c) ((c) & 0x1f)

/* coarse image format categorization.
 * to limit conversion variations
 */
typedef enum zbar_format_group_e {
    ZBAR_FMT_GRAY,
    ZBAR_FMT_YUV_PLANAR,
    ZBAR_FMT_YUV_PACKED,
    ZBAR_FMT_RGB_PACKED,
    ZBAR_FMT_YUV_NV,
    ZBAR_FMT_JPEG,

    /* enum size */
    ZBAR_FMT_NUM
} zbar_format_group_t;

/* Struct definition moved to zbar.h - this is just a forward declaration */

/* description of an image format */
typedef struct zbar_format_def_s {
    uint32_t format;	       /* fourcc */
    zbar_format_group_t group; /* coarse categorization */
    union {
	uint8_t gen[4]; /* raw bytes */
	struct {
	    uint8_t bpp;	      /* bits per pixel */
	    uint8_t red, green, blue; /* size/location a la RGB_BITS() */
	} rgb;
	struct {
	    uint8_t xsub2, ysub2; /* chroma subsampling in each axis */
	    uint8_t packorder;	  /* channel ordering flags
                             *   bit0: 0=UV, 1=VU
                             *   bit1: 0=Y/chroma, 1=chroma/Y
                             */
	} yuv;
	uint32_t cmp; /* quick compare equivalent formats */
    } p;
} zbar_format_def_t;

extern void _zbar_image_swap_symbols(zbar_image_t *a, zbar_image_t *b);
extern void _zbar_image_copy_size(zbar_image_t *dst, const zbar_image_t *src);
extern zbar_image_t *_zbar_image_copy(const zbar_image_t *src, int inverted);

#endif
