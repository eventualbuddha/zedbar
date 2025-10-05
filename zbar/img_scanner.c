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

#include "config.h"
#include <inttypes.h>
#include <unistd.h>

#include <assert.h>
#include <stdlib.h> /* malloc, free */
#include <string.h> /* memcmp, memset, memcpy */

#include "error.h"
#include "image.h"

#include "zbar.h"
#include "img_scanner.h"
#include "qrcode.h"
#include "sqcode.h"

/* FIXME cache setting configurability */

/* time interval for which two images are considered "nearby"
 */
#define CACHE_PROXIMITY 1000 /* ms */

/* time that a result must *not* be detected before
 * it will be reported again
 */
#define CACHE_HYSTERESIS 2000 /* ms */

/* time after which cache entries are invalidated
 */
#define CACHE_TIMEOUT (CACHE_HYSTERESIS * 2) /* ms */

#define NUM_SCN_CFGS (ZBAR_CFG_Y_DENSITY - ZBAR_CFG_X_DENSITY + 1)

#define CFG(iscn, cfg) ((iscn)->configs[(cfg) - ZBAR_CFG_X_DENSITY])
#define TEST_CFG(iscn, cfg) \
    (((iscn)->config >> ((cfg) - ZBAR_CFG_POSITION)) & 1)

#define RECYCLE_BUCKETS 5

typedef struct recycle_bucket_s {
    int nsyms;
    zbar_symbol_t *head;
} recycle_bucket_t;

/* image scanner state */
struct zbar_image_scanner_s {
    zbar_scanner_t *scn;   /* associated linear intensity scanner */
    zbar_decoder_t *dcode; /* associated symbol decoder */
    qr_reader *qr;	   /* QR Code 2D reader */
    sq_reader *sq;	   /* SQ Code 2D reader */

    const void *userdata; /* application data */
    /* user result callback */
    zbar_image_data_handler_t *handler;

    unsigned long time;	     /* scan start time */
    zbar_image_t *img;	     /* currently scanning image *root* */
    int dx, dy, du, umin, v; /* current scan direction */
    zbar_symbol_set_t *syms; /* previous decode results */
    /* recycled symbols in 4^n size buckets */
    recycle_bucket_t recycle[RECYCLE_BUCKETS];

    int enable_cache;	  /* current result cache state */
    zbar_symbol_t *cache; /* inter-image result cache entries */

    /* configuration settings */
    unsigned config; /* config flags */
    unsigned ean_config;
    int configs[NUM_SCN_CFGS];	  /* int valued configurations */
    int sym_configs[1][NUM_SYMS]; /* per-symbology configurations */
};

// Rust implementation - converted to src/img_scanner.rs
extern void _zbar_image_scanner_sq_handler(zbar_image_scanner_t *iscn);

static inline void sq_handler(zbar_image_scanner_t *iscn)
{
    _zbar_image_scanner_sq_handler(iscn);
}

void symbol_handler(zbar_decoder_t *dcode)
{
    zbar_image_scanner_t *iscn = zbar_decoder_get_userdata(dcode);
    zbar_symbol_type_t type    = zbar_decoder_get_type(dcode);
    int x = 0, y = 0, dir;
    const char *data;
    unsigned datalen;
    zbar_symbol_t *sym;

    if (type == ZBAR_QRCODE) {
	_zbar_image_scanner_qr_handler(iscn);
	return;
    }

    if (TEST_CFG(iscn, ZBAR_CFG_POSITION)) {
	/* tmp position fixup */
	int w = zbar_scanner_get_width(iscn->scn);
	int u = iscn->umin + iscn->du * zbar_scanner_get_edge(iscn->scn, w, 0);
	if (iscn->dx) {
	    x = u;
	    y = iscn->v;
	} else {
	    x = iscn->v;
	    y = u;
	}
    }

    /* FIXME debug flag to save/display all PARTIALs */
    if (type <= ZBAR_PARTIAL) {
	zprintf(256, "partial symbol @(%d,%d)\n", x, y);
	return;
    }

    data    = zbar_decoder_get_data(dcode);
    datalen = zbar_decoder_get_data_length(dcode);

    /* FIXME need better symbol matching */
    for (sym = iscn->syms->head; sym; sym = sym->next)
	if (sym->type == type && sym->datalen == datalen &&
	    !memcmp(sym->data, data, datalen)) {
	    sym->quality++;
	    zprintf(224, "dup symbol @(%d,%d): dup %s: %.20s\n", x, y,
		    zbar_get_symbol_name(type), data);
	    if (TEST_CFG(iscn, ZBAR_CFG_POSITION))
		/* add new point to existing set */
		/* FIXME should be polygon */
		_zbar_symbol_add_point(sym, x, y);
	    return;
	}

    sym		   = _zbar_image_scanner_alloc_sym(iscn, type, datalen + 1);
    sym->configs   = zbar_decoder_get_configs(dcode, type);
    sym->modifiers = zbar_decoder_get_modifiers(dcode);
    /* FIXME grab decoder buffer */
    memcpy(sym->data, data, datalen + 1);

    /* initialize first point */
    if (TEST_CFG(iscn, ZBAR_CFG_POSITION)) {
	zprintf(192, "new symbol @(%d,%d): %s: %.20s\n", x, y,
		zbar_get_symbol_name(type), data);
	_zbar_symbol_add_point(sym, x, y);
    }

    dir = zbar_decoder_get_direction(dcode);
    if (dir)
	sym->orient = (iscn->dy != 0) + ((iscn->du ^ dir) & 2);

    _zbar_image_scanner_add_sym(iscn, sym);
}

// Rust implementation
extern zbar_image_scanner_t *_zbar_image_scanner_create_rust(void);
extern void _zbar_image_scanner_destroy_rust(zbar_image_scanner_t *iscn);

inline zbar_image_scanner_t *zbar_image_scanner_create()
{
    return _zbar_image_scanner_create_rust();
}

inline void zbar_image_scanner_destroy(zbar_image_scanner_t *iscn)
{
    _zbar_image_scanner_destroy_rust(iscn);
}

int zbar_image_scanner_set_config(zbar_image_scanner_t *iscn,
				  zbar_symbol_type_t sym, zbar_config_t cfg,
				  int val)
{
    if ((sym == 0 || sym == ZBAR_COMPOSITE) && cfg == ZBAR_CFG_ENABLE) {
	iscn->ean_config = !!val;
	if (sym)
	    return (0);
    }

    if (cfg < ZBAR_CFG_UNCERTAINTY)
	return (zbar_decoder_set_config(iscn->dcode, sym, cfg, val));

    if (cfg < ZBAR_CFG_POSITION) {
	int c, i;
	if (cfg > ZBAR_CFG_UNCERTAINTY)
	    return (1);
	c = cfg - ZBAR_CFG_UNCERTAINTY;
	if (sym > ZBAR_PARTIAL) {
	    i			    = _zbar_get_symbol_hash(sym);
	    iscn->sym_configs[c][i] = val;
	} else
	    for (i = 0; i < NUM_SYMS; i++)
		iscn->sym_configs[c][i] = val;
	return (0);
    }

    /* Image scanner parameters apply only to ZBAR_PARTIAL */
    if (sym > ZBAR_PARTIAL)
	return (1);

    if (cfg >= ZBAR_CFG_X_DENSITY && cfg <= ZBAR_CFG_Y_DENSITY) {
	CFG(iscn, cfg) = val;
	return (0);
    }

    cfg -= ZBAR_CFG_POSITION;

    if (!val)
	iscn->config &= ~(1 << cfg);
    else if (val == 1)
	iscn->config |= (1 << cfg);
    else
	return (1);

    return (0);
}

#define movedelta(dx, dy)                  \
    do {                                   \
	x += (dx);                         \
	y += (dy);                         \
	p += (dx) + ((uintptr_t)(dy) * w); \
    } while (0);

static void *_zbar_scan_image(zbar_image_scanner_t *iscn, zbar_image_t *img)
{
    zbar_symbol_set_t *syms;
    const uint8_t *data;
    zbar_scanner_t *scn = iscn->scn;
    unsigned w, h;
    int density;
    char filter;
    int nean, naddon;

    /* timestamp image with simple counter
   * FIXME prefer video timestamp
   */
    static unsigned int scan_counter = 0;
    iscn->time			     = scan_counter++;

    _zbar_qr_reset(iscn->qr);

    _zbar_sq_reset(iscn->sq);

    /* image must be in grayscale format */
    if (img->format != fourcc('Y', '8', '0', '0') &&
	img->format != fourcc('G', 'R', 'E', 'Y'))
	return NULL;
    iscn->img = img;

    /* recycle previous scanner and image results */
    zbar_image_scanner_recycle_image(iscn, img);
    syms = iscn->syms;
    if (!syms) {
	syms = iscn->syms = _zbar_symbol_set_create();
	zbar_symbol_set_ref(syms, 1);
    } else
	zbar_symbol_set_ref(syms, 2);
    img->syms = syms;

    w	 = img->width;
    h	 = img->height;
    data = img->data;

    zbar_scanner_new_scan(scn);

    density = CFG(iscn, ZBAR_CFG_Y_DENSITY);
    if (density > 0) {
	const uint8_t *p = data;
	int x = 0, y = 0;

	int border = (((img->height - 1) % density) + 1) / 2;
	if ((unsigned)border > img->height / 2)
	    border = img->height / 2;
	assert((unsigned)border <= h);
	iscn->dy = 0;

	movedelta(0, border);
	iscn->v = y;

	while ((unsigned)y < h) {
	    iscn->dx = iscn->du = 1;
	    iscn->umin		= 0;
	    while ((unsigned)x < w) {
		uint8_t d = *p;
		movedelta(1, 0);
		zbar_scan_y(scn, d);
	    }
	    _zbar_image_scanner_quiet_border(iscn);

	    movedelta(-1, density);
	    iscn->v = y;
	    if ((unsigned)y >= h)
		break;

	    iscn->dx = iscn->du = -1;
	    iscn->umin		= w;
	    while (x >= 0) {
		uint8_t d = *p;
		movedelta(-1, 0);
		zbar_scan_y(scn, d);
	    }
	    _zbar_image_scanner_quiet_border(iscn);

	    movedelta(1, density);
	    iscn->v = y;
	}
    }
    iscn->dx = 0;

    density = CFG(iscn, ZBAR_CFG_X_DENSITY);
    if (density > 0) {
	const uint8_t *p = data;
	int x = 0, y = 0;

	int border = (((img->width - 1) % density) + 1) / 2;
	if ((unsigned)border > img->width / 2)
	    border = img->width / 2;
	assert((unsigned)border <= w);
	movedelta(border, 0);
	iscn->v = x;

	while ((unsigned)x < w) {
	    zprintf(128, "img_y+: %04d,%04d @%p\n", x, y, p);
	    iscn->dy = iscn->du = 1;
	    iscn->umin		= 0;
	    while ((unsigned)y < h) {
		uint8_t d = *p;
		movedelta(0, 1);
		zbar_scan_y(scn, d);
	    }
	    _zbar_image_scanner_quiet_border(iscn);

	    movedelta(density, -1);
	    iscn->v = x;
	    if ((unsigned)x >= w)
		break;

	    zprintf(128, "img_y-: %04d,%04d @%p\n", x, y, p);
	    iscn->dy = iscn->du = -1;
	    iscn->umin		= h;
	    while (y >= 0) {
		uint8_t d = *p;
		movedelta(0, -1);
		zbar_scan_y(scn, d);
	    }
	    _zbar_image_scanner_quiet_border(iscn);

	    movedelta(density, 1);
	    iscn->v = x;
	}
    }
    iscn->dy  = 0;
    iscn->img = NULL;

    _zbar_qr_decode(iscn->qr, iscn, img);

    sq_handler(iscn);
    _zbar_sq_decode(iscn->sq, iscn, img);

    /* FIXME tmp hack to filter bad EAN results */
    /* FIXME tmp hack to merge simple case EAN add-ons */
    filter = (!iscn->enable_cache &&
	      (density == 1 || CFG(iscn, ZBAR_CFG_Y_DENSITY) == 1));
    nean   = 0;
    naddon = 0;
    if (syms->nsyms) {
	zbar_symbol_t **symp;
	for (symp = &syms->head; *symp;) {
	    zbar_symbol_t *sym = *symp;
	    if (sym->cache_count <= 0 &&
		((sym->type < ZBAR_COMPOSITE && sym->type > ZBAR_PARTIAL) ||
		 sym->type == ZBAR_DATABAR || sym->type == ZBAR_DATABAR_EXP ||
		 sym->type == ZBAR_CODABAR)) {
		if ((sym->type == ZBAR_CODABAR || filter) && sym->quality < 4) {
		    if (iscn->enable_cache) {
			/* revert cache update */
			zbar_symbol_t *entry =
			    _zbar_image_scanner_cache_lookup(iscn, sym);
			if (entry)
			    entry->cache_count--;
			else
			    assert(0);
		    }

		    /* recycle */
		    *symp = sym->next;
		    syms->nsyms--;
		    sym->next = NULL;
		    _zbar_image_scanner_recycle_syms(iscn, sym);
		    continue;
		} else if (sym->type < ZBAR_COMPOSITE &&
			   sym->type != ZBAR_ISBN10) {
		    if (sym->type > ZBAR_EAN5)
			nean++;
		    else
			naddon++;
		}
	    }
	    symp = &sym->next;
	}

	if (nean == 1 && naddon == 1 && iscn->ean_config) {
	    int datalen;
	    zbar_symbol_t *ean_sym;
	    /* create container symbol for composite result */
	    zbar_symbol_t *ean = NULL, *addon = NULL;
	    for (symp = &syms->head; *symp;) {
		zbar_symbol_t *sym = *symp;
		if (sym->type < ZBAR_COMPOSITE && sym->type > ZBAR_PARTIAL) {
		    /* move to composite */
		    *symp = sym->next;
		    syms->nsyms--;
		    sym->next = NULL;
		    if (sym->type <= ZBAR_EAN5)
			addon = sym;
		    else
			ean = sym;
		} else
		    symp = &sym->next;
	    }
	    assert(ean);
	    assert(addon);

	    datalen = ean->datalen + addon->datalen + 1;
	    ean_sym =
		_zbar_image_scanner_alloc_sym(iscn, ZBAR_COMPOSITE, datalen);
	    ean_sym->orient = ean->orient;
	    ean_sym->syms   = _zbar_symbol_set_create();
	    memcpy(ean_sym->data, ean->data, ean->datalen);
	    memcpy(ean_sym->data + ean->datalen, addon->data,
		   addon->datalen + 1);
	    ean_sym->syms->head	 = ean;
	    ean->next		 = addon;
	    ean_sym->syms->nsyms = 2;
	    _zbar_image_scanner_add_sym(iscn, ean_sym);
	}
    }
    return syms;
}

int zbar_scan_image(zbar_image_scanner_t *iscn, zbar_image_t *img)
{
    zbar_symbol_set_t *syms;
    zbar_image_t *inv = NULL;

    syms = _zbar_scan_image(iscn, img);
    if (!syms)
	return -1;

    if (!syms->nsyms && TEST_CFG(iscn, ZBAR_CFG_TEST_INVERTED)) {
	inv = _zbar_image_copy(img, 1);
	if (inv) {
	    if (iscn->cache) {
		/* recycle all cached syms, if any */
		_zbar_image_scanner_recycle_syms(iscn, iscn->cache);
		iscn->cache = NULL;
	    }
	    syms = _zbar_scan_image(iscn, inv);
	    _zbar_image_swap_symbols(img, inv);
	}
    }

    if (syms->nsyms && iscn->handler)
	iscn->handler(img, iscn->userdata);

    if (inv)
	zbar_image_destroy(inv);

    return (syms->nsyms);
}
