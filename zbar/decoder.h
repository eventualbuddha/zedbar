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
#ifndef _DECODER_H_
#define _DECODER_H_

#include "config.h"
#include <limits.h>
#include <stdlib.h> /* realloc */

#include <zbar.h>

#include "debug.h"

#define NUM_CFGS (ZBAR_CFG_MAX_LEN - ZBAR_CFG_MIN_LEN + 1)

#include "decoder/ean.h"
#include "decoder/i25.h"
#include "decoder/databar.h"
#include "decoder/codabar.h"
#include "decoder/code39.h"
#include "decoder/code93.h"
#include "decoder/code128.h"
#include "decoder/qr_finder.h"
#include "decoder/sq_finder.h"

/* size of bar width history (implementation assumes power of two) */
#ifndef DECODE_WINDOW
#define DECODE_WINDOW 16
#endif

/* initial data buffer allocation */
#ifndef BUFFER_MIN
#define BUFFER_MIN 0x20
#endif

/* maximum data buffer allocation
 * (longer symbols are rejected)
 */
#ifndef BUFFER_MAX
#define BUFFER_MAX 0x100
#endif

/* buffer allocation increment */
#ifndef BUFFER_INCR
#define BUFFER_INCR 0x10
#endif

#define CFG(dcode, cfg)	      ((dcode).configs[(cfg) - ZBAR_CFG_MIN_LEN])
#define TEST_CFG(config, cfg) (((config) >> (cfg)) & 1)
#define MOD(mod)	      (1 << (mod))

/* symbology independent decoder state */
struct zbar_decoder_s {
    unsigned char idx;	       /* current width index */
    unsigned w[DECODE_WINDOW]; /* window of last N bar widths */
    zbar_symbol_type_t type;   /* type of last decoded data */
    zbar_symbol_type_t lock;   /* buffer lock */
    unsigned modifiers;	       /* symbology modifier */
    int direction;	       /* direction of last decoded data */
    unsigned s6;	       /* 6-element character width */

    /* everything above here is automatically reset */
    unsigned buf_alloc;		     /* dynamic buffer allocation */
    unsigned buflen;		     /* binary data length */
    unsigned char *buf;		     /* decoded characters */
    void *userdata;		     /* application data */
    zbar_decoder_handler_t *handler; /* application callback */

    /* symbology specific state */
    ean_decoder_t ean;		 /* EAN/UPC parallel decode attempts */
    i25_decoder_t i25;		 /* Interleaved 2 of 5 decode state */
    databar_decoder_t databar;	 /* DataBar decode state */
    codabar_decoder_t codabar;	 /* Codabar decode state */
    code39_decoder_t code39;	 /* Code 39 decode state */
    code93_decoder_t code93;	 /* Code 93 decode state */
    code128_decoder_t code128;	 /* Code 128 decode state */
    qr_finder_t qrf;		 /* QR Code finder state */
    sq_finder_t sqf;		 /* SQ Code finder state */
};

/* Helper functions for decoder - implementations in decoder.c */
extern char _zbar_decoder_get_color(const zbar_decoder_t *dcode);
extern unsigned _zbar_decoder_get_width(const zbar_decoder_t *dcode, unsigned char offset);
extern unsigned _zbar_decoder_pair_width(const zbar_decoder_t *dcode, unsigned char offset);
extern unsigned _zbar_decoder_calc_s(const zbar_decoder_t *dcode, unsigned char offset, unsigned char n);
extern int _zbar_decoder_decode_e(unsigned e, unsigned s, unsigned n);
extern unsigned _zbar_decoder_decode_sort3(zbar_decoder_t *dcode, int i0);
extern unsigned _zbar_decoder_decode_sortn(zbar_decoder_t *dcode, int n, int i0);
extern char _zbar_decoder_acquire_lock(zbar_decoder_t *dcode, zbar_symbol_type_t req);
extern char _zbar_decoder_release_lock(zbar_decoder_t *dcode, zbar_symbol_type_t req);
extern char _zbar_decoder_size_buf(zbar_decoder_t *dcode, unsigned len);

/* Compatibility macros for C code that still uses the old names */
#define get_color(dcode) _zbar_decoder_get_color(dcode)
#define get_width(dcode, offset) _zbar_decoder_get_width(dcode, offset)
#define pair_width(dcode, offset) _zbar_decoder_pair_width(dcode, offset)
#define calc_s(dcode, offset, n) _zbar_decoder_calc_s(dcode, offset, n)
#define decode_e(e, s, n) _zbar_decoder_decode_e(e, s, n)
#define decode_sort3(dcode, i0) _zbar_decoder_decode_sort3(dcode, i0)
#define decode_sortn(dcode, n, i0) _zbar_decoder_decode_sortn(dcode, n, i0)
#define acquire_lock(dcode, req) _zbar_decoder_acquire_lock(dcode, req)
#define release_lock(dcode, req) _zbar_decoder_release_lock(dcode, req)
#define size_buf(dcode, len) _zbar_decoder_size_buf(dcode, len)

extern const char *_zbar_decoder_buf_dump(unsigned char *buf,
					  unsigned int buflen);

#endif
