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
#ifndef _EAN_H_
#define _EAN_H_

/* state of each parallel decode attempt */
typedef struct ean_pass_s {
    signed char state;	  /* module position of w[idx] in symbol */
#define STATE_REV   0x80  /*   scan direction reversed */
#define STATE_ADDON 0x40  /*   scanning add-on */
#define STATE_IDX   0x3f  /*   element offset into symbol */
    unsigned width;	  /* width of last character */
    unsigned char raw[7]; /* decode in process */
} ean_pass_t;

/* EAN/UPC specific decode state */
typedef struct ean_decoder_s {
    ean_pass_t pass[4];	     /* state of each parallel decode attempt */
    zbar_symbol_type_t left; /* current holding buffer contents */
    zbar_symbol_type_t right;
    int direction;	 /* scan direction */
    unsigned s4, width;	 /* character width */
    signed char buf[18]; /* holding buffer */

    signed char enable;
    unsigned ean13_config;
    unsigned ean8_config;
    unsigned upca_config;
    unsigned upce_config;
    unsigned isbn10_config;
    unsigned isbn13_config;
    unsigned ean5_config;
    unsigned ean2_config;
} ean_decoder_t;

/* reset EAN/UPC pass specific state */
extern void _zbar_ean_new_scan(ean_decoder_t *ean);
#define ean_new_scan(ean) _zbar_ean_new_scan(ean)

/* reset all EAN/UPC state */
extern void _zbar_ean_reset(ean_decoder_t *ean);
#define ean_reset(ean) _zbar_ean_reset(ean)

extern unsigned _zbar_ean_get_config(ean_decoder_t *ean, zbar_symbol_type_t sym);
#define ean_get_config(ean, sym) _zbar_ean_get_config(ean, sym)

/* decode EAN/UPC symbols */
zbar_symbol_type_t _zbar_decode_ean(zbar_decoder_t *dcode);

#endif
