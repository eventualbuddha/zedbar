/*------------------------------------------------------------------------
 *  Copyright 2010 (c) Jeff Brown <spadix@users.sourceforge.net>
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

#include <stdio.h>
#include <stdlib.h>
#include "zbar.h"

#include "decoder.h"

#define GS ('\035')

/* DataBar helper function implementations */
/* Converted to Rust - see src/decoder.rs */

extern void _zbar_databar_new_scan(databar_decoder_t *db);
extern void _zbar_databar_reset(databar_decoder_t *db);

enum {
    SCH_NUM,
    SCH_ALNUM,
    SCH_ISO646
};

// Rust implementations - converted to src/databar_utils.rs
extern void _zbar_databar_append_check14(unsigned char *buf);
extern void _zbar_databar_decode10(unsigned char *buf, unsigned long n, int i);

#define VAR_MAX(l, i) ((((l) * 12 + (i)) * 2 + 6) / 7)

#define FEED_BITS(b)                         \
    while (i < (b) && len) {                 \
	d = (d << 12) | (*(data++) & 0xfff); \
	i += 12;                             \
	len--;                               \
    }

#define PUSH_CHAR(c) *(buf++) = (c)

#define PUSH_CHAR4(c0, c1, c2, c3) \
    do {                           \
	PUSH_CHAR(c0);             \
	PUSH_CHAR(c1);             \
	PUSH_CHAR(c2);             \
	PUSH_CHAR(c3);             \
    } while (0);

int databar_postprocess_exp(zbar_decoder_t *dcode, int *data)
{
    int i = 0, enc;
    unsigned n;
    unsigned char *buf;
    unsigned long d = *(data++);
    int len	    = d / 211 + 4, buflen;

    /* grok encodation method */
    d = *(data++);
    n = (d >> 4) & 0x7f;
    if (n >= 0x40) {
	i      = 10;
	enc    = 1;
	buflen = 2 + 14 + VAR_MAX(len, 10 - 2 - 44 + 6) + 2;
    } else if (n >= 0x38) {
	i      = 4;
	enc    = 6 + (n & 7);
	buflen = 2 + 14 + 4 + 6 + 2 + 6 + 2;
    } else if (n >= 0x30) {
	i      = 6;
	enc    = 2 + ((n >> 2) & 1);
	buflen = 2 + 14 + 4 + 3 + VAR_MAX(len, 6 - 2 - 44 - 2 - 10) + 2;
    } else if (n >= 0x20) {
	i      = 7;
	enc    = 4 + ((n >> 3) & 1);
	buflen = 2 + 14 + 4 + 6;
    } else {
	i      = 9;
	enc    = 0;
	buflen = VAR_MAX(len, 9 - 2) + 2;
    }
    zassert(buflen > 2, -1, "buflen=%d\n", buflen);

    if (enc < 4) {
	/* grok variable length symbol bit field */
	if ((len ^ (d >> (--i))) & 1)
	    /* even/odd length mismatch */
	    return (-1);
	if (((d >> (--i)) & 1) != (len > 14))
	    /* size group mismatch */
	    return (-1);
    }
    len -= 2;

    if (size_buf(dcode, buflen))
	return (-1);
    buf = dcode->buf;

    /* handle compressed fields */
    if (enc) {
	PUSH_CHAR('0');
	PUSH_CHAR('1');
    }

    if (enc == 1) {
	i -= 4;
	n = (d >> i) & 0xf;
	if (i >= 10)
	    return (-1);
	PUSH_CHAR('0' + n);
    } else if (enc)
	PUSH_CHAR('9');

    if (enc) {
	int j;
	for (j = 0; j < 4; j++) {
	    FEED_BITS(10);
	    i -= 10;
	    n = (d >> i) & 0x3ff;
	    if (n >= 1000)
		return (-1);
	    _zbar_databar_decode10(buf, n, 3);
	    buf += 3;
	}
	_zbar_databar_append_check14(buf - 13);
	buf++;
    }

    switch (enc) {
    case 2: /* 01100: AI 392x */
	FEED_BITS(2);
	i -= 2;
	n = (d >> i) & 0x3;
	PUSH_CHAR4('3', '9', '2', '0' + n);
	break;

    case 3: /* 01101: AI 393x */
	FEED_BITS(12);
	i -= 2;
	n = (d >> i) & 0x3;
	PUSH_CHAR4('3', '9', '3', '0' + n);
	i -= 10;
	n = (d >> i) & 0x3ff;
	if (n >= 1000)
	    return (-1);
	_zbar_databar_decode10(buf, n, 3);
	buf += 3;
	break;

    case 4: /* 0100: AI 3103 */
	FEED_BITS(15);
	i -= 15;
	n = (d >> i) & 0x7fff;
	PUSH_CHAR4('3', '1', '0', '3');
	_zbar_databar_decode10(buf, n, 6);
	buf += 6;
	break;

    case 5: /* 0101: AI 3202/3203 */
	FEED_BITS(15);
	i -= 15;
	n = (d >> i) & 0x7fff;
	PUSH_CHAR4('3', '2', '0', (n >= 10000) ? '3' : '2');
	if (n >= 10000)
	    n -= 10000;
	_zbar_databar_decode10(buf, n, 6);
	buf += 6;
	break;
    }
    if (enc >= 6) {
	/* 0111000 - 0111111: AI 310x/320x + AI 11/13/15/17 */
	PUSH_CHAR4('3', '1' + (enc & 1), '0', 'x');
	FEED_BITS(20);
	i -= 20;
	n = (d >> i) & 0xfffff;
	if (n >= 1000000)
	    return (-1);
	_zbar_databar_decode10(buf, n, 6);
	*(buf - 1) = *buf;
	*buf	   = '0';
	buf += 6;

	FEED_BITS(16);
	i -= 16;
	n = (d >> i) & 0xffff;
	if (n < 38400) {
	    int dd, mm, yy;
	    dd = n % 32;
	    n /= 32;
	    mm = n % 12 + 1;
	    n /= 12;
	    yy = n;
	    PUSH_CHAR('1');
	    PUSH_CHAR('0' + ((enc - 6) | 1));
	    _zbar_databar_decode10(buf, yy, 2);
	    buf += 2;
	    _zbar_databar_decode10(buf, mm, 2);
	    buf += 2;
	    _zbar_databar_decode10(buf, dd, 2);
	    buf += 2;
	} else if (n > 38400)
	    return (-1);
    }

    if (enc < 4) {
	/* remainder is general-purpose data compaction */
	int scheme = SCH_NUM;
	while (i > 0 || len > 0) {
	    FEED_BITS(8);

	    if (scheme == SCH_NUM) {
		int n1;
		i -= 4;
		if (i < 0)
		    break;
		if (!((d >> i) & 0xf)) {
		    scheme = SCH_ALNUM;
		    continue;
		}
		if (!len && i < 3) {
		    /* special case last digit */
		    n = ((d >> i) & 0xf) - 1;
		    if (n > 9)
			return (-1);
		    *(buf++) = '0' + n;
		    break;
		}
		i -= 3;
		zassert(i >= 0, -1, "\n");
		n	 = ((d >> i) & 0x7f) - 8;
		n1	 = n % 11;
		n	 = n / 11;
		*(buf++) = (n < 10) ? '0' + n : GS;
		*(buf++) = (n1 < 10) ? '0' + n1 : GS;
	    } else {
		unsigned c = 0;
		i -= 3;
		if (i < 0)
		    break;
		if (!((d >> i) & 0x7)) {
		    scheme = SCH_NUM;
		    continue;
		}
		i -= 2;
		if (i < 0)
		    break;
		n = (d >> i) & 0x1f;
		if (n == 0x04) {
		    scheme ^= 0x3;
		} else if (n == 0x0f)
		    c = GS;
		else if (n < 0x0f)
		    c = 43 + n;
		else if (scheme == SCH_ALNUM) {
		    i--;
		    if (i < 0)
			return (-1);
		    n = (d >> i) & 0x1f;
		    if (n < 0x1a)
			c = 'A' + n;
		    else if (n == 0x1a)
			c = '*';
		    else if (n < 0x1f)
			c = ',' + n - 0x1b;
		    else
			return (-1);
		} else if (scheme == SCH_ISO646 && n < 0x1d) {
		    i -= 2;
		    if (i < 0)
			return (-1);
		    n = (d >> i) & 0x3f;
		    if (n < 0x1a)
			c = 'A' + n;
		    else if (n < 0x34)
			c = 'a' + n - 0x1a;
		    else
			return (-1);
		} else if (scheme == SCH_ISO646) {
		    i -= 3;
		    if (i < 0)
			return (-1);
		    n = ((d >> i) & 0x1f);
		    if (n < 0xa)
			c = '!' + n - 8;
		    else if (n < 0x15)
			c = '%' + n - 0xa;
		    else if (n < 0x1b)
			c = ':' + n - 0x15;
		    else if (n == 0x1b)
			c = '_';
		    else if (n == 0x1c)
			c = ' ';
		    else
			return (-1);
		} else
		    return (-1);

		if (c) {
		    *(buf++) = c;
		}
	    }
	}
	/* FIXME check pad? */
    }

    i = buf - dcode->buf;
    zassert((int)i < (int)dcode->buf_alloc, -1, "i=%02x %s\n", i,
	    _zbar_decoder_buf_dump(dcode->buf, i));

    *buf	  = 0;
    dcode->buflen = i;
    if (i && *--buf == GS) {
	*buf = 0;
	dcode->buflen--;
    }

    return (0);
}
#undef FEED_BITS

// Rust implementation - converted to src/databar_utils.rs
extern int _zbar_databar_check_width(unsigned wf, unsigned wd, unsigned n);

extern signed lookup_sequence(databar_segment_t *seg, int fixed, int seq[22],
			      const size_t maxsize);

#define IDX(s) \
    (((s)->finder << 2) | ((s)->color << 1) | ((s)->color ^ (s)->side))

/* Temporarily exported for Rust */
zbar_symbol_type_t match_segment_exp(zbar_decoder_t *dcode,
				     databar_segment_t *seg, int dir)
{
    databar_decoder_t *db = &dcode->databar;
    int bestsegs[22], i = 0, segs[22], seq[22];
    int ifixed = seg - db->segs, fixed = IDX(seg), maxcnt = 0;
    int iseg[DATABAR_MAX_SEGMENTS];
    unsigned csegs = db->csegs, width = seg->width, maxage = 0x7fff;

    bestsegs[0] = segs[0] = seq[1] = -1;
    seq[0]			   = 0;
    for (i = csegs, seg = db->segs + csegs - 1; --i >= 0; seg--) {
	if (seg->exp && seg->finder >= 0 && (!seg->partial || seg->count >= 4))
	    iseg[i] = IDX(seg);
	else
	    iseg[i] = -1;
    }

    for (i = 0;; i--) {
	unsigned cnt, chk, age;
	unsigned data0, chk0;
	for (; i >= 0 && seq[i] >= 0; i--) {
	    int j;

	    if (seq[i] == fixed) {
		seg = db->segs + ifixed;
		if (segs[i] < 0 &&
		    _zbar_databar_check_width(width, seg->width, 14)) {
		    j = ifixed;
		} else
		    continue;
	    } else {
		for (j = segs[i] + 1; (int)j < (int)csegs; j++) {
		    if (iseg[j] == seq[i] &&
			(!i || _zbar_databar_check_width(
				   width, db->segs[j].width, 14))) {
			seg = db->segs + j;
			break;
		    }
		}
		if ((int)j == (int)csegs)
		    continue;
	    }

	    if (!i) {
		signed int lu = lookup_sequence(seg, fixed, seq,
						sizeof(seq) / sizeof(seq[0]));
		if (!lu) {
		    continue;
		}
		if (lu < 0) {
		    goto abort;
		}
		width = seg->width;
	    } else {
		width = (width + seg->width) / 2;
	    }
	    segs[i++] = j;
	    segs[i++] = -1;
	}
	if (i < 0)
	    break;

	seg = db->segs + segs[0];
	cnt = 0, chk = 0, age = (db->epoch - seg->epoch) & 0xff;
	for (i = 1; segs[i] >= 0; i++) {
	    seg = db->segs + segs[i];
	    chk += seg->check;
	    cnt += seg->count;
	    age += (db->epoch - seg->epoch) & 0xff;
	}

	data0 = db->segs[segs[0]].data;
	chk0  = data0 % 211;
	chk %= 211;
	if (chk != chk0)
	    continue;
	if (maxcnt > (int)cnt ||
	    (maxcnt == (int)cnt && (int)maxage <= (int)age))
	    continue;
	maxcnt = cnt;
	maxage = age;
	for (i = 0; segs[i] >= 0; i++)
	    bestsegs[i] = segs[i];
	bestsegs[i] = -1;
    }

    if (bestsegs[0] < 0)
	return (ZBAR_PARTIAL);

    if (acquire_lock(dcode, ZBAR_DATABAR_EXP))
	return (ZBAR_PARTIAL);

    for (i = 0; bestsegs[i] >= 0; i++)
	segs[i] = db->segs[bestsegs[i]].data;

    if (databar_postprocess_exp(dcode, segs)) {
	release_lock(dcode, ZBAR_DATABAR_EXP);
	return (ZBAR_PARTIAL);
    }

    for (i = 0; bestsegs[i] >= 0; i++)
	if (bestsegs[i] != ifixed) {
	    seg = db->segs + bestsegs[i];
	    if (!--seg->count)
		seg->finder = -1;
	}

    /* FIXME stacked rows are frequently reversed,
     * so direction is impossible to determine at this level
     */
    dcode->direction = (1 - 2 * (seg->side ^ seg->color)) * dir;
    dcode->modifiers = MOD(ZBAR_MOD_GS1);
    return (ZBAR_DATABAR_EXP);
abort:
    return (ZBAR_NONE);
}
#undef IDX
