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

/* DataBar character encoding groups */
struct group_s {
    unsigned short sum;
    unsigned char wmax;
    unsigned char todd;
    unsigned char teven;
} groups[] = {
    /* (17,4) DataBar Expanded character groups */
    { 0, 7, 87, 4 },
    { 348, 5, 52, 20 },
    { 1388, 4, 30, 52 },
    { 2948, 3, 10, 104 },
    { 3988, 1, 1, 204 },

    /* (16,4) DataBar outer character groups */
    { 0, 8, 161, 1 },
    { 161, 6, 80, 10 },
    { 961, 4, 31, 34 },
    { 2015, 3, 10, 70 },
    { 2715, 1, 1, 126 },

    /* (15,4) DataBar inner character groups */
    { 1516, 8, 81, 1 },
    { 1036, 6, 48, 10 },
    { 336, 4, 20, 35 },
    { 0, 2, 4, 84 },
};

/* DataBar expanded checksum multipliers */
static const unsigned char exp_checksums[] = { 1,   189, 62, 113, 46,  43,
					       109, 134, 6,  79,  161, 45 };

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

/* Converted to Rust - see src/decoder.rs */
extern void _zbar_databar_postprocess(zbar_decoder_t *dcode, unsigned d[4]);

// Rust implementation - converted to src/databar_utils.rs
extern int _zbar_databar_check_width(unsigned wf, unsigned wd, unsigned n);

// Compatibility wrapper
static inline int check_width(unsigned wf, unsigned wd, unsigned n)
{
    return _zbar_databar_check_width(wf, wd, n);
}

/* Converted to Rust - see src/decoder.rs */
extern void _zbar_databar_merge_segment(databar_decoder_t *db,
					databar_segment_t *seg);

static zbar_symbol_type_t match_segment(zbar_decoder_t *dcode,
					databar_segment_t *seg)
{
    databar_decoder_t *db = &dcode->databar;
    unsigned csegs = db->csegs, maxage = 0xfff;
    int i0, i1, i2, maxcnt = 0;
    databar_segment_t *smax[3] = {
	NULL,
    };
    unsigned d[4];

    if (seg->partial && seg->count < 4)
	return (ZBAR_PARTIAL);

    for (i0 = 0; (int)i0 < (int)csegs; i0++) {
	databar_segment_t *s0 = db->segs + i0;
	if (s0 == seg || s0->finder != seg->finder || s0->exp ||
	    s0->color != seg->color || s0->side == seg->side ||
	    (s0->partial && s0->count < 4) ||
	    !check_width(seg->width, s0->width, 14))
	    continue;

	for (i1 = 0; (int)i1 < (int)csegs; i1++) {
	    databar_segment_t *s1 = db->segs + i1;
	    int chkf, chks, chk;
	    unsigned age1;
	    if (i1 == i0 || s1->finder < 0 || s1->exp ||
		s1->color == seg->color || (s1->partial && s1->count < 4) ||
		!check_width(seg->width, s1->width, 14))
		continue;

	    if (seg->color)
		chkf = seg->finder + s1->finder * 9;
	    else
		chkf = s1->finder + seg->finder * 9;
	    if (chkf > 72)
		chkf--;
	    if (chkf > 8)
		chkf--;

	    chks = (seg->check + s0->check + s1->check) % 79;

	    if (chkf >= chks)
		chk = chkf - chks;
	    else
		chk = 79 + chkf - chks;

	    age1 = (((db->epoch - s0->epoch) & 0xff) +
		    ((db->epoch - s1->epoch) & 0xff));

	    for (i2 = i1 + 1; (int)i2 < (int)csegs; i2++) {
		databar_segment_t *s2 = db->segs + i2;
		unsigned cnt, age2, age;
		if (i2 == i0 || s2->finder != s1->finder || s2->exp ||
		    s2->color != s1->color || s2->side == s1->side ||
		    s2->check != chk || (s2->partial && s2->count < 4) ||
		    !check_width(seg->width, s2->width, 14))
		    continue;
		age2 = (db->epoch - s2->epoch) & 0xff;
		age  = age1 + age2;
		cnt  = s0->count + s1->count + s2->count;
		if (maxcnt < (int)cnt ||
		    (maxcnt == (int)cnt && (int)maxage > (int)age)) {
		    maxcnt  = cnt;
		    maxage  = age;
		    smax[0] = s0;
		    smax[1] = s1;
		    smax[2] = s2;
		}
	    }
	}
    }

    if (!smax[0])
	return (ZBAR_PARTIAL);

    d[(seg->color << 1) | seg->side] = seg->data;
    for (i0 = 0; i0 < 3; i0++) {
	d[(smax[i0]->color << 1) | smax[i0]->side] = smax[i0]->data;
	if (!--(smax[i0]->count))
	    smax[i0]->finder = -1;
    }
    seg->finder = -1;

    if (size_buf(dcode, 18))
	return (ZBAR_PARTIAL);

    if (acquire_lock(dcode, ZBAR_DATABAR))
	return (ZBAR_PARTIAL);

    _zbar_databar_postprocess(dcode, d);
    dcode->modifiers = MOD(ZBAR_MOD_GS1);
    dcode->direction = 1 - 2 * (seg->side ^ seg->color ^ 1);
    return (ZBAR_DATABAR);
}

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
		if (segs[i] < 0 && check_width(width, seg->width, 14)) {
		    j = ifixed;
		} else
		    continue;
	    } else {
		for (j = segs[i] + 1; (int)j < (int)csegs; j++) {
		    if (iseg[j] == seq[i] &&
			(!i || check_width(width, db->segs[j].width, 14))) {
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

// Rust implementation - converted to src/databar_utils.rs
extern unsigned _zbar_databar_calc_check(unsigned sig0, unsigned sig1,
					 unsigned side, unsigned mod);

// Compatibility wrapper
static inline unsigned calc_check(unsigned sig0, unsigned sig1, unsigned side,
				  unsigned mod)
{
    return _zbar_databar_calc_check(sig0, sig1, side, mod);
}

extern int calc_value4(unsigned sig, unsigned n, unsigned wmax,
		       unsigned nonarrow);

static zbar_symbol_type_t decode_char(zbar_decoder_t *dcode,
				      databar_segment_t *seg, int off, int dir)
{
    databar_decoder_t *db = &dcode->databar;
    unsigned s		  = calc_s(dcode, (dir > 0) ? off : off - 6, 8);
    int n, i, emin[2] = { 0, }, sum = 0;
    unsigned sig0 = 0, sig1 = 0;
    int diff, vodd, veven, v;
    unsigned sum0, sum1, chk;
    struct group_s *g;

    if (seg->exp)
	n = 17;
    else if (seg->side)
	n = 15;
    else
	n = 16;
    emin[1] = -n;

    if (s < 13 || !check_width(seg->width, s, n))
	return (ZBAR_NONE);

    for (i = 4; --i >= 0;) {
	int e = decode_e(pair_width(dcode, off), s, n);
	if (e < 0)
	    return (ZBAR_NONE);
	sum = e - sum;
	off += dir;
	sig1 <<= 4;
	if (emin[1] < -sum)
	    emin[1] = -sum;
	sig1 += sum;
	if (!i)
	    break;

	e = decode_e(pair_width(dcode, off), s, n);
	if (e < 0)
	    return (ZBAR_NONE);
	sum = e - sum;
	off += dir;
	sig0 <<= 4;
	if (emin[0] > sum)
	    emin[0] = sum;
	sig0 += sum;
    }

    diff = emin[~n & 1];
    diff = diff + (diff << 4);
    diff = diff + (diff << 8);

    sig0 -= diff;
    sig1 += diff;

    sum0 = sig0 + (sig0 >> 8);
    sum1 = sig1 + (sig1 >> 8);
    sum0 += sum0 >> 4;
    sum1 += sum1 >> 4;
    sum0 &= 0xf;
    sum1 &= 0xf;

    if ((int)(sum0 + sum1 + 8) != (int)n) {
	return (ZBAR_NONE);
    }

    if (((sum0 ^ (n >> 1)) | (sum1 ^ (n >> 1) ^ n)) & 1) {
	return (ZBAR_NONE);
    }

    i = ((n & 0x3) ^ 1) * 5 + (sum1 >> 1);
    zassert((size_t)i < sizeof(groups) / sizeof(*groups), -1,
	    "n=%d sum=%d/%d sig=%04x/%04x g=%d", n, sum0, sum1, sig0, sig1, i);
    g = groups + i;

    vodd = calc_value4(sig0 + 0x1111, sum0 + 4, g->wmax, ~n & 1);
    if (vodd < 0 || vodd > g->todd)
	return (ZBAR_NONE);

    veven = calc_value4(sig1 + 0x1111, sum1 + 4, 9 - g->wmax, n & 1);
    if (veven < 0 || veven > g->teven)
	return (ZBAR_NONE);

    v = g->sum;
    if (n & 2)
	v += vodd + veven * g->todd;
    else
	v += veven + vodd * g->teven;

    chk = 0;
    if (seg->exp) {
	unsigned side = seg->color ^ seg->side ^ 1;
	if (v >= 4096)
	    return (ZBAR_NONE);
	/* skip A1 left */
	chk = calc_check(sig0, sig1, side, 211);
	if (seg->finder || seg->color || seg->side) {
	    i = (seg->finder << 1) - side + seg->color;
	    zassert(i >= 0 && i < 12, ZBAR_NONE, "f=%d(%x%x%x) side=%d i=%d\n",
		    seg->finder, seg->exp, seg->color, seg->side, side, i);
	    chk = (chk * exp_checksums[i]) % 211;
	} else if (v >= 4009)
	    return (ZBAR_NONE);
	else
	    chk = 0;
    } else {
	chk = calc_check(sig0, sig1, seg->side, 79);
	if (seg->color)
	    chk = (chk * 16) % 79;
    }

    seg->check = chk;
    seg->data  = v;

    _zbar_databar_merge_segment(db, seg);

    if (seg->exp)
	return (match_segment_exp(dcode, seg, dir));
    else if (dir > 0)
	return (match_segment(dcode, seg));
    return (ZBAR_PARTIAL);
}

/* Converted to Rust - see src/decoder.rs */
extern int _zbar_databar_alloc_segment(databar_decoder_t *db);

static inline int alloc_segment(databar_decoder_t *db)
{
    return _zbar_databar_alloc_segment(db);
}

extern zbar_symbol_type_t decode_finder(zbar_decoder_t *dcode);

zbar_symbol_type_t _zbar_decode_databar(zbar_decoder_t *dcode)
{
    databar_decoder_t *db = &dcode->databar;
    databar_segment_t *seg, *pair;
    zbar_symbol_type_t sym;
    int iseg, i = dcode->idx & 0xf;

    sym = decode_finder(dcode);

    iseg = db->chars[i];
    if (iseg < 0)
	return (sym);

    db->chars[i] = -1;
    seg		 = db->segs + iseg;
    zassert(seg->finder >= 0, ZBAR_NONE, "i=%d f=%d(%x%x%x) part=%x\n", iseg,
	    seg->finder, seg->exp, seg->color, seg->side, seg->partial);

    if (seg->partial) {
	pair	  = NULL;
	seg->side = !seg->side;
    } else {
	int jseg     = alloc_segment(db);
	pair	     = db->segs + iseg;
	seg	     = db->segs + jseg;
	seg->finder  = pair->finder;
	seg->exp     = pair->exp;
	seg->color   = pair->color;
	seg->side    = !pair->side;
	seg->partial = 0;
	seg->count   = 1;
	seg->width   = pair->width;
	seg->epoch   = db->epoch;
    }

    sym = decode_char(dcode, seg, 1, 1);
    if (!sym) {
	seg->finder = -1;
	if (pair)
	    pair->partial = 1;
    } else
	db->epoch++;

    return (sym);
}
