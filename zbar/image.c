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

#include "image.h"
#include <string.h>

// _zbar_image_refcnt is implemented in Rust (src/image_ffi.rs)

void _zbar_image_swap_symbols(zbar_image_t *a, zbar_image_t *b)
{
    zbar_symbol_set_t *tmp = a->syms;
    a->syms                = b->syms;
    b->syms                = tmp;
}

void _zbar_image_copy_size(zbar_image_t *dst, const zbar_image_t *src)
{
    dst->width  = src->width;
    dst->height = src->height;
}

zbar_image_t *_zbar_image_copy(const zbar_image_t *src, int inverted)
{
    zbar_image_t *dst;

    if (inverted && (src->format != fourcc('Y', '8', '0', '0')) &&
        (src->format != fourcc('G', 'R', 'E', 'Y')))
        return NULL;

    dst         = zbar_image_create();
    dst->format = src->format;
    _zbar_image_copy_size(dst, src);
    dst->datalen = src->datalen;
    dst->data    = malloc(src->datalen);
    assert(dst->data);

    if (!inverted) {
        memcpy((void *)dst->data, src->data, src->datalen);
    } else {
        int i, len = src->datalen;
        long *sp = (void *)src->data, *dp = (void *)dst->data;
        char *spc, *dpc;

        /* Do it word per word, in order to speedup */
        for (i = 0; i < len; i += sizeof(long))
            *dp++ = ~(*sp++);

        /* Deal with non-aligned remains, if any */
        len -= i;
        spc = (char *)sp;
        dpc = (char *)dp;
        for (i = 0; i < len; i++)
            *dpc++ = ~(*spc++);
    }
    dst->cleanup = zbar_image_free_data;
    return (dst);
}
