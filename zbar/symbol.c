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

#include "symbol.h"

void _zbar_symbol_add_point(zbar_symbol_t *sym, int x, int y)
{
    int i = sym->npts;
    if (++sym->npts >= sym->pts_alloc)
        sym->pts = realloc(sym->pts, ++sym->pts_alloc * sizeof(point_t));
    sym->pts[i].x = x;
    sym->pts[i].y = y;
}

void _zbar_symbol_refcnt(zbar_symbol_t *sym, int delta)
{
    if (!_zbar_refcnt(&sym->refcnt, delta) && delta <= 0)
        _zbar_symbol_free(sym);
}

void _zbar_symbol_set_add(zbar_symbol_set_t *syms, zbar_symbol_t *sym)
{
    sym->next  = syms->head;
    syms->head = sym;
    syms->nsyms++;

    _zbar_symbol_refcnt(sym, 1);
}
