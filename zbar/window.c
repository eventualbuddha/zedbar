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

#include "window.h"
#include "image.h"
#include "timer.h"
#include <time.h> /* clock_gettime */

zbar_window_t *zbar_window_create() { return NULL; }

void zbar_window_destroy(zbar_window_t *w) {}

int zbar_window_attach(zbar_window_t *w, void *display,
                       unsigned long drawable) {
  return 0;
}

static inline int window_draw_overlay(zbar_window_t *w) { return 0; }

inline int zbar_window_redraw(zbar_window_t *w) { return 0; }

int zbar_window_draw(zbar_window_t *w, zbar_image_t *img) { return 0; }

void zbar_window_set_overlay(zbar_window_t *w, int lvl) {}

int zbar_window_get_overlay(const zbar_window_t *w) { return 0; }

int zbar_window_resize(zbar_window_t *w, unsigned width, unsigned height) {
  return 0;
}
