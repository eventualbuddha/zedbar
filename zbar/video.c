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

#include "video.h"
#include "image.h"

#ifdef HAVE_LIBJPEG
extern struct jpeg_decompress_struct *_zbar_jpeg_decomp_create(void);
extern void _zbar_jpeg_decomp_destroy(struct jpeg_decompress_struct *cinfo);
#endif

zbar_video_t *zbar_video_create() { return NULL; }

void zbar_video_destroy(zbar_video_t *vdo) {}

int zbar_video_open(zbar_video_t *vdo, const char *dev) { return 0; }

int zbar_video_get_fd(const zbar_video_t *vdo) { return 0; }

int zbar_video_request_size(zbar_video_t *vdo, unsigned width,
                            unsigned height) {
  return 0;
}

int zbar_video_request_interface(zbar_video_t *vdo, int ver) { return 0; }

int zbar_video_request_iomode(zbar_video_t *vdo, int iomode) { return 0; }

int zbar_video_get_width(const zbar_video_t *vdo) { return (vdo->width); }

int zbar_video_get_height(const zbar_video_t *vdo) { return (vdo->height); }

uint32_t zbar_video_get_format(const zbar_video_t *vdo) {
  return (vdo->format);
}

static inline int video_init_images(zbar_video_t *vdo) { return 0; }

int zbar_video_init(zbar_video_t *vdo, unsigned long fmt) { return 0; }

int zbar_video_enable(zbar_video_t *vdo, int enable) { return 0; }

zbar_image_t *zbar_video_next_image(zbar_video_t *vdo) { return NULL; }

/** @brief return if fun unsupported, otherwise continue */
#define return_if_not_supported(fun, name)                                     \
  {                                                                            \
    if (!(fun)) {                                                              \
      zprintf(1, "video driver does not implement %s\n", name);                \
      return ZBAR_ERR_UNSUPPORTED;                                             \
    }                                                                          \
  }
#define return_if_non_zero(a)                                                  \
  {                                                                            \
    int rv = a;                                                                \
    if (rv != 0)                                                               \
      return (rv);                                                             \
  }

int zbar_video_set_control(zbar_video_t *vdo, const char *control_name,
                           int value) {
  return 0;
}

int zbar_video_get_control(zbar_video_t *vdo, const char *control_name,
                           int *value) {
  return 0;
}

struct video_controls_s *zbar_video_get_controls(const zbar_video_t *vdo,
                                                 int index) {
  int i = 0;
  struct video_controls_s *p = vdo->controls;

  while (p && i != index) {
    i++;
    p = p->next;
  }

  if (!p)
    return NULL;

  return p;
}

struct video_resolution_s *zbar_video_get_resolutions(const zbar_video_t *vdo,
                                                      int index) {
  int i = 0;
  struct video_resolution_s *p = vdo->res;

  while (i != index) {
    if (!p->width || !p->height)
      return NULL;
    i++;
    p++;
  }

  if (!p->width || !p->height)
    return NULL;

  return p;
}
