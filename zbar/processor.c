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

#include "processor.h"
#include "image.h"
#include "img_scanner.h"

static inline int proc_enter(zbar_processor_t *proc) {
  _zbar_mutex_lock(&proc->mutex);
  return (_zbar_processor_lock(proc));
}

static inline int proc_leave(zbar_processor_t *proc) {
  int rc = _zbar_processor_unlock(proc, 0);
  _zbar_mutex_unlock(&proc->mutex);
  return (rc);
}

static inline int proc_open(zbar_processor_t *proc) {
  /* arbitrary default */
  int width = 640, height = 480;
  return (_zbar_processor_open(proc, "zbar barcode reader", width, height));
}

/* API lock is already held */
int _zbar_process_image(zbar_processor_t *proc, zbar_image_t *img) {
  if (img) {
    uint32_t format;
    zbar_image_t *tmp;
    int nsyms;

    format = zbar_image_get_format(img);
    zprintf(16, "processing: %.4s(%08x) %dx%d @%p\n", (char *)&format, format,
            zbar_image_get_width(img), zbar_image_get_height(img),
            zbar_image_get_data(img));

    /* FIXME locking all other interfaces while processing is conservative
     * but easier for now and we don't expect this to take long...
     */
    tmp = zbar_image_convert(img, fourcc('Y', '8', '0', '0'));
    if (!tmp)
      goto error;

    if (proc->syms) {
      zbar_symbol_set_ref(proc->syms, -1);
      proc->syms = NULL;
    }
    zbar_image_scanner_recycle_image(proc->scanner, img);
    nsyms = zbar_scan_image(proc->scanner, tmp);
    _zbar_image_swap_symbols(img, tmp);

    zbar_image_destroy(tmp);
    tmp = NULL;
    if (nsyms < 0)
      goto error;

    proc->syms = zbar_image_scanner_get_results(proc->scanner);
    if (proc->syms)
      zbar_symbol_set_ref(proc->syms, 1);

    if (_zbar_verbosity >= 8) {
      const zbar_symbol_t *sym = zbar_image_first_symbol(img);
      while (sym) {
        zbar_symbol_type_t type = zbar_symbol_get_type(sym);
        int count = zbar_symbol_get_count(sym);
        zprintf(8, "%s: %s (%d pts) (dir=%d) (q=%d) (%s)\n",
                zbar_get_symbol_name(type), zbar_symbol_get_data(sym),
                zbar_symbol_get_loc_size(sym), zbar_symbol_get_orientation(sym),
                zbar_symbol_get_quality(sym),
                (count < 0)   ? "uncertain"
                : (count > 0) ? "duplicate"
                              : "new");
        sym = zbar_symbol_next(sym);
      }
    }

    if (nsyms) {
      /* FIXME only call after filtering */
      _zbar_mutex_lock(&proc->mutex);
      _zbar_processor_notify(proc, EVENT_OUTPUT);
      _zbar_mutex_unlock(&proc->mutex);
    }
  }

  return 0;

error:
  return (err_capture(proc, SEV_ERROR, ZBAR_ERR_UNSUPPORTED, __func__,
                      "unknown image format"));
}

zbar_processor_t *zbar_processor_create(int threaded) {
  zbar_processor_t *proc = calloc(1, sizeof(zbar_processor_t));
  if (!proc)
    return (NULL);
  err_init(&proc->err, ZBAR_MOD_PROCESSOR);

  proc->scanner = zbar_image_scanner_create();
  if (!proc->scanner) {
    free(proc);
    return (NULL);
  }

  proc->threaded = !_zbar_mutex_init(&proc->mutex) && threaded;
  _zbar_processor_init(proc);
  return (proc);
}

void zbar_processor_destroy(zbar_processor_t *proc) {
  proc_waiter_t *w, *next;

  zbar_processor_init(proc, NULL, 0);

  if (proc->syms) {
    zbar_symbol_set_ref(proc->syms, -1);
    proc->syms = NULL;
  }
  if (proc->scanner) {
    zbar_image_scanner_destroy(proc->scanner);
    proc->scanner = NULL;
  }

  _zbar_mutex_destroy(&proc->mutex);
  _zbar_processor_cleanup(proc);

  assert(!proc->wait_head);
  assert(!proc->wait_tail);
  assert(!proc->wait_next);

  for (w = proc->free_waiter; w; w = next) {
    next = w->next;
    _zbar_event_destroy(&w->notify);
    free(w);
  }

  err_cleanup(&proc->err);
  free(proc);
}

int zbar_processor_init(zbar_processor_t *proc, const char *dev,
                        int enable_display) {
  int rc = 0;

  _zbar_mutex_lock(&proc->mutex);
  _zbar_thread_stop(&proc->input_thread, &proc->mutex);

  _zbar_processor_lock(proc);
  _zbar_mutex_unlock(&proc->mutex);

  rc = 0;

  _zbar_mutex_lock(&proc->mutex);
  proc_leave(proc);
  return (rc);
}

int zbar_processor_set_config(zbar_processor_t *proc, zbar_symbol_type_t sym,
                              zbar_config_t cfg, int val) {
  int rc;
  proc_enter(proc);
  rc = zbar_image_scanner_set_config(proc->scanner, sym, cfg, val);
  proc_leave(proc);
  return (rc);
}

int zbar_processor_set_control(zbar_processor_t *proc, const char *control_name,
                               int value) {
  return 0;
}

int zbar_processor_get_control(zbar_processor_t *proc, const char *control_name,
                               int *value) {
  return 0;
}

int zbar_processor_request_size(zbar_processor_t *proc, unsigned width,
                                unsigned height) {
  return 0;
}

int zbar_processor_is_visible(zbar_processor_t *proc) { return 0; }

int zbar_processor_set_visible(zbar_processor_t *proc, int visible) {
  return 0;
}

const zbar_symbol_set_t *
zbar_processor_get_results(const zbar_processor_t *proc) {
  const zbar_symbol_set_t *syms;
  zbar_processor_t *ncproc = (zbar_processor_t *)proc;
  proc_enter(ncproc);
  syms = proc->syms;
  if (syms)
    zbar_symbol_set_ref(syms, 1);
  proc_leave(ncproc);
  return (syms);
}

int zbar_processor_user_wait(zbar_processor_t *proc, int timeout) {
  int rc = -1;

  proc_enter(proc);
  _zbar_mutex_unlock(&proc->mutex);

  rc = err_capture(proc, SEV_WARNING, ZBAR_ERR_CLOSED, __func__,
                   "display window not available for input");

  _zbar_mutex_lock(&proc->mutex);
  proc_leave(proc);
  return (rc);
}

int zbar_processor_set_active(zbar_processor_t *proc, int active) {
  int rc;
  proc_enter(proc);
  rc = err_capture(proc, SEV_ERROR, ZBAR_ERR_INVALID, __func__,
                   "input not initialized");
  proc_leave(proc);
  return (rc);
}

int zbar_process_one(zbar_processor_t *proc, int timeout) {
  int rc;

  proc_enter(proc);
  _zbar_mutex_unlock(&proc->mutex);

  rc = err_capture(proc, SEV_ERROR, ZBAR_ERR_INVALID, __func__,
                   "input not initialized");
  _zbar_mutex_lock(&proc->mutex);
  proc_leave(proc);
  return (rc);
}

int zbar_process_image(zbar_processor_t *proc, zbar_image_t *img) {
  int rc = 0;

  proc_enter(proc);
  _zbar_mutex_unlock(&proc->mutex);

  zbar_image_scanner_enable_cache(proc->scanner, 0);
  rc = _zbar_process_image(proc, img);

  _zbar_mutex_lock(&proc->mutex);
  proc_leave(proc);
  return (rc);
}
