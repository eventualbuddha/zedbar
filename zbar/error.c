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

#include "error.h"
#include <errno.h>
#include <stdlib.h>
#include <string.h>

int _zbar_err_copy(void *dst_c, void *src_c)
{
    errinfo_t *dst = dst_c;
    errinfo_t *src = src_c;
    assert(dst->magic == ERRINFO_MAGIC);
    assert(src->magic == ERRINFO_MAGIC);

    dst->errnum  = src->errnum;
    dst->sev     = src->sev;
    dst->type    = src->type;
    dst->func    = src->func;
    dst->detail  = src->detail;
    dst->arg_str = src->arg_str;
    src->arg_str = NULL; /* unused at src, avoid double free */
    dst->arg_int = src->arg_int;
    return (-1);
}

int _zbar_err_capture(const void *container, errsev_t sev,
                      zbar_error_t type, const char *func,
                      const char *detail)
{
    errinfo_t *err = (errinfo_t *)container;
    assert(err->magic == ERRINFO_MAGIC);
    if (type == ZBAR_ERR_SYSTEM)
        err->errnum = errno;
    err->sev     = sev;
    err->type    = type;
    err->func    = func;
    err->detail  = detail;
    if (_zbar_verbosity >= 1)
        _zbar_error_spew(err, 0);
    return (-1);
}

int _zbar_err_capture_str(const void *container, errsev_t sev,
                          zbar_error_t type, const char *func,
                          const char *detail, const char *arg)
{
    errinfo_t *err = (errinfo_t *)container;
    assert(err->magic == ERRINFO_MAGIC);
    if (err->arg_str)
        free(err->arg_str);
    err->arg_str = strdup(arg);
    return (_zbar_err_capture(container, sev, type, func, detail));
}

int _zbar_err_capture_int(const void *container, errsev_t sev,
                          zbar_error_t type, const char *func,
                          const char *detail, int arg)
{
    errinfo_t *err = (errinfo_t *)container;
    assert(err->magic == ERRINFO_MAGIC);
    err->arg_int = arg;
    return (_zbar_err_capture(container, sev, type, func, detail));
}

int _zbar_err_capture_num(const void *container, errsev_t sev,
                          zbar_error_t type, const char *func,
                          const char *detail, int num)
{
    errinfo_t *err = (errinfo_t *)container;
    assert(err->magic == ERRINFO_MAGIC);
    err->errnum = num;
    return (_zbar_err_capture(container, sev, type, func, detail));
}

void _zbar_err_init(errinfo_t *err, errmodule_t module)
{
    err->magic  = ERRINFO_MAGIC;
    err->module = module;
}

void _zbar_err_cleanup(errinfo_t *err)
{
    assert(err->magic == ERRINFO_MAGIC);
    if (err->buf) {
        free(err->buf);
        err->buf = NULL;
    }
    if (err->arg_str) {
        free(err->arg_str);
        err->arg_str = NULL;
    }
}
