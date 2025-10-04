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
#ifndef _ERROR_H_
#define _ERROR_H_

#include "config.h"
#include <assert.h>
#include <errno.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <zbar.h>

#if __STDC_VERSION__ < 199901L
#if __GNUC__ >= 2
#define __func__ __FUNCTION__
#else
#define __func__ "<unknown>"
#endif
#endif

#define ERRINFO_MAGIC (0x5252457a) /* "zERR" (LE) */

typedef enum errsev_e {
    SEV_FATAL	= -2, /* application must terminate */
    SEV_ERROR	= -1, /* might be able to recover and continue */
    SEV_OK	= 0,
    SEV_WARNING = 1, /* unexpected condition */
    SEV_NOTE	= 2, /* fyi */
} errsev_t;

typedef enum errmodule_e {
    ZBAR_MOD_PROCESSOR,
    ZBAR_MOD_VIDEO,
    ZBAR_MOD_WINDOW,
    ZBAR_MOD_IMAGE_SCANNER,
    ZBAR_MOD_UNKNOWN,
} errmodule_t;

typedef struct errinfo_s {
    uint32_t magic;	/* just in case */
    errmodule_t module; /* reporting module */
    char *buf;		/* formatted and passed to application */
    int errnum;		/* errno for system errors */

    errsev_t sev;
    zbar_error_t type;
    const char *func;	/* reporting function */
    const char *detail; /* description */
    char *arg_str;	/* single string argument */
    int arg_int;	/* single integer argument */
} errinfo_t;

extern int _zbar_verbosity;

/* FIXME don't we need varargs hacks here? */

#define ZFLUSH

#define zprintf(level, format, args...)                       \
    do {                                                      \
	if (_zbar_verbosity >= level) {                       \
	    fprintf(stderr, "%s: " format, __func__, ##args); \
	    ZFLUSH                                            \
	}                                                     \
    } while (0)
#define zwprintf(level, format, args...)       \
    do {                                       \
	if (_zbar_verbosity >= level) {        \
	    fprintf(stderr, "%s: ", __func__); \
	    fwprintf(stderr, format, ##args);  \
	    ZFLUSH                             \
	}                                      \
    } while (0)

extern int _zbar_err_copy(void *dst_c, void *src_c);
extern int _zbar_err_capture(const void *container, errsev_t sev, zbar_error_t type, const char *func, const char *detail);
extern int _zbar_err_capture_str(const void *container, errsev_t sev, zbar_error_t type, const char *func, const char *detail, const char *arg);
extern int _zbar_err_capture_int(const void *container, errsev_t sev, zbar_error_t type, const char *func, const char *detail, int arg);
extern int _zbar_err_capture_num(const void *container, errsev_t sev, zbar_error_t type, const char *func, const char *detail, int num);
extern void _zbar_err_init(errinfo_t *err, errmodule_t module);
extern void _zbar_err_cleanup(errinfo_t *err);

#define err_copy(dst_c, src_c) _zbar_err_copy(dst_c, src_c)
#define err_capture(container, sev, type, func, detail) _zbar_err_capture(container, sev, type, func, detail)
#define err_capture_str(container, sev, type, func, detail, arg) _zbar_err_capture_str(container, sev, type, func, detail, arg)
#define err_capture_int(container, sev, type, func, detail, arg) _zbar_err_capture_int(container, sev, type, func, detail, arg)
#define err_capture_num(container, sev, type, func, detail, num) _zbar_err_capture_num(container, sev, type, func, detail, num)
#define err_init(err, module) _zbar_err_init(err, module)
#define err_cleanup(err) _zbar_err_cleanup(err)

#endif
