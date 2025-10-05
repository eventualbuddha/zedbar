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

#include "config.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "processor.h"

// Rust implementations - converted to src/processor.rs
// These are now just exported from Rust, no C implementation needed

int zbar_process_image(zbar_processor_t *proc, zbar_image_t *img)
{
    if (!proc || !img)
	return -1;

    // Clean up previous results
    if (proc->syms) {
	zbar_symbol_set_ref(proc->syms, -1);
	proc->syms = NULL;
    }

    // Process the image
    zbar_image_scanner_recycle_image(proc->scanner, img);
    int nsyms = zbar_scan_image(proc->scanner, img);

    if (nsyms < 0)
	return nsyms;

    // Store results
    proc->syms = zbar_image_scanner_get_results(proc->scanner);
    if (proc->syms)
	zbar_symbol_set_ref(proc->syms, 1);

    return nsyms;
}