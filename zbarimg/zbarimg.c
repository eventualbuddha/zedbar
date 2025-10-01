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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <jpeglib.h>
#include <png.h>
#include "zbar.h"

static const char *note_usage =
    ("usage: zbarimg [options] <image>...\n"
     "\n"
     "scan and decode bar codes from one or more image files\n"
     "\n"
     "options:\n"
     "    -h, --help      display this help text\n"
     "    --version       display version information and exit\n"
     "symbol data\n"
     "    -v, --verbose   increase debug output level\n"
     "    --verbose=N     set specific debug output level\n"
     // FIXME overlay level
     "\n");

static int notfound = 0, exit_code = 0;
static int num_images = 0;

static zbar_processor_t *processor = NULL;

int has_extension(const char *filename, const char *extension)
{
    char *dot = strrchr(filename, '.');
    if (dot && dot != filename) {
	return strcmp(dot, extension) == 0;
    }
    return 0;
}

static int load_image(const char *filename, zbar_image_t *zimage)
{
    if (exit_code == 3)
	return (-1);

    if (has_extension(filename, ".jpeg") || has_extension(filename, ".jpg")) {
	struct jpeg_decompress_struct cinfo;
	struct jpeg_error_mgr jerr;
	FILE *infile;
	int row_stride;

	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_decompress(&cinfo);

	if ((infile = fopen(filename, "rb")) == NULL) {
	    fprintf(stderr, "Can't open %s\n", filename);
	    return -1;
	}

	jpeg_stdio_src(&cinfo, infile);
	jpeg_read_header(&cinfo, TRUE);

	// Set the output color space to grayscale
	cinfo.out_color_space = JCS_GRAYSCALE;
	jpeg_start_decompress(&cinfo);

	// Calculate row stride for allocating buffer
	row_stride = cinfo.output_width * cinfo.output_components;

	JDIMENSION width    = cinfo.output_width;
	JDIMENSION height   = cinfo.output_height;
	int pixel_size	    = cinfo.output_components;
	size_t bloblen	    = width * height * pixel_size;
	unsigned char *blob = malloc(bloblen);
	zbar_image_set_format(zimage, zbar_fourcc('Y', '8', '0', '0'));
	zbar_image_set_size(zimage, width, height);
	zbar_image_set_data(zimage, blob, bloblen, zbar_image_free_data);

	while (cinfo.output_scanline < cinfo.output_height) {
	    unsigned char *buffer_array[1];
	    buffer_array[0] = blob + (cinfo.output_scanline) * row_stride;

	    jpeg_read_scanlines(&cinfo, buffer_array, 1);
	}

	// Finish decompression and release resources
	jpeg_finish_decompress(&cinfo);
	jpeg_destroy_decompress(&cinfo);
	fclose(infile);

	return 0;
    }

    if (has_extension(filename, ".png")) {
	png_structp png =
	    png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

	if (!png) {
	    return -1;
	}

	png_infop info = png_create_info_struct(png);

	if (!info) {
	    png_destroy_read_struct(&png, NULL, NULL);
	    return -1;
	}

	if (setjmp(png_jmpbuf(png))) {
	    png_destroy_read_struct(&png, &info, NULL);
	    return -1;
	}

	FILE *fp = fopen(filename, "rb");

	png_init_io(png, fp);
	png_read_info(png, info);

	int width	    = png_get_image_width(png, info);
	int height	    = png_get_image_height(png, info);
	png_byte color_type = png_get_color_type(png, info);
	png_byte bit_depth  = png_get_bit_depth(png, info);

	if (bit_depth == 16)
	    png_set_strip_16(png);

	png_set_strip_alpha(png);

	if (color_type == PNG_COLOR_TYPE_PALETTE)
	    png_set_palette_to_rgb(png);

	if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8)
	    png_set_expand_gray_1_2_4_to_8(png);

	// Convert RGB to grayscale
	if (color_type & PNG_COLOR_MASK_COLOR)
	    png_set_rgb_to_gray(png, 1, -1, -1);

	// Update info after transformations
	png_read_update_info(png, info);

	size_t bloblen	    = width * height;
	unsigned char *blob = malloc(bloblen);
	zbar_image_set_format(zimage, zbar_fourcc('Y', '8', '0', '0'));
	zbar_image_set_size(zimage, width, height);
	zbar_image_set_data(zimage, blob, bloblen, zbar_image_free_data);

	png_bytep *row_pointers = malloc(height * sizeof(png_bytep));
	for (int i = 0; i < height; i += 1) {
	    row_pointers[i] = blob + i * width;
	}

	png_read_image(png, row_pointers);
	free(row_pointers);

	fclose(fp);

	png_destroy_read_struct(&png, &info, NULL);

	return 0;
    }

    return -1;
}

static int scan_image(const char *filename)
{
    if (exit_code == 3)
	return (-1);

    zbar_image_t *zimage = zbar_image_create();
    load_image(filename, zimage);
    zbar_process_image(processor, zimage);
    // output result data
    const zbar_symbol_t *sym = zbar_image_first_symbol(zimage);
    for (; sym; sym = zbar_symbol_next(sym)) {
	zbar_symbol_type_t typ = zbar_symbol_get_type(sym);
	unsigned len	       = zbar_symbol_get_data_length(sym);
	if (typ == ZBAR_PARTIAL)
	    continue;
	else {
	    if (len && fwrite(zbar_symbol_get_data(sym), len, 1, stdout) != 1) {
		exit_code = 1;
		return (-1);
	    }
	}
    }
    zbar_image_destroy(zimage);

    return 0;
}

int usage(int rc, const char *msg, const char *arg)
{
    FILE *out = (rc) ? stderr : stdout;
    if (msg) {
	fprintf(out, "%s", msg);
	if (arg)
	    fprintf(out, "%s", arg);
	fprintf(out, "\n\n");
    }
    fprintf(out, "%s", (note_usage));
    return (rc);
}

int main(int argc, const char *argv[])
{
    // option pre-scan
    int i, j;

    for (i = 1; i < argc; i++) {
	const char *arg = argv[i];
	if (arg[0] != '-' || !arg[1])
	    // first pass, skip images
	    num_images++;
	else if (arg[1] != '-')
	    for (j = 1; arg[j]; j++) {
		switch (arg[j]) {
		case 'h':
		    return (usage(0, NULL, NULL));
		case 'v':
		    zbar_increase_verbosity();
		    break;
		default:
		    return (
			usage(1, "ERROR: unknown bundled option: -", arg + j));
		}
	    }
	else if (!strcmp(arg, "--help"))
	    return (usage(0, NULL, NULL));
	else if (!strcmp(arg, "--version")) {
	    printf("%s\n", PACKAGE_VERSION);
	    return (0);
	} else if (!strcmp(arg, "--verbose"))
	    zbar_increase_verbosity();
	else if (!strncmp(arg, "--verbose=", 10))
	    zbar_set_verbosity(strtol(argv[i] + 10, NULL, 0));
	else if (!strcmp(arg, "--")) {
	    num_images += argc - i - 1;
	    break;
	} else
	    return (usage(1, "ERROR: unknown option: ", arg));
    }

    if (!num_images)
	return (usage(1, "ERROR: specify image file(s) to scan", NULL));
    num_images = 0;

    processor = zbar_processor_create(0);
    assert(processor);

    for (i = 1; i < argc; i++) {
	const char *arg = argv[i];
	if (!arg)
	    continue;

	if (arg[0] != '-' || !arg[1]) {
	    if (scan_image(arg))
		return (exit_code);
	} else if (!strcmp(arg, "--"))
	    break;
    }
    for (i++; i < argc; i++)
	if (scan_image(argv[i]))
	    return (exit_code);

    /* ignore quit during last image */
    if (exit_code == 3)
	exit_code = 0;

    if (num_images && notfound && !exit_code)
	exit_code = 4;

    zbar_processor_destroy(processor);
    return (exit_code);
}
