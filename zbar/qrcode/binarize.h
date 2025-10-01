/*Copyright (C) 2008-2009  Timothy B. Terriberry (tterribe@xiph.org)
  You can redistribute this library and/or modify it under the terms of the
   GNU Lesser General Public License as published by the Free Software
   Foundation; either version 2.1 of the License, or (at your option) any later
   version.

  NOTE: This is a stub header for the Rust implementation of QR code binarization.
  The actual implementation is in src/qrcode/binarize.rs */
#if !defined(_qrcode_binarize_H)
#define _qrcode_binarize_H (1)

#ifdef __cplusplus
extern "C" {
#endif

/* These functions are not currently implemented */
void qr_image_cross_masking_median_filter(unsigned char *_img, int _width,
					  int _height);

void qr_wiener_filter(unsigned char *_img, int _width, int _height);

/* Binarizes a grayscale image using adaptive thresholding.
   Returns a binary mask where 0xFF represents foreground (black) and
   0x00 represents background (white).
   The returned pointer must be freed by the caller using free().
   Returns NULL if width or height is <= 0.
   Implemented in Rust (src/qrcode/binarize.rs) */
unsigned char *qr_binarize(const unsigned char *_img, int _width, int _height);

#ifdef __cplusplus
}
#endif

#endif
