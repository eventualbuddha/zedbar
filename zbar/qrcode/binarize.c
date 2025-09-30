/*Copyright (C) 2008-2009  Timothy B. Terriberry (tterribe@xiph.org)
  You can redistribute this library and/or modify it under the terms of the
   GNU Lesser General Public License as published by the Free Software
   Foundation; either version 2.1 of the License, or (at your option) any later
   version.*/
#include "binarize.h"
#include "image.h"
#include "util.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/*The above algorithms are computationally expensive, and do not work as well
   as the simple algorithm below.
  Sauvola by itself does an excellent job of classifying regions outside the
   QR code as background, which greatly reduces the chance of false alarms.
  However, it also tends to over-shrink isolated black dots inside the code,
   making them easy to miss with even slight mis-alignment.
  Since the Gatos method uses Sauvola as input to its background interpolation
   method, it cannot possibly mark any pixels as foreground which Sauvola
   classified as background, and thus suffers from the same problem.
  The following simple adaptive threshold method does not have this problem,
   though it produces essentially random noise outside the QR code region.
  QR codes are structured well enough that this does not seem to lead to any
   actual false alarms in practice, and it allows many more codes to be
   detected and decoded successfully than the Sauvola or Gatos binarization
   methods.*/

/*A simplified adaptive thresholder.
  This compares the current pixel value to the mean value of a (large) window
   surrounding it.*/
unsigned char *qr_binarize(const unsigned char *_img, int _width, int _height) {
  unsigned char *mask = NULL;
  if (_width > 0 && _height > 0) {
    unsigned *col_sums;
    int logwindw;
    int logwindh;
    int windw;
    int windh;
    int y0offs;
    int y1offs;
    unsigned g;
    int x;
    int y;
    mask = (unsigned char *)malloc(_width * _height * sizeof(*mask));
    /*We keep the window size fairly large to ensure it doesn't fit completely
   inside the center of a finder pattern of a version 1 QR code at full
   resolution.*/
    for (logwindw = 4; logwindw < 8 && (1 << logwindw) < (_width + 7 >> 3);
         logwindw++)
      ;
    for (logwindh = 4; logwindh < 8 && (1 << logwindh) < (_height + 7 >> 3);
         logwindh++)
      ;
    windw = 1 << logwindw;
    windh = 1 << logwindh;
    col_sums = (unsigned *)malloc(_width * sizeof(*col_sums));
    /*Initialize sums down each column.*/
    for (x = 0; x < _width; x++) {
      g = _img[x];
      col_sums[x] = (g << logwindh - 1) + g;
    }
    for (y = 1; y < (windh >> 1); y++) {
      y1offs = QR_MINI(y, _height - 1) * _width;
      for (x = 0; x < _width; x++) {
        g = _img[y1offs + x];
        col_sums[x] += g;
      }
    }
    for (y = 0; y < _height; y++) {
      unsigned m;
      int x0;
      int x1;
      /*Initialize the sum over the window.*/
      m = (col_sums[0] << logwindw - 1) + col_sums[0];
      for (x = 1; x < (windw >> 1); x++) {
        x1 = QR_MINI(x, _width - 1);
        m += col_sums[x1];
      }
      for (x = 0; x < _width; x++) {
        /*Perform the test against the threshold T = (m/n)-D,
   where n=windw*windh and D=3.*/
        g = _img[y * _width + x];
        mask[y * _width + x] = -(g + 3 << logwindw + logwindh < m) & 0xFF;
        /*Update the window sum.*/
        if (x + 1 < _width) {
          x0 = QR_MAXI(0, x - (windw >> 1));
          x1 = QR_MINI(x + (windw >> 1), _width - 1);
          m += col_sums[x1] - col_sums[x0];
        }
      }
      /*Update the column sums.*/
      if (y + 1 < _height) {
        y0offs = QR_MAXI(0, y - (windh >> 1)) * _width;
        y1offs = QR_MINI(y + (windh >> 1), _height - 1) * _width;
        for (x = 0; x < _width; x++) {
          col_sums[x] -= _img[y0offs + x];
          col_sums[x] += _img[y1offs + x];
        }
      }
    }
    free(col_sums);
  }
  return (mask);
}
