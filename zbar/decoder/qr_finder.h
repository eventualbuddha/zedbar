#ifndef _DECODER_QR_FINDER_H_
#define _DECODER_QR_FINDER_H_

#include "qrcode.h"

/* QR Code symbol finder state */
typedef struct qr_finder_s {
    unsigned s5;	 /* finder pattern width */
    qr_finder_line line; /* position info needed by decoder */

    unsigned config;
} qr_finder_t;

/* reset QR finder specific state */
extern void _zbar_qr_finder_reset(qr_finder_t *qrf);
#define qr_finder_reset(qrf) _zbar_qr_finder_reset(qrf)

/* find QR Code symbols */
zbar_symbol_type_t _zbar_find_qr(zbar_decoder_t *dcode);

#endif
