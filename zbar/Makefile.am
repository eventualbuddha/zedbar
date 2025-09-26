lib_LTLIBRARIES = libzbar.la
libzbar_la_CPPFLAGS = $(AM_CPPFLAGS)
libzbar_la_LDFLAGS = -no-undefined -version-info $(LIB_VERSION) \
    -export-symbols-regex "^(zbar|_zbar.*_error)_.*" $(AM_LDFLAGS)
libzbar_la_LIBADD = $(LTLIBICONV)

libzbar_la_SOURCES = debug.h config.c \
    error.h error.c symbol.h symbol.c \
    image.h image.c convert.c \
    processor.c processor.h processor/lock.c \
    refcnt.h refcnt.c timer.h mutex.h \
    event.h thread.h \
    img_scanner.h img_scanner.c scanner.c \
    decoder.h decoder.c misc.h misc.c

if ENABLE_EAN
libzbar_la_SOURCES += decoder/ean.h decoder/ean.c
endif
if ENABLE_DATABAR
libzbar_la_SOURCES += decoder/databar.h decoder/databar.c
endif
if ENABLE_CODE128
libzbar_la_SOURCES += decoder/code128.h decoder/code128.c
endif
if ENABLE_CODE93
libzbar_la_SOURCES += decoder/code93.h decoder/code93.c
endif
if ENABLE_CODE39
libzbar_la_SOURCES += decoder/code39.h decoder/code39.c
endif
if ENABLE_CODABAR
libzbar_la_SOURCES += decoder/codabar.h decoder/codabar.c
endif
if ENABLE_I25
libzbar_la_SOURCES += decoder/i25.h decoder/i25.c
endif
if ENABLE_QRCODE
libzbar_la_SOURCES += qrcode.h \
    decoder/qr_finder.h decoder/qr_finder.c \
    qrcode/qrdec.h qrcode/qrdec.c qrcode/qrdectxt.c \
    qrcode/rs.h qrcode/rs.c \
    qrcode/isaac.h qrcode/isaac.c \
    qrcode/bch15_5.h qrcode/bch15_5.c \
    qrcode/binarize.h qrcode/binarize.c \
    qrcode/util.h qrcode/util.c
endif
if ENABLE_SQCODE
libzbar_la_SOURCES += sqcode.h sqcode.c \
    decoder/sq_finder.h decoder/sq_finder.c
endif

libzbar_la_SOURCES += processor/posix.h processor/posix.c

if HAVE_JPEG
libzbar_la_SOURCES += jpeg.c
endif

libzbar_la_SOURCES += processor/null.c
libzbar_la_LDFLAGS += $(AM_LDFLAGS)
libzbar_la_LIBADD += $(AM_LIBADD)
