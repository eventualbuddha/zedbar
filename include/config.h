/* Simplified config.h for zbar - minimal configuration */
#ifndef CONFIG_H
#define CONFIG_H

/* Enable barcode symbologies */
#define ENABLE_QRCODE 1
#define ENABLE_CODE128 1
#define ENABLE_CODE39 1
#define ENABLE_CODE93 1
#define ENABLE_EAN 1
#define ENABLE_CODABAR 1
#define ENABLE_I25 1
#define ENABLE_DATABAR 1
#define ENABLE_SQCODE 1

/* System capabilities */
#define HAVE_STDLIB_H 1
#define HAVE_STRING_H 1
#define HAVE_UNISTD_H 1
#define HAVE_STDINT_H 1
#define HAVE_INTTYPES_H 1
#define HAVE_ERRNO_H 1
#define HAVE_FCNTL_H 1
#define HAVE_LIMITS_H 1
#define HAVE_POLL_H 1
#define HAVE_MALLOC 1
#define HAVE_MEMCHR 1
#define HAVE_MEMMOVE 1
#define HAVE_MEMSET 1
#define HAVE_STRCHR 1
#define HAVE_STRDUP 1
#define HAVE_STRRCHR 1
#define HAVE_STRSTR 1
#define HAVE_STRTOL 1
#define HAVE_STRTOUL 1
#define HAVE_GETTIMEOFDAY 1
#define HAVE_CLOCK_GETTIME 1
#define HAVE_PTHREAD 1

/* Image format support */
#define HAVE_LIBJPEG 1
#define HAVE_LIBPNG 1

/* Package information */
#define PACKAGE "zbar"
#define PACKAGE_NAME "zbar"
#define PACKAGE_VERSION "0.23.93"
#define VERSION "0.23.93"

/* Version components */
#define ZBAR_VERSION_MAJOR 0
#define ZBAR_VERSION_MINOR 23
#define ZBAR_VERSION_PATCH 93

#endif /* CONFIG_H */