bin_PROGRAMS += zbarimg/zbarimg
zbarimg_zbarimg_SOURCES = zbarimg/zbarimg.c
zbarimg_zbarimg_CPPFLAGS = $(AM_CPPFLAGS)
zbarimg_zbarimg_LDADD = zbar/libzbar.la
# automake bug in "monolithic mode"?
CLEANFILES += zbarimg/.libs/zbarimg
