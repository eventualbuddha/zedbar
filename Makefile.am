bin_PROGRAMS =
check_PROGRAMS =
EXTRA_PROGRAMS =
lib_LTLIBRARIES =
CLEANFILES =
MAINTAINERCLEANFILES =
BUILT_SOURCES =
EXTRA_DIST =
PHONY = $(SUBDIRS)

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = zbar.pc

dist_doc_DATA = COPYING HACKING.md INSTALL.md LICENSE.md NEWS.md README.md TODO.md

include $(srcdir)/include/Makefile.am.inc

SUBDIRS = zbar

zbar/libzbar.la:
	$(MAKE) -C @builddir@/zbar libzbar.la

include $(srcdir)/test/Makefile.am.inc
if HAVE_DOC
include $(srcdir)/doc/Makefile.am.inc
endif

EXTRA_DIST += zbar.ico

EXTRA_DIST += examples/*.png examples/sha1sum \
    examples/scan_image.c examples/scan_image.cpp examples/scan_image.vcproj

SUBDIRS += .

archive:
	git archive --format=tar --prefix=zbar-$(VERSION)/ -o zbar-$(VERSION).tar $(VERSION)
	bzip2 zbar-$(VERSION).tar

.PHONY : $(PHONY) archive

dist-hook:
	rm -f $(distdir)/debian $(distdir)/travis
