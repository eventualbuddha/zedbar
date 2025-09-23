ACLOCAL_AMFLAGS = -I config
AM_CFLAGS += -DLOCALEDIR=\"$(localedir)\"
AM_CXXFLAGS += -DLOCALEDIR=\"$(localedir)\"
bin_PROGRAMS =
check_PROGRAMS =
EXTRA_PROGRAMS =
lib_LTLIBRARIES =
CLEANFILES =
DISTCLEANFILES = \
	ABOUT-NLS \
	po/*.gmo \
	po/*.mo \
	po/Makefile.in.in \
	po/Makevars.template \
	po/POTFILES \
	po/Rules-quot \
	po/boldquot.sed \
	po/en@boldquot.header \
	po/en@quot.header \
	po/insert-header.sin \
	po/quot.sed \
	po/remove-potcdate.sed \
	po/remove-potcdate.sin \
	po/stamp-po \
	po/zbar.pot
MAINTAINERCLEANFILES =
BUILT_SOURCES =
EXTRA_DIST =
PHONY = $(SUBDIRS)

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = zbar.pc

dist_doc_DATA = COPYING HACKING.md INSTALL.md LICENSE.md NEWS.md README.md TODO.md

include $(srcdir)/include/Makefile.am.inc

SUBDIRS = zbar

if USE_NLS
  SUBDIRS += po
  dist_doc_DATA += ABOUT-NLS
endif

zbar/libzbar.la:
	$(MAKE) -C @builddir@/zbar libzbar.la

if HAVE_VIDEO
include $(srcdir)/zbarcam/Makefile.am.inc
endif
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
