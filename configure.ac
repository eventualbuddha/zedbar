changecom()dnl
dnl Process this file with autoconf to produce a configure script.
AC_PREREQ([2.68])
AC_INIT([zbar],[0.23.93],[mchehab+huawei@kernel.org])
m4_ifndef([AC_LANG_DEFINES_PROVIDED],
          [m4_define([AC_LANG_DEFINES_PROVIDED])])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR(config)
AM_INIT_AUTOMAKE([1.13 -Werror foreign subdir-objects std-options dist-bzip2])
m4_pattern_allow([AM_PROG_AR])
AC_CONFIG_HEADERS([include/config.h])
AC_CONFIG_SRCDIR(zbar/scanner.c)
LT_PREREQ([2.2])
LT_INIT([dlopen win32-dll])
LT_LANG([Windows Resource])
AM_SILENT_RULES([yes])

dnl update these just before each release (along w/package version above)
dnl   LIB_VERSION update instructions copied from libtool docs:
dnl   library version follows the form current:revision:age
dnl   - If the library source code has changed at all since the last update,
dnl     then increment revision (c:r:a becomes c:r+1:a).
dnl   - If any interfaces have been added, removed, or changed,
dnl     increment current, and set revision to 0.
dnl   - If any interfaces have been added since the last public release,
dnl     then increment age.
dnl   - If any interfaces have been removed since the last public release,
dnl     then set age to 0.
AC_SUBST([LIB_VERSION], [3:0:3])
AC_SUBST([RELDATE], [2017-04-11])

AC_DEFINE_UNQUOTED([ZBAR_VERSION_MAJOR],
  [[`echo "$PACKAGE_VERSION" | sed -e 's/\..*$//'`]],
  [Program major version (before the '.') as a number])
AC_DEFINE_UNQUOTED([ZBAR_VERSION_MINOR],
  [[`echo "$PACKAGE_VERSION" | sed -e 's/^[^\.]*\.\([^\.]*\).*/\1/'`]],
  [Program minor version (after '.') as a number])
AC_DEFINE_UNQUOTED([ZBAR_VERSION_PATCH],
  [[`echo "$PACKAGE_VERSION" | sed -e 's/^[^\.]*\.[^\.]*\.*//' | sed s,^$,0,`]],
  [Program patch version (after the second '.') as a number])

cur=`echo "$LIB_VERSION" | sed -e 's/:.*$//'`
age=`echo "$LIB_VERSION" | sed -e 's/^.*://'`
AC_DEFINE_UNQUOTED([LIB_VERSION_MAJOR], [[$(( $cur - $age ))]],
  [Library major version])
AC_DEFINE_UNQUOTED([LIB_VERSION_MINOR], [[$age]],
  [Library minor version])
AC_DEFINE_UNQUOTED([LIB_VERSION_REVISION],
  [[`echo "$LIB_VERSION" | sed -e 's/^[^:]*:\([^:]*\):.*$/\1/'`]],
  [Library revision])

AM_CPPFLAGS="-I\$(top_srcdir)/include"
AM_CFLAGS="-Wall -Wno-parentheses"
AM_CXXFLAGS="$AM_CFLAGS"
AC_SUBST([AM_CPPFLAGS])
AC_SUBST([AM_CFLAGS])
AC_SUBST([AM_CXXFLAGS])

dnl windows build

AC_CANONICAL_HOST
case $host_os in
  *cygwin* | *mingw* | *uwin* | *djgpp* | *ems* )
    win32="yes"
    AC_DEFINE([_WIN32_WINNT], [0x0500], [Minimum Windows API version])
    AC_FUNC_ALLOCA
    AC_FUNC_ERROR_AT_LINE
    AC_FUNC_FSEEKO
    AC_CHECK_HEADERS([arpa/inet.h malloc.h mntent.h netdb.h netinet/in.h shadow.h sys/file.h sys/mount.h sys/param.h sys/socket.h sys/statfs.h sys/statvfs.h sys/vfs.h unistd.h values.h])
  ;;
  * )
    win32="no"
  ;;
esac
AM_CONDITIONAL([WIN32], [test "x$win32" = "xyes"])

dnl programs

AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_CXX

PKG_PROG_PKG_CONFIG

dnl symbologies

AC_ARG_ENABLE([codes],
  [AS_HELP_STRING([--enable-codes=SYMS],
    [select symbologies to compile [default=ean,databar,code128,code93,code39,codabar,i25,qrcode,sqcode]])],
  [],
  [enable_codes="ean,databar,code128,code93,code39,codabar,i25,qrcode,sqcode"])

AC_DEFUN([AC_DEFINE_SUBST],
   [AC_DEFINE($1,$2,$3)
    AC_SUBST($1,$2)])

AC_DEFUN([ZBAR_CHK_CODE], [
  AC_MSG_CHECKING([whether to build $2])
  enable_$1="no"
  AH_TEMPLATE([ENABLE_]translit($1, a-z, A-Z),
    [whether to build support for $2])
  AS_CASE([$enable_codes],
    [*$1* | *all*],
    [enable_$1="yes"
     enabled_codes="$enabled_codes $1"
     AC_DEFINE_SUBST([ENABLE_]translit($1, a-z, A-Z), [1])
    ], [
     disabled_codes="$disabled_codes $1"
     AC_DEFINE_SUBST([ENABLE_]translit($1, a-z, A-Z), [0])
  ])
  AM_CONDITIONAL([ENABLE_]translit($1, a-z, A-Z),
    [test "x$enable_$1" = "xyes"])
  AC_MSG_RESULT([$enable_$1])
])dnl

ZBAR_CHK_CODE([ean], [EAN symbologies])
ZBAR_CHK_CODE([databar], [DataBar symbology])
ZBAR_CHK_CODE([code128], [Code 128 symbology])
ZBAR_CHK_CODE([code93], [Code 93 symbology])
ZBAR_CHK_CODE([code39], [Code 39 symbology])
ZBAR_CHK_CODE([codabar], [Codabar symbology])
ZBAR_CHK_CODE([i25], [Interleaved 2 of 5 symbology])
ZBAR_CHK_CODE([qrcode], [QR Code])
ZBAR_CHK_CODE([sqcode], [SQ Code])
ZBAR_CHK_CODE([pdf417], [PDF417 symbology (incomplete)])

dnl libraries

AC_SEARCH_LIBS([clock_gettime], [rt pthread])
AM_ICONV()
AM_GNU_GETTEXT([external])
AM_GNU_GETTEXT_VERSION(0.20)
AM_GNU_GETTEXT_REQUIRE_VERSION(0.18)
AM_CONDITIONAL([USE_NLS], [test "x${USE_NLS}" = "xyes"])

dnl poll support

AC_CHECK_HEADERS([poll.h], [have_poll="yes"], [have_poll="no"])
AM_CONDITIONAL([HAVE_POLL], [test "x$have_poll" = "xyes"])

dnl pthreads
dnl FIXME this doesn't port well, integrate something like this:
dnl     http://autoconf-archive.cryp.to/acx_pthread.html

AC_ARG_ENABLE([pthread],
  [AS_HELP_STRING([--disable-pthread],
    [omit support for threaded applications])],
  [],
  [AS_IF([test "x$win32" = "xno"],
   [enable_pthread="yes"],
   [enable_pthread="no"
])])

AS_IF([test "x$enable_pthread" != "xno"],
  [AC_CHECK_HEADERS([pthread.h], [],
     [AC_MSG_FAILURE([test for pthread support failed!
configure --disable-pthread to skip threaded support.])])
   AC_CHECK_LIB([pthread], [pthread_create], [],
     [AC_MSG_FAILURE([unable to link against -lpthread, although you
appear to have pthread.h? set LDFLAGS and/or LIBS to help the linker,
or configure --disable-pthread to skip threaded support.])])
   AC_DEFINE([__USE_UNIX98], [1], [used only for pthread debug attributes])
])

dnl doc
AC_ARG_ENABLE([doc],
  [AS_HELP_STRING([--disable-doc],
    [disable building docs])],
  [],
  [enable_doc="yes"])

AS_IF([test "x$enable_doc" != "xno"],
  [AC_ARG_VAR([XMLTO], [location of xmlto, used for optional documentation generation])
   AC_ARG_VAR([XMLTOFLAGS], [additional arguments for xmlto])
   AC_CHECK_PROGS([XMLTO], [xmlto])])

AS_IF([test "x$XMLTO" = "x"], enable_doc="no")

AM_CONDITIONAL([HAVE_DOC], [test "x$enable_doc" != "xno"])

dnl video
AC_ARG_ENABLE([video],
  [AS_HELP_STRING([--disable-video],
    [exclude video scanner features])],
  [],
  [enable_video="yes"])

dnl video, directshow
AC_ARG_WITH([directshow],
  [AS_HELP_STRING([--with-directshow],
     [compile directshow driver on windows instead of vfw,
available only when video support is enabled])],
  [],
  [with_directshow="no"])

have_v4l1="no"
have_v4l2="no"
have_libv4l="no"
AS_IF([test "x$enable_video" = "xno"],
  [],
  [test "x$win32" = "xno"],
  [AC_CHECK_HEADERS([linux/videodev.h], [have_v4l1="yes"])
   AC_CHECK_HEADERS([linux/videodev2.h], [have_v4l2="yes"])
   AC_CHECK_HEADERS([libv4l2.h], [have_libv4l="yes"])
   AS_IF([test "x$have_v4l2" = "xno" && test "x$have_v4l1" = "xno"],
         [AC_MSG_FAILURE([test for video support failed!
rebuild your kernel to include video4linux support or
configure --disable-video to skip building video support.])],
         [AS_IF([test "x$have_v4l2" = "xno"],
               [AC_MSG_WARN([v4l2 API not detected, upgrade your kernel!])])]
        )],
  [AS_IF([test "x$with_directshow" != "xno"],
    [with_video="directshow"],
    [AC_CHECK_HEADERS([vfw.h], [with_video="vfw"],
      [AC_MSG_FAILURE([test for VfW video support failed!
configure --disable-video to skip building video support.])],[#include <windows.h>])])])

AS_IF([test "x$have_libv4l" = "xyes"],
   [PKG_CHECK_MODULES([V4L2], [libv4l2], [],
     [AC_MSG_FAILURE([unable to find libv4l2.so])])],
  [AC_MSG_WARN([libv4l not detected. Install it to support more cameras!])])

AM_CONDITIONAL([HAVE_VIDEO], [test "x$enable_video" != "xno"])
AM_CONDITIONAL([HAVE_V4L1], [test "x$have_v4l1" != "xno"])
AM_CONDITIONAL([HAVE_V4L2], [test "x$have_v4l2" != "xno"])
AM_CONDITIONAL([HAVE_LIBV4L], [test "x$have_libv4l" != "xno"])
AM_CONDITIONAL([WITH_DIRECTSHOW], [test "x$with_directshow" != "xno"])

dnl X
AC_ARG_VAR([XSHM_LIBS], [linker flags for X shared memory extension])

AS_IF([test "x$win32" != "xno"],
  [have_x="no"],
  [AC_PATH_XTRA
   AH_BOTTOM([#ifndef X_DISPLAY_MISSING
# define HAVE_X
#endif
])])
AM_CONDITIONAL([HAVE_X], [test "x$have_x" = "xyes"])

AS_IF([test "x$XSHM_LIBS" = "x"], [XSHM_LIBS="-lXext"])
AC_ARG_WITH([xshm],
  [AS_HELP_STRING([--without-xshm],
    [disable support for X shared memory extension])],
  [],
  [with_xshm="check"])

AS_IF([test "x$with_xshm" != "xno"],
  [AC_CHECK_HEADERS([X11/extensions/XShm.h],
    [with_xshm="yes"],
    [AS_IF([test "x$with_xshm" = "xcheck"],
      [with_xshm="no"],
      [AC_MSG_FAILURE([test for X shared memory extension failed!
install the X shared memory extension, specify --x-includes or
configure --without-xshm to disable the extension])])],
    [[#include <X11/Xlib.h>
#include <sys/ipc.h>
#include <sys/shm.h>
]])
   AS_IF([test "x$with_xshm" != "xno"],
     [AC_CHECK_LIB([Xext], [XShmQueryVersion],
       [with_xshm="yes"],
       [AC_MSG_FAILURE([unable to find XShmQueryVersion in $XSHM_LIBS!
specify XSHM_LIBS or configure --without-xshm to disable the extension])],
       ["$X_LIBS" "$X_PRE_LIBS" -lX11 "$X_EXTRA_LIBS" "$XSHM_LIBS"])
   ])
])
AM_CONDITIONAL([HAVE_XSHM], [test "x$with_xshm" = "xyes"])

AC_ARG_VAR([XV_LIBS], [linker flags for XVideo extension])
AS_IF([test "x$XV_LIBS" = "x"], [XV_LIBS="-lXv"])
AC_ARG_WITH([xv],
  [AS_HELP_STRING([--without-xv],
    [disable support for XVideo extension])],
  [],
  [with_xv="check"])

AS_IF([test "x$with_xv" != "xno"],
  [AC_CHECK_HEADERS([X11/extensions/Xvlib.h],
    [with_xv="yes"],
    [AS_IF([test "x$with_xv" = "xcheck"],
      [with_xv="no"],
      [AC_MSG_FAILURE([test for XVideo extension failed!
install the XVideo extension, specify --x-includes or
configure --without-xv to disable the extension])])],
    [[#include <X11/Xlib.h>
]])
   AS_IF([test "x$with_xv" != "xno"],
     [AC_CHECK_LIB([Xv], [XvQueryExtension],
       [with_xv="yes"],
       [AC_MSG_FAILURE([unable to find XvQueryExtension in $XV_LIBS!
specify XV_LIBS or configure --without-xv to disable the extension])],
       ["$X_LIBS" "$X_PRE_LIBS" -lX11 "$X_EXTRA_LIBS" "$XV_LIBS"])
   ])
])
AM_CONDITIONAL([HAVE_XV], [test "x$with_xv" = "xyes"])

dnl libjpeg
AC_ARG_WITH([jpeg],
  [AS_HELP_STRING([--without-jpeg],
    [disable support for JPEG image conversions])],
  [],
  [with_jpeg="check"])

have_jpeg="maybe"
AS_IF([test "x$with_jpeg" != "xno"],
  [AC_CHECK_HEADERS([jpeglib.h jerror.h], [], [have_jpeg="no"])
   AC_CHECK_LIB([jpeg], [jpeg_read_header], [], [have_jpeg="no"])
   AS_IF([test "x$have_jpeg" != "xno"],
     [with_jpeg="yes"],
     [test "x$with_jpeg" = "xyes"],
     [AC_MSG_FAILURE([unable to find libjpeg! ensure CFLAGS/LDFLAGS are
set appropriately or configure --without-jpeg])],
     [with_jpeg="no"])
])
AM_CONDITIONAL([HAVE_JPEG], [test "x$with_jpeg" = "xyes"])

dnl libpng
AC_ARG_WITH([png],
  [AS_HELP_STRING([--without-png],
    [disable support for png image conversions])],
  [],
  [with_png="check"])

have_png="maybe"
AS_IF([test "x$with_png" != "xno"],
  [AC_CHECK_HEADERS([pnglib.h jerror.h], [], [have_png="no"])
   AC_CHECK_LIB([png], [png_access_version_number], [], [have_png="no"])
   AS_IF([test "x$have_png" != "xno"],
     [with_png="yes"],
     [test "x$with_png" = "xyes"],
     [AC_MSG_FAILURE([unable to find libpng! ensure CFLAGS/LDFLAGS are
set appropriately or configure --without-png])],
     [with_png="no"])
])
AM_CONDITIONAL([HAVE_PNG], [test "x$with_png" = "xyes"])

dnl Mozilla NPAPI
AC_ARG_WITH([npapi],
  [AS_HELP_STRING([--with-npapi],
    [enable support for Firefox/Mozilla/OpenOffice plugin])],
  [],
  [with_npapi="no"])

AS_IF([test "x$with_npapi" != "xno"],
  [PKG_CHECK_MODULES([NPAPI], [firefox-plugin])
   NPAPI_VERSION=`$PKG_CONFIG firefox-plugin --modversion`
   AC_MSG_NOTICE([using firefox-plugin version $NPAPI_VERSION])])

AM_CONDITIONAL([HAVE_NPAPI], [test "x$with_npapi" = "xyes"])

AC_ARG_VAR([GLIB_GENMARSHAL], [full path to glib-genmarshal])

dnl header files

dnl FIXME switches for shm, mmap
AC_HEADER_ASSERT
AC_CHECK_HEADERS([errno.h fcntl.h features.h inttypes.h float.h limits.h \
  locale.h stddef.h stdlib.h string.h unistd.h sys/types.h sys/stat.h \
  sys/ioctl.h sys/time.h sys/times.h sys/ipc.h sys/shm.h sys/mman.h])
AC_HEADER_MAJOR
AC_CHECK_HEADER_STDBOOL

dnl types

AC_TYPE_INT32_T
AC_TYPE_UINT32_T
AC_TYPE_UINT8_T
AC_TYPE_UINTPTR_T
AC_TYPE_UID_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T
AC_CHECK_MEMBERS([struct stat.st_rdev])

dnl compile characteristics

AC_C_CONST
AC_C_INLINE

dnl functions

AC_FUNC_MMAP
AC_CHECK_FUNCS([alarm clock_gettime floor getcwd gettimeofday localeconv memchr memmove memset modf munmap pow select setenv sqrt strcasecmp strchr strdup strerror strrchr strstr strtol strtoul malloc realloc])


dnl output generation

dnl avoid doc rebuilds unless revision info changes
AC_CONFIG_COMMANDS([doc/version.xml],
  [AS_IF([test -f doc/version.xml && \
          ! echo $VERSION | diff doc/version.xml - >/dev/null 2>&1 || \
          ! echo $VERSION | diff $srcdir/doc/version.xml - >/dev/null 2>&1 ],
    [echo "writing new doc/version.xml" ; echo $VERSION > $srcdir/doc/version.xml ])],
  [VERSION="$VERSION"]
)
AC_CONFIG_COMMANDS([doc/reldate.xml],
  [AS_IF([test -f doc/reldate.xml && \
          ! echo $RELDATE | diff doc/reldate.xml - >/dev/null 2>&1 || \
          ! echo $RELDATE | diff $srcdir/doc/reldate.xml - >/dev/null 2>&1 ],
    [echo "writing new doc/reldate.xml" ; echo $RELDATE > $srcdir/doc/reldate.xml ])],
  [RELDATE="$RELDATE"]
)

echo "Generating config files"

AC_CONFIG_FILES([
Makefile
po/Makefile.in
zbar/Makefile
zbar.pc
doc/doxygen.conf])

AC_CONFIG_FILES([test/test_examples.sh],[chmod 755 test/test_examples.sh])

AC_OUTPUT

dnl summary log

echo ""
echo "please verify that the detected configuration matches your expectations:"
echo "------------------------------------------------------------------------"
AS_IF([test "x$USE_NLS" = "xyes"],
  [echo    "gettext                $USE_NLS"]
)
AS_IF([test "x$win32" != "xno"],
  [AS_IF([test "x$with_directshow" != "xno"],
    [echo "DirectShow driver      --with-directshow=$with_directshow"],
    [echo "VfW driver             --with-directshow=$with_directshow"]
  )],
  [echo "X                      --with-x=$have_x"]
)

echo    "pthreads               --enable-pthread=$enable_pthread"
echo    "doc                    --enable-doc=$enable_doc"
echo    "v4l                    --enable-video=$enable_video"
echo    "jpeg                   --with-jpeg=$with_jpeg"
echo    "png                    --with-png=$with_png"

echo "Enabled codes:        $enabled_codes"
echo "Disabled codes:       $disabled_codes"

dnl Display "warnings" about disabled and experimental features

echo ""
AS_IF([test "x$enable_video" != "xyes"],
  [echo "        => zbarcam video scanner will *NOT* be built"])
AS_IF([test "x$have_libv4l" != "xyes"],
  [echo "        => libv4l will *NOT* be used"])
AS_IF([test "x$with_jpeg" != "xyes"],
  [echo "        => JPEG image conversions will *NOT* be supported"])
AS_IF([test "x$with_png" != "xyes"],
  [echo "        => PNG image conversions will *NOT* be supported"])
AS_IF([test "x$enable_pdf417" = "xyes"],
  [echo "        => the pdf417 code support is incomplete!"])
