ZBAR BAR CODE READER
====================

ZBar Bar Code Reader is an open source software suite for reading bar
codes from various sources, such as video streams, image files and raw
intensity sensors. It supports EAN-13/UPC-A, UPC-E, EAN-8, Code 128,
Code 93, Code 39, Codabar, Interleaved 2 of 5, QR Code and SQ Code.

Included with the library are basic applications for decoding captured bar
code images and using a video device (e.g. webcam) as a bar code scanner.
For application developers, language bindings are included for C, C++,
Python 2 and Perl as well as GUI widgets for Qt.0.

Zbar also supports sending the scanned codes via dbus, allowing its
integration with other applications.

Check the ZBar home page for the latest release, mailing lists, etc.:

- <https://github.com/mchehab/zbar>

Tarballs with ZBar can be obtained from:

- <https://linuxtv.org/downloads/zbar/>

Since ZBar version 0.23.90, binaries auto-generated from Github's
Actions workflows are auto-generated for each release:

- <https://linuxtv.org/downloads/zbar/binaries/>

They contain binaries for:

- Ubuntu SID, generated via pbuilder;
- Mac OS;
- Windows, for 4 different configurations:
  - 32 bits/64 bits;
  - Video for Windows (VfW) or DirectShow (DShow).

License information can be found in `COPYING`.

You may find some outdated documentation at the original ZBar's
site at Sourceforge, but please notice that the content there is not
updated for ages:
 <http://zbar.sourceforge.net/>

BUILDING
========

See `INSTALL.md` for generic configuration and build instructions.

Please notice that at least autotools related packages and a
C compiler are needed, in order to generate the configure script.

So, on Debian, at least those packages are needed:
 autoconf autopoint pkg-config libtool gcc make

If you have installed all needed dependencies, all you need to do is to run:

```
autoreconf -vfi
./configure
make
```

- NOTES

  1) Currently, we maintain a Continuous Integration build test at
     TravisCI:

        <https://travis-ci.org/github/mchehab/zbar/>

     Due to that, there are scripts meant to test ZBar build on
     Linux, Windows and MacOS, that could be helpful. Please see
     the `.travis.yml` file, and the corresponding scripts under `travis/`.

The scanner/decoder library itself only requires a few standard
library functions which should be available almost anywhere.

The zbarcam program uses the video4linux API (v4l1 or v4l2) to access
the video device.  This interface is part of the linux kernel, a 3.16
kernel or upper is recommended for full support.  More information is
available at:

- <http://www.linuxtv.org/wiki/>

`pkg-config` is used to locate installed libraries.  You should have
installed `pkg-config` if you need any of the remaining components.
pkg-config may be obtained from:

- <http://pkg-config.freedesktop.org/>

The `zbarimg` program uses `ImageMagick` to read image files in many
different formats.  You will need at least `ImageMagick` version 6.2.6
if you want to scan image files. You may also use `GraphicsMagick`
package instead.

`ImageMagick` may be obtained from:

- <http://www.imagemagick.org/>

Qt Widget
---------

The Qt widget requires Qt4 or Qt5. You will need Qt if you would like to
use or develop a Qt GUI application with an integrated bar code
scanning widget. Qt4 may be obtained from:

- <https://www.qt.io/>

Python widgets
**Python bindings**

The Python bindings require Python 2 or 3 and provide only non-GUI functions.
You will need Python and PIL or Pillow if you would like to scan images or
video directly using Python. Python is available from:

- <http://python.org/>

Java Widget
-----------

The Java ZBar widget uses Java Native Interface (JNI), meaning that the
widget will contain machine-dependent code. It works with Java version
7 and above.  Java open JDK is available from:

- <https://openjdk.java.net/>

RUNNING
=======

`make install` will install the library and application programs.  Run
`zbarcam-qt` or `zbarcam` to start the video scanner. Use `zbarimg <file>`
to decode a saved image file.

Check the manual to find specific options for each program.

DBUS TESTING
============

In order to test if dbus is working, you could use:

 dbus-monitor --system interface=org.linuxtv.Zbar1.Code

or build the test programs with:

 make test_progs

And run:
 $ ./test/test_dbus

With that, running this command on a separate shell:

 $ ./zbarimg/zbarimg examples/code-128.png
 CODE-128:<https://github.com/mchehab/zbar>
 scanned 1 barcode symbols from 1 images in 0.01 seconds

Will produce this output at test_dbus shell window:

 Waiting for Zbar events
 Type = CODE-128
 Value = <https://github.com/mchehab/zbar>

REPORTING BUGS
==============

Bugs can be reported on the project page:

- <https://github.com/mchehab/zbar>

Please include the ZBar version number and a detailed description of
the problem.  You'll probably have better luck if you're also familiar
with the concepts from:

- <http://www.catb.org/~esr/faqs/smart-questions.html>
