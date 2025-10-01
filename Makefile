# Simple Makefile for ZBar - Barcode Reader Library
# This is a simplified build system for preparing the codebase for Rust migration

CC = gcc
CFLAGS = -Wall -Wextra -O2 -fPIC -Izbar
LDFLAGS = -ljpeg -lpng -lm

# Library source files
LIB_SOURCES = \
	zbar/decoder.c \
	zbar/error.c \
	zbar/image.c \
	zbar/img_scanner.c \
	zbar/processor.c \
	zbar/scanner.c \
	zbar/symbol.c \
	zbar/decoder/codabar.c \
	zbar/decoder/code128.c \
	zbar/decoder/code39.c \
	zbar/decoder/code93.c \
	zbar/decoder/databar.c \
	zbar/decoder/ean.c \
	zbar/decoder/i25.c \
	zbar/decoder/qr_finder.c \
	zbar/decoder/sq_finder.c \
	zbar/qrcode/qrdec.c \
	zbar/qrcode/qrdectxt.c \
	zbar/qrcode/rs.c

# Application source files
ZBARIMG_SOURCES = \
	zbarimg/zbarimg.c

# Object files
LIB_OBJECTS = $(LIB_SOURCES:.c=.o)
ZBARIMG_OBJECTS = $(ZBARIMG_SOURCES:.c=.o)

# Targets
all: libzbar.a zbarimg-bin

# Static library
libzbar.a: $(LIB_OBJECTS)
	ar rcs $@ $^

# zbarimg application
zbarimg-bin: $(ZBARIMG_OBJECTS) libzbar.a
	$(CC) -o zbarimg/zbarimg $^ $(LDFLAGS)

# Object file compilation
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Test with sample images
test: zbarimg-bin
	@echo "Testing zbarimg..."
	@if [ -f examples/test-qr.png ]; then \
		./zbarimg/zbarimg examples/test-qr.png; \
	else \
		echo "No test images found in examples/"; \
	fi

# Clean up build artifacts
clean:
	rm -f $(LIB_OBJECTS) $(ZBARIMG_OBJECTS) libzbar.a zbarimg/zbarimg

# Install (basic)
install: libzbar.a zbarimg-bin
	@echo "Installing to /usr/local..."
	install -d /usr/local/lib /usr/local/bin /usr/local/include
	install -m 644 libzbar.a /usr/local/lib/
	install -m 755 zbarimg/zbarimg /usr/local/bin/zbarimg
	install -m 644 zbar/zbar.h /usr/local/include/

# Show source files (useful for Rust migration planning)
list-sources:
	@echo "Library sources:"
	@echo $(LIB_SOURCES) | tr ' ' '\n' | sort
	@echo
	@echo "Application sources:"
	@echo $(ZBARIMG_SOURCES) | tr ' ' '\n' | sort

.PHONY: all test clean install list-sources zbarimg-bin
