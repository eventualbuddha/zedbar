# zbar - Simplified for Rust Migration

This is a simplified fork of the [ZBar bar code reader](https://github.com/mchehab/zbar) library, focused on QR code scanning and being migrated to Rust.

## Project Goal

This project aims to:

1. Simplify the original ZBar C codebase by removing unnecessary complexity
2. Gradually migrate functionality to Rust
3. Maintain QR code scanning capabilities throughout the transition
4. Eventually create a pure Rust implementation

## Current Status

The codebase has been significantly simplified from the original ZBar:

- Removed autotools build system in favor of simple Makefile
- Removed threading, video capture, and GUI components
- Removed C++ bindings and language wrappers
- Focused on core QR code decoding functionality
- Simple command-line tool (`zbarimg`) for testing

## Building

### C Implementation (current)

```bash
make clean && make
```

This builds:

- `libzbar.a` - Static library with barcode decoding
- `zbarimg/zbarimg` - Command-line tool for decoding images

### Testing

```bash
./zbarimg/zbarimg examples/test-qr.png
./zbarimg/zbarimg examples/test-qr.jpg
```

## Dependencies

- Standard C library and math library

## Rust Migration

The `rust-wrapper` branch contains initial work on Rust FFI bindings and a gradual migration path. The goal is to eventually replace all C code with safe Rust implementations.

## License

LGPL 2.1 - See LICENSE.md for details.

This is a fork for educational/experimental purposes focused on Rust migration, not intended as a replacement for the upstream ZBar project.
