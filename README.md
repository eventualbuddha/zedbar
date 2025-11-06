# zbar-rust

A pure Rust implementation of barcode scanning, based on the [ZBar bar code reader](https://github.com/mchehab/zbar) C library.

This is a port of the ZBar library to Rust, providing barcode detection and decoding capabilities with a type-safe, idiomatic Rust API.

## Features

- **Multiple Barcode Formats**: QR Code, EAN-13, EAN-8, UPC-A, UPC-E, ISBN-10, ISBN-13, Code 128, Code 93, Code 39, Codabar, Interleaved 2 of 5, DataBar (RSS), SQCode
- **Pure Rust**: No C dependencies, fully memory-safe implementation
- **Type-Safe Configuration**: Compile-time validated configuration API
- **Command-line Tool**: `zbarimg` utility for scanning images from the command line
- **Position Tracking**: Optional tracking of barcode positions in images
- **Inverted Image Support**: Can detect barcodes in inverted images

## Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
zbar-rust = "0.1"
```

### Cargo Features

By default, all symbologies are enabled. You can selectively enable only the ones you need to reduce compile time and binary size:

```toml
[dependencies]
zbar-rust = { version = "0.1", default-features = false, features = ["qrcode", "ean"] }
```

#### Symbology Features

- `qrcode` - QR Code 2D barcode
- `sqcode` - SQ Code 2D barcode
- `ean` - EAN-8, EAN-13, UPC-A, UPC-E, ISBN-10, ISBN-13
- `code128` - Code 128
- `code39` - Code 39
- `code93` - Code 93
- `codabar` - Codabar
- `databar` - GS1 DataBar (RSS)
- `i25` - Interleaved 2 of 5

#### Optional Dependencies

**Heavy dependencies (tied to features):**

- `encoding_rs`, `reed-solomon`, `rand`, `rand_chacha` - Required for QR code decoding (enabled with `qrcode` feature)
- `image` - Image format loading (PNG, JPEG, etc.) - needed for tests and the `zbarimg` binary
- `clap` - Command-line parsing (needed for the `zbarimg` binary)

**Note:** 1D barcode decoders (EAN, Code39, Code128, etc.) have **zero external dependencies**!

The `default` feature enables all symbologies plus optional dependencies:

```toml
default = ["qrcode", "sqcode", "ean", "code128", "code39", "code93", "codabar", "databar", "i25", "image", "clap"]
```

#### Minimal Library

For the absolute minimal build with zero external dependencies (1D barcodes only):

```toml
[dependencies]
zbar-rust = { version = "0.1", default-features = false, features = ["ean"] }
```

For QR codes only (with necessary dependencies):

```toml
[dependencies]
zbar-rust = { version = "0.1", default-features = false, features = ["qrcode"] }
```

Note: Disabling a feature at compile-time means that symbology will not be compiled into the binary at all, which is different from disabling it via runtime configuration.

## Usage

### Library

```rust
use zbar::{Image, Scanner};

// Load and convert image to grayscale
let img = image::open("barcode.png")?;
let gray = img.to_luma8();
let (width, height) = gray.dimensions();

// Create scanner and scan image
let mut zbar_img = Image::from_gray(gray.as_raw(), width, height)?;
let mut scanner = Scanner::new();
let symbols = scanner.scan(&mut zbar_img);

// Process decoded symbols
for symbol in symbols {
    println!("{:?}: {}", symbol.symbol_type(), symbol.data_string().unwrap_or(""));
}
```

### Advanced Configuration

```rust
use zbar::config::*;
use zbar::{DecoderConfig, Scanner};

let config = DecoderConfig::new()
    .enable(QrCode)
    .enable(Ean13)
    .set_binary(QrCode, true)          // Preserve binary data in QR codes
    .set_length_limits(Code39, 4, 20)  // Code39 must be 4-20 chars
    .test_inverted(true)               // Try inverted image if no symbols found
    .scan_density(2, 2);               // Scan every 2nd line (faster)

let mut scanner = Scanner::with_config(config);
```

### Command-line Tool

```bash
# Build the tool
cargo build --release --bin zbarimg

# Scan an image
cargo run --bin zbarimg examples/test-qr.png
cargo run --bin zbarimg examples/test-ean13.png
```

## Testing

```bash
cargo test
```

## Benchmarks

Comprehensive benchmarks comparing this library with rqrr, rxing, and the original C zbar library are available:

```bash
# Compare with rqrr (default)
cargo bench

# Compare with all libraries (requires optional dependencies)
cargo bench --features bench_zbar_c,bench_rxing
```

See [benches/README.md](benches/README.md) for detailed benchmark documentation and results.

## Credits

### Original ZBar Library

This project is based on the [ZBar bar code reader](https://github.com/mchehab/zbar) library:

- Original C implementation: Copyright (C) 2007-2010 Jeff Brown <spadix@users.sourceforge.net>
- QR code decoder components: Copyright (C) 2008-2009 Timothy B. Terriberry (<tterribe@xiph.org>)
- SQCode decoder: Copyright (C) 2018 Javier Serrano Polo <javier@jasp.net>
- Current C library maintenance: Mauro Carvalho Chehab and contributors

The original ZBar library is licensed under LGPL 2.1 or later.

### Rust Port

This Rust implementation preserves the algorithms and structure of the original C library while providing a safe, idiomatic Rust API.

## Alternatives

If this library doesn't meet your needs, consider these alternatives:

### Rust Libraries

- **[rqrr](https://crates.io/crates/rqrr)** - Pure Rust QR code reader with a different algorithm. Fast and reliable for QR codes specifically, but only supports QR codes.
- **[bardecoder](https://crates.io/crates/bardecoder)** - Another Rust barcode decoder supporting various 1D formats.
- **[rxing](https://crates.io/crates/rxing)** - Rust port of ZXing (Zebra Crossing) library, supports many formats.
- **[quircs](https://crates.io/crates/quircs)** - Rust bindings to the quirc QR code library.

### C/C++ Libraries (via FFI)

- **[ZBar](https://github.com/mchehab/zbar)** - The original C library this project is based on. Mature and well-tested.
- **[ZXing](https://github.com/zxing/zxing)** - Popular Java library with C++ port available.
- **[ZBar-lite](https://github.com/xantares/zbar-lite)** - Lightweight fork of ZBar.

### Choosing an Alternative

- **For QR codes only**: Consider `rqrr` - it's fast, pure Rust, and has a simpler API.
- **For maximum format support**: The original ZBar C library or ZXing are very mature.
- **For pure Rust with broad format support**: This library (`zbar-rust`) or `rxing`.

## License

LGPL 3.0 or later - See [LICENSE](LICENSE) for details.

This library is licensed under the GNU Lesser General Public License v3.0 or later, consistent with the original ZBar library's licensing.
