# zedbar

A pure Rust implementation of barcode scanning, based on the [ZBar bar code reader](https://github.com/mchehab/zbar) C library.

This is a port of the ZBar library to Rust, providing barcode detection and decoding capabilities with a type-safe, idiomatic Rust API.

## Features

- **Multiple Barcode Formats**: QR Code, EAN-13, EAN-8, EAN-2, EAN-5, UPC-A, UPC-E, ISBN-10, ISBN-13, Code 128, Code 93, Code 39, Codabar, Interleaved 2 of 5, DataBar (RSS), DataBar Expanded, SQCode
- **Pure Rust**: No C dependencies, fully memory-safe implementation
- **Type-Safe Configuration**: Compile-time validated configuration API
- **Command-line Tool**: `zedbarimg` utility for scanning images from the command line
- **Position Tracking**: Optional tracking of barcode positions in images
- **Inverted Image Support**: Can detect barcodes in inverted images

## Installation

```bash
cargo add zedbar
```

### Cargo Features

By default, all symbologies are enabled. You can selectively enable only the ones you need to reduce compile time and binary size:

```bash
cargo add zedbar --no-default-features --features qrcode,ean
```

#### Symbology Features

- `qrcode` - QR Code 2D barcode
- `sqcode` - SQ Code 2D barcode
- `ean` - EAN-8, EAN-13, EAN-2, EAN-5, UPC-A, UPC-E, ISBN-10, ISBN-13
- `code128` - Code 128
- `code39` - Code 39
- `code93` - Code 93
- `codabar` - Codabar
- `databar` - GS1 DataBar (RSS) and DataBar Expanded
- `i25` - Interleaved 2 of 5

#### Optional Dependencies

**Heavy dependencies (tied to features):**

- `encoding_rs`, `reed-solomon`, `rand`, `rand_chacha` - Required for QR code decoding (enabled with `qrcode` feature)
- `image` - Image format loading (PNG, JPEG, etc.) - needed for tests and the `zedbarimg` binary
- `clap` - Command-line parsing (needed for the `zedbarimg` binary)
- `wasm-bindgen`, `js-sys`, `getrandom` - JavaScript bindings (enabled with the `wasm` feature, which is off by default and builds the [npm package](npm/README.md))

**Note:** 1D barcode decoders (EAN, Code39, Code128, etc.) have **zero external dependencies**!

The `default` feature enables all symbologies plus optional dependencies:

```toml
default = ["qrcode", "sqcode", "ean", "code128", "code39", "code93", "codabar", "databar", "i25", "image", "clap"]
```

#### Minimal Library

For the absolute minimal build with zero external dependencies (1D barcodes only):

```bash
cargo add zedbar --no-default-features --features ean
```

For QR codes only (with necessary dependencies):

```bash
cargo add zedbar --no-default-features --features qrcode
```

Note: Disabling a feature at compile-time means that symbology will not be compiled into the binary at all, which is different from disabling it via runtime configuration.

## Usage

### Library

```rust
use zedbar::{Image, Scanner};

// Load an image and convert it for scanning
let img = image::open("barcode.png")?;
let mut img = Image::from_dynamic(&img)?;

// Create scanner and scan image
let mut scanner = Scanner::new();
let symbols = scanner.scan(&mut img);

// Process decoded symbols
for symbol in symbols {
    println!("{:?}: {}", symbol.symbol_type(), symbol.data_string().unwrap_or(""));
}
```

`Image::from_dynamic` (available with the `image` feature) composites any
transparency over white before converting to grayscale. Prefer it over
`img.to_luma8()`: a barcode drawn in black on a transparent background — the
usual shape of one rendered from SVG — loses its alpha in a direct conversion
and becomes black on black. If you already hold grayscale pixels, pass them
straight to `Image::from_gray(data, width, height)`.

### Locating a Symbol in the Image

Each symbol records the image-coordinate points where it was detected,
which lets you draw a box around the decode or crop the source image.

```rust
for symbol in symbols {
    // QR codes record the four corners of their bounding rectangle;
    // linear barcodes record per-scan touchpoints.
    for point in symbol.points() {
        println!("  point at ({}, {})", point.x, point.y);
    }

    // Or, the axis-aligned bounding box of all those points:
    if let Some(b) = symbol.bounds() {
        println!("  bounds: {}×{} at ({}, {})", b.width, b.height, b.x, b.y);
    }
}
```

`bounds()` returns the AABB of `points()`, with `width` and `height`
reported as `max - min` of the recorded points (the extent between the
outermost points), not as a pixel count. Both return empty / `None` for
symbols that did not record any points.

### Choosing Symbologies

`DecoderConfig::new()` starts **empty** — opt in to each symbology you need:

```rust
use zedbar::config::*;
use zedbar::{DecoderConfig, Scanner};

let config = DecoderConfig::new()
    .enable(QrCode)
    .enable(Ean13);

let mut scanner = Scanner::with_config(config);
```

For exploratory use ("just scan whatever's there"), `DecoderConfig::all()`
or `Scanner::new()` enables every supported symbology in one call.

### Advanced Configuration

```rust
use zedbar::config::*;
use zedbar::{DecoderConfig, Scanner};

let config = DecoderConfig::new()
    .enable(QrCode)
    .enable(Ean13)
    .set_length_limits(Code39, 4, 20)  // also enables Code39
    .test_inverted(true)               // Try inverted image if no symbols found
    .retry_undecoded_regions(true)     // Crop+upscale small QR codes automatically
    .scan_density(2, 2);               // Scan every 2nd line (faster)

let mut scanner = Scanner::with_config(config);
```

The per-symbology setters (`set_length_limits`, `set_checksum`,
`set_uncertainty`) auto-enable the symbology they target — if you've
already mentioned a symbology by name, it's on.

### Small QR Codes in Large Images

When a QR code is too small relative to the image (e.g. a QR code on a scanned page),
the scanner reports undecoded finder regions. You can handle these manually, or enable
automatic retry:

```rust
use zedbar::config::*;
use zedbar::{DecoderConfig, Scanner, Image};

// Option 1: Automatic retry. Each undecoded region is cropped, upscaled and
// rescanned, and any QR or SQ code recovered from it joins the results with
// its coordinates mapped back to the original image.
let config = DecoderConfig::new()
    .enable(QrCode)
    .retry_undecoded_regions(true);
let mut scanner = Scanner::with_config(config);
let result = scanner.scan(&mut img);

// Option 2: Manual control via finder_regions(), which reports one entry per
// cluster of finder patterns that did not yield a symbol.
let mut scanner = Scanner::new();
let result = scanner.scan(&mut img);
for region in result.finder_regions() {
    let pad = region.width.max(region.height) / 2;
    let x = region.x.saturating_sub(pad);
    let y = region.y.saturating_sub(pad);
    let w = (region.width + 2 * pad).min(img.width() - x);
    let h = (region.height + 2 * pad).min(img.height() - y);
    if let Some(mut upscaled) = img.crop(x, y, w, h).and_then(|c| c.upscale(4)) {
        let retry = scanner.scan(&mut upscaled);
        // process retry.symbols()...
    }
}
```

With automatic retry, the regions left in `finder_regions()` are the ones the
retry could not resolve.

### Command-line Tool

`zedbarimg` scans one or more image files and prints what it finds, in the
manner of ZBar's `zbarimg`:

```bash
cargo install zedbar
zedbarimg barcode.png qrcode.jpg
```

`--quiet` prints only the decoded data, `--raw` leaves it unconverted by any
charset, and `--disable-all` plus the `--enable-*` flags narrow the scan to
chosen symbologies. `zedbarimg --help` lists them all.

From a checkout of this repository, where the test images live:

```bash
# Install the binary from the working tree
cargo install --path .
zedbarimg examples/test-qr.png

# Or run it without installing. Use --release: a debug build carries
# overflow checks and no optimization, and scans far slower.
cargo run --release --bin zedbarimg -- examples/test-ean13.png
```

### JavaScript

The same scanner is published to npm as [`zedbar`](https://www.npmjs.com/package/zedbar),
compiled to WebAssembly, with its own `zedbarimg` command. See
[npm/README.md](npm/README.md).

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
- **[quircs](https://crates.io/crates/quircs)** - Pure Rust port of the quirc QR code library.

### C/C++ Libraries (via FFI)

- **[ZBar](https://github.com/mchehab/zbar)** - The original C library this project is based on. Mature and well-tested.
- **[ZXing](https://github.com/zxing/zxing)** - Popular Java library with C++ port available.
- **[ZBar-lite](https://github.com/xantares/zbar-lite)** - Lightweight fork of ZBar.

### Choosing an Alternative

- **For QR codes only**: Consider `rqrr` - it's fast, pure Rust, and has a simpler API.
- **For maximum format support**: The original ZBar C library or ZXing are very mature.
- **For pure Rust with broad format support**: This library (`zedbar`) or `rxing`.
- **For C bindings to ZBar**: Use `zbar-rust` which provides FFI bindings to the original C library.

## License

LGPL 3.0 or later - See [LICENSE](LICENSE) for details.

This library is licensed under the GNU Lesser General Public License v3.0 or later, consistent with the original ZBar library's licensing.
