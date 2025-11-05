# zbar

A Rust implementation of barcode scanning, based on the [ZBar bar code reader](https://github.com/mchehab/zbar) library.

## Features

- QR code decoding
- Support for multiple barcode formats (EAN, UPC, Code128, Code39, etc.)
- Pure Rust implementation
- Command-line tool (`zbarimg`) for scanning images

## Building

```bash
cargo build --release
```

## Usage

### Command-line tool

```bash
cargo run --bin zbarimg examples/test-qr.png
cargo run --bin zbarimg examples/test-qr.jpg
```

### Library

```rust
use zbar::{Scanner, Image};

let img = image::open("examples/test-qr.png")?;
let gray = img.to_luma8();
let (width, height) = gray.dimensions();

let mut zbar_img = Image::from_gray(gray.as_raw(), width, height)?;
let mut scanner = Scanner::new();
let num_symbols = scanner.scan(&mut zbar_img)?;

for symbol in zbar_img.symbols() {
    println!("Decoded {}: {}", symbol.symbol_type(), symbol.data_string()?);
}
```

## Testing

```bash
cargo test
```

## License

LGPL 3.0 or later - See LICENSE for details.
