# zedbar

Fast QR code and barcode scanner for Node.js, powered by WebAssembly.

A port of the [ZBar](http://zbar.sourceforge.net/) barcode scanning library from C to Rust, compiled to WebAssembly for use in Node.js.

## Features

- **Fast**: Native-speed scanning via WebAssembly
- **No native dependencies**: Works on any platform Node.js supports
- **Multiple formats**: QR Code, EAN-13, EAN-8, UPC-A, UPC-E, Code 128, Code 93, Code 39, Codabar, Interleaved 2 of 5, DataBar, SQ Code
- **CLI tool**: Command-line `zedbarimg` binary for scanning images (similar to `zbarimg`)

## Installation

```bash
npm install zedbar
```

## Command-line Usage

The package includes a `zedbarimg` command-line tool:

```bash
# Scan a single image
zedbarimg barcode.png

# Scan multiple images
zedbarimg image1.png image2.jpg image3.gif

# Quiet mode (only output barcode data)
zedbarimg --quiet qrcode.png

# Show help
zedbarimg --help
```

After installing globally, you can use it directly:

```bash
npm install -g zedbar
zedbarimg barcode.png
```

## Library Usage

### Scan from image file bytes

```javascript
import { scanImageBytes } from 'zedbar';
import { readFileSync } from 'fs';

// Read image file (PNG, JPEG, WebP, BMP)
const imageBytes = readFileSync('barcode.png');

// Scan for barcodes
for (const { symbolType, text } of scanImageBytes(imageBytes)) {
  console.log(`${symbolType}: ${text}`);
}
```

### Scan from grayscale data

If you already have grayscale image data:

```javascript
import { scanGrayscale } from 'zedbar';

// data is Uint8Array of grayscale pixels (1 byte per pixel, row-major)
for (const { symbolType, text } of scanGrayscale(data, width, height)) {
  console.log(`${symbolType}: ${text}`);
}
```

## API

### `scanImageBytes(bytes)`

Scans an encoded image (PNG, JPEG, WebP, BMP) for barcodes and QR codes.

**Parameters:**

- `bytes` (`Uint8Array` or `Buffer`) - Raw bytes of an image file

**Returns:** `DecodeResult[]` - Array of decoded barcodes

### `scanGrayscale(data, width, height)`

Scans grayscale image data for barcodes and QR codes.

**Parameters:**

- `data` (`Uint8Array`) - Grayscale pixel data, 1 byte per pixel, row-major order
- `width` (`number`) - Image width in pixels
- `height` (`number`) - Image height in pixels

**Returns:** `DecodeResult[]` - Array of decoded barcodes

### `DecodeResult`

- `symbolType` (`string`) - Barcode format (e.g., `"QR-Code"`, `"EAN-13"`)
- `data` (`Uint8Array`) - Raw decoded bytes
- `text` (`string | undefined`) - Decoded data as UTF-8 string, or `undefined` if not valid UTF-8

## Supported Formats

| Format | `symbolType` |
|--------|-------------|
| QR Code | `QR-Code` |
| EAN-13 | `EAN-13` |
| EAN-8 | `EAN-8` |
| UPC-A | `UPC-A` |
| UPC-E | `UPC-E` |
| Code 128 | `CODE-128` |
| Code 93 | `CODE-93` |
| Code 39 | `CODE-39` |
| Codabar | `Codabar` |
| Interleaved 2 of 5 | `I2/5` |
| DataBar | `DataBar` |
| SQ Code | `SQ-Code` |

## License

LGPL-3.0-or-later

Based on the [ZBar](http://zbar.sourceforge.net/) library. See [GitHub](https://github.com/eventualbuddha/zedbar) for source code.
