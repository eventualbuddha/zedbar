# zedbar

Fast QR code and barcode scanner for Node.js, powered by WebAssembly.

A port of the [ZBar](http://zbar.sourceforge.net/) barcode scanning library from C to Rust, compiled to WebAssembly for use in Node.js.

## Features

- **Fast**: Native-speed scanning via WebAssembly
- **No native dependencies**: Works on any platform Node.js supports
- **Multiple formats**: QR Code, EAN-13, EAN-8, UPC-A, UPC-E, Code 128, Code 93, Code 39, Codabar, Interleaved 2 of 5, DataBar, SQ Code

## Installation

```bash
npm install zedbar
```

## Usage

```javascript
import { scanGrayscale } from 'zedbar';
import sharp from 'sharp';

// Load image and convert to grayscale
const { data, info } = await sharp('barcode.png')
  .grayscale()
  .raw()
  .toBuffer({ resolveWithObject: true });

// Scan for barcodes
for (const { symbolType, text } of scanGrayscale(data, info.width, info.height)) {
  console.log(`${symbolType}: ${text}`);
}
```

## API

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
