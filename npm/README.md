# zedbar

Fast QR code and barcode scanner for Node.js, powered by WebAssembly.

A port of the [ZBar](https://github.com/mchehab/zbar) barcode scanning library from C to Rust, compiled to WebAssembly for use in Node.js.

## Features

- **Fast**: Native-speed scanning via WebAssembly
- **No native dependencies**: Works on any platform Node.js supports
- **Multiple formats**: QR Code, EAN-13, EAN-8, EAN-2, EAN-5, UPC-A, UPC-E, ISBN-10, ISBN-13, Code 128, Code 93, Code 39, Codabar, Interleaved 2 of 5, DataBar, DataBar Expanded, SQ Code
- **CLI tool**: Command-line `zedbarimg` binary for scanning images (similar to `zbarimg`)

## Installation

```bash
npm install zedbar
```

## Command-line Usage

The package includes a `zedbarimg` command-line tool. Install the package
globally to put it on your `PATH`:

```bash
npm install -g zedbar
```

```bash
# Scan a single image
zedbarimg barcode.png

# Scan multiple images
zedbarimg image1.png image2.jpg image3.webp

# Quiet mode (only output barcode data)
zedbarimg --quiet qrcode.png

# Raw mode (decoded bytes, no charset conversion)
zedbarimg --raw qrcode.png

# Show help
zedbarimg --help
```

Without a global install, run it as `npx --package=zedbar zedbarimg barcode.png`.

It exits non-zero when a file cannot be read or no barcode is found, and prints
a `scanned N barcode symbols` summary to stderr so the decoded data on stdout
stays pipeable.

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

### Restrict what is scanned

```javascript
import { scanImageBytes, ScanOptions } from 'zedbar';

const options = new ScanOptions();
options.symbologies = ['QR-Code'];

const results = scanImageBytes(imageBytes, options);
```

`ScanOptions` is a class: construct it with `new` and assign its properties.
Passing a plain object literal throws `expected instance of ScanOptions`.

The package is CommonJS, so `require` works too:

```javascript
const { scanImageBytes, ScanOptions } = require('zedbar');
```

## API

### `scanImageBytes(bytes, options?)`

Scans an encoded image (PNG, JPEG, WebP, BMP) for barcodes and QR codes.

**Parameters:**

- `bytes` (`Uint8Array` or `Buffer`) - Raw bytes of an image file
- `options` (`ScanOptions`, optional) - Scanning options

**Returns:** `DecodeResult[]` - Array of decoded barcodes

### `scanGrayscale(data, width, height, options?)`

Scans grayscale image data for barcodes and QR codes.

**Parameters:**

- `data` (`Uint8Array`) - Grayscale pixel data, 1 byte per pixel, row-major order
- `width` (`number`) - Image width in pixels
- `height` (`number`) - Image height in pixels
- `options` (`ScanOptions`, optional) - Scanning options

**Returns:** `DecodeResult[]` - Array of decoded barcodes

### `ScanOptions`

A class, not a plain object: construct it with `new ScanOptions()` and assign
the properties you want. Both are write-only setters.

- `retryUndecodedRegions` (`boolean`, default: `true`) - Automatically retry undecoded QR finder regions by cropping and upscaling. Disable to skip the retry and reduce processing time for images that are known to have sufficient resolution.
- `symbologies` (`string[]`, optional) - Restrict scanning to the listed symbologies. When omitted, every supported symbology is enabled. Names match the `symbolType` field on results — see [Supported Formats](#supported-formats) — and are case-sensitive; an unrecognized name throws `unknown symbology: "..."`.

```javascript
// Scan QR codes only
const options = new ScanOptions();
options.symbologies = ['QR-Code'];
const results = scanImageBytes(imageBytes, options);
```

### `DecodeResult`

- `symbolType` (`string`) - Barcode format (e.g., `"QR-Code"`, `"EAN-13"`), or `"Partial"` for something read but incomplete — see [Supported Formats](#supported-formats)
- `data` (`Uint8Array`) - Raw decoded bytes. For QR codes these are the bytes as encoded in the barcode, before any charset conversion; for SQ codes they are the sampled bit matrix, which is all that symbology decodes to.
- `text` (`string | undefined`) - The data as a string. For QR codes the encoding is detected (UTF-8, Shift-JIS, Windows-1252, or whatever an ECI designator names) and converted; for linear barcodes the data is always ASCII. SQ codes carry arbitrary binary, so this is `undefined` for them unless those bytes happen to be valid UTF-8.
- `points` (`Point[]`) - Position points in image coordinates. For QR codes, the four corners of the QR's bounding rectangle (in implementation-defined order). For linear barcodes, one or more touchpoints accumulated as the symbol was scanned. Empty if no points were recorded.
- `bounds` (`Bounds | undefined`) - Axis-aligned bounding rectangle of `points`, or `undefined` if no points were recorded.

### `Point`

- `x` (`number`)
- `y` (`number`)

### `Bounds`

- `x`, `y` (`number`) - Top-left corner in image coordinates.
- `width`, `height` (`number`) - Reported as `max - min` of the recorded points (the horizontal and vertical extent between the outermost points), not as a pixel count.

```javascript
for (const result of scanImageBytes(imageBytes)) {
  for (const { x, y } of result.points) {
    console.log(`  point at (${x}, ${y})`);
  }
  if (result.bounds) {
    const { x, y, width, height } = result.bounds;
    console.log(`  bounds: ${width}×${height} at (${x}, ${y})`);
  }
}
```

## Supported Formats

| Format | `symbolType` |
|--------|-------------|
| QR Code | `QR-Code` |
| EAN-13 | `EAN-13` |
| EAN-8 | `EAN-8` |
| EAN-2 add-on | `EAN-2` |
| EAN-5 add-on | `EAN-5` |
| UPC-A | `UPC-A` |
| UPC-E | `UPC-E` |
| ISBN-10 | `ISBN-10` |
| ISBN-13 | `ISBN-13` |
| Code 128 | `CODE-128` |
| Code 93 | `CODE-93` |
| Code 39 | `CODE-39` |
| Codabar | `Codabar` |
| Interleaved 2 of 5 | `I2/5` |
| DataBar | `DataBar` |
| DataBar Expanded | `DataBar-Exp` |
| SQ Code | `SQ-Code` |

A result may also carry the `symbolType` `Partial`, which is not a format and
cannot be passed to `symbologies`. It means something was read but not a whole
symbol: a structured-append QR group with at least one of its codes missing
from the image. Its `data` is the parts that were found, joined with a NUL byte
wherever an absent part belongs.

## License

LGPL-3.0-or-later

Based on the [ZBar](http://zbar.sourceforge.net/) library. See [GitHub](https://github.com/eventualbuddha/zedbar) for source code.
