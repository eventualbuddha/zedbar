#!/usr/bin/env node

/**
 * zedbarimg - Scan and decode barcodes from image files
 *
 * A command-line tool for scanning barcodes and QR codes from images.
 * Powered by the zedbar Rust library compiled to WebAssembly.
 */

import { readFileSync } from "fs";
import { dirname, join } from "path";
import { fileURLToPath } from "url";
import { createRequire } from "module";

const __filename = fileURLToPath(import.meta.url);
const __dirname = dirname(__filename);
const require = createRequire(import.meta.url);

// Import the CommonJS WASM module
const zedbar = require("../zedbar.js");

// Parse command-line arguments
function parseArgs() {
  const args = process.argv.slice(2);
  const options = {
    quiet: false,
    raw: false,
    files: [],
  };

  for (let i = 0; i < args.length; i++) {
    const arg = args[i];

    if (arg === "-q" || arg === "--quiet") {
      options.quiet = true;
    } else if (arg === "--raw") {
      options.raw = true;
    } else if (arg === "-h" || arg === "--help") {
      showHelp();
      process.exit(0);
    } else if (arg === "-v" || arg === "--version") {
      const pkgPath = join(__dirname, "..", "package.json");
      const pkg = JSON.parse(readFileSync(pkgPath, "utf8"));
      console.log(`zedbarimg ${pkg.version}`);
      process.exit(0);
    } else if (!arg.startsWith("-")) {
      options.files.push(arg);
    } else {
      console.error(`Unknown option: ${arg}`);
      showHelp();
      process.exit(1);
    }
  }

  if (options.files.length === 0) {
    console.error("Error: No image files specified\n");
    showHelp();
    process.exit(1);
  }

  return options;
}

function showHelp() {
  console.log(`
Usage: zedbarimg [OPTIONS] <IMAGE_FILES...>

Scan and decode barcodes from one or more image files

Options:
  -q, --quiet      Minimal output, only print decoded symbol data
  --raw            Output decoded symbol data without converting charsets
  -h, --help       Print help information
  -v, --version    Print version information

Supported formats:
  - Images: PNG, JPEG, BMP, WebP
  - Barcodes: QR Code, EAN-13, EAN-8, UPC-A, UPC-E, Code 128, Code 93,
              Code 39, Codabar, Interleaved 2/5, DataBar, SQ Code

Examples:
  zedbarimg barcode.png
  zedbarimg --quiet qrcode.jpg
  zedbarimg image1.png image2.jpg image3.gif
`);
}

function main() {
  const options = parseArgs();
  let totalSymbols = 0;

  for (const filename of options.files) {
    // Read the image file as raw bytes
    let imageBytes;
    try {
      imageBytes = readFileSync(filename);
    } catch (err) {
      if (!options.quiet) {
        console.error(`Failed to read file '${filename}': ${err.message}`);
      }
      process.exit(1);
    }

    // Scan the image using the WASM library
    let results;
    try {
      results = zedbar.scanImageBytes(imageBytes);
    } catch (err) {
      if (!options.quiet) {
        console.error(`Failed to scan image '${filename}': ${err.message}`);
      }
      process.exit(1);
    }

    totalSymbols += results.length;

    // Print results
    for (const symbol of results) {
      if (options.raw) {
        // Raw mode: just output the data bytes with newline
        process.stdout.write(Buffer.from(symbol.data));
        process.stdout.write("\n");
      } else {
        // Normal mode: include symbol type
        const text = symbol.text || Buffer.from(symbol.data).toString("utf8");
        console.log(`${symbol.symbolType}:${text}`);
      }
    }
  }

  // Print statistics unless quiet mode
  if (!options.quiet) {
    if (totalSymbols === 0) {
      console.error("No barcodes found");
      process.exit(1);
    } else {
      console.error(
        `scanned ${totalSymbols} barcode symbol${totalSymbols === 1 ? "" : "s"} from ${options.files.length} image${options.files.length === 1 ? "" : "s"}`,
      );
    }
  } else if (totalSymbols === 0) {
    // In quiet mode, still exit with error if no barcodes found
    process.exit(1);
  }
}

main();
