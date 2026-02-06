#!/usr/bin/env node
/**
 * Test script to validate the zedbar WASM package works correctly.
 * Run with: node test.mjs
 */

import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import { test } from 'node:test'
import zedbar from './pkg/zedbar.js';
import sharp from 'sharp';
import assert from 'node:assert';

const __dirname = dirname(fileURLToPath(import.meta.url));
const examplesDir = join(__dirname, '..', 'examples');

// Test cases: [filename, expectedType, expectedData]
// Note: EAN-8 and UPC-A are reported as EAN-13 (with leading zeros) since they're subsets
const testCases = [
  ['test-qr.png', 'QR-Code', 'Hello, simplified zbar!\n'],
  ['test-ean13.png', 'EAN-13', '5901234123457'],
  ['test-ean8.png', 'EAN-13', '0000963850742'],  // EAN-8 padded to EAN-13
  ['test-upca.png', 'EAN-13', '0012345678905'],  // UPC-A padded to EAN-13
  ['test-code128.png', 'CODE-128', 'HELLO123\0'],
  ['test-code39.png', 'CODE-39', 'TEST123'],
  ['test-code93.png', 'CODE-93', 'CODE93'],
  ['test-codabar.png', 'Codabar', 'A40156B'],
  ['test-i25.png', 'I2/5', '1234567890'],
];

for (const [filename, expectedType, expectedData] of testCases) {
  const imagePath = join(examplesDir, filename);

  test(filename, async () => {
    // Load and convert image to grayscale
    const { data, info } = await sharp(imagePath)
      .grayscale()
      .raw()
      .toBuffer({ resolveWithObject: true });

    // Scan the image
    const results = zedbar.scanGrayscale(data, info.width, info.height);
    assert(results.length > 0, 'No barcodes found');

    const result = results[0];
    const actualType = result.symbolType;
    const actualData = result.text ?? Buffer.from(result.data).toString('utf-8');

    assert.equal(actualType, expectedType);
    assert.equal(actualData, expectedData);
  });

}
