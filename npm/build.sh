#!/usr/bin/env bash
set -euo pipefail

# Build the WASM package for npm publishing.
# Prerequisites: wasm-pack (cargo install wasm-pack)

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

echo "Building zedbar WASM for npm..."
wasm-pack build "$PROJECT_DIR" \
  --target nodejs \
  --out-dir "$SCRIPT_DIR/pkg" \
  --no-default-features \
  --features wasm

# Customize the generated package.json
node -e "
const fs = require('fs');
const pkg = JSON.parse(fs.readFileSync('$SCRIPT_DIR/pkg/package.json', 'utf8'));

// Update package metadata
pkg.name = 'zedbar';
pkg.description = 'Fast QR code and barcode scanner for Node.js, powered by WebAssembly';
pkg.repository = {
  type: 'git',
  url: 'https://github.com/eventualbuddha/zedbar.git'
};
pkg.keywords = ['qrcode', 'qr-code', 'barcode', 'scanner', 'decoder', 'wasm', 'webassembly', 'ean', 'code128', 'code39'];
pkg.author = 'eventualbuddha';
pkg.license = 'LGPL-3.0-or-later';
pkg.homepage = 'https://github.com/eventualbuddha/zedbar';
pkg.bugs = { url: 'https://github.com/eventualbuddha/zedbar/issues' };

// Add engines requirement
pkg.engines = { node: '>=16.0.0' };

fs.writeFileSync('$SCRIPT_DIR/pkg/package.json', JSON.stringify(pkg, null, 2) + '\n');
"

# Replace auto-generated README with npm-specific one
cp "$SCRIPT_DIR/README.md" "$SCRIPT_DIR/pkg/README.md"

# Remove unnecessary files
rm -f "$SCRIPT_DIR/pkg/.gitignore"

echo ""
echo "Build complete. Package ready in: $SCRIPT_DIR/pkg/"
echo ""
echo "To test: cd $SCRIPT_DIR && npm test"
echo "To publish: cd $SCRIPT_DIR/pkg && npm publish"
