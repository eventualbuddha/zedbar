#!/usr/bin/env bash
set -euo pipefail

# Build the WASM package and copy it into the demo site.
# Prerequisites: wasm-pack (cargo install wasm-pack)

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

echo "Building zedbar WASM..."
wasm-pack build "$PROJECT_DIR" \
  --target web \
  --out-dir "$SCRIPT_DIR/pkg" \
  --no-default-features \
  --features wasm

# wasm-pack generates files we don't need for the demo
rm -f "$SCRIPT_DIR/pkg/.gitignore" "$SCRIPT_DIR/pkg/package.json" "$SCRIPT_DIR/pkg/README.md"

echo ""
echo "Build complete. To preview locally:"
echo "  cd $SCRIPT_DIR && python3 -m http.server 8080"
echo ""
echo "For Cloudflare Pages, set:"
echo "  Build command: demo/build.sh"
echo "  Build output:  demo/"
