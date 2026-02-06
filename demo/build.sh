#!/usr/bin/env bash
set -euo pipefail

# Build the WASM package and HTML for the demo site.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

# Install Rust if not available (needed for Cloudflare Pages)
if ! command -v cargo &>/dev/null; then
  echo "Installing Rust..."
  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal
  # shellcheck source=/dev/null
  source "$HOME/.cargo/env"
fi

# Install wasm-pack if not available
if ! command -v wasm-pack &>/dev/null; then
  echo "Installing wasm-pack..."
  cargo install wasm-pack
fi

echo "Installing demo dependencies..."
npm ci --prefix "$SCRIPT_DIR"

echo ""
echo "Building HTML with syntax highlighting..."
npm run --prefix "$SCRIPT_DIR" build:html

echo ""
echo "Building zedbar WASM for web..."
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
