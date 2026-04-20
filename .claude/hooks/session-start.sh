#!/bin/bash
# Installs the system zbar library and zbarimg tool. These are needed by:
#   - tests/decode_examples.rs (compares zedbar's output against zbarimg)
#   - benches/comparison.rs     (links against the zbar C library)
#   - cargo clippy --all-targets --all-features (pulls in the bench)
set -euo pipefail

# Only run in the Claude Code on the web sandbox. Locally the user manages
# their own system packages.
if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

if command -v zbarimg >/dev/null 2>&1 && dpkg -s libzbar-dev >/dev/null 2>&1; then
  # Already installed; nothing to do.
  exit 0
fi

# Run the install in the background so session startup isn't blocked.
echo '{"async": true, "asyncTimeout": 300000}'

SUDO=""
if [ "$(id -u)" -ne 0 ]; then
  SUDO="sudo"
fi

export DEBIAN_FRONTEND=noninteractive
# `apt-get update` can fail transiently on unrelated third-party repositories
# in the sandbox image. A failure there shouldn't block installing our packages
# — apt will use whatever cached index is available.
$SUDO apt-get update -qq || true
$SUDO apt-get install -y --no-install-recommends \
  libzbar-dev \
  libzbar0 \
  zbar-tools
