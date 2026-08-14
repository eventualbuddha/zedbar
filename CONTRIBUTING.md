# Contributing to zedbar

## Testing

```sh
cargo test
```

Most of the decoders are a close port of the [ZBar][zbar] C library, so the
tests lean on two external references rather than on this crate's own idea of
what is correct:

- **`zbarimg`** (Debian `zbar-tools`, `brew install zbar`) is the behavioral
  oracle. `tests/decode_examples.rs` decodes each fixture with both and
  asserts they agree. The cross-checks are skipped when `zbarimg` is not
  installed, so install it to get the full suite.
- **`rqrr`** gives an independent second opinion on QR codes, and is a normal
  dev-dependency.

Every symbology has at least one end-to-end fixture in `examples/`. Fixtures
are generated with [`zint`][zint] where possible; the two that aren't have
generators alongside them (`cargo run --example generate_...`).

## Coverage

Line and region coverage comes from [`cargo-llvm-cov`][llvm-cov], which drives
LLVM's source-based instrumentation and so counts the same code the compiler
sees:

```sh
cargo install cargo-llvm-cov
rustup component add llvm-tools-preview

cargo llvm-cov --summary-only          # per-file table
cargo llvm-cov --open                  # annotated HTML, opens a browser
cargo llvm-cov --html --output-dir cov # annotated HTML, no browser
```

`--summary-only` is the quick check; the HTML report is what to use when
adding tests, since it marks the individual uncovered lines.

Coverage is a guide, not a target — a decoder can sit at 90% while the
handful of uncovered lines are exactly the interesting ones. It is most
useful for spotting whole functions or error paths that no test reaches.

[zbar]: https://github.com/mchehab/zbar
[zint]: https://www.zint.org.uk/
[llvm-cov]: https://github.com/taiki-e/cargo-llvm-cov

## Commit Message Format

This project uses [Conventional Commits](https://www.conventionalcommits.org/) for automated releases via release-please.

**Format**: `<type>: <description>`

### Types:
- `feat:` - New feature (triggers **minor** version bump: 0.2.0 → 0.3.0)
- `fix:` - Bug fix (triggers **patch** version bump: 0.2.0 → 0.2.1)
- `feat!:` or `fix!:` - Breaking change (triggers **major** version bump: 0.2.0 → 1.0.0)
- `chore:` - Maintenance (no release)
- `docs:` - Documentation only (no release)
- `test:` - Tests only (no release)
- `refactor:` - Code refactoring (no release)

### Examples:
```
feat: add WebP support to CLI
fix: handle empty image files correctly
feat!: remove deprecated scanImage function
chore: update dependencies
docs: add CLI usage examples
```

### Breaking Changes:
Add `!` after type OR include `BREAKING CHANGE:` in commit body:
```
feat!: redesign API interface

BREAKING CHANGE: scanImage() renamed to scanImageBytes()
```

## Release Process

1. Commit with the conventional format in a feature branch and open a PR targeting `main`
2. When PRs are merged into `main`, release-please automatically opens/updates a single release PR that includes all unreleased conventional commits since the last release
3. Merge the release PR → auto-publishes to crates.io and npm

**Note**: `feat:` and `fix:` commits are batched together into the next release based on what has been merged into `main`, rather than creating a separate release for each individual commit. Keep commits focused and atomic so the generated changelog remains clear.
