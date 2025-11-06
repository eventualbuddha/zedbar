# Barcode Library Benchmarks

Comprehensive benchmarks comparing `zedbar` with other Rust barcode decoding libraries.

## Libraries Compared

- **zedbar** (this library) - Pure Rust port of ZBar
- **rqrr** - Pure Rust QR code decoder (always included)
- **zbars** - Rust bindings to the C ZBar library (optional)
- **rxing** - Rust port of ZXing (optional)

## Quick Start

### Basic Benchmark (zedbar vs rqrr)

```bash
# Full benchmark run (takes a few minutes)
cargo bench

# Quick run for testing (less accurate but faster)
cargo bench -- --quick
```

This runs the default benchmark comparing zedbar with rqrr on both QR codes and 1D barcodes.

**Note:** rqrr only supports QR codes, so only QR tests will run for it.

### Compare All Libraries

```bash
cargo bench --features bench_zbar_c,bench_rxing
```

### Compare Specific Libraries

```bash
# Compare with C library
cargo bench --features bench_zbar_c

# Compare with rxing
cargo bench --features bench_rxing
```

## Test Images

The benchmarks use various test images from the `examples/` directory:

### QR Codes

- `test-qr.png` - Simple QR code
- `qr-code-low-contrast.png` - Challenging low contrast
- `qr-code-capstone-interference.png` - QR with pattern interference
- `pixel-wifi-sharing-qr-code.png` - Large, high-resolution QR

### 1D Barcodes

- `test-ean13.png` - EAN-13
- `test-ean8.png` - EAN-8
- `test-code128.png` - Code 128
- `test-code39.png` - Code 39

**Note:** rqrr only supports QR codes, so it skip 1D barcode tests.

## Understanding Results

### Output Location

- Console: Summary statistics
- `target/criterion/` - Detailed HTML reports with graphs

### Viewing HTML Reports

```bash
# Open the report index
open target/criterion/report/index.html  # macOS
xdg-open target/criterion/report/index.html  # Linux
start target/criterion/report/index.html  # Windows
```

### Key Metrics

- **Time** - Lower is better (median decode time)
- **Throughput** - Higher is better (images per second)
- **Variance** - Lower is better (consistency)

### Initial Findings

Based on quick benchmarks, here are some observations:

**QR Codes:**

- **Simple QR**: zedbar (~730µs) is significantly faster than rqrr (~5.5ms)
- **Low Contrast**: rqrr (~520µs) is faster than zedbar (~2.6ms)
- **Capstone Interference**: rqrr (~690µs) is faster than zedbar (~2.4ms)
- **Large QR**: zedbar (~25ms) is faster than rqrr (~8.8ms)

**Trade-offs:**

- zedbar excels at simple, clean QR codes and large images
- rqrr is more robust for difficult/challenging QR codes
- zedbar supports many 1D formats that rqrr doesn't support

**1D Barcodes (zedbar only):**

- EAN-13: ~12.4ms
- EAN-8: ~12.6ms
- Code 128: ~13.1ms
- Code 39: ~14.7ms

These are baseline numbers - run your own benchmarks with your specific images!

### Interpreting Comparisons

Criterion automatically compares against previous runs:

- `No change in performance detected` - Within noise threshold
- `Performance has improved` / `regressed` - Statistically significant change
- `+X%` or `-X%` - Percentage change from baseline

## Adding More Test Images

To add your own test images:

1. Place images in `examples/` directory
2. Edit `benches/comparison.rs`
3. Add to `TEST_IMAGES` array:

```rust
TestImage {
    name: "my_test",
    path: "examples/my-barcode.png",
    expected_type: "QR",  // or "EAN13", "CODE128", etc.
},
```

## Customizing Benchmarks

### Measurement Time

Edit `comparison.rs` to add before `group.finish()`:

```rust
group.measurement_time(std::time::Duration::from_secs(10));
```

### Sample Size

```rust
group.sample_size(100);
```

### Warm-up Time

```rust
group.warm_up_time(std::time::Duration::from_secs(3));
```

## Benchmark Features

The benchmark code is feature-gated to only compile comparisons when requested:

- `bench_zbar_c` - Include C library via zbars
- `bench_rxing` - Include rxing

This keeps the default build fast and doesn't require installing optional dependencies.

## Comparing Against Your Own Changes

```bash
# Run baseline
cargo bench

# Make changes to code
# ...

# Run comparison (Criterion automatically compares)
cargo bench
```

Criterion saves baselines in `target/criterion/` and automatically compares new runs.

### Saving Named Baselines

```bash
# Save current performance as baseline
cargo bench -- --save-baseline main

# After changes, compare against it
cargo bench -- --baseline main
```

## Troubleshooting

### zbars fails to compile

The C library may not be installed. Install it:

```bash
# Ubuntu/Debian
sudo apt-get install libzbar-dev

# macOS
brew install zbar

# Fedora
sudo dnf install zbar-devel
```

Or skip it:

```bash
cargo bench  # Don't use --features bench_zbar_c
```

### Image not found errors

Ensure you're running from the repository root where `examples/` exists.

## Performance Tips

These are measured by the benchmarks:

1. **Image preprocessing** - Converting to grayscale is measured
2. **Scanner configuration** - Default settings are used
3. **Multiple passes** - Some libraries try multiple strategies

To isolate pure decode time, you'd need to modify the benchmark to exclude image loading/conversion.
