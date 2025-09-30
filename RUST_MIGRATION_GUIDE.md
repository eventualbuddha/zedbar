# ZBar Rust Migration Guide

This guide outlines the process for converting the ZBar C library to Rust using `c2rust`.

## Project Setup Complete ✓

The basic Rust project structure is now in place:

```
├── Cargo.toml          # Rust project configuration
├── build.rs            # Build script (compiles C code for now)  
├── src/
│   ├── lib.rs          # Main library entry point
│   ├── ffi.rs          # FFI bindings to C code
│   ├── error.rs        # Rust error types
│   ├── image.rs        # Image handling API
│   ├── symbol.rs       # Symbol/barcode result types
│   ├── scanner.rs      # Image scanner API
│   ├── decoder.rs      # Low-level decoder (placeholder)
│   └── bin/
│       └── zbarimg.rs  # Command-line tool
```

## Next Steps for C2Rust Conversion

### Phase 1: Install and Setup c2rust

```bash
# Install c2rust (requires Rust nightly)
rustup toolchain install nightly
cargo +nightly install c2rust

# Or install from source
git clone https://github.com/immunant/c2rust
cd c2rust
cargo install --path c2rust-translate
```

### Phase 2: Generate Rust Code Module by Module

Start with the smallest, most self-contained modules:

```bash
# 1. Start with utility modules (39 lines)
c2rust transpile compile_commands.json -- --translate-const-macros \
  --emit-modules --output-dir rust-modules zbar/refcnt.c

# 2. Error handling (152 lines)  
c2rust transpile compile_commands.json -- --translate-const-macros \
  --emit-modules --output-dir rust-modules zbar/error.c

# 3. Configuration (168 lines)
c2rust transpile compile_commands.json -- --translate-const-macros \
  --emit-modules --output-dir rust-modules zbar/config.c

# 4. Core data structures
c2rust transpile compile_commands.json -- --translate-const-macros \
  --emit-modules --output-dir rust-modules zbar/image.c zbar/symbol.c

# 5. Scanners and decoders
c2rust transpile compile_commands.json -- --translate-const-macros \
  --emit-modules --output-dir rust-modules zbar/scanner.c zbar/decoder.c
```

### Phase 3: Gradual Integration Strategy

1. **Replace FFI bindings**: As each module is converted, replace the corresponding FFI bindings in `src/ffi.rs`

2. **Refactor unsafe code**: c2rust generates unsafe Rust. Gradually make it safe:
   - Replace raw pointers with Box, Rc, Arc as appropriate
   - Replace malloc/free with Rust's allocator
   - Add proper lifetimes and borrowing

3. **Improve APIs**: Transform C-style APIs into idiomatic Rust:
   - Use Result<T, E> instead of error codes
   - Use Iterator traits instead of manual loops
   - Use Option<T> instead of null pointers

### Phase 4: Module Replacement Order (Recommended)

Based on dependencies and complexity:

1. **zbar/refcnt.c** → Custom reference counting (39 lines)
2. **zbar/error.c** → Error handling (152 lines) 
3. **zbar/misc.c** → Utility functions (105 lines)
4. **zbar/config.c** → Configuration (168 lines)
5. **zbar/symbol.c** → Symbol data structures (478 lines)
6. **zbar/image.c** → Image handling (195 lines)
7. **zbar/scanner.c** → Line scanner (291 lines)
8. **zbar/decoder.c** → Main decoder (522 lines)
9. **Simple decoders**:
   - zbar/decoder/sq_finder.c (7 lines)
   - zbar/decoder/qr_finder.c (100 lines)
   - zbar/decoder/i25.c (265 lines)
   - zbar/decoder/code39.c (336 lines)
10. **Complex decoders**:
    - zbar/decoder/ean.c (735 lines)
    - zbar/decoder/code128.c (582 lines)  
    - zbar/decoder/databar.c (1240 lines)
11. **QR Code implementation** (most complex):
    - zbar/qrcode/*.c (6000+ lines total)

## Testing Strategy

1. **Unit tests**: Each converted module should have comprehensive tests
2. **Integration tests**: Test against known barcode images
3. **Benchmarks**: Ensure performance is maintained
4. **Compatibility**: Ensure C API compatibility is maintained during transition

## C2Rust Commands Reference

To start the conversion, you'll need a `compile_commands.json` file:

```bash
# Generate compile_commands.json using bear
bear -- make clean all

# Or manually create it based on your Makefile
# Then run c2rust on specific files
c2rust transpile compile_commands.json -- \
  --translate-const-macros \
  --emit-modules \
  --output-dir rust-converted \
  zbar/refcnt.c
```

This approach allows for incremental migration while maintaining a working library throughout the process.