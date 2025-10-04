fn main() {
    // For now, just build the C library as a fallback
    // Later we'll remove this as modules are converted to Rust

    let mut build = cc::Build::new();

    build
        .files([
            "zbar/decoder.c",
            "zbar/error.c",
            "zbar/image.c",
            "zbar/img_scanner.c",
            "zbar/processor.c",
            "zbar/scanner.c",
            // "zbar/sqcode.c", // Converted to Rust
            // "zbar/symbol.c", // Converted to Rust
            // Decoder modules
            "zbar/decoder/codabar.c",
            "zbar/decoder/code128.c",
            "zbar/decoder/code39.c",
            "zbar/decoder/code93.c",
            "zbar/decoder/databar.c",
            // "zbar/decoder/ean.c", // Converted to Rust
            // "zbar/decoder/i25.c", // Converted to Rust
            // "zbar/decoder/qr_finder.c", // Converted to Rust
            // "zbar/decoder/sq_finder.c", // Converted to Rust
            // QR code modules
            // "zbar/qrcode/bch15_5.c", // Converted to Rust
            // "zbar/qrcode/binarize.c", // Converted to Rust
            // "zbar/qrcode/isaac.c", // Converted to Rust
            "zbar/qrcode/qrdec.c",
            // "zbar/qrcode/qrdectxt.c", // Converted to Rust
            "zbar/qrcode/rs.c",
            // "zbar/qrcode/util.c", // Converted to Rust
        ])
        .include("zbar")
        .compile("zbar_c");

    println!("cargo:rustc-link-lib=static=zbar_c");
    println!("cargo:rerun-if-changed=zbar/");
}
