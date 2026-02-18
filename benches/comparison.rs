use std::hint::black_box;

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use image::DynamicImage;

// Test image structure
struct TestImage {
    name: &'static str,
    path: &'static str,
    expected_type: &'static str,
}

const TEST_IMAGES: &[TestImage] = &[
    // QR codes - all libraries can attempt these
    TestImage {
        name: "qr_simple",
        path: "examples/test-qr.png",
        expected_type: "QR",
    },
    TestImage {
        name: "qr_low_contrast",
        path: "examples/qr-code-low-contrast.png",
        expected_type: "QR",
    },
    TestImage {
        name: "qr_capstone_interference",
        path: "examples/qr-code-capstone-interference.png",
        expected_type: "QR",
    },
    TestImage {
        name: "qr_large",
        path: "examples/pixel-wifi-sharing-qr-code.png",
        expected_type: "QR",
    },
    // 1D barcodes - only zedbar, zbar-c, and rxing support these
    TestImage {
        name: "ean13",
        path: "examples/test-ean13.png",
        expected_type: "EAN13",
    },
    TestImage {
        name: "ean8",
        path: "examples/test-ean8.png",
        expected_type: "EAN8",
    },
    TestImage {
        name: "code128",
        path: "examples/test-code128.png",
        expected_type: "CODE128",
    },
    TestImage {
        name: "code39",
        path: "examples/test-code39.png",
        expected_type: "CODE39",
    },
];

fn load_image(path: &str) -> Option<DynamicImage> {
    image::open(path).ok()
}

fn benchmark_zedbar(c: &mut Criterion) {
    let mut group = c.benchmark_group("zedbar");

    for test_img in TEST_IMAGES {
        if let Some(img) = load_image(test_img.path) {
            let gray = img.to_luma8();
            let (width, height) = gray.dimensions();

            group.bench_with_input(
                BenchmarkId::from_parameter(test_img.name),
                &(gray, width, height),
                |b, (gray, width, height)| {
                    b.iter(|| {
                        let mut zbar_img = zedbar::Image::from_gray(
                            black_box(gray.as_raw()),
                            black_box(*width),
                            black_box(*height),
                        )
                        .unwrap();
                        let mut scanner = zedbar::Scanner::new();
                        let symbols = scanner.scan(&mut zbar_img);
                        black_box(symbols.len())
                    });
                },
            );
        }
    }

    group.finish();
}

fn benchmark_rqrr(c: &mut Criterion) {
    let mut group = c.benchmark_group("rqrr");

    for test_img in TEST_IMAGES {
        // rqrr only supports QR codes
        if !test_img.expected_type.starts_with("QR") {
            continue;
        }

        if let Some(img) = load_image(test_img.path) {
            let gray = img.to_luma8();

            group.bench_with_input(
                BenchmarkId::from_parameter(test_img.name),
                &gray,
                |b, gray| {
                    b.iter(|| {
                        let mut img = rqrr::PreparedImage::prepare(black_box(gray.clone()));
                        let grids = img.detect_grids();
                        black_box(grids.len())
                    });
                },
            );
        }
    }

    group.finish();
}

#[cfg(feature = "bench_zbar_c")]
fn benchmark_zbar_c(c: &mut Criterion) {
    let mut group = c.benchmark_group("zbar-c");

    for test_img in TEST_IMAGES {
        if let Some(img) = load_image(test_img.path) {
            let gray = img.to_luma8();
            let (width, height) = gray.dimensions();

            group.bench_with_input(
                BenchmarkId::from_parameter(test_img.name),
                &(gray, width, height),
                |b, (gray, width, height)| {
                    b.iter(|| {
                        let mut scanner = zbar_rust::ZBarImageScanner::new();
                        let results = scanner.scan_y800(black_box(gray.as_raw()), *width, *height);
                        black_box(results.map(|syms| syms.len()).unwrap_or(0))
                    });
                },
            );
        }
    }

    group.finish();
}

#[cfg(feature = "bench_rxing")]
fn benchmark_rxing(c: &mut Criterion) {
    use rxing::{BinaryBitmap, MultiFormatReader, Reader};

    let mut group = c.benchmark_group("rxing");

    for test_img in TEST_IMAGES {
        if let Some(img) = load_image(test_img.path) {
            group.bench_with_input(
                BenchmarkId::from_parameter(test_img.name),
                &img,
                |b, img| {
                    b.iter(|| {
                        let mut reader = MultiFormatReader::default();

                        // rxing expects DynamicImage
                        let lum_source =
                            rxing::BufferedImageLuminanceSource::new(black_box(img.clone()));
                        let mut bitmap =
                            BinaryBitmap::new(rxing::common::HybridBinarizer::new(lum_source));

                        let result = reader.decode(&mut bitmap);
                        black_box(result.is_ok())
                    });
                },
            );
        }
    }

    group.finish();
}

// Define benchmark groups based on enabled features
#[cfg(all(not(feature = "bench_zbar_c"), not(feature = "bench_rxing")))]
criterion_group!(benches, benchmark_zedbar, benchmark_rqrr);

#[cfg(all(feature = "bench_zbar_c", not(feature = "bench_rxing")))]
criterion_group!(benches, benchmark_zedbar, benchmark_rqrr, benchmark_zbar_c);

#[cfg(all(not(feature = "bench_zbar_c"), feature = "bench_rxing"))]
criterion_group!(benches, benchmark_zedbar, benchmark_rqrr, benchmark_rxing);

#[cfg(all(feature = "bench_zbar_c", feature = "bench_rxing"))]
criterion_group!(
    benches,
    benchmark_zedbar,
    benchmark_rqrr,
    benchmark_zbar_c,
    benchmark_rxing
);

criterion_main!(benches);
