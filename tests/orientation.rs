//! Integration test for orientation detection
//!
//! This test verifies that all four orientation variants can be produced
//! by scanning barcode images rotated to different orientations.

use std::collections::HashSet;
use zbar::{Image, Orientation, Scanner};

#[test]
fn test_orientation_all_variants() {
    // Load a simple 1D barcode (Code 128)
    let img = image::open("examples/test-code128.png").expect("Failed to load test image");

    // Test all four rotations using image crate's built-in rotation
    struct TestCase {
        name: &'static str,
        rotation: &'static str,
        expected: Orientation,
    }

    let test_cases = vec![
        TestCase {
            name: "0 degrees",
            rotation: "none",
            expected: Orientation::Up,
        },
        TestCase {
            name: "90 degrees counter-clockwise",
            rotation: "90",
            expected: Orientation::Down,
        },
        TestCase {
            name: "180 degrees",
            rotation: "180",
            expected: Orientation::Right,
        },
        TestCase {
            name: "270 degrees counter-clockwise (90 degrees clockwise)",
            rotation: "270",
            expected: Orientation::Left,
        },
    ];

    let mut orientations_found = HashSet::new();

    for test_case in test_cases {
        // Rotate the image
        let rotated = match test_case.rotation {
            "none" => img.to_luma8(),
            "90" => img.rotate90().to_luma8(),
            "180" => img.rotate180().to_luma8(),
            "270" => img.rotate270().to_luma8(),
            _ => panic!("Unknown rotation"),
        };

        let mut scanner = Scanner::new();

        let mut zbar_image = Image::from_gray(rotated.as_raw(), rotated.width(), rotated.height())
            .expect("Failed to create zbar image");

        let symbols = scanner.scan(&mut zbar_image);

        assert!(
            !symbols.is_empty(),
            "Failed to decode barcode at orientation: {}",
            test_case.name
        );

        let symbol = &symbols[0];
        let orientation = symbol.orientation();

        orientations_found.insert(orientation);

        assert_eq!(
            orientation, test_case.expected,
            "Wrong orientation for {}: got {:?}, expected {:?}",
            test_case.name, orientation, test_case.expected
        );
    }

    // Verify all four orientations were produced
    assert_eq!(
        orientations_found.len(),
        4,
        "Not all four orientations were produced: {:?}",
        orientations_found
    );

    // Verify we have each specific orientation
    assert!(orientations_found.contains(&Orientation::Up));
    assert!(orientations_found.contains(&Orientation::Right));
    assert!(orientations_found.contains(&Orientation::Down));
    assert!(orientations_found.contains(&Orientation::Left));
}
