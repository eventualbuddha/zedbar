use std::mem::{offset_of, size_of};
use zbar::decoder_types::*;

#[test]
fn test_decoder_struct_sizes() {
    // Verify struct sizes match C
    assert_eq!(
        size_of::<ean_decoder_t>(),
        136,
        "ean_decoder_t size mismatch"
    );
    assert_eq!(
        size_of::<i25_decoder_t>(),
        28,
        "i25_decoder_t size mismatch"
    );
    assert_eq!(
        size_of::<databar_decoder_t>(),
        40,
        "databar_decoder_t size mismatch"
    );
    assert_eq!(
        size_of::<codabar_decoder_t>(),
        32,
        "codabar_decoder_t size mismatch"
    );
    assert_eq!(
        size_of::<code39_decoder_t>(),
        24,
        "code39_decoder_t size mismatch"
    );
    assert_eq!(
        size_of::<code93_decoder_t>(),
        24,
        "code93_decoder_t size mismatch"
    );
    assert_eq!(
        size_of::<code128_decoder_t>(),
        24,
        "code128_decoder_t size mismatch"
    );
    assert_eq!(size_of::<qr_finder_t>(), 28, "qr_finder_t size mismatch");
    assert_eq!(size_of::<sq_finder_t>(), 4, "sq_finder_t size mismatch");
    assert_eq!(
        size_of::<zbar_decoder_t>(),
        464,
        "zbar_decoder_t size mismatch"
    );
}

#[test]
fn test_decoder_field_offsets() {
    // Verify field offsets match C
    assert_eq!(offset_of!(zbar_decoder_t, idx), 0);
    assert_eq!(offset_of!(zbar_decoder_t, w), 4);
    assert_eq!(offset_of!(zbar_decoder_t, type_), 68);
    assert_eq!(offset_of!(zbar_decoder_t, lock), 72);
    assert_eq!(offset_of!(zbar_decoder_t, modifiers), 76);
    assert_eq!(offset_of!(zbar_decoder_t, direction), 80);
    assert_eq!(offset_of!(zbar_decoder_t, s6), 84);
    assert_eq!(offset_of!(zbar_decoder_t, buf_alloc), 88);
    assert_eq!(offset_of!(zbar_decoder_t, buflen), 92);
    assert_eq!(offset_of!(zbar_decoder_t, buf), 96);
    assert_eq!(offset_of!(zbar_decoder_t, userdata), 104);
    assert_eq!(offset_of!(zbar_decoder_t, handler), 112);
    assert_eq!(offset_of!(zbar_decoder_t, ean), 120);
    assert_eq!(offset_of!(zbar_decoder_t, i25), 256);
    assert_eq!(offset_of!(zbar_decoder_t, databar), 288);
    assert_eq!(offset_of!(zbar_decoder_t, codabar), 328);
    assert_eq!(offset_of!(zbar_decoder_t, code39), 360);
    assert_eq!(offset_of!(zbar_decoder_t, code93), 384);
    assert_eq!(offset_of!(zbar_decoder_t, code128), 408);
    assert_eq!(offset_of!(zbar_decoder_t, qrf), 432);
    assert_eq!(offset_of!(zbar_decoder_t, sqf), 460);
}
