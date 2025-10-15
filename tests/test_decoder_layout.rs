use std::mem::size_of;
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
