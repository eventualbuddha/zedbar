use std::mem::size_of;
use zbar::{
    decoder::{zbar_decoder_create, zbar_decoder_destroy},
    decoder_types::*,
};

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
fn test_decoder_reset_preserves_heap_allocations() {
    unsafe {
        let decoder = zbar_decoder_create();
        assert!(!decoder.is_null(), "decoder allocation failed");

        let buf_ptr = (*decoder).buf;
        let segs_ptr = (*decoder).databar.segs;

        assert!(!buf_ptr.is_null(), "decoder buffer not allocated");
        assert!(!segs_ptr.is_null(), "databar segment array not allocated");

        (*decoder).reset();

        assert_eq!((*decoder).buf, buf_ptr, "reset should not free buffer");
        assert_eq!(
            (*decoder).databar.segs,
            segs_ptr,
            "reset should not free databar segments"
        );

        zbar_decoder_destroy(decoder);
    }
}
