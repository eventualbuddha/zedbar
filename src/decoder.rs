//! Decoder type definitions
//!
//! This module contains type definitions for all the barcode decoder types.

/// Bounding box (x, y, width, height) in pixel coordinates.
pub(crate) type BBox = (u32, u32, u32, u32);

#[cfg(feature = "databar")]
use crate::color::Color;

/// Window size for bar width history (must be power of 2)
pub(crate) const DECODE_WINDOW: usize = 16;

/// Fixed character width decode assist
///
/// Bar+space width are compared as a fraction of the reference dimension "x"
/// - +/- 1/2 x tolerance
/// - measured total character width (s) compared to symbology baseline (n)
/// - bar+space *pair width* "e" is used to factor out bad "exposures"
///
/// Returns encoded number of units - 2 (for use as zero based index)
/// or -1 if invalid
pub(crate) fn decode_e(e: u32, s: u32, n: u32) -> i32 {
    // A reference character of zero width is not a character, so nothing can
    // be classified against it. zbar divides by `s` unguarded and every caller
    // is expected to have ruled this out — but `ean::aux_end` computes its own
    // reference width and does not, so a run of zero-width elements reaches
    // here and takes the C process down with SIGFPE. No image produces that
    // run (the scanner never emits four zero widths in a row), but the
    // decoder accepts it, and "invalid" is the answer the classifier already
    // has for input it cannot use.
    if s == 0 {
        return -1;
    }
    let e_val = ((e * n * 2 + 1) / s).wrapping_sub(3) / 2;
    if e_val >= n - 3 { -1 } else { e_val as i32 }
}

/// Barcode data modifiers (GS1 and AIM identifiers)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Modifier {
    Gs1 = 0,
    Aim = 1,
}

impl Modifier {
    /// Convert to bit flag for use in bitfield
    pub(crate) fn bit(self) -> u32 {
        1 << (self as u32)
    }
}

// ============================================================================
// Simple decoder types
// ============================================================================

/// Interleaved 2 of 5 decoder state
#[cfg(feature = "i25")]
#[derive(Default)]
pub(crate) struct I25Decoder {
    // Bitfields packed into first 32 bits:
    // direction: 1 bit, element: 4 bits, character: 12 bits = 17 bits used
    // We'll use a u32 and provide accessor methods
    bitfields: u32,
    pub(crate) s10: u32,
    pub(crate) width: u32,
    pub(crate) buffer: Vec<u8>,
}

#[cfg(feature = "i25")]
impl I25Decoder {
    pub(crate) fn direction(&self) -> bool {
        (self.bitfields & 0x1) != 0
    }

    pub(crate) fn set_direction(&mut self, val: bool) {
        self.bitfields = (self.bitfields & !0x1) | (val as u32);
    }

    pub(crate) fn element(&self) -> u8 {
        ((self.bitfields >> 1) & 0xF) as u8
    }

    pub(crate) fn set_element(&mut self, val: u8) {
        self.bitfields = (self.bitfields & !(0xF << 1)) | ((val as u32 & 0xF) << 1);
    }

    pub(crate) fn character(&self) -> i16 {
        // Sign extend the 12-bit value
        let val = ((self.bitfields >> 5) & 0xFFF) as i16;
        // If the sign bit (bit 11) is set, extend it
        if val & 0x800 != 0 { val | !0xFFF } else { val }
    }

    pub(crate) fn set_character(&mut self, val: i16) {
        self.bitfields = (self.bitfields & !(0xFFF << 5)) | (((val as u32) & 0xFFF) << 5);
    }

    /// Reset i25 decoder state
    pub(crate) fn reset(&mut self) {
        self.set_direction(false);
        self.set_element(0);
        self.set_character(-1);
        self.s10 = 0;
    }

    pub(crate) fn set_byte(&mut self, index: usize, value: u8) {
        if self.buffer.len() <= index {
            self.buffer.resize(index + 1, 0);
        }
        self.buffer[index] = value;
    }
}

/// Code 39 decoder state
#[cfg(feature = "code39")]
#[derive(Default)]
pub(crate) struct Code39Decoder {
    // Bitfields: direction: 1, element: 4, character: 12
    bitfields: u32,
    pub s9: u32,
    pub width: u32,
    pub buffer: Vec<u8>,
}

#[cfg(feature = "code39")]
impl Code39Decoder {
    pub(crate) fn direction(&self) -> bool {
        (self.bitfields & 0x1) != 0
    }

    pub(crate) fn set_direction(&mut self, val: bool) {
        self.bitfields = (self.bitfields & !0x1) | (val as u32);
    }

    pub(crate) fn element(&self) -> u8 {
        ((self.bitfields >> 1) & 0xF) as u8
    }

    pub(crate) fn set_element(&mut self, val: u8) {
        self.bitfields = (self.bitfields & !(0xF << 1)) | ((val as u32 & 0xF) << 1);
    }

    pub(crate) fn character(&self) -> i16 {
        // Sign extend the 12-bit value
        let val = ((self.bitfields >> 5) & 0xFFF) as i16;
        // If the sign bit (bit 11) is set, extend it
        if val & 0x800 != 0 { val | !0xFFF } else { val }
    }

    pub(crate) fn set_character(&mut self, val: i16) {
        self.bitfields = (self.bitfields & !(0xFFF << 5)) | (((val as u32) & 0xFFF) << 5);
    }

    /// Reset code39 decoder state
    pub(crate) fn reset(&mut self) {
        self.set_direction(false);
        self.set_element(0);
        self.set_character(-1);
        self.s9 = 0;
    }

    pub(crate) fn set_byte(&mut self, index: usize, value: u8) {
        if self.buffer.len() <= index {
            self.buffer.resize(index + 1, 0);
        }
        self.buffer[index] = value;
    }
}

/// Code 93 decoder state
#[cfg(feature = "code93")]
#[derive(Default)]
pub(crate) struct Code93Decoder {
    // Bitfields: direction: 1, element: 3, character: 12
    bitfields: u32,
    pub(crate) width: u32,
    pub(crate) buf: u8,
}

#[cfg(feature = "code93")]
impl Code93Decoder {
    pub(crate) fn direction(&self) -> bool {
        (self.bitfields & 0x1) != 0
    }

    pub(crate) fn set_direction(&mut self, val: bool) {
        self.bitfields = (self.bitfields & !0x1) | (val as u32);
    }

    pub(crate) fn element(&self) -> u8 {
        ((self.bitfields >> 1) & 0x7) as u8
    }

    pub(crate) fn set_element(&mut self, val: u8) {
        self.bitfields = (self.bitfields & !(0x7 << 1)) | ((val as u32 & 0x7) << 1);
    }

    pub(crate) fn character(&self) -> i16 {
        // Sign extend the 12-bit value
        let val = ((self.bitfields >> 4) & 0xFFF) as i16;
        // If the sign bit (bit 11) is set, extend it
        if val & 0x800 != 0 { val | !0xFFF } else { val }
    }

    pub(crate) fn set_character(&mut self, val: i16) {
        self.bitfields = (self.bitfields & !(0xFFF << 4)) | (((val as u32) & 0xFFF) << 4);
    }

    /// Reset code93 decoder state
    pub(crate) fn reset(&mut self) {
        self.set_direction(false);
        self.set_element(0);
        self.set_character(-1);
    }
}

/// Codabar decoder state
#[cfg(feature = "codabar")]
#[derive(Default)]
pub(crate) struct CodabarDecoder {
    // Bitfields: direction: 1, element: 4, character: 12
    bitfields: u32,
    pub(crate) s7: u32,
    pub(crate) width: u32,
    pub(crate) buf: [u8; 6],
}

#[cfg(feature = "codabar")]
impl CodabarDecoder {
    pub(crate) fn direction(&self) -> bool {
        (self.bitfields & 0x1) != 0
    }

    pub(crate) fn set_direction(&mut self, val: bool) {
        self.bitfields = (self.bitfields & !0x1) | (val as u32);
    }

    pub(crate) fn element(&self) -> u8 {
        ((self.bitfields >> 1) & 0xF) as u8
    }

    pub(crate) fn set_element(&mut self, val: u8) {
        self.bitfields = (self.bitfields & !(0xF << 1)) | ((val as u32 & 0xF) << 1);
    }

    pub(crate) fn character(&self) -> i16 {
        // Sign extend the 12-bit value
        let val = ((self.bitfields >> 5) & 0xFFF) as i16;
        // If the sign bit (bit 11) is set, extend it
        if val & 0x800 != 0 { val | !0xFFF } else { val }
    }

    pub(crate) fn set_character(&mut self, val: i16) {
        self.bitfields = (self.bitfields & !(0xFFF << 5)) | (((val as u32) & 0xFFF) << 5);
    }

    /// Reset codabar decoder state
    pub(crate) fn reset(&mut self) {
        self.set_direction(false);
        self.set_element(0);
        self.set_character(-1);
        self.s7 = 0;
    }
}

/// Code 128 decoder state
#[cfg(feature = "code128")]
#[derive(Default)]
pub(crate) struct Code128Decoder {
    // Bitfields: direction: 1, element: 3, character: 12 (16 bits)
    // start: 8 bits - packed into same u32
    // Total: 24 bits used in first u32
    bitfields_and_start: u32,
    pub(crate) s6: u32,
    pub(crate) width: u32,
}

#[cfg(feature = "code128")]
impl Code128Decoder {
    pub(crate) fn direction(&self) -> u8 {
        (self.bitfields_and_start & 0x1) as u8
    }

    pub(crate) fn set_direction(&mut self, val: u8) {
        self.bitfields_and_start = (self.bitfields_and_start & !0x1) | (val as u32 & 0x1);
    }

    pub(crate) fn element(&self) -> u8 {
        ((self.bitfields_and_start >> 1) & 0x7) as u8
    }

    pub(crate) fn set_element(&mut self, val: u8) {
        self.bitfields_and_start =
            (self.bitfields_and_start & !(0x7 << 1)) | ((val as u32 & 0x7) << 1);
    }

    pub(crate) fn character(&self) -> i16 {
        // Sign extend the 12-bit value
        let val = ((self.bitfields_and_start >> 4) & 0xFFF) as i16;
        // If the sign bit (bit 11) is set, extend it
        if val & 0x800 != 0 { val | !0xFFF } else { val }
    }

    pub(crate) fn set_character(&mut self, val: i16) {
        self.bitfields_and_start =
            (self.bitfields_and_start & !(0xFFF << 4)) | (((val as u32) & 0xFFF) << 4);
    }

    pub(crate) fn start(&self) -> u8 {
        ((self.bitfields_and_start >> 16) & 0xFF) as u8
    }

    pub(crate) fn set_start(&mut self, val: u8) {
        self.bitfields_and_start =
            (self.bitfields_and_start & !(0xFF << 16)) | ((val as u32) << 16);
    }

    /// Reset code128 decoder state
    pub(crate) fn reset(&mut self) {
        self.set_direction(0);
        self.set_element(0);
        self.set_character(-1);
        self.s6 = 0;
    }
}

// ============================================================================
// Complex decoder types
// ============================================================================

/// DataBar segment (partial)
#[cfg(feature = "databar")]
#[derive(Clone)]
pub(crate) struct DatabarSegment {
    // First 32 bits of bitfields:
    // finder: 5, exp: 1, color: 1, side: 1,
    // partial: 1, count: 7, epoch: 8, check: 8 = 32 bits
    bitfields: u32,
    /// Decoded character data (C: `signed short data`).
    pub(crate) data: i16,
    /// Measured width of the finder, 14 modules (C: `unsigned short width`).
    pub(crate) width: u16,
}

#[cfg(feature = "databar")]
impl DatabarSegment {
    pub(crate) fn finder(&self) -> i8 {
        // finder is a signed 5-bit field (bits 0-4)
        let val = (self.bitfields & 0x1F) as i8;
        // Sign extend from 5 bits to 8 bits
        if val & 0x10 != 0 {
            val | 0xE0_u8 as i8
        } else {
            val
        }
    }

    pub(crate) fn set_finder(&mut self, val: i8) {
        self.bitfields = (self.bitfields & !0x1F) | ((val as u32) & 0x1F);
    }

    pub(crate) fn partial(&self) -> bool {
        // partial is bit 8
        (self.bitfields & (1 << 8)) != 0
    }

    pub(crate) fn set_partial(&mut self, val: bool) {
        if val {
            self.bitfields |= 1 << 8;
        } else {
            self.bitfields &= !(1 << 8);
        }
    }

    pub(crate) fn exp(&self) -> bool {
        // exp is bit 5
        (self.bitfields & (1 << 5)) != 0
    }

    pub(crate) fn set_exp(&mut self, val: bool) {
        if val {
            self.bitfields |= 1 << 5;
        } else {
            self.bitfields &= !(1 << 5);
        }
    }

    pub(crate) fn color(&self) -> Color {
        // color is bit 6
        (((self.bitfields >> 6) & 1) as u8).into()
    }

    pub(crate) fn set_color(&mut self, val: Color) {
        self.bitfields = (self.bitfields & !(1 << 6)) | ((val as u32) << 6);
    }

    pub(crate) fn side(&self) -> u8 {
        // side is bit 7
        ((self.bitfields >> 7) & 1) as u8
    }

    pub(crate) fn set_side(&mut self, val: u8) {
        self.bitfields = (self.bitfields & !(1 << 7)) | (((val & 1) as u32) << 7);
    }

    pub(crate) fn count(&self) -> u8 {
        // count is bits 9-15 (7 bits)
        ((self.bitfields >> 9) & 0x7F) as u8
    }

    pub(crate) fn set_count(&mut self, val: u8) {
        self.bitfields = (self.bitfields & !(0x7F << 9)) | (((val & 0x7F) as u32) << 9);
    }

    pub(crate) fn epoch(&self) -> u8 {
        // epoch is bits 16-23 (8 bits)
        ((self.bitfields >> 16) & 0xFF) as u8
    }

    pub(crate) fn set_epoch(&mut self, val: u8) {
        self.bitfields = (self.bitfields & !(0xFF << 16)) | ((val as u32) << 16);
    }

    pub(crate) fn check(&self) -> u8 {
        // check is bits 24-31 (8 bits)
        ((self.bitfields >> 24) & 0xFF) as u8
    }

    pub(crate) fn set_check(&mut self, val: u8) {
        self.bitfields = (self.bitfields & !(0xFF << 24)) | ((val as u32) << 24);
    }

    pub(crate) fn segment_index(&self) -> i32 {
        ((self.finder() as i32) << 2)
            | ((self.color() as i32) << 1)
            | (((self.color() as u8 ^ self.side()) as i32) & 1)
    }
}

#[cfg(feature = "databar")]
impl Default for DatabarSegment {
    fn default() -> Self {
        let mut seg = Self {
            bitfields: 0,
            data: 0,
            width: 0,
        };
        seg.set_finder(-1);
        seg
    }
}

/// DataBar decoder state
#[cfg(feature = "databar")]
pub(crate) struct DatabarDecoder {
    epoch: u8,
    pub(crate) segs: Vec<DatabarSegment>,
    chars: [i8; 16],
}

#[cfg(feature = "databar")]
impl Default for DatabarDecoder {
    fn default() -> Self {
        Self {
            epoch: 0,
            segs: vec![Default::default(); 4],
            chars: [-1; 16],
        }
    }
}

#[cfg(feature = "databar")]
impl DatabarDecoder {
    pub(crate) fn csegs(&self) -> usize {
        self.segs.len()
    }

    pub(crate) fn seg(&self, index: usize) -> &DatabarSegment {
        &self.segs[index]
    }

    pub(crate) fn seg_mut(&mut self, index: usize) -> &mut DatabarSegment {
        &mut self.segs[index]
    }

    pub(crate) fn char(&self, index: usize) -> i8 {
        self.chars[index]
    }

    pub(crate) fn set_char(&mut self, index: usize, value: i8) {
        self.chars[index] = value;
    }

    pub(crate) fn resize_segs(&mut self, size: usize) {
        self.segs.resize_with(size, DatabarSegment::default);
    }

    pub(crate) fn epoch(&self) -> u8 {
        self.epoch
    }

    pub(crate) fn set_epoch(&mut self, val: u8) {
        self.epoch = val;
    }

    /// Reset DataBar decoder state
    pub(crate) fn reset(&mut self) {
        let n = self.segs.len();
        self.new_scan();
        for i in 0..n {
            let seg = self.seg_mut(i);
            seg.set_finder(-1);
        }
    }

    /// Prepare DataBar decoder for new scan
    pub(crate) fn new_scan(&mut self) {
        for i in 0..16 {
            if self.chars[i] >= 0 {
                let seg = &mut self.segs[self.chars[i] as usize];
                if seg.partial() {
                    seg.set_finder(-1);
                }
                self.chars[i] = -1;
            }
        }
    }
}

/// QR finder line (from qrcode.h)
#[cfg(feature = "qrcode")]
#[derive(Default, Copy, Clone)]
pub(crate) struct QrFinderLine {
    pub(crate) pos: [i32; 2], // qr_point
    pub(crate) len: i32,
    pub(crate) boffs: i32,
    pub(crate) eoffs: i32,
}

/// QR Code finder state
#[cfg(feature = "qrcode")]
#[derive(Default)]
pub(crate) struct QrFinder {
    pub(crate) s5: u32,
    pub(crate) line: QrFinderLine,
}

#[cfg(feature = "qrcode")]
impl QrFinder {
    /// Reset QR finder state
    pub(crate) fn reset(&mut self) {
        self.s5 = 0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `decode_e` classifies a bar+space pair width against the character
    /// width. zbar computes it as
    ///
    /// ```c
    /// unsigned char E = ((e * n * 2 + 1) / s - 3) / 2;
    /// return E >= n - 3 ? -1 : E;
    /// ```
    ///
    /// The values below were taken from that C expression directly; the port
    /// was checked against it exhaustively for every `n` the decoders use and
    /// every `e <= s <= 4096`.
    #[test]
    fn decode_e_matches_reference_classification() {
        // For a character of exactly `n` modules, a pair spanning `e` modules
        // classifies as `e - 2`, valid for 2 <= e <= n - 2. These rows come
        // straight from the C expression; the port was also checked against it
        // exhaustively for every `n` the decoders use and every e <= s <= 4096.
        let rows: [(u32, &[i32]); 3] = [
            // n = 7 (EAN)
            (7, &[-1, -1, 0, 1, 2, 3, -1, -1]),
            // n = 9 (Code 93)
            (9, &[-1, -1, 0, 1, 2, 3, 4, 5, -1, -1]),
            // n = 11 (Code 128)
            (11, &[-1, -1, 0, 1, 2, 3, 4, 5, 6, 7, -1, -1]),
        ];

        for (n, expected) in rows {
            for (e, &want) in expected.iter().enumerate() {
                assert_eq!(decode_e(e as u32, n, n), want, "n={n} e={e}");
            }
        }
    }

    /// Widths below the smallest valid pair make the C expression underflow
    /// through zero; every such case must land outside the valid range rather
    /// than wrapping back into it.
    #[test]
    fn decode_e_rejects_undersized_pairs() {
        for n in [7u32, 9, 11, 14, 45, 72] {
            for e in 0..=1 {
                assert_eq!(decode_e(e, 100, n), -1, "e={e} n={n}");
            }
        }
    }

    /// `ean::aux_end` builds its own reference width by summing four elements
    /// and does not check the result, so a run of zero-width elements reaches
    /// the classifier with nothing to classify against. zbar divides by it
    /// regardless and dies with SIGFPE; there is no width for which that is a
    /// useful answer, so it is rejected like any other unusable input.
    #[test]
    fn decode_e_rejects_a_zero_width_character() {
        for n in [7u32, 9, 11, 14, 45, 72] {
            for e in [0u32, 1, 7, 100, u32::MAX / (2 * n)] {
                assert_eq!(decode_e(e, 0, n), -1, "e={e} n={n}");
            }
        }
    }

    /// A pair can never be wider than the character it sits in, but the
    /// classifier is still expected to reject the boundary cleanly.
    #[test]
    fn decode_e_rejects_oversized_pairs() {
        for n in [7u32, 9, 11, 14, 45, 72] {
            assert_eq!(decode_e(100, 100, n), -1, "n={n}");
        }
    }
}

#[cfg(test)]
mod bitfield_tests {
    use super::*;

    /// The decoders mirror C bitfields, where `character` is a 12-bit signed
    /// field. Every decoder holds -1 in it to mean "no symbol in progress", so
    /// the sign extension on read has to work — a plain mask would return 4095
    /// and the `character < 0` guards would never fire.
    macro_rules! assert_character_roundtrip {
        ($decoder:expr) => {{
            let mut d = $decoder;

            // The idle sentinel every decoder resets to.
            d.set_character(-1);
            assert_eq!(d.character(), -1, "-1 must survive the round trip");

            // Full range of a 12-bit signed field.
            for v in [-2048i16, -2047, -100, -2, 0, 1, 100, 2046, 2047] {
                d.set_character(v);
                assert_eq!(d.character(), v, "character {v}");
            }

            // Writing `character` must not disturb the neighboring fields.
            d.set_element(5);
            d.set_character(-1);
            assert_eq!(d.element(), 5, "element clobbered by set_character");
            d.set_character(1234);
            assert_eq!(d.element(), 5, "element clobbered by set_character");
        }};
    }

    #[cfg(feature = "i25")]
    #[test]
    fn i25_character_field_sign_extends() {
        assert_character_roundtrip!(I25Decoder::default());
    }

    #[cfg(feature = "code39")]
    #[test]
    fn code39_character_field_sign_extends() {
        assert_character_roundtrip!(Code39Decoder::default());
    }

    #[cfg(feature = "code93")]
    #[test]
    fn code93_character_field_sign_extends() {
        assert_character_roundtrip!(Code93Decoder::default());
    }

    #[cfg(feature = "codabar")]
    #[test]
    fn codabar_character_field_sign_extends() {
        assert_character_roundtrip!(CodabarDecoder::default());
    }

    /// Code 128 packs `start` into the same word as the bitfields, so the two
    /// must not tread on each other.
    #[cfg(feature = "code128")]
    #[test]
    fn code128_character_and_start_are_independent() {
        let mut d = Code128Decoder::default();

        d.set_character(-1);
        assert_eq!(d.character(), -1);
        for v in [-2048i16, -1, 0, 1, 2047] {
            d.set_character(v);
            assert_eq!(d.character(), v, "character {v}");
        }

        d.set_start(0xA5);
        d.set_character(-1);
        assert_eq!(d.start(), 0xA5, "start clobbered by set_character");
        assert_eq!(d.character(), -1, "character clobbered by set_start");

        d.set_character(1000);
        d.set_start(0x5A);
        assert_eq!(d.character(), 1000, "character clobbered by set_start");
        assert_eq!(d.start(), 0x5A);
    }

    /// DataBar's `finder` is a 5-bit signed field holding -1 for "unused".
    #[cfg(feature = "databar")]
    #[test]
    fn databar_segment_fields_are_independent() {
        let mut seg = DatabarSegment::default();
        assert_eq!(seg.finder(), -1, "a fresh segment is unused");

        for v in [-16i8, -1, 0, 1, 15] {
            seg.set_finder(v);
            assert_eq!(seg.finder(), v, "finder {v}");
        }

        // Set every field, then check none of them disturbed another.
        seg.set_finder(-1);
        seg.set_exp(true);
        seg.set_color(Color::Bar);
        seg.set_side(1);
        seg.set_partial(true);
        seg.set_count(0x7F);
        seg.set_epoch(0xAB);
        seg.set_check(0xCD);

        assert_eq!(seg.finder(), -1);
        assert!(seg.exp());
        assert_eq!(seg.color(), Color::Bar);
        assert_eq!(seg.side(), 1);
        assert!(seg.partial());
        assert_eq!(seg.count(), 0x7F);
        assert_eq!(seg.epoch(), 0xAB);
        assert_eq!(seg.check(), 0xCD);
    }
}
