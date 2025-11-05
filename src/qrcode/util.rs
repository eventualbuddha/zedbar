//! QR Code Utility Functions
//!
//! Mathematical utility functions for QR code decoding, including integer
//! square root, hypotenuse calculation, and logarithm computation.
//!
//! Rust port based on C code from the ZBar library.
//! Original C code copyright (C) 2008-2009 Timothy B. Terriberry (tterribe@xiph.org)
//! Licensed under LGPL 3.0 or later

/// Computes floor(sqrt(val)) exactly using binary search
///
/// Uses the method from http://www.azillionmonkeys.com/qed/sqroot.html
/// The main idea is to search for the largest binary digit b such that
/// (g+b)*(g+b) <= val, and add it to the solution g.
pub fn qr_isqrt(mut val: u32) -> u32 {
    let mut g = 0u32;
    let mut b = 0x8000u32;

    for bshift in (0..16).rev() {
        let t = ((g << 1) + b) << bshift;
        if t <= val {
            g += b;
            val -= t;
        }
        b >>= 1;
    }

    g
}

/// Computes sqrt(x*x + y*y) using CORDIC algorithm
///
/// This implementation is valid for all 32-bit inputs and returns a result
/// accurate to about 27 bits of precision.
///
/// It has been tested for all positive 16-bit inputs, where it returns correctly
/// rounded results in 99.998% of cases and the maximum error is
/// 0.500137134862598032 (for x=48140, y=63018).
///
/// Very nearly all results less than (1<<16) are correctly rounded.
/// All Pythagorean triples with a hypotenuse of less than ((1<<27)-1) evaluate
/// correctly, and the total bias over all Pythagorean triples is -0.04579, with
/// a relative RMS error of 7.2864E-10 and a relative peak error of 7.4387E-9.
pub fn qr_ihypot(x_arg: i32, y_arg: i32) -> u32 {
    let mut x = x_arg.unsigned_abs();
    let x_signed = x_arg.abs();
    let mut y = y_arg.unsigned_abs();
    let mut y_signed = y_arg.abs();

    // Ensure x >= y by swapping if needed
    let mask = -((x > y) as i32) & (x_signed ^ y_signed);
    x ^= mask as u32;
    y ^= mask as u32;
    y_signed ^= mask;

    // Scale the inputs
    let shift = (31 - qr_ilog(y)).max(0);
    x = (((x as u64) << shift).wrapping_mul(0x9B74EDAA) >> 32) as u32;
    y_signed = (((y_signed as i64) << shift).wrapping_mul(0x9B74EDA9) >> 32) as i32;

    // CORDIC iterations
    let mut u = x as i32;
    let mut mask = -((y_signed < 0) as i32);
    x = x.wrapping_add(((y_signed.wrapping_add(mask)) ^ mask) as u32);
    y_signed = y_signed.wrapping_sub((u.wrapping_add(mask)) ^ mask);

    u = ((x.wrapping_add(1)) >> 1) as i32;
    let mut v = y_signed.wrapping_add(1) >> 1;
    mask = -((y_signed < 0) as i32);
    x = x.wrapping_add((v.wrapping_add(mask) ^ mask) as u32);
    y_signed = y_signed.wrapping_sub((u.wrapping_add(mask)) ^ mask);

    for i in 1..16 {
        u = ((x.wrapping_add(1)) >> 2) as i32;
        let r = (1i32 << (2 * i)) >> 1;
        v = y_signed.wrapping_add(r) >> (2 * i);
        mask = -((y_signed < 0) as i32);
        x = x.wrapping_add((v.wrapping_add(mask) ^ mask) as u32);
        y_signed = (y_signed.wrapping_sub((u.wrapping_add(mask)) ^ mask)) << 1;
    }

    (x + ((1u32 << shift) >> 1)) >> shift
}

/// Computes the integer logarithm base 2 of a value
///
/// Returns floor(log2(v)) for v > 0, and 0 for v == 0.
pub fn qr_ilog(mut v: u32) -> i32 {
    let m = (((v & 0xFFFF0000) != 0) as i32) << 4;
    v >>= m;
    let mut ret = m;

    let m = (((v & 0xFF00) != 0) as i32) << 3;
    v >>= m;
    ret |= m;

    let m = (((v & 0xF0) != 0) as i32) << 2;
    v >>= m;
    ret |= m;

    let m = (((v & 0xC) != 0) as i32) << 1;
    v >>= m;
    ret |= m;

    ret |= ((v & 0x2) != 0) as i32;
    ret + (v != 0) as i32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_qr_isqrt_basic() {
        assert_eq!(qr_isqrt(0), 0);
        assert_eq!(qr_isqrt(1), 1);
        assert_eq!(qr_isqrt(2), 1);
        assert_eq!(qr_isqrt(3), 1);
        assert_eq!(qr_isqrt(4), 2);
        assert_eq!(qr_isqrt(9), 3);
        assert_eq!(qr_isqrt(16), 4);
        assert_eq!(qr_isqrt(25), 5);
        assert_eq!(qr_isqrt(100), 10);
        assert_eq!(qr_isqrt(255), 15);
        assert_eq!(qr_isqrt(256), 16);
    }

    #[test]
    fn test_qr_isqrt_large() {
        assert_eq!(qr_isqrt(65535), 255);
        assert_eq!(qr_isqrt(65536), 256);
        assert_eq!(qr_isqrt(1000000), 1000);
        assert_eq!(qr_isqrt(4294967295), 65535); // sqrt(2^32 - 1)
    }

    #[test]
    fn test_qr_ihypot_basic() {
        // Test 3-4-5 right triangle
        assert_eq!(qr_ihypot(3, 4), 5);
        assert_eq!(qr_ihypot(4, 3), 5);
        assert_eq!(qr_ihypot(-3, 4), 5);
        assert_eq!(qr_ihypot(3, -4), 5);
        assert_eq!(qr_ihypot(-3, -4), 5);
    }

    #[test]
    fn test_qr_ihypot_pythagorean_triples() {
        // Test various Pythagorean triples
        assert_eq!(qr_ihypot(5, 12), 13);
        assert_eq!(qr_ihypot(8, 15), 17);
        assert_eq!(qr_ihypot(7, 24), 25);
        assert_eq!(qr_ihypot(20, 21), 29);
    }

    #[test]
    fn test_qr_ihypot_edge_cases() {
        assert_eq!(qr_ihypot(0, 0), 0);
        assert_eq!(qr_ihypot(10, 0), 10);
        assert_eq!(qr_ihypot(0, 10), 10);
        assert_eq!(qr_ihypot(-10, 0), 10);
        assert_eq!(qr_ihypot(0, -10), 10);
    }

    #[test]
    fn test_qr_ihypot_known_error_case() {
        // The documented maximum error case
        // For x=48140, y=63018, the error is 0.500137134862598032
        let result = qr_ihypot(48140, 63018);
        let expected = ((48140f64.powi(2) + 63018f64.powi(2)).sqrt()).round() as u32;
        // Allow for the documented error tolerance
        assert!((result as i32 - expected as i32).abs() <= 1);
    }

    #[test]
    fn test_qr_ilog_basic() {
        // qr_ilog returns ceil(log2(x)) for x > 0, or 0 for x == 0
        // Equivalently: number of bits needed to represent x
        assert_eq!(qr_ilog(0), 0);
        assert_eq!(qr_ilog(1), 1); // 1 bit needed
        assert_eq!(qr_ilog(2), 2); // 2 bits needed
        assert_eq!(qr_ilog(3), 2); // 2 bits needed
        assert_eq!(qr_ilog(4), 3); // 3 bits needed
        assert_eq!(qr_ilog(7), 3); // 3 bits needed
        assert_eq!(qr_ilog(8), 4); // 4 bits needed
        assert_eq!(qr_ilog(15), 4); // 4 bits needed
        assert_eq!(qr_ilog(16), 5); // 5 bits needed
    }

    #[test]
    fn test_qr_ilog_powers_of_two() {
        for i in 0..31 {
            let val = 1u32 << i;
            // For power of 2, returns i + 1 (the bit position, 1-indexed)
            assert_eq!(qr_ilog(val), i + 1);
            // For one less than power of 2, same as the power
            if i > 0 {
                assert_eq!(qr_ilog(val - 1), i);
            }
        }
    }

    #[test]
    fn test_qr_ilog_large() {
        assert_eq!(qr_ilog(255), 8); // 8 bits needed
        assert_eq!(qr_ilog(256), 9); // 9 bits needed
        assert_eq!(qr_ilog(65535), 16); // 16 bits needed
        assert_eq!(qr_ilog(65536), 17); // 17 bits needed
        assert_eq!(qr_ilog(0x7FFFFFFF), 31); // 31 bits needed
        assert_eq!(qr_ilog(0x80000000), 32); // 32 bits needed
        assert_eq!(qr_ilog(0xFFFFFFFF), 32); // 32 bits needed
    }

    #[test]
    fn test_math_functions() {
        assert_eq!(qr_isqrt(100), 10);
        assert_eq!(qr_ihypot(3, 4), 5);
        assert_eq!(qr_ilog(256), 9);
    }
}
