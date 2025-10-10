//! DataBar utility functions
//!
//! This module provides utility functions for DataBar barcode decoding.

use libc::c_int;

/// Check if two widths are compatible within tolerance
///
/// Validates that two measured widths (wf and wd) are within an acceptable
/// tolerance range for n modules. Used to match DataBar segments.
///
/// # Parameters
/// - `wf`: First width measurement
/// - `wd`: Second width measurement
/// - `n`: Number of modules
///
/// # Returns
/// 1 if widths match within tolerance, 0 otherwise
#[no_mangle]
pub extern "C" fn _zbar_databar_check_width(wf: u32, wd: u32, n: u32) -> c_int {
    let dwf = wf * 3;
    let wd = wd * 14;
    let wf = wf * n;

    // Check: wf - dwf <= wd && wd <= wf + dwf
    // In C, this relies on unsigned wraparound if wf < dwf
    // For unsigned subtraction: wf - dwf will wrap if wf < dwf,
    // resulting in a very large number, making the condition false
    if wf.wrapping_sub(dwf) <= wd && wd <= wf + dwf {
        1
    } else {
        0
    }
}

/// Decode a number into decimal digits
///
/// Converts an unsigned long value into ASCII decimal digits and stores
/// them in the provided buffer.
///
/// # Parameters
/// - `buf`: Buffer to write digits to
/// - `n`: Number to decode
/// - `i`: Number of digits to write
///
/// # Safety
/// Buffer must have at least `i` bytes available
#[no_mangle]
pub unsafe extern "C" fn _zbar_databar_decode10(buf: *mut u8, mut n: u64, i: c_int) {
    if buf.is_null() || i <= 0 {
        return;
    }

    let mut pos = buf.add(i as usize);
    let mut remaining = i;

    while remaining > 0 {
        let d = (n % 10) as u8;
        n /= 10;
        pos = pos.sub(1);
        *pos = b'0' + d;
        remaining -= 1;
    }
}

/// Append checksum digit to 13-digit buffer
///
/// Calculates and appends a check digit to a 13-digit numeric string
/// using weighted sum modulo 10 algorithm.
///
/// # Parameters
/// - `buf`: Buffer containing 13 ASCII digits, with space for 14th digit
///
/// # Safety
/// Buffer must contain at least 14 bytes, with first 13 being ASCII digits '0'-'9'
#[no_mangle]
pub unsafe extern "C" fn _zbar_databar_append_check14(buf: *mut u8) {
    if buf.is_null() {
        return;
    }

    let mut chk: u8 = 0;
    let mut ptr = buf;

    for i in (0..13).rev() {
        let d = *ptr - b'0';
        chk = chk.wrapping_add(d);
        if (i & 1) == 0 {
            chk = chk.wrapping_add(d << 1);
        }
        ptr = ptr.add(1);
    }

    chk %= 10;
    if chk != 0 {
        chk = 10 - chk;
    }
    *ptr = chk + b'0';
}

/// Calculate DataBar checksum
///
/// Computes a checksum value for DataBar symbols based on signature values.
///
/// # Parameters
/// - `sig0`: First signature value (4 nibbles)
/// - `sig1`: Second signature value (4 nibbles)
/// - `side`: Side indicator (0 or 1)
/// - `mod_val`: Modulus value for checksum calculation
///
/// # Returns
/// Calculated checksum value
#[no_mangle]
pub extern "C" fn _zbar_databar_calc_check(
    mut sig0: u32,
    mut sig1: u32,
    side: u32,
    mod_val: u32,
) -> u32 {
    let mut chk: u32 = 0;

    for i in (0..4).rev() {
        chk = (chk * 3 + (sig1 & 0xf) + 1) * 3 + (sig0 & 0xf) + 1;
        sig1 >>= 4;
        sig0 >>= 4;
        if (i & 1) == 0 {
            chk %= mod_val;
        }
    }

    if side != 0 {
        chk = (chk * (6561 % mod_val)) % mod_val;
    }

    chk
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_check_width() {
        // Test exact match: wf=100, wd=100, n=14
        // After scaling: dwf=300, wd=1400, wf=1400
        // Check: (1400-300 <= 1400) && (1400 <= 1400+300) => (1100<=1400) && (1400<=1700) => true
        assert_eq!(_zbar_databar_check_width(100, 100, 14), 1);

        // Test within tolerance: wf=100, wd=105, n=14
        // After scaling: dwf=300, wd=1470, wf=1400
        // Check: (1400-300 <= 1470) && (1470 <= 1400+300) => (1100<=1470) && (1470<=1700) => true
        assert_eq!(_zbar_databar_check_width(100, 105, 14), 1);

        // Test within tolerance: wf=100, wd=95, n=14
        // After scaling: dwf=300, wd=1330, wf=1400
        // Check: (1400-300 <= 1330) && (1330 <= 1400+300) => (1100<=1330) && (1330<=1700) => true
        assert_eq!(_zbar_databar_check_width(100, 95, 14), 1);

        // Test outside tolerance: wf=100, wd=130, n=14
        // After scaling: dwf=300, wd=1820, wf=1400
        // Check: (1400-300 <= 1820) && (1820 <= 1400+300) => (1100<=1820) && (1820<=1700) => false
        assert_eq!(_zbar_databar_check_width(100, 130, 14), 0);

        // Test outside tolerance: wf=100, wd=70, n=14
        // After scaling: dwf=300, wd=980, wf=1400
        // Check: (1400-300 <= 980) && (980 <= 1400+300) => (1100<=980) && (980<=1700) => false
        assert_eq!(_zbar_databar_check_width(100, 70, 14), 0);
    }

    #[test]
    fn test_decode10() {
        unsafe {
            let mut buf = [0u8; 10];

            // Test simple number
            _zbar_databar_decode10(buf.as_mut_ptr(), 123, 3);
            assert_eq!(&buf[0..3], b"123");

            // Test with leading zeros
            _zbar_databar_decode10(buf.as_mut_ptr(), 45, 6);
            assert_eq!(&buf[0..6], b"000045");

            // Test zero
            _zbar_databar_decode10(buf.as_mut_ptr(), 0, 3);
            assert_eq!(&buf[0..3], b"000");
        }
    }

    #[test]
    fn test_append_check14() {
        unsafe {
            // Test with 13 ASCII digits
            let mut buf = [
                b'9', b'7', b'8', b'0', b'1', b'4', b'3', b'0', b'0', b'7', b'2', b'3', b'0', 0,
            ];
            _zbar_databar_append_check14(buf.as_mut_ptr());
            // Check digit: (9+7+8+0+1+4+3+0+0+7+2+3) + (9+8+1+3+0+2)*2 = 44 + 46 = 90, 90%10=0, check=0
            assert_eq!(buf[13], b'0');

            // Test another example: 1234567890120
            let mut buf2 = [
                b'1', b'2', b'3', b'4', b'5', b'6', b'7', b'8', b'9', b'0', b'1', b'2', b'0', 0,
            ];
            _zbar_databar_append_check14(buf2.as_mut_ptr());
            // Check: (1+2+3+4+5+6+7+8+9+0+1+2) + (1+3+5+7+9+1)*2 = 48 + 52 = 100, 100%10=0, check=0
            assert_eq!(buf2[13], b'0');
        }
    }

    #[test]
    fn test_calc_check() {
        // Test basic checksum calculation
        // These values are based on understanding the algorithm
        let chk1 = _zbar_databar_calc_check(0x1234, 0x5678, 0, 211);
        assert!(chk1 < 211);

        let chk2 = _zbar_databar_calc_check(0x1234, 0x5678, 1, 211);
        assert!(chk2 < 211);

        // Side should affect the result
        assert_ne!(chk1, chk2);

        // Test with different modulus
        let chk3 = _zbar_databar_calc_check(0x1234, 0x5678, 0, 79);
        assert!(chk3 < 79);

        // Same inputs should give same output
        let chk4 = _zbar_databar_calc_check(0x1234, 0x5678, 0, 211);
        assert_eq!(chk1, chk4);
    }
}
