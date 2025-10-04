//! Reed-Solomon Error Correction
//!
//! Implementation of Reed-Solomon encoder and decoder for QR codes.
//!
//! Original implementation (C) Henry Minsky (hqm@ua.com, hqm@ai.mit.edu),
//! Universal Access 1991-1995.
//!
//! Updates by Timothy B. Terriberry (C) 2008-2009:
//! - Properly reject codes when error-locator polynomial has repeated roots or
//!   non-trivial irreducible factors.
//! - Removed the hard-coded parity size and performed general API cleanup.
//! - Allow multiple representations of GF(2**8), since different standards use
//!   different irreducible polynomials.
//! - Allow different starting indices for the generator polynomial, since
//!   different standards use different values.
//! - Greatly reduced the computation by eliminating unnecessary operations.
//! - Explicitly solve for the roots of low-degree polynomials instead of using
//!   an exhaustive search.
//!   This is another major speed boost when there are few errors.

use libc::{c_int, c_uint};

/// Galois Field GF(2**8) for Reed-Solomon operations
#[repr(C)]
pub struct rs_gf256 {
    /// Logarithm table in GF(2**8)
    pub log: [u8; 256],
    /// Exponential table in GF(2**8): exp[i] contains x^i reduced modulo the
    /// irreducible primitive polynomial used to define the field.
    /// The extra 256 entries are used to do arithmetic mod 255, since some extra
    /// table lookups are generally faster than doing the modulus.
    pub exp: [u8; 511],
}

// ============================================================================
// Galois Field Arithmetic
// ============================================================================

/// Initialize discrete logarithm tables for GF(2**8) using a given primitive
/// irreducible polynomial.
#[no_mangle]
pub unsafe extern "C" fn rs_gf256_init(gf: *mut rs_gf256, ppoly: c_uint) {
    let gf = &mut *gf;
    let mut p: c_uint = 1;

    // Initialize the table of powers of a primitive root, alpha=0x02
    for i in 0..256 {
        gf.exp[i] = p as u8;
        gf.exp[i + 255] = p as u8;
        let mask = if (p >> 7) != 0 { ppoly } else { 0 };
        p = ((p << 1) ^ mask) & 0xFF;
    }

    // Invert the table to recover the logs
    for i in 0..255 {
        gf.log[gf.exp[i] as usize] = i as u8;
    }
    // Note that we rely on the fact that gf.log[0]=0 below
    gf.log[0] = 0;
}

/// Multiplication in GF(2**8) using logarithms
#[inline]
fn rs_gmul(gf: &rs_gf256, a: u8, b: u8) -> u8 {
    if a == 0 || b == 0 {
        0
    } else {
        gf.exp[(gf.log[a as usize] as usize) + (gf.log[b as usize] as usize)]
    }
}

/// Division in GF(2**8) using logarithms
/// The result of division by zero is undefined
#[inline]
fn rs_gdiv(gf: &rs_gf256, a: u8, b: u8) -> u8 {
    if a == 0 {
        0
    } else {
        gf.exp[(gf.log[a as usize] as usize) + 255 - (gf.log[b as usize] as usize)]
    }
}

/// Multiplication in GF(2**8) when one of the numbers is known to be non-zero
/// (proven by representing it by its logarithm)
#[inline]
fn rs_hgmul(gf: &rs_gf256, a: u8, logb: usize) -> u8 {
    if a == 0 {
        0
    } else {
        gf.exp[(gf.log[a as usize] as usize) + logb]
    }
}

/// Square root in GF(2**8) using logarithms
#[inline]
fn rs_gsqrt(gf: &rs_gf256, a: u8) -> u8 {
    if a == 0 {
        0
    } else {
        let loga = gf.log[a as usize] as usize;
        gf.exp[(loga + (255 & (-(loga as i32 & 1) as usize))) >> 1]
    }
}

// ============================================================================
// Polynomial Root Finding
// ============================================================================

/// Solve a quadratic equation x**2 + b*x + c in GF(2**8)
/// Returns the number of distinct roots
fn rs_quadratic_solve(gf: &rs_gf256, b: u8, c: u8, x: &mut [u8]) -> c_int {
    // If b is zero, all we need is a square root
    if b == 0 {
        x[0] = rs_gsqrt(gf, c);
        return 1;
    }

    // If c is zero, 0 and b are the roots
    if c == 0 {
        x[0] = 0;
        x[1] = b;
        return 2;
    }

    let mut logb = gf.log[b as usize] as usize;
    let mut logc = gf.log[c as usize] as usize;

    // If b lies in GF(2**4), scale x to move it out
    let inc = if logb % (255 / 15) == 0 { 1 } else { 0 };
    let b = if inc != 0 {
        let b = gf.exp[logb + 254];
        logb = gf.log[b as usize] as usize;
        let c_val = gf.exp[logc + 253];
        logc = gf.log[c_val as usize] as usize;
        b
    } else {
        b
    };

    let logb2 = gf.log[gf.exp[logb << 1] as usize] as usize;
    let logb4 = gf.log[gf.exp[logb2 << 1] as usize] as usize;
    let logb8 = gf.log[gf.exp[logb4 << 1] as usize] as usize;
    let logb12 = gf.log[gf.exp[logb4 + logb8] as usize] as usize;
    let logb14 = gf.log[gf.exp[logb2 + logb12] as usize] as usize;
    let logc2 = gf.log[gf.exp[logc << 1] as usize] as usize;
    let logc4 = gf.log[gf.exp[logc2 << 1] as usize] as usize;
    let c8 = gf.exp[logc4 << 1];

    let g3 = rs_hgmul(gf,
        gf.exp[logb14 + logc] ^ gf.exp[logb12 + logc2] ^ gf.exp[logb8 + logc4] ^ c8,
        logb);

    // If g3 doesn't lie in GF(2**4), then our roots lie in an extension field
    if gf.log[g3 as usize] as usize % (255 / 15) != 0 {
        return 0;
    }

    let z3 = rs_gdiv(gf, g3, gf.exp[logb8 << 1] ^ b);
    let l3 = rs_hgmul(gf,
        rs_gmul(gf, z3, z3) ^ rs_hgmul(gf, z3, logb) ^ c,
        255 - logb2);
    let c0 = rs_hgmul(gf, l3, 255 - 2 * (255 / 15));

    let g2 = rs_hgmul(gf,
        rs_hgmul(gf, c0, 255 - 2 * (255 / 15)) ^ rs_gmul(gf, c0, c0),
        255 - 255 / 15);
    let z2 = rs_gdiv(gf, g2, gf.exp[255 - (255 / 15) * 4] ^ gf.exp[255 - (255 / 15)]);
    let l2 = rs_hgmul(gf,
        rs_gmul(gf, z2, z2) ^ rs_hgmul(gf, z2, 255 - (255 / 15)) ^ c0,
        2 * (255 / 15));

    x[0] = gf.exp[gf.log[
        (z3 ^ rs_hgmul(gf,
            rs_hgmul(gf, l2, 255 / 3) ^ rs_hgmul(gf, z2, 255 / 15),
            logb))as usize] as usize + inc];
    x[1] = x[0] ^ b;
    2
}

/// Solve a cubic equation x**3 + a*x**2 + b*x + c in GF(2**8)
/// Returns the number of distinct roots
fn rs_cubic_solve(gf: &rs_gf256, a: u8, b: u8, c: u8, x: &mut [u8]) -> c_int {
    // If c is zero, factor out the 0 root
    if c == 0 {
        let mut nroots = rs_quadratic_solve(gf, a, b, x);
        if b != 0 {
            x[nroots as usize] = 0;
            nroots += 1;
        }
        return nroots;
    }

    let k = rs_gmul(gf, a, b) ^ c;
    let d2 = rs_gmul(gf, a, a) ^ b;

    if d2 == 0 {
        if k == 0 {
            // We have a triple root
            x[0] = a;
            return 1;
        }
        let mut logx = gf.log[k as usize] as c_int;
        if logx % 3 != 0 {
            return 0;
        }
        logx /= 3;
        x[0] = a ^ gf.exp[logx as usize];
        x[1] = a ^ gf.exp[(logx + 255 / 3) as usize];
        x[2] = a ^ x[0] ^ x[1];
        return 3;
    }

    let logd2 = gf.log[d2 as usize] as usize;
    let logd = (logd2 + (255 & (-(logd2 as i32 & 1) as usize))) >> 1;
    let k = rs_gdiv(gf, k, gf.exp[logd + logd2]);

    let nroots = rs_quadratic_solve(gf, k, 1, x);
    if nroots < 1 {
        return 0;
    }

    let logw = gf.log[x[0] as usize] as usize;
    if logw != 0 {
        if logw % 3 != 0 {
            return 0;
        }
        let logw = logw / 3;
        x[0] = gf.exp[gf.log[(gf.exp[logw] ^ gf.exp[255 - logw]) as usize] as usize + logd] ^ a;
        let logw = logw + 255 / 3;
        x[1] = gf.exp[gf.log[(gf.exp[logw] ^ gf.exp[255 - logw]) as usize] as usize + logd] ^ a;
        x[2] = x[0] ^ x[1] ^ a;
        return 3;
    } else {
        x[0] = a;
        return 1;
    }
}

/// Solve a quartic equation x**4 + a*x**3 + b*x**2 + c*x + d in GF(2**8)
/// Returns the number of distinct roots
fn rs_quartic_solve(gf: &rs_gf256, a: u8, b: u8, c: u8, d: u8, x: &mut [u8]) -> c_int {
    // If d is zero, factor out the 0 root
    if d == 0 {
        let mut nroots = rs_cubic_solve(gf, a, b, c, x);
        if c != 0 {
            x[nroots as usize] = 0;
            nroots += 1;
        }
        return nroots;
    }

    if a != 0 {
        let loga = gf.log[a as usize] as usize;
        let r = rs_hgmul(gf, c, 255 - loga);
        let s = rs_gsqrt(gf, r);
        let t = d ^ rs_gmul(gf, b, r) ^ rs_gmul(gf, r, r);

        if t != 0 {
            let logti = 255 - gf.log[t as usize] as usize;
            let nroots = rs_quartic_solve(gf, 0,
                rs_hgmul(gf, b ^ rs_hgmul(gf, s, loga), logti),
                gf.exp[loga + logti],
                gf.exp[logti],
                x);
            for i in 0..nroots as usize {
                x[i] = gf.exp[255 - gf.log[x[i] as usize] as usize] ^ s;
            }
            return nroots;
        } else {
            let mut nroots = rs_quadratic_solve(gf, a, b ^ r, x);
            if nroots != 2 || (x[0] != s && x[1] != s) {
                x[nroots as usize] = s;
                nroots += 1;
            }
            return nroots;
        }
    }

    // If there are no odd powers, it's really just a quadratic in disguise
    if c == 0 {
        return rs_quadratic_solve(gf, rs_gsqrt(gf, b), rs_gsqrt(gf, d), x);
    }

    let nroots = rs_cubic_solve(gf, 0, b, c, x);
    if nroots < 1 {
        return 0;
    }

    let r = x[0];
    let b = rs_gdiv(gf, c, r);
    let nroots = rs_quadratic_solve(gf, b, d, x);
    if nroots < 2 {
        return 0;
    }

    let s = x[0];
    let t = x[1];
    let nroots = rs_quadratic_solve(gf, r, s, x);
    nroots + rs_quadratic_solve(gf, r, t, &mut x[nroots as usize..])
}

// ============================================================================
// Polynomial Arithmetic
// ============================================================================

#[inline]
fn rs_poly_zero(p: &mut [u8]) {
    for val in p.iter_mut() {
        *val = 0;
    }
}

#[inline]
fn rs_poly_copy(p: &mut [u8], q: &[u8]) {
    p[..q.len()].copy_from_slice(q);
}

/// Multiply the polynomial by the free variable, x (shift the coefficients)
fn rs_poly_mul_x(p: &mut [u8], q: &[u8], dp1: usize) {
    p[1..dp1].copy_from_slice(&q[0..dp1-1]);
    p[0] = 0;
}

/// Compute the first (d+1) coefficients of the product of a degree e and a
/// degree f polynomial
fn rs_poly_mult(gf: &rs_gf256, p: &mut [u8], dp1: usize,
                q: &[u8], ep1: usize,
                r: &[u8], fp1: usize) {
    rs_poly_zero(&mut p[..dp1]);
    let m = if ep1 < dp1 { ep1 } else { dp1 };

    for i in 0..m {
        if q[i] != 0 {
            let logqi = gf.log[q[i] as usize] as usize;
            let n = if dp1 - i < fp1 { dp1 - i } else { fp1 };
            for j in 0..n {
                p[i + j] ^= rs_hgmul(gf, r[j], logqi);
            }
        }
    }
}

// ============================================================================
// Decoding
// ============================================================================

/// Computes the syndrome of a codeword
fn rs_calc_syndrome(gf: &rs_gf256, m0: c_int, s: &mut [u8], npar: c_int,
                    data: &[u8], ndata: c_int) {
    for j in 0..npar as usize {
        let mut sj: u8 = 0;
        let alphaj = gf.log[gf.exp[j + m0 as usize] as usize] as usize;
        for i in 0..ndata as usize {
            sj = data[i] ^ rs_hgmul(gf, sj, alphaj);
        }
        s[j] = sj;
    }
}

/// Initialize lambda to the product of (1-x*alpha**e[i]) for erasure locations e[i]
fn rs_init_lambda(gf: &rs_gf256, lambda: &mut [u8], npar: c_int,
                  erasures: &[u8], nerasures: c_int, ndata: c_int) {
    let size = if npar < 4 { 4 } else { npar } as usize + 1;
    rs_poly_zero(&mut lambda[..size]);
    lambda[0] = 1;

    for i in 0..nerasures as usize {
        for j in (1..=i+1).rev() {
            lambda[j] ^= rs_hgmul(gf, lambda[j - 1],
                (ndata - 1 - erasures[i] as c_int) as usize);
        }
    }
}

/// Modified Berlekamp-Massey algorithm
/// Returns the number of errors detected (degree of lambda)
fn rs_modified_berlekamp_massey(gf: &rs_gf256, lambda: &mut [u8],
                                s: &[u8], omega: &mut [u8], npar: c_int,
                                erasures: &[u8], nerasures: c_int,
                                ndata: c_int) -> c_int {
    let mut tt = [0u8; 256];

    rs_init_lambda(gf, lambda, npar, erasures, nerasures, ndata);
    rs_poly_copy(&mut tt[..npar as usize + 1], &lambda[..npar as usize + 1]);

    let mut l = nerasures;
    let mut k = 0;

    for n in (nerasures + 1)..=npar {
        let tt_copy = tt.clone();
        rs_poly_mul_x(&mut tt, &tt_copy, (n - k + 1) as usize);
        let mut d: u8 = 0;
        for i in 0..=l as usize {
            d ^= rs_gmul(gf, lambda[i], s[(n - 1 - i as c_int) as usize]);
        }

        if d != 0 {
            let logd = gf.log[d as usize] as usize;
            if l < n - k {
                for i in 0..=(n - k) as usize {
                    let tti = tt[i];
                    tt[i] = rs_hgmul(gf, lambda[i], 255 - logd);
                    lambda[i] = lambda[i] ^ rs_hgmul(gf, tti, logd);
                }
                let t = n - k;
                k = n - l;
                l = t;
            } else {
                for i in 0..=l as usize {
                    lambda[i] = lambda[i] ^ rs_hgmul(gf, tt[i], logd);
                }
            }
        }
    }

    rs_poly_mult(gf, omega, npar as usize, lambda, (l + 1) as usize, s, npar as usize);
    l
}

/// Finds all the roots of an error-locator polynomial lambda
/// Returns the number of valid roots identified
fn rs_find_roots(gf: &rs_gf256, epos: &mut [u8], lambda: &[u8],
                 nerrors: c_int, ndata: c_int) -> c_int {
    let mut nroots = 0;

    if nerrors <= 4 {
        let nerrors_found = rs_quartic_solve(gf, lambda[1], lambda[2],
            lambda[3], lambda[4], epos);
        for i in 0..nerrors_found as usize {
            if epos[i] != 0 {
                let alpha = gf.log[epos[i] as usize] as c_int;
                if alpha < ndata {
                    epos[nroots as usize] = alpha as u8;
                    nroots += 1;
                }
            }
        }
        return nroots;
    } else {
        for alpha in 0..ndata as usize {
            let mut sum: u8 = 0;
            let mut alphai: usize = 0;
            for i in 0..=nerrors as usize {
                sum ^= rs_hgmul(gf, lambda[nerrors as usize - i], alphai);
                alphai = gf.log[gf.exp[alphai + alpha] as usize] as usize;
            }
            if sum == 0 {
                epos[nroots as usize] = alpha as u8;
                nroots += 1;
            }
        }
    }

    nroots
}

/// Corrects a codeword with ndata<256 bytes, of which the last npar are parity bytes
///
/// Known locations of errors can be passed in the erasures array.
/// Twice as many (up to npar) errors with a known location can be corrected
/// compared to errors with an unknown location.
///
/// Returns the number of errors corrected if successful, or a negative number if
/// the message could not be corrected because too many errors were detected.
#[no_mangle]
pub unsafe extern "C" fn rs_correct(
    gf: *const rs_gf256,
    m0: c_int,
    data: *mut u8,
    ndata: c_int,
    npar: c_int,
    erasures: *const u8,
    nerasures: c_int
) -> c_int {
    let gf = &*gf;
    let data_slice = std::slice::from_raw_parts_mut(data, ndata as usize);
    let erasures_slice = if erasures.is_null() || nerasures == 0 {
        &[]
    } else {
        std::slice::from_raw_parts(erasures, nerasures as usize)
    };

    let mut lambda = [0u8; 256];
    let mut omega = [0u8; 256];
    let mut epos = [0u8; 256];
    let mut s = [0u8; 256];

    // If we already have too many erasures, we can't possibly succeed
    if nerasures > npar {
        return -1;
    }

    // Compute the syndrome values
    rs_calc_syndrome(gf, m0, &mut s, npar, data_slice, ndata);

    // Check for a non-zero value
    for i in 0..npar as usize {
        if s[i] != 0 {
            // Construct the error locator polynomial
            let nerrors = rs_modified_berlekamp_massey(gf, &mut lambda, &s, &mut omega,
                npar, erasures_slice, nerasures, ndata);

            // If we can't locate any errors, or have too many errors, fail
            if nerrors <= 0 || nerrors - nerasures > ((npar - nerasures) >> 1) {
                return -1;
            }

            // Compute the locations of the errors
            if rs_find_roots(gf, &mut epos, &lambda, nerrors, ndata) < nerrors {
                return -1;
            }

            // Now compute the error magnitudes
            for i in 0..nerrors as usize {
                let alpha = epos[i] as usize;
                let alphan1 = 255 - alpha;

                // Evaluate omega at alpha**-1
                let mut a: u8 = 0;
                let mut alphanj = 0;
                for j in 0..npar as usize {
                    a ^= rs_hgmul(gf, omega[j], alphanj);
                    alphanj = gf.log[gf.exp[alphanj + alphan1] as usize] as usize;
                }

                // Evaluate the derivative of lambda at alpha**-1
                let mut b: u8 = 0;
                let alphan2 = gf.log[gf.exp[alphan1 << 1] as usize] as usize;
                let mut alphanj = (alphan1 + m0 as usize * alpha % 255) as usize;
                for j in (1..=npar as usize).step_by(2) {
                    b ^= rs_hgmul(gf, lambda[j], alphanj);
                    alphanj = gf.log[gf.exp[alphanj + alphan2] as usize] as usize;
                }

                // Apply the correction
                data_slice[ndata as usize - 1 - alpha] ^= rs_gdiv(gf, a, b);
            }

            return nerrors;
        }
    }

    0
}
