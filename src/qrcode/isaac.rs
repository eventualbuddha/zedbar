//! ISAAC - Indirection, Shift, Accumulate, Add, and Count
//!
//! ISAAC is a cryptographically secure pseudo-random number generator designed by
//! Robert J. Jenkins Jr. in 1996. This is one of the most advanced PRNGs available.
//!
//! According to the designer:
//! - No efficient method is known for deducing their internal states
//! - ISAAC requires an amortized 18.75 instructions to produce a 32-bit value
//! - There are no cycles in ISAAC shorter than 2^40 values
//! - The expected cycle length is 2^8295 values
//!
//! Written by Timothy B. Terriberry (tterribe@xiph.org) 1999-2009
//! Public domain implementation based on Robert J. Jenkins Jr.'s work
//! http://www.burtleburtle.net/bob/rand/isaac.html

const ISAAC_SZ_LOG: usize = 8;
const ISAAC_SZ: usize = 1 << ISAAC_SZ_LOG; // 256
const ISAAC_SEED_SZ_MAX: usize = ISAAC_SZ << 2; // 1024
const ISAAC_MASK: u32 = 0xFFFFFFFF;

/// ISAAC random number generator context
#[repr(C)]
pub struct IsaacCtx {
    /// Number of values remaining in r before needing to generate more
    n: u32,
    /// The results buffer
    r: [u32; ISAAC_SZ],
    /// The internal state
    m: [u32; ISAAC_SZ],
    /// Accumulator
    a: u32,
    /// The last result
    b: u32,
    /// Counter (incremented on each update)
    c: u32,
}

impl IsaacCtx {
    /// Create a new ISAAC context with the given seed
    pub fn new(seed: &[u8]) -> Self {
        let mut ctx = IsaacCtx {
            n: 0,
            r: [0; ISAAC_SZ],
            m: [0; ISAAC_SZ],
            a: 0,
            b: 0,
            c: 0,
        };
        ctx.init(seed);
        ctx
    }

    /// Initialize the ISAAC context with a seed
    fn init(&mut self, seed: &[u8]) {
        self.a = 0;
        self.b = 0;
        self.c = 0;

        // Initialize the mixing state with a golden ratio constant
        let mut x = [0x9E3779B9u32; 8];

        // Mix the state 4 times
        for _ in 0..4 {
            Self::mix(&mut x);
        }

        // Process the seed
        let nseed = seed.len().min(ISAAC_SEED_SZ_MAX);

        // Convert seed bytes to u32 values in little-endian format
        let full_words = nseed / 4;
        for i in 0..full_words {
            let base = i * 4;
            self.r[i] =
                u32::from_le_bytes([seed[base], seed[base + 1], seed[base + 2], seed[base + 3]]);
        }

        // Handle remaining bytes
        let remainder = nseed % 4;
        if remainder > 0 {
            let base = full_words * 4;
            let mut bytes = [0u8; 4];
            bytes[..remainder].copy_from_slice(&seed[base..base + remainder]);
            self.r[full_words] = u32::from_le_bytes(bytes);
        }

        // Zero out the rest
        for i in nseed.div_ceil(4)..ISAAC_SZ {
            self.r[i] = 0;
        }

        // Mix in the seed data
        for i in (0..ISAAC_SZ).step_by(8) {
            for (j, value) in x.iter_mut().enumerate() {
                *value = value.wrapping_add(self.r[i + j]);
            }
            Self::mix(&mut x);
            self.m[i..i + 8].copy_from_slice(&x);
        }

        // Mix in the state
        for i in (0..ISAAC_SZ).step_by(8) {
            for (j, value) in x.iter_mut().enumerate() {
                *value = value.wrapping_add(self.m[i + j]);
            }
            Self::mix(&mut x);
            self.m[i..i + 8].copy_from_slice(&x);
        }

        // Generate the first batch of results
        self.update();
    }

    /// Mix function for initialization
    fn mix(x: &mut [u32; 8]) {
        const SHIFTS: [usize; 8] = [11, 2, 8, 16, 10, 4, 8, 9];

        for i in (0..8).step_by(2) {
            x[i] ^= x[(i + 1) & 7] << SHIFTS[i];
            x[(i + 3) & 7] = x[(i + 3) & 7].wrapping_add(x[i]);
            x[(i + 1) & 7] = x[(i + 1) & 7].wrapping_add(x[(i + 2) & 7]);

            let i = i + 1;
            x[i] ^= x[(i + 1) & 7] >> SHIFTS[i];
            x[(i + 3) & 7] = x[(i + 3) & 7].wrapping_add(x[i]);
            x[(i + 1) & 7] = x[(i + 1) & 7].wrapping_add(x[(i + 2) & 7]);
        }
    }

    /// Generate the next batch of random values
    fn update(&mut self) {
        self.c = self.c.wrapping_add(1);
        self.b = self.b.wrapping_add(self.c);

        let mut a = self.a;
        let mut b = self.b;

        // First half
        for i in (0..ISAAC_SZ / 2).step_by(4) {
            let x = self.m[i];
            a = (a ^ (a << 13)).wrapping_add(self.m[i + ISAAC_SZ / 2]);
            self.m[i] = self.m[((x as usize) & ((ISAAC_SZ - 1) << 2)) >> 2]
                .wrapping_add(a)
                .wrapping_add(b);
            b = self.m[(self.m[i] as usize >> (ISAAC_SZ_LOG + 2)) & (ISAAC_SZ - 1)].wrapping_add(x);
            self.r[i] = b;

            let x = self.m[i + 1];
            a = (a ^ (a >> 6)).wrapping_add(self.m[i + ISAAC_SZ / 2 + 1]);
            self.m[i + 1] = self.m[((x as usize) & ((ISAAC_SZ - 1) << 2)) >> 2]
                .wrapping_add(a)
                .wrapping_add(b);
            b = self.m[(self.m[i + 1] as usize >> (ISAAC_SZ_LOG + 2)) & (ISAAC_SZ - 1)]
                .wrapping_add(x);
            self.r[i + 1] = b;

            let x = self.m[i + 2];
            a = (a ^ (a << 2)).wrapping_add(self.m[i + ISAAC_SZ / 2 + 2]);
            self.m[i + 2] = self.m[((x as usize) & ((ISAAC_SZ - 1) << 2)) >> 2]
                .wrapping_add(a)
                .wrapping_add(b);
            b = self.m[(self.m[i + 2] as usize >> (ISAAC_SZ_LOG + 2)) & (ISAAC_SZ - 1)]
                .wrapping_add(x);
            self.r[i + 2] = b;

            let x = self.m[i + 3];
            a = (a ^ (a >> 16)).wrapping_add(self.m[i + ISAAC_SZ / 2 + 3]);
            self.m[i + 3] = self.m[((x as usize) & ((ISAAC_SZ - 1) << 2)) >> 2]
                .wrapping_add(a)
                .wrapping_add(b);
            b = self.m[(self.m[i + 3] as usize >> (ISAAC_SZ_LOG + 2)) & (ISAAC_SZ - 1)]
                .wrapping_add(x);
            self.r[i + 3] = b;
        }

        // Second half
        for i in (ISAAC_SZ / 2..ISAAC_SZ).step_by(4) {
            let x = self.m[i];
            a = (a ^ (a << 13)).wrapping_add(self.m[i - ISAAC_SZ / 2]);
            self.m[i] = self.m[((x as usize) & ((ISAAC_SZ - 1) << 2)) >> 2]
                .wrapping_add(a)
                .wrapping_add(b);
            b = self.m[(self.m[i] as usize >> (ISAAC_SZ_LOG + 2)) & (ISAAC_SZ - 1)].wrapping_add(x);
            self.r[i] = b;

            let x = self.m[i + 1];
            a = (a ^ (a >> 6)).wrapping_add(self.m[i - ISAAC_SZ / 2 + 1]);
            self.m[i + 1] = self.m[((x as usize) & ((ISAAC_SZ - 1) << 2)) >> 2]
                .wrapping_add(a)
                .wrapping_add(b);
            b = self.m[(self.m[i + 1] as usize >> (ISAAC_SZ_LOG + 2)) & (ISAAC_SZ - 1)]
                .wrapping_add(x);
            self.r[i + 1] = b;

            let x = self.m[i + 2];
            a = (a ^ (a << 2)).wrapping_add(self.m[i - ISAAC_SZ / 2 + 2]);
            self.m[i + 2] = self.m[((x as usize) & ((ISAAC_SZ - 1) << 2)) >> 2]
                .wrapping_add(a)
                .wrapping_add(b);
            b = self.m[(self.m[i + 2] as usize >> (ISAAC_SZ_LOG + 2)) & (ISAAC_SZ - 1)]
                .wrapping_add(x);
            self.r[i + 2] = b;

            let x = self.m[i + 3];
            a = (a ^ (a >> 16)).wrapping_add(self.m[i - ISAAC_SZ / 2 + 3]);
            self.m[i + 3] = self.m[((x as usize) & ((ISAAC_SZ - 1) << 2)) >> 2]
                .wrapping_add(a)
                .wrapping_add(b);
            b = self.m[(self.m[i + 3] as usize >> (ISAAC_SZ_LOG + 2)) & (ISAAC_SZ - 1)]
                .wrapping_add(x);
            self.r[i + 3] = b;
        }

        self.a = a;
        self.b = b;
        self.n = ISAAC_SZ as u32;
    }

    /// Get the next random 32-bit unsigned integer
    pub fn next_u32(&mut self) -> u32 {
        if self.n == 0 {
            self.update();
        }
        self.n -= 1;
        self.r[self.n as usize]
    }

    /// Get a uniformly distributed random integer in the range [0, n)
    ///
    /// Returns a random integer uniformly distributed between 0 (inclusive)
    /// and n (exclusive). n must be strictly less than 2^32.
    pub fn next_uint(&mut self, n: u32) -> u32 {
        loop {
            let r = self.next_u32();
            let v = r % n;
            let d = r.wrapping_sub(v);
            // Check if this sample is unbiased
            if d.wrapping_add(n).wrapping_sub(1) & ISAAC_MASK >= d {
                return v;
            }
        }
    }
}

// C FFI exports

/// Initialize an ISAAC context with a seed
///
/// # Safety
///
/// `ctx` must be a valid pointer to an `IsaacCtx` structure
/// `seed` must be a valid pointer to at least `nseed` bytes
#[no_mangle]
pub unsafe extern "C" fn isaac_init(ctx: *mut IsaacCtx, seed: *const u8, nseed: i32) {
    let ctx = &mut *ctx;

    let seed_slice = if seed.is_null() || nseed <= 0 {
        &[]
    } else {
        std::slice::from_raw_parts(seed, nseed as usize)
    };

    ctx.init(seed_slice);
}

/// Get the next random 32-bit unsigned integer
///
/// # Safety
///
/// `ctx` must be a valid pointer to an initialized `IsaacCtx` structure
#[no_mangle]
pub unsafe extern "C" fn isaac_next_uint32(ctx: *mut IsaacCtx) -> u32 {
    (*ctx).next_u32()
}

/// Get a uniformly distributed random integer less than the given maximum
///
/// # Safety
///
/// `ctx` must be a valid pointer to an initialized `IsaacCtx` structure
#[no_mangle]
pub unsafe extern "C" fn isaac_next_uint(ctx: *mut IsaacCtx, n: u32) -> u32 {
    (*ctx).next_uint(n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_isaac_basic() {
        let seed = b"test seed";
        let mut ctx = IsaacCtx::new(seed);

        // Just verify we can generate numbers without panicking
        let _val1 = ctx.next_u32();
        let _val2 = ctx.next_u32();
        let _val3 = ctx.next_u32();
    }

    #[test]
    fn test_isaac_deterministic() {
        let seed = b"deterministic seed";

        let mut ctx1 = IsaacCtx::new(seed);
        let mut ctx2 = IsaacCtx::new(seed);

        // Same seed should produce same sequence
        for _ in 0..100 {
            assert_eq!(ctx1.next_u32(), ctx2.next_u32());
        }
    }

    #[test]
    fn test_isaac_next_uint() {
        let seed = b"bounded test";
        let mut ctx = IsaacCtx::new(seed);

        // Test that bounded values are in range
        for _ in 0..100 {
            let val = ctx.next_uint(10);
            assert!(val < 10);
        }
    }

    #[test]
    fn test_isaac_empty_seed() {
        let mut ctx = IsaacCtx::new(&[]);

        // Should work with empty seed
        let _val = ctx.next_u32();
    }

    #[test]
    fn test_isaac_generates_full_batch() {
        let seed = b"batch test";
        let mut ctx = IsaacCtx::new(seed);

        // Generate more than ISAAC_SZ values to test update
        for _ in 0..ISAAC_SZ * 2 {
            let _val = ctx.next_u32();
        }
    }

    #[test]
    fn test_isaac_c_ffi() {
        unsafe {
            let seed = b"ffi test";
            let mut ctx = IsaacCtx {
                n: 0,
                r: [0; ISAAC_SZ],
                m: [0; ISAAC_SZ],
                a: 0,
                b: 0,
                c: 0,
            };

            isaac_init(&mut ctx, seed.as_ptr(), seed.len() as i32);

            let val1 = isaac_next_uint32(&mut ctx);
            let val2 = isaac_next_uint32(&mut ctx);

            // Values should be different (with extremely high probability)
            assert_ne!(val1, val2);

            // Test bounded generation
            let bounded = isaac_next_uint(&mut ctx, 100);
            assert!(bounded < 100);
        }
    }
}
