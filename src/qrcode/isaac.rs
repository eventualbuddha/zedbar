//! Random number generator for RANSAC
//!
//! This module provides a simple RNG wrapper using the `rand` crate's ChaCha8
//! algorithm, which is a modern, fast, and cryptographically secure PRNG.
//!
//! This replaces the original ISAAC implementation with a more maintainable
//! solution from the Rust ecosystem.

use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;

/// Random number generator context for RANSAC
pub struct IsaacCtx {
    rng: ChaCha8Rng,
}

impl IsaacCtx {
    /// Create a new RNG with the given seed
    pub fn new(seed: &[u8]) -> Self {
        let rng = if seed.is_empty() {
            // Use default seed (all zeros) for deterministic behavior
            ChaCha8Rng::from_seed([0u8; 32])
        } else {
            // Hash the seed to get a 32-byte seed for ChaCha8
            let mut seed_array = [0u8; 32];
            for (i, &byte) in seed.iter().enumerate() {
                seed_array[i % 32] ^= byte;
            }
            ChaCha8Rng::from_seed(seed_array)
        };

        IsaacCtx { rng }
    }

    /// Initialize the RNG with a seed
    fn init(&mut self, seed: &[u8]) {
        *self = Self::new(seed);
    }

    /// Get the next random 32-bit unsigned integer
    pub fn next_u32(&mut self) -> u32 {
        self.rng.gen()
    }

    /// Get a uniformly distributed random integer in the range [0, n)
    ///
    /// Returns a random integer uniformly distributed between 0 (inclusive)
    /// and n (exclusive).
    pub fn next_uint(&mut self, n: u32) -> u32 {
        self.rng.gen_range(0..n)
    }
}

// C FFI exports

/// Initialize an RNG context with a seed
///
/// # Safety
///
/// `ctx` must be a valid pointer to an `IsaacCtx` structure
/// `seed` must be a valid pointer to at least `nseed` bytes
pub unsafe fn isaac_init(ctx: *mut IsaacCtx, seed: *const u8, nseed: i32) {
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
pub unsafe fn isaac_next_uint32(ctx: *mut IsaacCtx) -> u32 {
    let ctx = &mut *ctx;
    ctx.next_u32()
}

/// Get a uniformly distributed random integer less than the given maximum
///
/// # Safety
///
/// `ctx` must be a valid pointer to an initialized `IsaacCtx` structure
pub unsafe fn isaac_next_uint(ctx: *mut IsaacCtx, n: u32) -> u32 {
    let ctx = &mut *ctx;
    ctx.next_uint(n)
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
    fn test_isaac_c_ffi() {
        unsafe {
            let seed = b"ffi test";
            let mut ctx = IsaacCtx::new(&[]);

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
