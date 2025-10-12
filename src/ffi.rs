//! FFI bindings to the C ZBar library
//!
//! This module provides direct bindings to the original C implementation.
//! As modules are converted to Rust, these bindings will be gradually removed.

use libc::c_int;

// Reference counting helper
pub(crate) unsafe fn refcnt(cnt: *mut c_int, delta: c_int) -> c_int {
    let rc = *cnt + delta;
    *cnt = rc;
    debug_assert!(rc >= 0);
    rc
}
