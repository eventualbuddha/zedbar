//! Error types and handling

use libc::{c_int, c_uint};
use std::fmt;

// Version constants - must match zbar/config.h
const ZBAR_VERSION_MAJOR: c_uint = 0;
const ZBAR_VERSION_MINOR: c_uint = 23;
const ZBAR_VERSION_PATCH: c_uint = 93;

/// Global verbosity level
///
/// This is accessed by C macros in error.h for conditional debug output.
/// Must be pub and no_mangle so C code can link to it.
#[no_mangle]
pub static mut _zbar_verbosity: c_int = 0;

#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    OutOfMemory,
    Internal,
    Unsupported,
    Invalid,
    System,
    Locking,
    Busy,
    XDisplay,
    XProto,
    Closed,
    WinApi,
    Unknown(i32),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::OutOfMemory => write!(f, "out of memory"),
            Error::Internal => write!(f, "internal library error"),
            Error::Unsupported => write!(f, "unsupported request"),
            Error::Invalid => write!(f, "invalid request"),
            Error::System => write!(f, "system error"),
            Error::Locking => write!(f, "locking error"),
            Error::Busy => write!(f, "all resources busy"),
            Error::XDisplay => write!(f, "X11 display error"),
            Error::XProto => write!(f, "X11 protocol error"),
            Error::Closed => write!(f, "output window is closed"),
            Error::WinApi => write!(f, "windows system error"),
            Error::Unknown(code) => write!(f, "unknown error code: {code}"),
        }
    }
}

impl std::error::Error for Error {}

pub type Result<T> = std::result::Result<T, Error>;

impl From<i32> for Error {
    fn from(code: i32) -> Self {
        match code {
            1 => Error::OutOfMemory,
            2 => Error::Internal,
            3 => Error::Unsupported,
            4 => Error::Invalid,
            5 => Error::System,
            6 => Error::Locking,
            7 => Error::Busy,
            8 => Error::XDisplay,
            9 => Error::XProto,
            10 => Error::Closed,
            11 => Error::WinApi,
            _ => Error::Unknown(code),
        }
    }
}

/// Get the version information
///
/// Fills in the major, minor, and patch version numbers if the pointers are non-null.
/// Returns 0 on success.
#[no_mangle]
pub unsafe extern "C" fn zbar_version(
    major: *mut c_uint,
    minor: *mut c_uint,
    patch: *mut c_uint,
) -> libc::c_int {
    if !major.is_null() {
        *major = ZBAR_VERSION_MAJOR;
    }
    if !minor.is_null() {
        *minor = ZBAR_VERSION_MINOR;
    }
    if !patch.is_null() {
        *patch = ZBAR_VERSION_PATCH;
    }
    0
}

/// Set the verbosity level for debug output
///
/// Sets the global verbosity level used by debug macros in error.h.
#[no_mangle]
pub extern "C" fn zbar_set_verbosity(level: c_int) {
    unsafe {
        _zbar_verbosity = level;
    }
}

/// Increase the verbosity level
///
/// Increments verbosity from 0 to 1, then doubles it on each subsequent call
/// (creating the pattern: 0 -> 1 -> 2 -> 4 -> 8 -> 16...).
#[no_mangle]
pub extern "C" fn zbar_increase_verbosity() {
    unsafe {
        if _zbar_verbosity == 0 {
            _zbar_verbosity = 1;
        } else {
            _zbar_verbosity <<= 1;
        }
    }
}
