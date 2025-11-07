//! Error types and result handling
//!
//! This module provides error types used throughout the library:
//! - [`Error`] - The main error type
//! - [`Result`] - Type alias for `Result<T, Error>`
//!
//! Most ZBar operations return [`Result`] to indicate success or failure.

use std::fmt;

pub const ZBAR_VERSION_MAJOR: u32 = 0;
pub const ZBAR_VERSION_MINOR: u32 = 23;
pub const ZBAR_VERSION_PATCH: u32 = 93;

/// Global verbosity level
pub static mut ZBAR_VERBOSITY: i32 = 0;

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
            Self::OutOfMemory => write!(f, "out of memory"),
            Self::Internal => write!(f, "internal library error"),
            Self::Unsupported => write!(f, "unsupported request"),
            Self::Invalid => write!(f, "invalid request"),
            Self::System => write!(f, "system error"),
            Self::Locking => write!(f, "locking error"),
            Self::Busy => write!(f, "all resources busy"),
            Self::XDisplay => write!(f, "X11 display error"),
            Self::XProto => write!(f, "X11 protocol error"),
            Self::Closed => write!(f, "output window is closed"),
            Self::WinApi => write!(f, "windows system error"),
            Self::Unknown(code) => write!(f, "unknown error code: {code}"),
        }
    }
}

impl std::error::Error for Error {}

pub type Result<T, E = Error> = std::result::Result<T, E>;
