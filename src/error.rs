//! Error types and handling

use std::fmt;

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

