//! Error types and handling

use libc::{c_char, c_int, c_uint};
use std::ffi::CStr;
use std::fmt;

pub const ZBAR_VERSION_MAJOR: c_uint = 0;
pub const ZBAR_VERSION_MINOR: c_uint = 23;
pub const ZBAR_VERSION_PATCH: c_uint = 93;

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

pub type Result<T> = std::result::Result<T, Error>;

impl From<i32> for Error {
    fn from(code: i32) -> Self {
        match code {
            1 => Self::OutOfMemory,
            2 => Self::Internal,
            3 => Self::Unsupported,
            4 => Self::Invalid,
            5 => Self::System,
            6 => Self::Locking,
            7 => Self::Busy,
            8 => Self::XDisplay,
            9 => Self::XProto,
            10 => Self::Closed,
            11 => Self::WinApi,
            _ => Self::Unknown(code),
        }
    }
}

/// Set the verbosity level for debug output
///
/// Sets the global verbosity level used by debug macros in error.h.
pub fn zbar_set_verbosity(level: c_int) {
    unsafe {
        _zbar_verbosity = level;
    }
}

// Error severity levels (must match error.h)
const SEV_FATAL: c_int = -2;
const SEV_NOTE: c_int = 2;

// Error module types (must match error.h)
pub const ZBAR_MOD_PROCESSOR: c_int = 0;
const ZBAR_MOD_UNKNOWN: c_int = 4;

// Error type codes (must match zbar.h)
const ZBAR_ERR_NUM: c_int = 12;

// Static error strings
static SEV_STR: [&str; 5] = ["FATAL ERROR", "ERROR", "OK", "WARNING", "NOTE"];
static MOD_STR: [&str; 5] = ["processor", "video", "window", "image scanner", "<unknown>"];
static ERR_STR: [&str; 13] = [
    "no error",
    "out of memory",
    "internal library error",
    "unsupported request",
    "invalid request",
    "system error",
    "locking error",
    "all resources busy",
    "X11 display error",
    "X11 protocol error",
    "output window is closed",
    "windows system error",
    "unknown error",
];

pub struct ErrInfo {
    magic: u32,
    module: c_int,
    errnum: c_int,
    sev: c_int,
    type_: c_int,
    func: *const c_char,
    detail: *const c_char,
    arg_str: *mut c_char,
    arg_int: c_int,
}

/// Print error to stderr
///
/// Formats and prints the error information to stderr.
/// Returns the negative severity value.
pub unsafe fn _zbar_error_spew(err: *mut ErrInfo, verbosity: c_int) -> c_int {
    let err_str = _zbar_error_string(err, verbosity);
    eprint!("{err_str}");

    -(*err).sev
}

/// Format error message as string
///
/// Allocates and formats a detailed error message string.
/// The string is stored in err->buf and returned.
pub unsafe fn _zbar_error_string(err: *mut ErrInfo, _verbosity: c_int) -> String {
    // Get severity string
    let sev = if (*err).sev >= SEV_FATAL && (*err).sev <= SEV_NOTE {
        SEV_STR[((*err).sev + 2) as usize]
    } else {
        SEV_STR[1]
    };

    // Get module string
    let module = if (*err).module >= ZBAR_MOD_PROCESSOR && (*err).module < ZBAR_MOD_UNKNOWN {
        MOD_STR[(*err).module as usize]
    } else {
        MOD_STR[ZBAR_MOD_UNKNOWN as usize]
    };

    // Get function name
    let func = if !(*err).func.is_null() {
        CStr::from_ptr((*err).func).to_str().unwrap_or("<unknown>")
    } else {
        "<unknown>"
    };

    // Get error type string
    let err_type = if (*err).type_ >= 0 && (*err).type_ < ZBAR_ERR_NUM {
        ERR_STR[(*err).type_ as usize]
    } else {
        ERR_STR[ZBAR_ERR_NUM as usize]
    };

    // Build base error message
    let base = format!(
        "{}: zbar {} in {}():\n    {}: ",
        sev, module, func, err_type
    );

    // Add detail if present
    let mut msg = base;
    if !(*err).detail.is_null() {
        let detail = CStr::from_ptr((*err).detail).to_str().unwrap_or("");

        if detail.contains("%s") {
            let arg = if !(*err).arg_str.is_null() {
                CStr::from_ptr((*err).arg_str).to_str().unwrap_or("<?>")
            } else {
                "<?>"
            };
            // Simple replacement - C uses sprintf which is more complex
            msg.push_str(&detail.replace("%s", arg));
        } else if detail.contains("%d") || detail.contains("%x") {
            // Format integer argument
            let formatted = if detail.contains("%x") {
                detail.replace("%x", &format!("{:x}", (*err).arg_int))
            } else {
                detail.replace("%d", &format!("{}", (*err).arg_int))
            };
            msg.push_str(&formatted);
        } else {
            msg.push_str(detail);
        }
    }

    // Add system error if applicable
    if (*err).type_ == 5 {
        // ZBAR_ERR_SYSTEM
        let syserr = libc::strerror((*err).errnum);
        if !syserr.is_null() {
            let syserr_str = CStr::from_ptr(syserr).to_str().unwrap_or("unknown");
            msg.push_str(&format!(": {} ({})\n", syserr_str, (*err).errnum));
        } else {
            msg.push_str(&format!(": ({})\n", (*err).errnum));
        }
    } else {
        msg.push('\n');
    }

    msg
}

// Error handling helper functions
const ERRINFO_MAGIC: u32 = 0x5252457a; // "zERR" (LE)
const ZBAR_ERR_SYSTEM: c_int = 2; // from zbar.h

#[no_mangle]
pub unsafe extern "C" fn _zbar_err_copy(
    dst_c: *mut libc::c_void,
    src_c: *mut libc::c_void,
) -> c_int {
    let dst = dst_c as *mut ErrInfo;
    let src = src_c as *mut ErrInfo;
    debug_assert!((*dst).magic == ERRINFO_MAGIC);
    debug_assert!((*src).magic == ERRINFO_MAGIC);

    (*dst).errnum = (*src).errnum;
    (*dst).sev = (*src).sev;
    (*dst).type_ = (*src).type_;
    (*dst).func = (*src).func;
    (*dst).detail = (*src).detail;
    (*dst).arg_str = (*src).arg_str;
    (*src).arg_str = std::ptr::null_mut(); // unused at src, avoid double free
    (*dst).arg_int = (*src).arg_int;
    -1
}

#[no_mangle]
pub unsafe extern "C" fn _zbar_err_capture(
    err: *mut ErrInfo,
    sev: c_int,
    type_: c_int,
    func: *const c_char,
    detail: *const c_char,
) -> c_int {
    debug_assert!((*err).magic == ERRINFO_MAGIC);
    if type_ == ZBAR_ERR_SYSTEM {
        (*err).errnum = *libc::__errno_location();
    }
    (*err).sev = sev;
    (*err).type_ = type_;
    (*err).func = func;
    (*err).detail = detail;
    if _zbar_verbosity >= 1 {
        _zbar_error_spew(err, 0);
    }
    -1
}

#[no_mangle]
pub unsafe extern "C" fn _zbar_err_capture_str(
    err: *mut ErrInfo,
    sev: c_int,
    type_: c_int,
    func: *const c_char,
    detail: *const c_char,
    arg: *const c_char,
) -> c_int {
    debug_assert!((*err).magic == ERRINFO_MAGIC);
    if !(*err).arg_str.is_null() {
        libc::free((*err).arg_str as *mut libc::c_void);
    }
    (*err).arg_str = libc::strdup(arg);
    _zbar_err_capture(err, sev, type_, func, detail)
}

#[no_mangle]
pub unsafe extern "C" fn _zbar_err_capture_int(
    err: *mut ErrInfo,
    sev: c_int,
    type_: c_int,
    func: *const c_char,
    detail: *const c_char,
    arg: c_int,
) -> c_int {
    debug_assert!((*err).magic == ERRINFO_MAGIC);
    (*err).arg_int = arg;
    _zbar_err_capture(err, sev, type_, func, detail)
}

#[no_mangle]
pub unsafe extern "C" fn _zbar_err_capture_num(
    err: *mut ErrInfo,
    sev: c_int,
    type_: c_int,
    func: *const c_char,
    detail: *const c_char,
    num: c_int,
) -> c_int {
    debug_assert!((*err).magic == ERRINFO_MAGIC);
    (*err).errnum = num;
    _zbar_err_capture(err, sev, type_, func, detail)
}

#[no_mangle]
pub unsafe extern "C" fn _zbar_err_init(err: *mut ErrInfo, module: c_int) {
    (*err).magic = ERRINFO_MAGIC;
    (*err).module = module;
}

#[no_mangle]
pub unsafe extern "C" fn _zbar_err_cleanup(err: *mut ErrInfo) {
    debug_assert!((*err).magic == ERRINFO_MAGIC);
    if !(*err).arg_str.is_null() {
        libc::free((*err).arg_str as *mut libc::c_void);
        (*err).arg_str = std::ptr::null_mut();
    }
}
