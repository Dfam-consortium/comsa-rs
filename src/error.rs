use std::fmt;
use std::io;

/// Errors produced while reading, compressing or decompressing alignments.
#[derive(Debug)]
pub enum Error {
    Io(io::Error),
    /// The input alignment is not acceptable (bad symbol, unequal lengths, interleaved format, ...).
    InvalidInput(String),
    /// A compressed archive or block is damaged or not a CoMSA archive.
    Corrupt(String),
    Lzma(liblzma::stream::Error),
}

pub type Result<T> = std::result::Result<T, Error>;

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Io(e) => write!(f, "I/O error: {e}"),
            Error::InvalidInput(s) => write!(f, "invalid input: {s}"),
            Error::Corrupt(s) => write!(f, "corrupt archive: {s}"),
            Error::Lzma(e) => write!(f, "LZMA error: {e}"),
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Error::Io(e) => Some(e),
            Error::Lzma(e) => Some(e),
            _ => None,
        }
    }
}

impl From<io::Error> for Error {
    fn from(e: io::Error) -> Self {
        Error::Io(e)
    }
}

impl From<liblzma::stream::Error> for Error {
    fn from(e: liblzma::stream::Error) -> Self {
        Error::Lzma(e)
    }
}

pub(crate) fn corrupt<T>(msg: impl Into<String>) -> Result<T> {
    Err(Error::Corrupt(msg.into()))
}
