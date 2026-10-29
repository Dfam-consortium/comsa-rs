//! Compression of multiple sequence alignments (MSA) in Stockholm and aligned FASTA format.
//!
//! This is a Rust port of CoMSA 1.1 (Deorowicz et al., "CoMSA: compression of protein
//! multiple sequence alignment files", Bioinformatics 2019) and is format compatible with
//! it: archives (`.smsac` for Stockholm, `.msac` for FASTA) written by either
//! implementation can be read by the other.
//!
//! The alignment is coded column-wise: generalized positional BWT (PBWT), then weighted
//! frequency count (WFC) or move-to-front (MTF, "fast" variant), then zero-run-length
//! coding and a context-modelled range coder. Sequence names and Stockholm annotation
//! (`#=GF`, `#=GS`, `#=GR`, `#=GC` lines) are compressed with LZMA.
//!
//! # Example
//! ```no_run
//! use std::fs::File;
//! use std::io::{BufReader, BufWriter};
//! use comsa::{CompressOptions, archive, stockholm};
//!
//! # fn main() -> comsa::Result<()> {
//! // Compress a Stockholm file into an archive.
//! let reader = stockholm::Reader::new(BufReader::new(File::open("pfam.sto")?));
//! let out = BufWriter::new(File::create("pfam.smsac")?);
//! archive::compress_stockholm(reader, out, &CompressOptions::default(), |_| {})?;
//!
//! // Random access to a single family.
//! let mut arc = archive::ArchiveReader::open(File::open("pfam.smsac")?)?;
//! if let Some(i) = arc.families().iter().position(|d| d.ac == "PF00005.26") {
//!     let family = arc.read_family(i, 1)?;
//!     let mut text = Vec::new();
//!     stockholm::write_family(&mut text, &family)?;
//! }
//! # Ok(())
//! # }
//! ```

pub mod archive;
pub mod codec;
mod entropy;
mod error;
pub mod fasta;
mod lzma;
mod mtf;
mod parallel;
mod pbwt;
mod rangecoder;
mod rle;
pub mod stockholm;
mod varint;
mod wfc;

pub use codec::CompressOptions;
pub use error::{Error, Result};
pub use stockholm::Family;

/// Symbols accepted in aligned sequence data: gaps (`-`, `.`), residues (`A`-`Z`, `a`-`z`),
/// stop (`*`) and Easel's missing-data symbol (`~`).
pub const fn is_valid_sequence_symbol(c: u8) -> bool {
    matches!(c, b'-' | b'.' | b'A'..=b'Z' | b'a'..=b'z' | b'*' | b'~')
}

/// Initial symbol list of the WFC/MTF stage: valid symbols first (so their ranks stay
/// below 64), then all remaining 7-bit values. Identical to CoMSA's list with `~` added
/// after `*`, which leaves every other valid symbol's rank unchanged.
pub(crate) fn initial_symbol_order() -> [u8; 128] {
    let mut order = [0u8; 128];
    let mut n = 0;
    let first = b"-.ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz*~";
    for &c in first {
        order[n] = c;
        n += 1;
    }
    for c in 0u8..128 {
        if !first.contains(&c) {
            order[n] = c;
            n += 1;
        }
    }
    order
}
