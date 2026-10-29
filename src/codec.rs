//! Compression of a single alignment (one Stockholm family or one FASTA file) into a
//! CoMSA block.
//!
//! Block layout (all integers are [`varint`](crate::varint)s):
//! `flags` byte (context size + 64 if fast variant), n_sequences, n_columns,
//! text_size, seq_size, n_tokens, LZMA text (`text_size` bytes), range-coded alignment
//! (`seq_size` bytes).

use std::collections::HashSet;
use std::thread;

use crate::entropy::{CtxLength, EntropyDecoder, EntropyEncoder};
use crate::error::{Error, Result, corrupt};
use crate::mtf::Mtf;
use crate::parallel::ordered_map;
use crate::pbwt::Pbwt;
use crate::stockholm::Family;
use crate::wfc::Wfc;
use crate::{lzma, rle, varint};

const LZMA_PRESET: u32 = 9;
/// Alignments smaller than this (rows * columns) are coded on the calling thread.
const PARALLEL_MIN_CELLS: u64 = 1 << 22;
const BLOCK_COLUMNS: usize = 64;

#[derive(Clone, Debug)]
pub struct CompressOptions {
    /// Use MTF instead of WFC: faster, slightly worse compression.
    pub fast: bool,
    /// Worker threads (for archives: families compressed concurrently; for large
    /// alignments: threads used inside one alignment).
    pub threads: usize,
}

impl Default for CompressOptions {
    fn default() -> Self {
        CompressOptions { fast: false, threads: default_threads() }
    }
}

pub fn default_threads() -> usize {
    thread::available_parallelism().map_or(4, |n| n.get()).min(16)
}

/// A compressed alignment block.
#[derive(Clone, Debug)]
pub struct EncodedBlock {
    pub data: Vec<u8>,
    /// Size of the LZMA-compressed names/annotation.
    pub text_size: usize,
    /// Size of the range-coded alignment.
    pub seq_size: usize,
}

// ---------------------------------------------------------------------------------------
// Public API

/// Compress one Stockholm family. Fails on invalid alignment data.
pub fn compress_family(f: &Family, opts: &CompressOptions) -> Result<EncodedBlock> {
    validate_alignment(&f.names, &f.sequences, true).map_err(|e| with_family(e, f))?;
    let mut text = Vec::new();
    put_lines(&mut text, &f.meta);
    put_lines(&mut text, &f.names);
    varint::put(&mut text, f.offsets.len() as u64);
    for &o in &f.offsets {
        varint::put(&mut text, o as u64);
    }
    compress_block(&text, &f.sequences, false, opts)
}

/// Decompress one Stockholm family block.
pub fn decompress_family(data: &[u8], threads: usize) -> Result<Family> {
    let (text, sequences) = decompress_block(data, threads)?;
    let mut pos = 0;
    let meta = get_lines(&text, &mut pos)?;
    let names = get_lines(&text, &mut pos)?;
    let n_offsets = varint::get(&text, &mut pos)? as usize;
    if n_offsets > meta.len() {
        return corrupt("bad annotation offsets");
    }
    let offsets = (0..n_offsets).map(|_| varint::get(&text, &mut pos).map(|x| x as u32)).collect::<Result<Vec<_>>>()?;
    if names.len() != sequences.len() {
        return corrupt("number of names does not match number of sequences");
    }
    Ok(Family { meta, offsets, names, sequences })
}

/// Compress an aligned FASTA file given as names (including the leading `>`) and sequences.
pub fn compress_fasta(names: &[Vec<u8>], sequences: &[Vec<u8>], opts: &CompressOptions) -> Result<EncodedBlock> {
    validate_alignment(names, sequences, false)?;
    let mut text = Vec::new();
    put_lines(&mut text, names);
    compress_block(&text, sequences, true, opts)
}

/// Decompress a FASTA block into names (with leading `>`) and sequences.
pub fn decompress_fasta(data: &[u8], threads: usize) -> Result<crate::fasta::Records> {
    let (text, sequences) = decompress_block(data, threads)?;
    let names = get_lines(&text, &mut 0)?;
    if names.len() != sequences.len() {
        return corrupt("number of names does not match number of sequences");
    }
    Ok((names, sequences))
}

/// Check that all sequences have the same length and use only valid symbols; with
/// `unique_names`, also that no name repeats (as in interleaved Stockholm).
pub fn validate_alignment(names: &[Vec<u8>], sequences: &[Vec<u8>], unique_names: bool) -> Result<()> {
    let Some(first) = sequences.first() else { return Ok(()) };
    let name = |i: usize| String::from_utf8_lossy(name_token(names.get(i).map_or(&[][..], |n| n))).into_owned();
    if unique_names {
        let mut seen = HashSet::with_capacity(names.len());
        for n in names {
            if !seen.insert(name_token(n)) {
                return Err(Error::InvalidInput(format!(
                    "sequence name '{}' occurs more than once (interleaved/multi-block Stockholm is not supported)",
                    String::from_utf8_lossy(name_token(n))
                )));
            }
        }
    }
    for (i, s) in sequences.iter().enumerate() {
        if s.len() != first.len() {
            return Err(Error::InvalidInput(format!(
                "sequence '{}' has {} columns but '{}' has {} (all sequences must have the same length; \
                 interleaved/multi-block Stockholm is not supported)",
                name(i),
                s.len(),
                name(0),
                first.len()
            )));
        }
        if let Some(col) = s.iter().position(|&c| !crate::is_valid_sequence_symbol(c)) {
            let c = s[col];
            let shown = if c.is_ascii_graphic() { format!("'{}'", c as char) } else { format!("byte 0x{c:02X}") };
            return Err(Error::InvalidInput(format!(
                "invalid symbol {shown} in sequence '{}' at column {} (allowed: A-Z a-z - . * ~)",
                name(i),
                col + 1
            )));
        }
    }
    Ok(())
}

fn name_token(n: &[u8]) -> &[u8] {
    let end = n.iter().position(|&c| c == b' ' || c == b'\t').unwrap_or(n.len());
    &n[..end]
}

fn with_family(e: Error, f: &Family) -> Error {
    match e {
        Error::InvalidInput(msg) => {
            let id = f.id();
            Error::InvalidInput(if id.is_empty() { msg } else { format!("family {id}: {msg}") })
        }
        e => e,
    }
}

// ---------------------------------------------------------------------------------------
// Block assembly

fn compress_block(text: &[u8], sequences: &[Vec<u8>], extreme: bool, opts: &CompressOptions) -> Result<EncodedBlock> {
    let n = sequences.len();
    let n_cols = sequences.first().map_or(0, |s| s.len());
    let cells = (n * n_cols) as u64;
    let ctx = CtxLength::for_cells(cells);

    let (text_c, (seq_c, n_tokens)) = if opts.threads > 1 && cells >= PARALLEL_MIN_CELLS {
        thread::scope(|s| {
            let t = s.spawn(|| lzma::compress(text, LZMA_PRESET, extreme));
            let a = encode_alignment(sequences, ctx, opts.fast, opts.threads);
            (t.join().expect("LZMA thread panicked"), a)
        })
    } else {
        (lzma::compress(text, LZMA_PRESET, extreme), encode_alignment(sequences, ctx, opts.fast, 1))
    };
    let text_c = text_c?;

    let mut data = Vec::with_capacity(text_c.len() + seq_c.len() + 32);
    data.push(ctx as u8 + if opts.fast { 64 } else { 0 });
    varint::put(&mut data, n as u64);
    varint::put(&mut data, if n_tokens > 0 { n_cols as u64 } else { 0 });
    varint::put(&mut data, text_c.len() as u64);
    varint::put(&mut data, seq_c.len() as u64);
    varint::put(&mut data, n_tokens);
    data.extend_from_slice(&text_c);
    data.extend_from_slice(&seq_c);
    Ok(EncodedBlock { data, text_size: text_c.len(), seq_size: seq_c.len() })
}

fn decompress_block(data: &[u8], threads: usize) -> Result<(Vec<u8>, Vec<Vec<u8>>)> {
    let Some(&flags) = data.first() else { return corrupt("empty block") };
    let fast = flags & 64 != 0;
    let Some(ctx) = CtxLength::from_u8(flags & !64) else { return corrupt("bad block header") };
    let mut pos = 1;
    let n = varint::get(data, &mut pos)?;
    let n_cols = varint::get(data, &mut pos)?;
    let text_size = varint::get(data, &mut pos)?;
    let seq_size = varint::get(data, &mut pos)?;
    let n_tokens = varint::get(data, &mut pos)?;
    if text_size.saturating_add(seq_size) > (data.len() - pos) as u64 {
        return corrupt("block shorter than its header says");
    }
    // Every cell needs at least one coded bit, so this bounds memory on damaged input.
    let cells = n.saturating_mul(n_cols);
    if n > u32::MAX as u64 || n_cols > u32::MAX as u64 || cells > seq_size.saturating_mul(8 * 64) + 64 {
        return corrupt("implausible alignment size");
    }
    let (n, n_cols) = (n as usize, n_cols as usize);
    let text_c = &data[pos..pos + text_size as usize];
    let seq_c = &data[pos + text_size as usize..pos + (text_size + seq_size) as usize];

    if threads > 1 && cells >= PARALLEL_MIN_CELLS {
        thread::scope(|s| {
            let t = s.spawn(|| lzma::decompress(text_c));
            let a = decode_alignment(seq_c, ctx, fast, n, n_cols, n_tokens, threads);
            Ok((t.join().expect("LZMA thread panicked")?, a?))
        })
    } else {
        Ok((lzma::decompress(text_c)?, decode_alignment(seq_c, ctx, fast, n, n_cols, n_tokens, 1)?))
    }
}

fn put_lines(v: &mut Vec<u8>, lines: &[Vec<u8>]) {
    let len: usize = lines.iter().map(|l| l.len() + 1).sum();
    varint::put(v, len as u64);
    v.reserve(len);
    for l in lines {
        v.extend_from_slice(l);
        v.push(b'\n');
    }
}

fn get_lines(text: &[u8], pos: &mut usize) -> Result<Vec<Vec<u8>>> {
    let len = varint::get(text, pos)? as usize;
    if len > text.len() - *pos {
        return corrupt("text section truncated");
    }
    let mut lines: Vec<Vec<u8>> = text[*pos..*pos + len].split(|&c| c == b'\n').map(|l| l.to_vec()).collect();
    if lines.last().is_some_and(|l| l.is_empty()) {
        lines.pop();
    }
    *pos += len;
    Ok(lines)
}

// ---------------------------------------------------------------------------------------
// Alignment coding: columns right to left; PBWT -> WFC/MTF -> RLE-0 -> range coder

enum SecondStage {
    Wfc(Box<Wfc>),
    Mtf(Box<Mtf>),
}

impl SecondStage {
    fn new(fast: bool) -> Self {
        if fast { SecondStage::Mtf(Box::new(Mtf::new())) } else { SecondStage::Wfc(Box::new(Wfc::new())) }
    }

    fn forward(&mut self, src: &[u8], out: &mut Vec<u8>) {
        match self {
            SecondStage::Wfc(w) => w.forward(src, out),
            SecondStage::Mtf(m) => m.forward(src, out),
        }
    }

    fn reverse(&mut self, src: &[u8], out: &mut Vec<u8>) {
        match self {
            SecondStage::Wfc(w) => w.reverse(src, out),
            SecondStage::Mtf(m) => m.reverse(src, out),
        }
    }
}

/// Column blocks [lo, hi) from the right end of the alignment.
fn column_blocks(n_cols: usize) -> impl Iterator<Item = (usize, usize)> {
    (0..n_cols.div_ceil(BLOCK_COLUMNS)).map(move |b| {
        let hi = n_cols - b * BLOCK_COLUMNS;
        (hi.saturating_sub(BLOCK_COLUMNS), hi)
    })
}

/// Copy columns [lo, hi) of the rows into `cols` (column-major).
fn gather_columns(rows: &[Vec<u8>], lo: usize, hi: usize, cols: &mut Vec<Vec<u8>>) {
    cols.resize_with(hi - lo, Vec::new);
    for c in cols.iter_mut() {
        c.clear();
        c.reserve(rows.len());
    }
    for row in rows {
        for (c, &x) in cols.iter_mut().zip(&row[lo..hi]) {
            c.push(x);
        }
    }
}

/// Returns (range coder output, number of RLE tokens).
fn encode_alignment(rows: &[Vec<u8>], ctx: CtxLength, fast: bool, threads: usize) -> (Vec<u8>, u64) {
    let n = rows.len();
    if n == 0 {
        return (Vec::new(), 0);
    }
    let n_cols = rows[0].len();
    let mut ent = EntropyEncoder::new(ctx);

    if threads <= 1 || ((n * n_cols) as u64) < PARALLEL_MIN_CELLS {
        let mut pbwt = Pbwt::new(n);
        let mut ss = SecondStage::new(fast);
        let (mut cols, mut p, mut r, mut t) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        for (lo, hi) in column_blocks(n_cols) {
            gather_columns(rows, lo, hi, &mut cols);
            for col in cols.iter().rev() {
                pbwt.forward(col, &mut p);
                ss.forward(&p, &mut r);
                rle::encode(&r, &mut t);
                ent.encode_column(&t);
            }
        }
    } else {
        let _ = ordered_map(
            threads,
            threads * 2,
            |send| {
                let mut pbwt = Pbwt::new(n);
                let mut cols = Vec::new();
                for (lo, hi) in column_blocks(n_cols) {
                    gather_columns(rows, lo, hi, &mut cols);
                    let batch: Vec<Vec<u8>> = cols
                        .iter()
                        .rev()
                        .map(|col| {
                            let mut p = Vec::new();
                            pbwt.forward(col, &mut p);
                            p
                        })
                        .collect();
                    if !send(batch) {
                        break;
                    }
                }
                Ok(())
            },
            |batch: Vec<Vec<u8>>| {
                let mut ss = SecondStage::new(fast);
                let mut r = Vec::new();
                batch
                    .into_iter()
                    .map(|p| {
                        let mut t = Vec::new();
                        ss.forward(&p, &mut r);
                        rle::encode(&r, &mut t);
                        t
                    })
                    .collect::<Vec<_>>()
            },
            |tokens: Vec<Vec<u8>>| {
                for t in &tokens {
                    ent.encode_column(t);
                }
                Ok(())
            },
        );
    }
    let n_tokens = ent.n_tokens;
    (ent.finish(), n_tokens)
}

fn decode_alignment(
    data: &[u8],
    ctx: CtxLength,
    fast: bool,
    n: usize,
    n_cols: usize,
    n_tokens: u64,
    threads: usize,
) -> Result<Vec<Vec<u8>>> {
    if n_tokens == 0 {
        return Ok(vec![Vec::new(); n]);
    }
    if n == 0 || n_cols == 0 {
        return corrupt("alignment data present but no rows or columns");
    }
    let mut rows = vec![vec![0u8; n_cols]; n];
    let mut dec = EntropyDecoder::new(ctx, data);

    // Columns of one block arrive right to left; scatter them back into the rows.
    let scatter = |rows: &mut [Vec<u8>], lo: usize, cols: &[Vec<u8>]| {
        for (r, row) in rows.iter_mut().enumerate() {
            for (k, x) in row[lo..lo + cols.len()].iter_mut().rev().enumerate() {
                *x = cols[k][r];
            }
        }
    };

    if threads <= 1 || ((n * n_cols) as u64) < PARALLEL_MIN_CELLS {
        let mut pbwt = Pbwt::new(n);
        let mut ss = SecondStage::new(fast);
        let (mut r, mut p) = (Vec::new(), Vec::new());
        let mut cols: Vec<Vec<u8>> = Vec::new();
        for (lo, hi) in column_blocks(n_cols) {
            cols.resize_with(hi - lo, Vec::new);
            for col in cols.iter_mut() {
                dec.decode_column(n, &mut r)?;
                ss.reverse(&r, &mut p);
                pbwt.reverse(&p, col);
            }
            scatter(&mut rows, lo, &cols);
        }
    } else {
        let mut pbwt = Pbwt::new(n);
        let mut blocks = column_blocks(n_cols);
        let dec = &mut dec;
        ordered_map(
            threads,
            threads * 2,
            |send| {
                for (lo, hi) in column_blocks(n_cols) {
                    let mut batch = Vec::with_capacity(hi - lo);
                    for _ in lo..hi {
                        let mut r = Vec::new();
                        dec.decode_column(n, &mut r)?;
                        batch.push(r);
                    }
                    if !send(batch) {
                        break;
                    }
                }
                Ok(())
            },
            |batch: Vec<Vec<u8>>| {
                let mut ss = SecondStage::new(fast);
                batch
                    .iter()
                    .map(|r| {
                        let mut p = Vec::new();
                        ss.reverse(r, &mut p);
                        p
                    })
                    .collect::<Vec<_>>()
            },
            |batch: Vec<Vec<u8>>| {
                let (lo, _) = blocks.next().expect("more column blocks than expected");
                let cols: Vec<Vec<u8>> = batch
                    .iter()
                    .map(|p| {
                        let mut col = Vec::new();
                        pbwt.reverse(p, &mut col);
                        col
                    })
                    .collect();
                scatter(&mut rows, lo, &cols);
                Ok(())
            },
        )?;
    }
    if dec.n_tokens != n_tokens {
        return corrupt("alignment stream length mismatch");
    }
    Ok(rows)
}
