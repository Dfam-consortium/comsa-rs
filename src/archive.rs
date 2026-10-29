//! The multi-family Stockholm archive (`.smsac`) and the FASTA container (`.msac`).
//!
//! `.smsac` layout: a sequence of blocks, each an 8-byte little-endian length followed by
//! the block; then a footer with one record per family (varints n_sequences, n_columns,
//! raw_size, compressed_size, block offset; NUL-terminated ID and AC); then the footer
//! length as 8-byte little-endian. `.msac` is a single block with no framing.

use std::io::{self, BufReader, Read, Seek, SeekFrom, Write};

use crate::codec::{self, CompressOptions};
use crate::error::{Error, Result, corrupt};
use crate::parallel::ordered_map;
use crate::stockholm::{self, Family};
use crate::varint;

/// Footer record describing one family.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct FamilyDesc {
    pub n_sequences: u64,
    pub n_columns: u64,
    /// Size of the family in the original Stockholm input.
    pub raw_size: u64,
    /// Compressed size of names/annotation plus alignment.
    pub compressed_size: u64,
    /// Offset of the family's block in the archive.
    pub offset: u64,
    pub id: String,
    pub ac: String,
}

/// Writes blocks and, on [`finish`](Self::finish), the footer.
pub struct ArchiveWriter<W: Write> {
    out: W,
    pos: u64,
    families: Vec<FamilyDesc>,
}

impl<W: Write> ArchiveWriter<W> {
    pub fn new(out: W) -> Self {
        ArchiveWriter { out, pos: 0, families: Vec::new() }
    }

    /// Append a compressed family block.
    pub fn add(&mut self, block: &codec::EncodedBlock, mut desc: FamilyDesc) -> io::Result<()> {
        desc.offset = self.pos;
        desc.compressed_size = (block.text_size + block.seq_size) as u64;
        self.out.write_all(&(block.data.len() as u64).to_le_bytes())?;
        self.out.write_all(&block.data)?;
        self.pos += 8 + block.data.len() as u64;
        self.families.push(desc);
        Ok(())
    }

    pub fn families(&self) -> &[FamilyDesc] {
        &self.families
    }

    /// Write the footer and flush; returns the underlying writer.
    pub fn finish(mut self) -> io::Result<W> {
        let mut footer = Vec::new();
        for d in &self.families {
            for x in [d.n_sequences, d.n_columns, d.raw_size, d.compressed_size, d.offset] {
                varint::put(&mut footer, x);
            }
            for s in [&d.id, &d.ac] {
                footer.extend(s.bytes().filter(|&c| c != 0));
                footer.push(0);
            }
        }
        footer.extend_from_slice(&(footer.len() as u64).to_le_bytes());
        self.out.write_all(&footer)?;
        self.out.flush()?;
        Ok(self.out)
    }
}

/// Random-access reader of a `.smsac` archive.
pub struct ArchiveReader<R: Read + Seek> {
    input: BufReader<R>,
    data_size: u64,
    families: Vec<FamilyDesc>,
}

impl<R: Read + Seek> ArchiveReader<R> {
    pub fn open(mut input: R) -> Result<Self> {
        let size = input.seek(SeekFrom::End(0))?;
        if size < 8 {
            return corrupt("file too small to be a CoMSA archive");
        }
        input.seek(SeekFrom::Start(size - 8))?;
        let mut b = [0u8; 8];
        input.read_exact(&mut b)?;
        let footer_size = u64::from_le_bytes(b);
        if footer_size > size - 8 {
            return corrupt("bad footer size (truncated or not a CoMSA archive)");
        }
        let data_size = size - 8 - footer_size;
        input.seek(SeekFrom::Start(data_size))?;
        let mut footer = vec![0u8; footer_size as usize];
        input.read_exact(&mut footer)?;

        let mut families = Vec::new();
        let mut pos = 0;
        let cstr = |pos: &mut usize| -> Result<String> {
            let Some(len) = footer[*pos..].iter().position(|&c| c == 0) else {
                return corrupt("bad footer");
            };
            let s = String::from_utf8_lossy(&footer[*pos..*pos + len]).into_owned();
            *pos += len + 1;
            Ok(s)
        };
        while pos < footer.len() {
            let mut v = [0u64; 5];
            for x in &mut v {
                *x = varint::get(&footer, &mut pos)?;
            }
            let id = cstr(&mut pos)?;
            let ac = cstr(&mut pos)?;
            if v[4] >= data_size {
                return corrupt("bad block offset in footer");
            }
            families.push(FamilyDesc {
                n_sequences: v[0],
                n_columns: v[1],
                raw_size: v[2],
                compressed_size: v[3],
                offset: v[4],
                id,
                ac,
            });
        }
        input.seek(SeekFrom::Start(0))?;
        Ok(ArchiveReader { input: BufReader::with_capacity(1 << 20, input), data_size, families })
    }

    pub fn families(&self) -> &[FamilyDesc] {
        &self.families
    }

    /// Read the raw block starting at `offset`; returns the block and the next offset.
    fn read_block_at(&mut self, offset: u64) -> Result<(Vec<u8>, u64)> {
        self.input.seek(SeekFrom::Start(offset))?;
        let mut b = [0u8; 8];
        self.input.read_exact(&mut b)?;
        let len = u64::from_le_bytes(b);
        if len > self.data_size.saturating_sub(offset + 8) {
            return corrupt("block extends past end of data");
        }
        let mut block = vec![0u8; len as usize];
        self.input.read_exact(&mut block)?;
        Ok((block, offset + 8 + len))
    }

    /// Decompress the family with index `i` (in archive order).
    pub fn read_family(&mut self, i: usize, threads: usize) -> Result<Family> {
        let offset = self.families.get(i).ok_or_else(|| Error::InvalidInput("no such family".into()))?.offset;
        let (block, _) = self.read_block_at(offset)?;
        codec::decompress_family(&block, threads)
    }

    /// Iterate over all raw blocks in file order.
    pub fn blocks(&mut self) -> impl Iterator<Item = Result<Vec<u8>>> + '_ {
        let mut offset = 0u64;
        std::iter::from_fn(move || {
            if offset >= self.data_size {
                return None;
            }
            match self.read_block_at(offset) {
                Ok((b, next)) => {
                    offset = next;
                    Some(Ok(b))
                }
                Err(e) => {
                    offset = self.data_size;
                    Some(Err(e))
                }
            }
        })
    }
}

/// Progress information passed to callbacks.
#[derive(Clone, Copy, Debug, Default)]
pub struct Progress {
    pub families: u64,
    pub raw_bytes: u64,
    pub text_bytes: u64,
    pub seq_bytes: u64,
}

/// Compress Stockholm families into an archive, `opts.threads` families at a time.
/// Families are written in input order.
pub fn compress_stockholm<I, W>(families: I, out: W, opts: &CompressOptions, mut progress: impl FnMut(&Progress)) -> Result<Progress>
where
    I: Iterator<Item = Result<(Family, u64)>> + Send,
    W: Write,
{
    let mut writer = ArchiveWriter::new(out);
    let mut stats = Progress::default();
    // Small families are coded single-threaded inside their worker; large ones also
    // parallelize internally.
    let inner = opts.clone();
    ordered_map(
        opts.threads,
        opts.threads * 2,
        |send| {
            for f in families {
                if !send(f) {
                    break;
                }
            }
            Ok(())
        },
        |f: Result<(Family, u64)>| -> Result<(codec::EncodedBlock, FamilyDesc)> {
            let (f, raw_size) = f?;
            let block = codec::compress_family(&f, &inner)?;
            let desc = FamilyDesc {
                n_sequences: f.sequences.len() as u64,
                n_columns: f.n_columns() as u64,
                raw_size,
                id: f.id(),
                ac: f.ac(),
                ..Default::default()
            };
            Ok((block, desc))
        },
        |r| {
            let (block, desc) = r?;
            stats.families += 1;
            stats.raw_bytes += desc.raw_size;
            stats.text_bytes += block.text_size as u64;
            stats.seq_bytes += block.seq_size as u64;
            writer.add(&block, desc)?;
            progress(&stats);
            Ok(())
        },
    )?;
    writer.finish()?;
    Ok(stats)
}

/// Decompress every family of an archive to Stockholm text, in archive order.
pub fn decompress_stockholm<R: Read + Seek + Send, W: Write>(
    archive: &mut ArchiveReader<R>,
    mut out: W,
    threads: usize,
    mut progress: impl FnMut(u64),
) -> Result<u64> {
    let inner = threads;
    let mut n = 0u64;
    ordered_map(
        threads,
        threads * 2,
        |send| {
            for b in archive.blocks() {
                let failed = b.is_err();
                if !send(b) || failed {
                    break;
                }
            }
            Ok(())
        },
        |b: Result<Vec<u8>>| -> Result<Vec<u8>> {
            let f = codec::decompress_family(&b?, inner)?;
            let mut text = Vec::new();
            stockholm::write_family(&mut text, &f)?;
            Ok(text)
        },
        |text| {
            out.write_all(&text?)?;
            n += 1;
            progress(n);
            Ok(())
        },
    )?;
    out.flush()?;
    Ok(n)
}
