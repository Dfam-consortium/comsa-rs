//! Aligned FASTA reading and writing.

use std::io::{self, Read, Write};

use crate::error::{Error, Result};

/// Sequence names (with leading `>`) and sequences.
pub type Records = (Vec<Vec<u8>>, Vec<Vec<u8>>);

/// Read an aligned FASTA file. Names keep their leading `>`; sequence line breaks
/// (LF or CRLF) are removed.
pub fn read<R: Read>(mut input: R) -> Result<Records> {
    let mut data = Vec::new();
    input.read_to_end(&mut data)?;
    let start = data.iter().position(|&c| c != b'\n' && c != b'\r').unwrap_or(data.len());
    if start == data.len() {
        return Ok((Vec::new(), Vec::new()));
    }
    if data[start] != b'>' {
        return Err(Error::InvalidInput("FASTA input must start with '>'".into()));
    }
    let (mut names, mut seqs) = (Vec::new(), Vec::new());
    for record in data[start + 1..].split(|&c| c == b'>') {
        let name_end = record.iter().position(|&c| c == b'\n' || c == b'\r').unwrap_or(record.len());
        let mut name = Vec::with_capacity(name_end + 1);
        name.push(b'>');
        name.extend_from_slice(&record[..name_end]);
        names.push(name);
        seqs.push(record[name_end..].iter().copied().filter(|&c| c != b'\n' && c != b'\r').collect());
    }
    Ok((names, seqs))
}

/// Write FASTA; `width` = 0 writes each sequence on one line.
pub fn write<W: Write>(w: &mut W, names: &[Vec<u8>], seqs: &[Vec<u8>], width: usize) -> io::Result<()> {
    for (name, seq) in names.iter().zip(seqs) {
        w.write_all(name)?;
        w.write_all(b"\n")?;
        if width == 0 {
            w.write_all(seq)?;
            w.write_all(b"\n")?;
        } else {
            for chunk in seq.chunks(width) {
                w.write_all(chunk)?;
                w.write_all(b"\n")?;
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    #[test]
    fn read_without_trailing_newline() {
        let (n, s) = super::read(&b">a x\nAC-\nD\r\n>b\nAC-E"[..]).unwrap();
        assert_eq!(n, vec![b">a x".to_vec(), b">b".to_vec()]);
        assert_eq!(s, vec![b"AC-D".to_vec(), b"AC-E".to_vec()]);
    }
}
