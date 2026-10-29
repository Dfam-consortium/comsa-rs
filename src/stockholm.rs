//! Stockholm reading and writing.
//!
//! Only single-block (Pfam-style, non-interleaved) Stockholm is supported. Blank lines
//! are dropped, CR/CRLF line ends become LF, and a family always ends with `//`, so the
//! round trip is exact for files already in that normal form.

use std::io::{self, BufRead, Write};

use crate::error::Result;

/// One Stockholm family (alignment).
///
/// Annotation lines (`#...`) are kept verbatim in `meta`. The first
/// `meta.len() - offsets.len()` of them come before the first sequence; each later
/// annotation line `k` follows `offsets[k]` sequence lines after the previous
/// annotation line. `names[i]` includes the whitespace that separates it from
/// `sequences[i]`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Family {
    pub meta: Vec<Vec<u8>>,
    pub offsets: Vec<u32>,
    pub names: Vec<Vec<u8>>,
    pub sequences: Vec<Vec<u8>>,
}

impl Family {
    /// Value of the `#=GF ID` line, or "" if there is none.
    pub fn id(&self) -> String {
        self.gf_value(b"#=GF ID")
    }

    /// Value of the `#=GF AC` line, or "" if there is none.
    pub fn ac(&self) -> String {
        self.gf_value(b"#=GF AC")
    }

    pub fn n_columns(&self) -> usize {
        self.sequences.first().map_or(0, |s| s.len())
    }

    fn gf_value(&self, tag: &[u8]) -> String {
        let trim = |c: &u8| !matches!(c, b' ' | b'\n' | b'\r');
        self.meta
            .iter()
            .rev()
            .find(|l| l.starts_with(tag))
            .map(|l| {
                let v = &l[tag.len()..];
                let s = v.iter().position(trim).unwrap_or(v.len());
                let e = v.iter().rposition(trim).map_or(s, |e| e + 1);
                String::from_utf8_lossy(&v[s..e]).into_owned()
            })
            .unwrap_or_default()
    }
}

/// Streaming reader yielding `(family, bytes of input consumed)`.
pub struct Reader<R: BufRead> {
    input: R,
    line: Vec<u8>,
    eof: bool,
}

impl<R: BufRead> Reader<R> {
    pub fn new(input: R) -> Self {
        Reader { input, line: Vec::new(), eof: false }
    }

    /// Next non-empty line (without terminator) into `self.line`; returns bytes consumed.
    fn read_line(&mut self) -> io::Result<(bool, u64)> {
        self.line.clear();
        let mut consumed = 0u64;
        loop {
            let buf = match self.input.fill_buf() {
                Ok(b) => b,
                Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
                Err(e) => return Err(e),
            };
            if buf.is_empty() {
                self.eof = true;
                return Ok((!self.line.is_empty(), consumed));
            }
            match buf.iter().position(|&c| c == b'\n' || c == b'\r') {
                Some(i) if i > 0 || !self.line.is_empty() => {
                    self.line.extend_from_slice(&buf[..i]);
                    self.input.consume(i + 1);
                    return Ok((true, consumed + i as u64 + 1));
                }
                Some(_) => {
                    // Line terminator before any content: skip it.
                    self.input.consume(1);
                    consumed += 1;
                }
                None => {
                    let len = buf.len();
                    self.line.extend_from_slice(buf);
                    self.input.consume(len);
                    consumed += len as u64;
                }
            }
        }
    }

    /// Read the next family; `Ok(None)` at end of input.
    pub fn read_family(&mut self) -> Result<Option<(Family, u64)>> {
        loop {
            if self.eof {
                return Ok(None);
            }
            let mut consumed = 0u64;
            let mut f = Family::default();
            let mut line_no = 0u32;
            let mut last_meta_line_no = 0u32;
            loop {
                let (got, n) = self.read_line()?;
                consumed += n;
                if !got || self.line == b"//" {
                    break;
                }
                line_no += 1;
                if self.line[0] == b'#' {
                    if line_no - last_meta_line_no > 1 || !f.sequences.is_empty() {
                        f.offsets.push(line_no - last_meta_line_no - 1);
                    }
                    last_meta_line_no = line_no;
                    f.meta.push(self.line.clone());
                } else {
                    // Name runs up to and including the whitespace before the sequence.
                    let l = &self.line;
                    let ws = l.iter().position(|&c| c == b' ' || c == b'\t').unwrap_or(l.len());
                    let seq = l[ws..].iter().position(|&c| c != b' ' && c != b'\t').map_or(l.len(), |p| ws + p);
                    f.names.push(l[..seq].to_vec());
                    f.sequences.push(l[seq..].to_vec());
                }
            }
            if !f.meta.is_empty() || !f.sequences.is_empty() {
                return Ok(Some((f, consumed)));
            }
        }
    }
}

impl<R: BufRead> Iterator for Reader<R> {
    type Item = Result<(Family, u64)>;

    fn next(&mut self) -> Option<Self::Item> {
        self.read_family().transpose()
    }
}

/// Write a family in Stockholm format, terminated by `//`.
pub fn write_family<W: Write>(w: &mut W, f: &Family) -> io::Result<()> {
    let n_leading = f.meta.len().saturating_sub(f.offsets.len());
    let mut meta = f.meta.iter();
    for m in meta.by_ref().take(n_leading) {
        w.write_all(m)?;
        w.write_all(b"\n")?;
    }
    let mut offsets = f.offsets.iter();
    let mut until_meta = offsets.next().copied();
    let mut flush_meta = |w: &mut W, until_meta: &mut Option<u32>| -> io::Result<()> {
        while *until_meta == Some(0) {
            if let Some(m) = meta.next() {
                w.write_all(m)?;
                w.write_all(b"\n")?;
            }
            *until_meta = offsets.next().copied();
        }
        Ok(())
    };
    for (name, seq) in f.names.iter().zip(&f.sequences) {
        flush_meta(w, &mut until_meta)?;
        w.write_all(name)?;
        w.write_all(seq)?;
        w.write_all(b"\n")?;
        if let Some(u) = until_meta.as_mut() {
            *u -= 1;
        }
    }
    flush_meta(w, &mut until_meta)?;
    w.write_all(b"//\n")
}

#[cfg(test)]
mod tests {
    use super::*;

    const TEXT: &[u8] = b"# STOCKHOLM 1.0\n#=GF ID   Fam1 \n#=GF AC   PF00001.1\n\nseq1/1-3   AC-D\n#=GR seq1/1-3 SS ..HH\nseq2/2-4\tA-CD\n#=GC RF xxxx\n//\n# STOCKHOLM 1.0\r\n#=GF ID Fam2\r\ns a.\r\n//\r\n";

    #[test]
    fn parse_and_write() {
        let fams: Vec<_> = Reader::new(TEXT).collect::<Result<_>>().unwrap();
        assert_eq!(fams.len(), 2);
        let f = &fams[0].0;
        assert_eq!(f.id(), "Fam1");
        assert_eq!(f.ac(), "PF00001.1");
        assert_eq!(f.names, vec![b"seq1/1-3   ".to_vec(), b"seq2/2-4\t".to_vec()]);
        assert_eq!(f.offsets, vec![1, 1]);
        let mut out = Vec::new();
        for (f, _) in &fams {
            write_family(&mut out, f).unwrap();
        }
        let expect = String::from_utf8_lossy(TEXT).replace("\r\n", "\n").replace("\n\n", "\n");
        assert_eq!(String::from_utf8_lossy(&out), expect);
        // As in CoMSA, the LF of the final CRLF is consumed after the last family ends.
        assert_eq!(fams[0].1 + fams[1].1, TEXT.len() as u64 - 1);
    }
}
