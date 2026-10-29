use std::io::Cursor;

use comsa::archive::{self, ArchiveReader};
use comsa::{CompressOptions, Error, Family, codec, stockholm};

/// Small deterministic PRNG (xorshift64*).
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
}

const SYMBOLS: &[u8] = b"-.ACDEFGHIKLMNPQRSTVWYXBZUOJacdefghiklmnpqrstvwy*~";

/// Random family whose rows are mutated copies of a few ancestors, with annotation
/// lines before, between and after the sequences.
fn random_family(rng: &mut Rng, id: usize, rows: usize, cols: usize) -> Family {
    let ancestors: Vec<Vec<u8>> =
        (0..4).map(|_| (0..cols).map(|_| SYMBOLS[rng.below(SYMBOLS.len())]).collect()).collect();
    let mut f = Family::default();
    f.meta.push(b"# STOCKHOLM 1.0".to_vec());
    f.meta.push(format!("#=GF ID   fam{id}").into_bytes());
    f.meta.push(format!("#=GF AC   PF{id:05}.1").into_bytes());
    let mut since_meta = 0;
    for r in 0..rows {
        let mut s = ancestors[rng.below(4)].clone();
        for _ in 0..cols / 20 {
            let c = rng.below(cols.max(1));
            s[c] = SYMBOLS[rng.below(SYMBOLS.len())];
        }
        f.names.push(format!("seq{r}/1-{cols}{}", " ".repeat(1 + rng.below(3))).into_bytes());
        f.sequences.push(s);
        since_meta += 1;
        if rng.below(5) == 0 {
            f.meta.push(format!("#=GR seq{r}/1-{cols} SS {}", "H".repeat(cols)).into_bytes());
            f.offsets.push(since_meta);
            since_meta = 0;
        }
    }
    f.meta.push(format!("#=GC seq_cons {}", "x".repeat(cols)).into_bytes());
    f.offsets.push(since_meta);
    f
}

fn to_text(fams: &[Family]) -> Vec<u8> {
    let mut v = Vec::new();
    for f in fams {
        stockholm::write_family(&mut v, f).unwrap();
    }
    v
}

#[test]
fn family_roundtrip_all_variants() {
    let mut rng = Rng(42);
    for (rows, cols) in [(1, 1), (1, 300), (7, 325), (300, 50), (2000, 3)] {
        let f = random_family(&mut rng, rows, rows, cols);
        for fast in [false, true] {
            let opts = CompressOptions { fast, threads: 1 };
            let block = codec::compress_family(&f, &opts).unwrap();
            assert_eq!(codec::decompress_family(&block.data, 1).unwrap(), f, "rows={rows} cols={cols} fast={fast}");
        }
    }
}

#[test]
fn large_family_parallel_matches_sequential() {
    // Above the 4M-cell threshold, so the threaded pipeline is used.
    let mut rng = Rng(7);
    let f = random_family(&mut rng, 1, 2100, 2000);
    for fast in [false, true] {
        let seq = codec::compress_family(&f, &CompressOptions { fast, threads: 1 }).unwrap();
        let par = codec::compress_family(&f, &CompressOptions { fast, threads: 8 }).unwrap();
        assert_eq!(seq.data, par.data);
        assert_eq!(codec::decompress_family(&par.data, 8).unwrap(), f);
        assert_eq!(codec::decompress_family(&par.data, 1).unwrap(), f);
    }
}

#[test]
fn edge_case_families() {
    let opts = CompressOptions { fast: false, threads: 1 };
    let meta_only = Family { meta: vec![b"#=GF ID x".to_vec()], ..Default::default() };
    let empty_seqs = Family {
        meta: vec![b"#=GF ID y".to_vec()],
        names: vec![b"a ".to_vec(), b"b ".to_vec()],
        sequences: vec![Vec::new(), Vec::new()],
        ..Default::default()
    };
    for f in [meta_only, empty_seqs] {
        let block = codec::compress_family(&f, &opts).unwrap();
        assert_eq!(codec::decompress_family(&block.data, 1).unwrap(), f);
    }
}

#[test]
fn rejects_invalid_input() {
    let opts = CompressOptions::default();
    let mut f = Family {
        names: vec![b"a ".to_vec(), b"b ".to_vec()],
        sequences: vec![b"AC-D".to_vec(), b"AC-E".to_vec()],
        ..Default::default()
    };
    f.sequences[1][2] = b'}';
    assert!(matches!(codec::compress_family(&f, &opts), Err(Error::InvalidInput(_))));
    f.sequences[1] = b"AC-".to_vec();
    assert!(matches!(codec::compress_family(&f, &opts), Err(Error::InvalidInput(_))));
    f.sequences[1] = b"AC-E".to_vec();
    f.names[1] = b"a   ".to_vec();
    assert!(matches!(codec::compress_family(&f, &opts), Err(Error::InvalidInput(_))));
}

#[test]
fn archive_roundtrip_and_random_access() {
    let mut rng = Rng(3);
    let fams: Vec<Family> = (0..40)
        .map(|i| {
            let (rows, cols) = (1 + rng.below(60), 1 + rng.below(200));
            random_family(&mut rng, i, rows, cols)
        })
        .collect();
    let text = to_text(&fams);

    let reader = stockholm::Reader::new(&text[..]);
    let mut arc = Vec::new();
    let stats = archive::compress_stockholm(reader, &mut arc, &CompressOptions { fast: false, threads: 4 }, |_| {}).unwrap();
    assert_eq!(stats.families, 40);
    assert_eq!(stats.raw_bytes, text.len() as u64);

    // Single-threaded compression gives the same archive.
    let mut arc1 = Vec::new();
    archive::compress_stockholm(stockholm::Reader::new(&text[..]), &mut arc1, &CompressOptions { fast: false, threads: 1 }, |_| {})
        .unwrap();
    assert_eq!(arc, arc1);

    let mut ar = ArchiveReader::open(Cursor::new(&arc)).unwrap();
    let mut out = Vec::new();
    assert_eq!(archive::decompress_stockholm(&mut ar, &mut out, 4, |_| {}).unwrap(), 40);
    assert_eq!(out, text);

    let descs = ar.families().to_vec();
    assert_eq!(descs[17].id, "fam17");
    assert_eq!(descs[17].ac, "PF00017.1");
    assert_eq!(ar.read_family(17, 1).unwrap(), fams[17]);
}

#[test]
fn damaged_archives_fail_cleanly() {
    let mut rng = Rng(5);
    let fams: Vec<Family> = (0..5).map(|i| random_family(&mut rng, i, 30, 80)).collect();
    let text = to_text(&fams);
    let mut arc = Vec::new();
    archive::compress_stockholm(stockholm::Reader::new(&text[..]), &mut arc, &CompressOptions::default(), |_| {}).unwrap();

    let decode = |data: &[u8]| -> comsa::Result<Vec<u8>> {
        let mut ar = ArchiveReader::open(Cursor::new(data))?;
        let mut out = Vec::new();
        archive::decompress_stockholm(&mut ar, &mut out, 2, |_| {})?;
        Ok(out)
    };
    for len in [0, 5, 8, arc.len() / 2, arc.len() - 1] {
        assert!(decode(&arc[..len]).is_err(), "truncated to {len}");
    }
    for i in (0..arc.len()).step_by(37) {
        let mut d = arc.clone();
        d[i] ^= 0x10;
        // Must not panic or hang. (The alignment stream has no checksum, so a flip there
        // is not guaranteed to be detected.)
        let _ = decode(&d);
    }
}
