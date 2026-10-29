//! Command-line interface, argument-compatible with CoMSA 1.1.

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::process::ExitCode;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Instant;

use comsa::archive::{self, ArchiveReader};
use comsa::{CompressOptions, Error, Result, codec, fasta, stockholm};

const IO_BUFFER: usize = 4 << 20;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Mode {
    FastaCompress,
    FastaDecompress,
    StockholmCompress,
    StockholmDecompress,
    StockholmExtract,
    StockholmList,
}

struct Args {
    mode: Mode,
    inputs: Vec<String>,
    output: String,
    wrap: usize,
    opts: CompressOptions,
    extract_id: Option<String>,
    extract_ac: Option<String>,
}

fn usage() {
    eprintln!(
        "Usage: comsa-rs <mode> [options] <in_file> <out_file>
       comsa-rs <mode> [options] @<in_file_list> <out_file>
Modes:
   Fc - compress FASTA file           Fd - decompress FASTA file
   Sc - compress Stockholm file(s)    Sd - decompress Stockholm archive
   Se - extract families from a Stockholm archive
   Sl - list families in a Stockholm archive (no out_file)
Options:
   -w <width>  wrap sequences in FASTA output (Fd only); default 0 (no wrapping)
   -f          fast variant (MTF instead of WFC)
   -t <n>      number of threads; default: number of cores (max 16)
   -eID <id>   extract family with this ID (Se only)
   -eAC <ac>   extract family with this accession (Se only)
   @file       file with one Stockholm file name per line (Sc only)
   --version   print the version and exit
Examples:
   comsa-rs Sc pfam.sto pfam.smsac
   comsa-rs Sd pfam.smsac pfam.sto
   comsa-rs Se -eAC PF00005.26 pfam.smsac PF00005.sto
   comsa-rs Fd -w 60 PF00005.msac PF00005.fasta"
    );
}

fn parse_args() -> std::result::Result<Args, String> {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    if matches!(argv.first().map(String::as_str), Some("-version" | "--version" | "-V")) {
        println!("comsa-rs {}", env!("CARGO_PKG_VERSION"));
        std::process::exit(0);
    }
    let mode = match argv.first().map(String::as_str) {
        Some("Fc") => Mode::FastaCompress,
        Some("Fd") => Mode::FastaDecompress,
        Some("Sc") => Mode::StockholmCompress,
        Some("Sd") => Mode::StockholmDecompress,
        Some("Se") => Mode::StockholmExtract,
        Some("Sl") => Mode::StockholmList,
        Some(m) => return Err(format!("invalid mode: {m}")),
        None => return Err(String::new()),
    };
    let n_files = if mode == Mode::StockholmList { 1 } else { 2 };
    let mut args = Args {
        mode,
        inputs: Vec::new(),
        output: String::new(),
        wrap: 0,
        opts: CompressOptions::default(),
        extract_id: None,
        extract_ac: None,
    };
    let mut i = 1;
    while i + n_files < argv.len() {
        let value = |i: usize| argv.get(i + 1).cloned().ok_or_else(|| format!("option {} needs a value", argv[i]));
        match argv[i].as_str() {
            "-w" => args.wrap = value(i)?.parse().map_err(|_| "bad -w value".to_string())?,
            "-t" => args.opts.threads = value(i)?.parse::<usize>().map_err(|_| "bad -t value".to_string())?.max(1),
            "-eID" => args.extract_id = Some(value(i)?),
            "-eAC" => args.extract_ac = Some(value(i)?),
            "-f" => {
                args.opts.fast = true;
                i += 1;
                continue;
            }
            o => return Err(format!("invalid option: {o}")),
        }
        i += 2;
    }
    if argv.len() != i + n_files {
        return Err(String::new());
    }
    let input = &argv[i];
    if let Some(list) = input.strip_prefix('@').filter(|_| mode == Mode::StockholmCompress) {
        let f = File::open(list).map_err(|e| format!("cannot open {list}: {e}"))?;
        for line in BufReader::new(f).lines() {
            let line = line.map_err(|e| format!("cannot read {list}: {e}"))?;
            let name = line.trim();
            if !name.is_empty() {
                args.inputs.push(name.to_string());
            }
        }
    } else {
        args.inputs.push(input.clone());
    }
    if n_files == 2 {
        args.output = argv[i + 1].clone();
    }
    Ok(args)
}

fn open_in(path: &str) -> Result<File> {
    File::open(path).map_err(|e| Error::InvalidInput(format!("cannot open {path}: {e}")))
}

/// Set once the output file has been created, so a failed run only removes its own output.
static OUTPUT_CREATED: AtomicBool = AtomicBool::new(false);

fn create_out(path: &str) -> Result<BufWriter<File>> {
    OUTPUT_CREATED.store(true, Ordering::Relaxed);
    let f = File::create(path).map_err(|e| Error::InvalidInput(format!("cannot create {path}: {e}")))?;
    Ok(BufWriter::with_capacity(IO_BUFFER, f))
}

fn stockholm_compress(a: &Args) -> Result<()> {
    let mut readers = Vec::new();
    for name in &a.inputs {
        readers.push((name, BufReader::with_capacity(IO_BUFFER, open_in(name)?)));
    }
    let families = readers.into_iter().flat_map(|(_, r)| stockholm::Reader::new(r));
    let started = Instant::now();
    let stats = archive::compress_stockholm(families, create_out(&a.output)?, &a.opts, |p| {
        if p.families % 256 == 0 {
            eprint!("\rfamilies: {}", p.families);
        }
    })
    .map_err(|e| match (e, a.inputs.len()) {
        (Error::InvalidInput(m), 1) => Error::InvalidInput(format!("{}: {m}", a.inputs[0])),
        (e, _) => e,
    })?;
    let total = stats.text_bytes + stats.seq_bytes;
    eprint!("\r");
    println!("Families               : {}", stats.families);
    println!("Total size of metadata : {}B", stats.text_bytes);
    println!("Total size of alignment: {}B", stats.seq_bytes);
    if total > 0 {
        println!("Compression ratio      : {:.3}", stats.raw_bytes as f64 / total as f64);
    }
    println!("Compression time       : {:.3} s", started.elapsed().as_secs_f64());
    Ok(())
}

fn stockholm_decompress(a: &Args) -> Result<()> {
    let started = Instant::now();
    let mut arc = ArchiveReader::open(open_in(&a.inputs[0])?)?;
    let n = archive::decompress_stockholm(&mut arc, create_out(&a.output)?, a.opts.threads, |_| {})?;
    println!("Families           : {n}");
    println!("Decompression time : {:.3} s", started.elapsed().as_secs_f64());
    Ok(())
}

fn stockholm_extract(a: &Args) -> Result<()> {
    let started = Instant::now();
    let mut arc = ArchiveReader::open(open_in(&a.inputs[0])?)?;
    let wanted: Vec<usize> = (0..arc.families().len())
        .filter(|&i| {
            let d = &arc.families()[i];
            a.extract_id.as_ref().is_none_or(|x| *x == d.id) && a.extract_ac.as_ref().is_none_or(|x| *x == d.ac)
        })
        .collect();
    if wanted.is_empty() {
        return Err(Error::InvalidInput("no family matches the given ID/AC".into()));
    }
    let mut out = create_out(&a.output)?;
    for &i in &wanted {
        let f = arc.read_family(i, a.opts.threads)?;
        stockholm::write_family(&mut out, &f)?;
    }
    out.flush()?;
    println!("Families extracted : {}", wanted.len());
    println!("Extraction time    : {:.3} s", started.elapsed().as_secs_f64());
    Ok(())
}

fn stockholm_list(a: &Args) -> Result<()> {
    let arc = ArchiveReader::open(open_in(&a.inputs[0])?)?;
    let stdout = std::io::stdout();
    let mut out = stdout.lock();
    writeln!(out, "ID\t AC\t no. sequences\t no. columns\t uncompressed size\t compressed size")?;
    for d in arc.families() {
        writeln!(out, "{}\t {}\t {}\t {}\t {}\t {}", d.id, d.ac, d.n_sequences, d.n_columns, d.raw_size, d.compressed_size)?;
    }
    Ok(())
}

fn fasta_compress(a: &Args) -> Result<()> {
    let started = Instant::now();
    let (names, seqs) = fasta::read(BufReader::with_capacity(IO_BUFFER, open_in(&a.inputs[0])?))
        .map_err(|e| match e {
            Error::InvalidInput(m) => Error::InvalidInput(format!("{}: {m}", a.inputs[0])),
            e => e,
        })?;
    let block = codec::compress_fasta(&names, &seqs, &a.opts).map_err(|e| match e {
        Error::InvalidInput(m) => Error::InvalidInput(format!("{}: {m}", a.inputs[0])),
        e => e,
    })?;
    let mut out = create_out(&a.output)?;
    out.write_all(&block.data)?;
    out.flush()?;
    println!("Sequences compressed to: {} bytes", block.seq_size);
    println!("Names compressed to    : {} bytes", block.text_size);
    println!("Total size             : {} bytes", block.seq_size + block.text_size);
    println!("Compression time       : {:.3} s", started.elapsed().as_secs_f64());
    Ok(())
}

fn fasta_decompress(a: &Args) -> Result<()> {
    let started = Instant::now();
    let data = std::fs::read(&a.inputs[0])?;
    let (names, seqs) = codec::decompress_fasta(&data, a.opts.threads)?;
    let mut out = create_out(&a.output)?;
    fasta::write(&mut out, &names, &seqs, a.wrap)?;
    out.flush()?;
    println!("Decompression time : {:.3} s", started.elapsed().as_secs_f64());
    Ok(())
}

fn main() -> ExitCode {
    let args = match parse_args() {
        Ok(a) => a,
        Err(msg) => {
            if !msg.is_empty() {
                eprintln!("comsa-rs: {msg}");
            }
            usage();
            return ExitCode::from(2);
        }
    };
    let result = match args.mode {
        Mode::FastaCompress => fasta_compress(&args),
        Mode::FastaDecompress => fasta_decompress(&args),
        Mode::StockholmCompress => stockholm_compress(&args),
        Mode::StockholmDecompress => stockholm_decompress(&args),
        Mode::StockholmExtract => stockholm_extract(&args),
        Mode::StockholmList => stockholm_list(&args),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("\ncomsa-rs: {e}");
            // Don't leave a partial output file behind.
            if OUTPUT_CREATED.load(Ordering::Relaxed) && Path::new(&args.output).is_file() {
                let _ = std::fs::remove_file(&args.output);
            }
            ExitCode::FAILURE
        }
    }
}
