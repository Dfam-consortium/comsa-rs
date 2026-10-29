# comsa-rs

Fast, lossless compression of multiple sequence alignments in Stockholm and aligned
FASTA format: a Rust port of **CoMSA**, with bug fixes and parallel processing.

## Credit: CoMSA

The compression method, the archive format and the original implementation are the
work of **Sebastian Deorowicz, Joanna Walczyszyn and Agnieszka Debudaj-Grabysz**
([REFRESH Bioinformatics Group](https://github.com/refresh-bio), Silesian University
of Technology). This project is a translation of their CoMSA 1.1 into Rust. The
algorithms are entirely theirs. If you use this software, please cite their paper:

> Deorowicz S, Walczyszyn J, Debudaj-Grabysz A. **CoMSA: compression of protein
> multiple sequence alignment files.** *Bioinformatics* 35(2):227–234 (2019).
> [doi:10.1093/bioinformatics/bty619](https://doi.org/10.1093/bioinformatics/bty619)

Original software: <https://github.com/refresh-bio/CoMSA>

CoMSA codes the alignment column by column. Each column goes through a generalized
positional Burrows–Wheeler transform (PBWT), then a weighted-frequency-count (or
move-to-front) transform, zero-run-length coding and a context-modelled range coder.
Names and annotation are compressed with LZMA. The range coder comes from the same
group's FaStore and ORCOM tools (after Subbotin's carryless range coder).

## Why a port?

* **Fits the Dfam/RepeatMasker Rust ecosystem.** It is usable as a library as well
  as a command-line tool.
* **Fixes bugs found in review.** CoMSA 1.1 silently corrupts sequences containing
  `}` and makes undecodable archives from `~`. It crashes on digits, `?` or
  non-ASCII bytes, loops forever on inputs whose size is a multiple of 128 MiB, and
  hangs on truncated archives. See [CHANGELOG.md](CHANGELOG.md).
* **Much faster on archives of many families.** Families are processed in parallel,
  and per-family start-up costs were removed. See [Performance](#performance).
* **Format compatible.** Archives are interchangeable with CoMSA 1.1 in both
  directions; see [Compatibility](#compatibility-with-comsa).

The same fixes are also available for the C++ code as a patch; see
[Patching CoMSA 1.1 (C++)](#patching-comsa-11-c).

## Installation

### Prebuilt binaries

Download the archive for your platform (Linux x86_64 or macOS arm64) from the
[releases page](https://github.com/Dfam-consortium/comsa-rs/releases) and unpack it:

```sh
tar xzf comsa-<version>-<target>.tar.gz
comsa-<version>/bin/comsa --version
```

### From source

Requires a Rust toolchain (<https://rustup.rs>, Rust 1.85 or newer) and a C compiler
(liblzma is built from source and linked statically).

```sh
git clone https://github.com/Dfam-consortium/comsa-rs.git
cd comsa-rs
make                         # target/release/comsa
make test                    # optional
make install PREFIX=/opt     # installs /opt/comsa-<version>/bin/comsa and docs
```

or simply `cargo install --git https://github.com/Dfam-consortium/comsa-rs`.

### As a library

```toml
[dependencies]
comsa = { git = "https://github.com/Dfam-consortium/comsa-rs" }
```

## Usage

The command line is the same as CoMSA's, with `-t` and `--version` added.

```text
comsa <mode> [options] <in_file> <out_file>

Modes:
  Sc   compress Stockholm file(s) into an archive (.smsac)
  Sd   decompress all families of an archive
  Se   extract families by ID or accession (-eID / -eAC)
  Sl   list families in an archive (no out_file)
  Fc   compress an aligned FASTA file (.msac)
  Fd   decompress an aligned FASTA file

Options:
  -t <n>      threads (default: number of cores, max 16)
  -f          fast variant (MTF instead of WFC; about 1% larger)
  -w <width>  line width for FASTA output (Fd); 0 = one line per sequence
  -eID <id>   select family by #=GF ID (Se)
  -eAC <ac>   select family by #=GF AC (Se)
  @file       (Sc) file listing Stockholm files to compress, one per line
```

Examples:

```sh
comsa Sc Dfam.seed.sto Dfam.seed.smsac
comsa Sc @stockholm_files.txt families.smsac
comsa Sl Dfam.seed.smsac | head
comsa Se -eAC DF000000001.4 Dfam.seed.smsac MIR.sto
comsa Sd Dfam.seed.smsac Dfam.seed.sto
comsa Fc PF00005.fasta PF00005.msac && comsa Fd -w 60 PF00005.msac PF00005.fasta
```

Any failure exits with status 1 and removes the partial output file.

### Library

```rust
use std::fs::File;
use std::io::{BufReader, BufWriter};
use comsa::{archive, stockholm, CompressOptions};

// Compress a Stockholm file (families are compressed in parallel)
let reader = stockholm::Reader::new(BufReader::new(File::open("in.sto")?));
archive::compress_stockholm(reader, BufWriter::new(File::create("out.smsac")?),
                            &CompressOptions::default(), |_progress| {})?;

// Random access to one family
let mut arc = archive::ArchiveReader::open(File::open("out.smsac")?)?;
let i = arc.families().iter().position(|d| d.ac == "DF000000001.4").unwrap();
let family = arc.read_family(i, 4)?;          // comsa::Family
stockholm::write_family(&mut std::io::stdout(), &family)?;
```

`codec::compress_family` / `codec::decompress_family` (Stockholm) and
`codec::compress_fasta` / `codec::decompress_fasta` work on single in-memory blocks.

## Input rules

* **Sequence symbols:** `A`–`Z`, `a`–`z`, `-`, `.`, `*`, and `~` (Easel's
  missing-data symbol). Anything else is rejected with the family, sequence and
  column.
* **Shape:** all sequences in a family must have the same length and unique names.
  Interleaved (multi-block) Stockholm is not supported and is rejected.
* **Normalization:** blank lines are dropped, CR/CRLF line ends become LF, and every
  family ends with `//`. Files already in this form, such as Pfam and Dfam
  downloads, round-trip byte for byte.

## Performance

These are crude measurements: wall-clock times from single runs on a shared
2 × Intel Xeon Silver 4314 server (64 hardware threads), with a load average of
45–72 from other jobs. Treat differences under about 30% as noise. "Patched C++" is
CoMSA 1.1 with `patches/comsa-1.1-fixes.patch`. Times are compress / decompress, in
seconds, with default thread counts.

| Data set | CoMSA 1.1 | Patched C++ | comsa-rs |
|---|---|---|---|
| Pfam sample: 6 families, 34 MB (CoMSA's sample file) | 2.06 / 0.65 | 1.3 / 0.54 | 1.3 / 0.6 |
| 3,000 small Pfam families, 12 MB | 15.3 / 7.8 | 0.34 / 0.25 | 0.32 / 0.27 |
| One alignment, 19,205 × 1,437, as FASTA (28 MB) | 0.68 / 0.55 | – | 0.65 / 0.58 |
| Input of exactly 128 MiB | never finishes | 1.82 / – | 1.28 / – |
| **Dfam: 10,016 random seed alignments, 3.43 GB** | **179 / 91** | **47 / 19** | **13.5 / 12.3** |

On the Dfam sample, single-threaded comsa-rs (`-t 1`) decompresses in 131 s.

Compression ratios are identical for all three implementations:

| Data set | comsa | `xz -9` | `gzip -9` |
|---|---|---|---|
| Pfam sample | 15.6× | 13.9× | – |
| Dfam sample (3.43 GB → 67.6 MB) | 50.7× (13.5 s) | 44.3× (83.5 s, 16 threads) | 20.8× (564 s) |

Most of the speed-up over CoMSA 1.1 on many-family archives comes from three changes,
all of which are also in the C++ patch:

* The weight table of the WFC transform is built once instead of per family and
  thread. It was about 60% of CPU time on small families.
* The LZMA dictionary is sized to the input instead of always 64 MiB.
* Whole families are processed in parallel.

Within a single large alignment, both implementations are limited by the sequential
PBWT and range-coding stages and perform about the same.

## Validation

Correctness was checked in five ways:

1. **Bit-exactness with CoMSA.**
   * The WFC weight table is compared against one dumped from the C++ code (a unit
     test).
   * For every input tested, comsa-rs writes archives byte-identical to the patched
     C++ build, in default and fast (`-f`) modes, Stockholm and FASTA.
   * Its alignment streams are byte-identical to unmodified CoMSA 1.1.
   * Archives cross-decode in every direction between comsa-rs, patched C++ and
     unmodified CoMSA 1.1.
2. **Dfam round trip.** 10,016 seed alignments, downloaded at random from the Dfam
   API (505 curated `DF` and 9,511 uncurated `DR` families, 3.43 GB, up to 2,071
   sequences and 55,910 columns):
   * Every file round-tripped byte for byte through comsa-rs (default and `-f`) and
     through the FASTA path.
   * Every archive was byte-identical to the patched C++ archive and decoded
     correctly with patched and unmodified CoMSA 1.1.
   * All 10,016 families in a single 3.43 GB archive decoded exactly with all three
     decoders, including the `@list` input path.
   * 200 random `Se` extractions matched.
   * The scripts are in [`validation/`](validation/).
3. **Edge cases.**
   * Each rejected symbol class (`}`, `?`, digits, non-ASCII bytes) and interleaved
     input gives an error, exit status 1, and no output file.
   * Blank lines, CRLF and a missing final `//` are normalized; `~` and lowercase
     round-trip.
   * Families with no sequences, and families whose sequences are all empty,
     round-trip.
   * FASTA files without a final newline round-trip.
   * An input of exactly 128 MiB round-trips.
4. **Damaged archives.** Truncations, empty and 5-byte files, and 80 random
   single-bit flips of a Pfam archive all gave a clean error, with no crash or hang.
   Only the LZMA part carries a checksum, so damage to the alignment stream is
   *usually* but not *always* detected.
5. **Tests.** `cargo test` runs unit tests (range coder, WFC/MTF, varints, Stockholm
   and FASTA parsing) and integration tests:
   * random families in every size class, including ones large enough for the
     parallel path, which must give byte-identical output to the sequential path;
   * archive round trip and random access;
   * input rejection;
   * damaged-archive handling.

   CI runs clippy, the tests, and a check that the C++ patch still applies to
   CoMSA v1.1.

## Compatibility with CoMSA

* Archives (`.smsac` for Stockholm, `.msac` for FASTA) are readable by CoMSA 1.1, and
  CoMSA 1.1 archives are readable by comsa-rs.
* The one exception is alignments containing `~`. CoMSA 1.1 never handled `~`
  correctly; comsa-rs and the patched C++ give it a new symbol rank, so such
  archives need one of those two to decode.
* The LZMA dictionary size differs from CoMSA 1.1, so whole archives are not
  byte-identical to 1.1's. The alignment data is identical, and archives are
  byte-identical to the patched C++ build.
* CoMSA's `master` branch has an untagged version 1.2 that adds gzip input and a
  sequences-only extraction mode. Its codec is unchanged from 1.1. Those features
  are not in comsa-rs yet.

The block and archive layouts are documented in `src/codec.rs` and `src/archive.rs`.

## Patching CoMSA 1.1 (C++)

If you would rather use the original C++ code, `patches/comsa-1.1-fixes.patch`
applies the same fixes and speed-ups to the CoMSA v1.1 release. The format stays
compatible. The patch also adds a `-t <threads>` option.

```sh
curl -LO https://github.com/refresh-bio/CoMSA/archive/refs/tags/v1.1.tar.gz
tar xzf v1.1.tar.gz
cd CoMSA-1.1
patch -p1 < /path/to/comsa-rs/patches/comsa-1.1-fixes.patch
make                      # Linux/g++; produces ./CoMSA
./CoMSA Sc ../input.sto input.smsac
```

The patch applies only to the tagged 1.1 release, not to the 1.2 code on CoMSA's
`master` branch. The Visual Studio project files are not updated, but the patch adds
no new source files. The top of the patch file lists every change.

## License

CoMSA's source files state that it is "distributed under GNU GPL 3 licence". As a
derivative work, comsa-rs is released under the **GNU General Public License,
version 3** (`GPL-3.0-only`); see [LICENSE](LICENSE). The C++ patch is likewise
GPL-3.

The range coder in CoMSA came from the same group's FaStore (GPL-3) and ORCOM
(GPL-2) projects. liblzma, which is statically linked, is public domain / 0BSD.
