# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## 0.1.0 - 2026-10-29
Initial release: a Rust port of CoMSA 1.1 (Deorowicz, Walczyszyn &
Debudaj-Grabysz), format compatible with it.
### Breaking Changes
- None (archives are interchangeable with CoMSA 1.1)
### Added
- `comsa` command-line tool with the same modes and options as CoMSA 1.1
  (`Sc`, `Sd`, `Se`, `Sl`, `Fc`, `Fd`, `-f`, `-w`, `-eID`, `-eAC`, `@list`)
- `-t <threads>`: families are compressed and decompressed in parallel, and large
  alignments are also pipelined internally
- `--version`
- Library API: streaming Stockholm reader and writer, per-family codec, and an
  archive reader with random access to families
- Support for `~` (Easel's missing-data symbol) in sequences
- `patches/comsa-1.1-fixes.patch`: the same fixes for the CoMSA 1.1 C++ sources
- `validation/`: scripts for the Dfam round-trip validation
### Fixed (relative to CoMSA 1.1)
- Sequence symbols outside `A-Z a-z - . * ~` are rejected with a clear error.
  CoMSA 1.1 silently corrupted `}`, made undecodable archives from `~`, and
  aborted on digits, `?` and non-ASCII bytes.
- Unequal sequence lengths and duplicate names (interleaved Stockholm) are
  reported clearly
- Inputs whose size is a multiple of 128 MiB no longer loop until memory runs out
- Damaged or truncated archives fail with an error instead of hanging or crashing
- Non-zero exit status on failure, and the partial output file is removed
- FASTA input without a trailing newline no longer loses its last byte
- Families whose sequences are all empty round-trip
### Changed
- The LZMA dictionary is sized to the input instead of always 64 MiB. Compression
  is the same, it is much faster on small families, and older CoMSA versions still
  read the result.
### Known Issues
- Interleaved (multi-block) Stockholm is not supported
- Only the names/annotation part of a block carries a checksum (LZMA CRC64);
  damage to the alignment stream is usually, but not always, detected
