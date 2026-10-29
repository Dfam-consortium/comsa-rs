# Validation against Dfam

These scripts reproduce the round-trip validation described in the top-level README.
Downloaded data and working files go in `validation/data/` and `validation/work/`,
which are git-ignored. The full run needs about 4 GB of disk.

You need three binaries:

* `RS`: this project (`target/release/comsa-rs`)
* `CF`: CoMSA 1.1 with `patches/comsa-1.1-fixes.patch` applied
* `CO`: unmodified CoMSA 1.1

```sh
cd validation
python3 fetch.py data/files 10000 500        # 500 curated (DF) + 9,500 uncurated (DR) seeds

python3 verify.py data/files work/files "$RS" "$CF" "$CO" 16
python3 archive_test.py data/files work/archive work/files/results.json "$RS" "$CF" "$CO"
```

`verify.py` checks each file separately:

* `comsa` compress, then decompress, gives back the input, after only the documented
  normalization (blank lines dropped, CRLF to LF, `//` appended).
* The archive is byte-identical to the patched C++ build's archive, for both the
  default and `-f` variants.
* Rust archives decode correctly with the patched C++ build and with unmodified
  CoMSA 1.1 (skipped for files containing `~`, which 1.1 cannot handle).
* `Sl` works.
* The first family round-trips as aligned FASTA (`Fc`/`Fd`), and the patched C++
  build decodes Rust's FASTA archive.
* When a file is rejected, both tools must reject it, and neither may leave output
  behind.

`archive_test.py` puts all accepted families into one archive. It checks:

* Rust, patched C++, unmodified C++, and Rust with `@list` input: archive sizes,
  byte-identity, and times.
* Every decoder against every archive.
* `Sl` lists every family.
* 200 random `Se -eAC` extractions.
* xz and gzip sizes, for reference.

Note: unmodified CoMSA 1.1 exits with status 0 even when it fails (for example on
an unknown option), so the scripts compare its output instead of trusting its exit
status.
