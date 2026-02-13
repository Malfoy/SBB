# biobloom-rs

A Rust reimplementation of the core `biobloom` workflow with a focus on throughput and low overhead.

Implemented binaries:
- `biobloommaker`
- `biobloomcategorizer`
- `biobloomrecruit`

## Build

```bash
cargo build --release
```

## Quick Start

Create a filter:

```bash
./target/release/biobloommaker -p ecoli -o filters -k 25 -f 0.0075 ref1.fa ref2.fa.gz
```

Categorize reads:

```bash
./target/release/biobloomcategorizer \
  -f "filters/ecoli.bf filters/other.bf" \
  --fq -p run1 \
  reads_1.fq.gz reads_2.fq.gz
```

Recruit reads against one filter:

```bash
./target/release/biobloomrecruit \
  -f filters/ecoli.bf \
  --fq -p recruited \
  reads.fq.gz
```

## Performance Design

- Canonical 2-bit k-mer rolling iterator (`A/C/G/T`, reset on ambiguous bases).
- Double hashing (`h1 + i*h2`) to minimize hashing cost per k-mer.
- Lock-free bit setting with `AtomicU64::fetch_or` during filter construction.
- Batched classification with `rayon` parallel scoring.
- Streaming FASTA/FASTQ(+gz) parsing via `needletail`.

## CLI Notes

`biobloommaker` supports:
- `-p/--file_prefix`
- `-o/--output_dir`
- `-f/--fal_pos_rate`
- `-g/--hash_num`
- `-k/--kmer_size`
- `-n/--num_ele`
- `--bit_len`
- blocked layout is default (`--block_words` tunes block size)
- `--classic` to force classic (non-blocked) layout
- `-r/--progressive`
- `-s/--subtract`
- `-e/--iterations`
- `--seed_files`
- `-t/--threads`

`biobloomcategorizer` supports:
- `-p/--prefix`
- `-f/--filter_files`
- `-e/--paired_mode` (exactly two files)
- `-i/--inclusive`
- `-s/--score`
- `-b/--best_hit`
- `-w/--with_score`
- `-t/--threads`
- `-g/--gz_output`
- `--fa` or `--fq`

`biobloomrecruit` supports single-filter recruiting with analogous output and scoring options.

## Output Behavior

`biobloomcategorizer` can emit:
- per-filter files: `<prefix>_<filter_id>.fa|fq[.gz]`
- `<prefix>_nomatch.fa|fq[.gz]`
- `<prefix>_multimatch.fa|fq[.gz]`

It also prints tab-separated summary counts to stdout.

## Format Compatibility

This project currently uses an internal Rust Bloom filter format (`.bf`) and is **not** guaranteed to be binary-compatible with upstream C++ `biobloom` filter files.

## Tests

```bash
cargo test
```
