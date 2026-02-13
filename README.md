# SBB

A Rust reimplementation of the core `biobloom` workflow with a focus on throughput and low overhead.

Implemented CLI:
- `sbb maker`
- `sbb categorizer`
- `sbb recruit`

## Build

```bash
cargo build --release
```

## Quick Start

Create a filter:

```bash
./target/release/sbb maker -p ecoli -o filters -k 25 -f 0.0075 ref1.fa ref2.fa.gz
```

Categorize reads:

```bash
./target/release/sbb categorizer \
  -f "filters/ecoli.bf filters/other.bf" \
  --fq -p run1 \
  reads_1.fq.gz reads_2.fq.gz
```

Recruit reads against one filter:

```bash
./target/release/sbb recruit \
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

`sbb maker` supports:
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

`sbb categorizer` supports:
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

`sbb recruit` supports single-filter recruiting with analogous output and scoring options.

## Output Behavior

`sbb categorizer` can emit:
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

Expanded coverage now includes:
- Bloom layout behavior and false-negative checks
- k-mer iterator edge cases (`N` reset, `k=1`, `k=32`, short reads)
- Bloom sizing/hash math sanity checks
- `sbb maker` CLI validation failures
- Progressive mode behavior (with and without subtract filter)
- Format roundtrip checks for classic and blocked layouts

## Benchmarks (Rust Only)

Dataset generation is done by the multithreaded Rust binary `sbbbenchgen` (release build), not Python.

Bloom-size sweep (classic and blocked for each point):

```bash
bash bench/benchmark_bloom_size.sh
```

Hash-count sweep (classic and blocked for each point):

```bash
bash bench/benchmark_hash_count.sh
```

Useful environment overrides for both scripts:
- `BENCH_ROOT` (default: `/tmp/sbb-bench-rust`)
- `THREADS` (default: `nproc`)
- `REF_LEN` (default: `100000000`)
- `READ_LEN` (default: `20000`)
- `READ_COUNT` (default: `150000`)
- `SEED` (default: `1337`)
- `KMER_SIZE` (default: `25`)
- `FPR` (default: `0.000001`)
- `CAT_SCORE` (default: `1.0`)
- `BLOCK_WORDS` (default: `1`)

Bloom-size script specific:
- `BITS_PER_ELEMENT_LIST` (default: `8,12,16,24,32`)
- `HASH_NUM` (default: `4`)

Hash-count script specific:
- `BITS_PER_ELEMENT` (default: `16`)
- `HASH_LIST` (default: `1,2,3,4,6,8`)
