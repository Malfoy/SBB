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
  -f "filters/ecoli.bf.zst filters/other.bf.zst" \
  --fq -p run1 \
  reads_1.fq.gz reads_2.fq.gz
```

Recruit reads against one filter:

```bash
./target/release/sbb recruit \
  -f filters/ecoli.bf.zst \
  --fq -p recruited \
  reads.fq.gz
```

## Performance Design

- Canonical 2-bit k-mer rolling iterator (`A/C/G/T`, reset on ambiguous bases).
- Double hashing (`h1 + i*h2`) to minimize hashing cost per k-mer.
- Lock-free bit setting with `AtomicU64::fetch_or` during filter construction.
- Batched classification with `rayon` parallel scoring.
- Streaming FASTA/FASTQ parsing via SIMD-accelerated `helicase` (plain, `.gz`, `.zst`, `.xz`, `.bz2`).

## CLI Notes

`sbb maker` supports:
- `-p/--file_prefix`
- `-o/--output_dir`
- `-f/--fal_pos_rate`
- `-g/--hash_num`
- `-k/--kmer_size`
- `-n/--num_ele`
- `--sample_alpha`
- `--bit_len`
- blocked layout is default (`--block_words` tunes block size)
- `--classic` to force classic (non-blocked) layout
- `--bloomybloom` to use the vendored minimizer-blocked implementation from `bloomybloom`
- `-r/--progressive`
- `-s/--subtract`
- `-e/--iterations`
- `--seed_files`
- `-t/--threads`

`--bloomybloom` notes:
- requires an odd `-k/--kmer_size`
- supports k-mer sizes above `32`
- does not support `--sample_alpha < 1.0`
- does not support `--subtract` during progressive insertion

Input reads/references can be FASTA or FASTQ, plain or compressed as:
- `.gz`
- `.zst`
- `.xz`
- `.bz2`

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

`sbb maker` writes one compressed index file (`.bf.zst`) by default.
The internal Bloom payload format is Rust-specific and is **not** guaranteed to be binary-compatible with upstream C++ `biobloom` filter files.

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
- Format roundtrip checks for classic, blocked, and bloomybloom layouts

## Benchmarks

Single-reference comparison benchmarks (C++ vs Rust classic vs Rust blocked):

- `bench/benchmark_single_ref_baseline.sh`
- `bench/benchmark_single_ref_accuracy.sh`
- `bench/benchmark_single_ref_sweep.sh`
- `bench/benchmark_single_ref_scaling.sh`
- `bench/benchmark_fixed_geometry.sh`
- `bench/benchmark_random_queries.sh` (indexes a provided FASTA, uses provided queries or generates random ones in Rust, reports positive hits per tool)

Each benchmark is independent and writes a TSV plus plots (via its paired plot script):

- `bench/plot_single_ref_baseline.py`
- `bench/plot_single_ref_accuracy.py`
- `bench/plot_single_ref_sweep.py`
- `bench/plot_single_ref_scaling.py`

Plot scripts require `python3` and `matplotlib`.

All benchmark scripts accept `REPEATS` (default `3`) as an environment parameter.

Minimal example:

```bash
REPEATS=5 bash bench/benchmark_single_ref_baseline.sh
```

Common environment variables:
- `BENCH_ROOT` benchmark output root
- `REPEATS` number of repeats per benchmark point
- `THREADS` worker threads (or `THREAD_LIST` for scaling)
- `REF_LEN`, `READ_LEN`, `READ_COUNT`, `SEED`
- `KMER_SIZE`, `FPR`, `CAT_SCORE`, `BLOCK_WORDS`
- `CPP_BIN_DIR` repo-local directory for upstream C++ `biobloommaker` and `biobloomcategorizer` (default `bench/bin`)
- `CPP_LIB_DIR` optional repo-local shared-library directory for upstream C++ runtime deps (default `bench/lib` if present)
- `CPP_MAKER`, `CPP_CAT`
- `SBB_BIN`, `BENCHGEN`
- `PLOT` (`1` to auto-generate plots, `0` to skip)

Fixed-geometry benchmark specific:
- `INDEX_DATASET`, `QUERY_DATASET`
- `BIT_LEN` power-of-two Bloom size in bits
- `HASH_NUM` fixed hash count
- `CPP_NUM_ELE` required when `RUN_CPP=1`; the script derives BioBloomMaker FPR from `BIT_LEN`, `HASH_NUM`, and `CPP_NUM_ELE`
- `RUN_CPP`, `RUN_RUST_CLASSIC`, `RUN_RUST_BLOCKED`, `RUN_RUST_BLOOMY`

Random-query benchmark specific:
- `INDEX_FASTA`
- `QUERY_FASTA` optional prebuilt FASTA query dataset; when set, random-query generation is skipped
- `QUERY_COUNT`, `QUERY_LEN`, `SEED`
- `HASH_NUM` for fixed hash count
- `BIT_EXP` for fixed bloom size as a power-of-two exponent (`BIT_EXP=31` means `2^31` bits)
- `CAT_SCORE` (default `0.001`)
- comparable-memory mode is default; to bypass it set `ALLOW_NON_COMPARABLE=1`

Accuracy benchmark specific:
- `POS_COUNT`, `NEAR_COUNT`, `NEG_COUNT`, `NEAR_MUTATION_RATE`

Sweep benchmark specific:
- `SWEEP_MODES` (`fpr,hash,block_words`)
- `FPR_LIST`, `HASH_LIST`, `BLOCK_WORDS_LIST`
- `FIXED_FPR`, `FIXED_HASH`, `FIXED_BLOCK_WORDS`

Scaling benchmark specific:
- `THREAD_LIST` (for example `1,2,4,8,16,32`)
