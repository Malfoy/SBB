# Benchmarks

Each script is standalone and accepts:

```bash
python3 benchmarks/<script>.py <index_file> <query_file>
```

Outputs:
- one TSV in `benchmarks/results/`
- one PNG plot in `benchmarks/results/`

The scripts parse these lines from `bloomybloom` output:
- `total indexing time (s)`
- `index cpu time (s)`
- `total query time (s)`
- `query cpu time (s)`

Requirements:
- `python3`
- `matplotlib`
- a built `./target/release/bloomybloom` binary, or let the script build it

Scripts:
- `bench_k.py`: sweep `k`
- `bench_hashes.py`: sweep `n_hashes`
- `bench_threads.py`: sweep thread count
- `bench_s.py`: sweep `s`
- `bench_block_size.py`: sweep block-size exponent

All tunable parameters are defined near the top of each file.
