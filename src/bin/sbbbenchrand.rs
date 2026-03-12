use anyhow::{Context, Result, bail};
use clap::Parser;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
const IO_BUF_CAP: usize = 16 * 1024 * 1024;
const CHUNK_RECORDS: usize = 16_384;

#[derive(Parser, Debug)]
#[command(
    name = "sbbbenchrand",
    about = "Generate deterministic random FASTA queries with parallel Rust generation"
)]
struct Cli {
    #[arg(long = "out_fa")]
    out_fa: PathBuf,

    #[arg(long = "count")]
    count: usize,

    #[arg(long = "len", default_value_t = 1000)]
    len: usize,

    #[arg(long = "seed", default_value_t = 1337)]
    seed: u64,

    #[arg(short = 't', long = "threads")]
    threads: Option<usize>,
}

#[inline]
fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E3779B97F4A7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D049BB133111EB);
    x ^ (x >> 31)
}

#[inline]
fn digits(mut n: usize) -> usize {
    let mut d = 1;
    while n >= 10 {
        n /= 10;
        d += 1;
    }
    d
}

#[inline]
fn seq_seed(seed: u64, idx: usize) -> u64 {
    splitmix64(seed ^ (idx as u64).wrapping_mul(0xD6E8FEB86659FD93))
}

#[inline]
fn fill_random_sequence(mut state: u64, seq: &mut [u8]) {
    for base in seq.iter_mut() {
        state = splitmix64(state);
        *base = BASES[(state & 0b11) as usize];
    }
}

fn generate_chunk_bases(start_idx: usize, count: usize, len: usize, seed: u64) -> Vec<u8> {
    let mut bases = vec![0_u8; count * len];
    bases
        .par_chunks_mut(len)
        .enumerate()
        .for_each(|(local_idx, seq)| {
            let global_idx = start_idx + local_idx;
            fill_random_sequence(seq_seed(seed, global_idx), seq);
        });
    bases
}

fn write_chunk_fasta(
    w: &mut BufWriter<File>,
    bases: &[u8],
    start_idx: usize,
    count: usize,
    len: usize,
    width: usize,
) -> Result<()> {
    for local_idx in 0..count {
        let global_idx = start_idx + local_idx;
        writeln!(w, ">q{global_idx:0width$}")?;
        let begin = local_idx * len;
        let end = begin + len;
        w.write_all(&bases[begin..end])?;
        w.write_all(b"\n")?;
    }
    Ok(())
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    if cli.count == 0 {
        bail!("--count must be > 0");
    }
    if cli.len == 0 {
        bail!("--len must be > 0");
    }
    if let Some(t) = cli.threads {
        if t == 0 {
            bail!("--threads must be > 0");
        }
        let _ = rayon::ThreadPoolBuilder::new().num_threads(t).build_global();
    }

    if let Some(parent) = cli.out_fa.parent() {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("failed to create {}", parent.display()))?;
    }

    eprintln!(
        "generating random queries: count={} len={} seed={} threads={}",
        cli.count,
        cli.len,
        cli.seed,
        rayon::current_num_threads()
    );

    let out = File::create(&cli.out_fa)
        .with_context(|| format!("failed to create {}", cli.out_fa.display()))?;
    let mut w = BufWriter::with_capacity(IO_BUF_CAP, out);

    let width = digits(cli.count - 1);
    let mut start = 0usize;
    while start < cli.count {
        let remaining = cli.count - start;
        let chunk_count = remaining.min(CHUNK_RECORDS);
        let bases = generate_chunk_bases(start, chunk_count, cli.len, cli.seed);
        write_chunk_fasta(&mut w, &bases, start, chunk_count, cli.len, width)?;
        start += chunk_count;
    }

    w.flush()?;

    println!("queries_path={}", cli.out_fa.display());
    Ok(())
}
