use anyhow::{Context, Result, bail};
use clap::Parser;
use rayon::prelude::*;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
const WRAP: usize = 80;
const IO_BUF_CAP: usize = 16 * 1024 * 1024;

#[derive(Parser, Debug)]
#[command(
    name = "sbbbenchgen",
    about = "Generate deterministic FASTA benchmark datasets with parallel Rust generation"
)]
struct Cli {
    #[arg(long = "out_dir")]
    out_dir: PathBuf,

    #[arg(long = "ref_len", default_value_t = 60_000_000)]
    ref_len: usize,

    #[arg(long = "read_len", default_value_t = 20_000)]
    read_len: usize,

    #[arg(long = "read_count", default_value_t = 300_000)]
    read_count: usize,

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
fn base_for(seed: u64, idx: usize) -> u8 {
    let x = splitmix64(seed ^ (idx as u64).wrapping_mul(0xD6E8FEB86659FD93));
    BASES[(x & 0b11) as usize]
}

fn generate_reference(ref_len: usize, seed: u64) -> Vec<u8> {
    let mut ref_seq = vec![0_u8; ref_len];
    ref_seq
        .par_chunks_mut(1 << 20)
        .enumerate()
        .for_each(|(chunk_idx, chunk)| {
            let start = chunk_idx * (1 << 20);
            for (offset, base) in chunk.iter_mut().enumerate() {
                *base = base_for(seed, start + offset);
            }
        });
    ref_seq
}

fn generate_read_starts(read_count: usize, span: usize, seed: u64) -> Vec<usize> {
    (0..read_count)
        .into_par_iter()
        .map(|i| {
            let x =
                splitmix64(seed ^ 0xA24BAED4963EE407 ^ (i as u64).wrapping_mul(0x9E3779B97F4A7C15));
            (x % span as u64) as usize
        })
        .collect()
}

fn write_reference_fasta(path: &Path, seq: &[u8]) -> Result<()> {
    let file =
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let mut w = BufWriter::with_capacity(IO_BUF_CAP, file);
    w.write_all(b">ref1\n")?;
    for chunk in seq.chunks(WRAP) {
        w.write_all(chunk)?;
        w.write_all(b"\n")?;
    }
    w.flush()?;
    Ok(())
}

fn write_reads_fasta(path: &Path, ref_seq: &[u8], starts: &[usize], read_len: usize) -> Result<()> {
    let file =
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let mut w = BufWriter::with_capacity(IO_BUF_CAP, file);
    for (i, &start) in starts.iter().enumerate() {
        w.write_all(b">r")?;
        write!(w, "{i}")?;
        w.write_all(b"\n")?;
        w.write_all(&ref_seq[start..start + read_len])?;
        w.write_all(b"\n")?;
    }
    w.flush()?;
    Ok(())
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    if cli.ref_len == 0 {
        bail!("--ref_len must be > 0");
    }
    if cli.read_len == 0 {
        bail!("--read_len must be > 0");
    }
    if cli.read_len > cli.ref_len {
        bail!("--read_len must be <= --ref_len");
    }
    if let Some(t) = cli.threads {
        if t == 0 {
            bail!("--threads must be > 0");
        }
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(t)
            .build_global();
    }

    fs::create_dir_all(&cli.out_dir)
        .with_context(|| format!("failed to create {}", cli.out_dir.display()))?;

    let ref_path = cli.out_dir.join("ref.fa");
    let reads_path = cli.out_dir.join("reads.fa");
    let meta_path = cli.out_dir.join("dataset.meta");

    eprintln!(
        "generating reference: len={} seed={} threads={}",
        cli.ref_len,
        cli.seed,
        rayon::current_num_threads()
    );
    let ref_seq = generate_reference(cli.ref_len, cli.seed);
    write_reference_fasta(&ref_path, &ref_seq)?;

    let span = cli.ref_len - cli.read_len + 1;
    eprintln!(
        "generating reads: count={} read_len={} span={}",
        cli.read_count, cli.read_len, span
    );
    let starts = generate_read_starts(cli.read_count, span, cli.seed);
    write_reads_fasta(&reads_path, &ref_seq, &starts, cli.read_len)?;

    let meta = format!(
        "ref_len={}\nread_len={}\nread_count={}\nseed={}\nreads_format=fasta\n",
        cli.ref_len, cli.read_len, cli.read_count, cli.seed
    );
    fs::write(&meta_path, meta)
        .with_context(|| format!("failed to write {}", meta_path.display()))?;

    println!("ref_path={}", ref_path.display());
    println!("reads_path={}", reads_path.display());
    println!("meta_path={}", meta_path.display());
    Ok(())
}
