use anyhow::{Context, Result, bail};
use biobloom_rs::bloom::{
    ConcurrentBloomFilter, bit_len_for_fpr_with_hash_count, optimal_bit_len, optimal_hash_count,
};
use biobloom_rs::fastx::{count_kmers_in_file, insert_file_into_filter};
use clap::{ArgAction, Parser};
use rayon::prelude::*;
use std::fs;
use std::path::PathBuf;
use std::time::Instant;

const AUTO_MAX_HASHES: u8 = 3;

#[derive(Parser, Debug)]
#[command(
    name = "biobloommaker",
    about = "Build Bloom filters from FASTA/FASTQ (optionally gz)"
)]
struct Cli {
    #[arg(short = 'p', long = "file_prefix")]
    file_prefix: String,

    #[arg(short = 'o', long = "output_dir", default_value = ".")]
    output_dir: PathBuf,

    #[arg(short = 'f', long = "fal_pos_rate", default_value_t = 0.0075)]
    false_positive_rate: f64,

    #[arg(short = 'g', long = "hash_num")]
    hash_num: Option<u8>,

    #[arg(short = 'k', long = "kmer_size", default_value_t = 25)]
    kmer_size: usize,

    #[arg(short = 'n', long = "num_ele", default_value_t = 0)]
    num_ele: u64,

    #[arg(short = 't', long = "threads")]
    threads: Option<usize>,

    #[arg(short = 'v', long = "verbose", action = ArgAction::SetTrue)]
    verbose: bool,

    #[arg(required = true)]
    files: Vec<PathBuf>,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    if !(1..=32).contains(&cli.kmer_size) {
        bail!("k-mer size must be in [1, 32]");
    }
    if !(0.0..1.0).contains(&cli.false_positive_rate) {
        bail!("false positive rate must be in (0, 1)");
    }

    if let Some(t) = cli.threads {
        if t == 0 {
            bail!("threads must be > 0");
        }
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(t)
            .build_global();
    }

    fs::create_dir_all(&cli.output_dir)
        .with_context(|| format!("failed to create {}", cli.output_dir.display()))?;

    let expected_elements = if cli.num_ele > 0 {
        cli.num_ele
    } else {
        if cli.verbose {
            eprintln!("estimating element count from input files...");
        }
        cli.files
            .par_iter()
            .map(|p| count_kmers_in_file(p, cli.kmer_size))
            .try_reduce(|| 0_u64, |a, b| Ok(a + b))?
            .max(1)
    };

    let base_bit_len = optimal_bit_len(expected_elements, cli.false_positive_rate).max(64);
    let (hash_count, bit_len) = if let Some(k) = cli.hash_num {
        (k, base_bit_len)
    } else {
        let chosen_hashes =
            optimal_hash_count(base_bit_len, expected_elements).min(AUTO_MAX_HASHES);
        let tuned_bit_len = bit_len_for_fpr_with_hash_count(
            expected_elements,
            cli.false_positive_rate,
            chosen_hashes,
        )
        .max(64);
        (chosen_hashes, tuned_bit_len)
    };
    let bit_len = bit_len
        .checked_next_power_of_two()
        .unwrap_or(bit_len)
        .max(64);

    if cli.verbose {
        eprintln!(
            "building filter: k={} hashes={} bits={} expected_elements={}",
            cli.kmer_size, hash_count, bit_len, expected_elements
        );
    }

    let filter = ConcurrentBloomFilter::new(
        cli.kmer_size as u8,
        hash_count,
        bit_len,
        expected_elements,
        cli.false_positive_rate,
    )?;

    let insert_start = Instant::now();
    let inserted: u64 = cli
        .files
        .par_iter()
        .map(|p| insert_file_into_filter(p, &filter))
        .try_reduce(|| 0_u64, |a, b| Ok(a + b))?;
    let insert_elapsed = insert_start.elapsed();

    let finalized = filter.finalize();

    let bf_path = cli.output_dir.join(format!("{}.bf", cli.file_prefix));
    finalized.save(&bf_path)?;

    let info_path = cli.output_dir.join(format!("{}.txt", cli.file_prefix));
    let mut info = String::new();
    info.push_str(&format!("filter_id\t{}\n", cli.file_prefix));
    info.push_str(&format!("kmer_size\t{}\n", finalized.kmer_size));
    info.push_str(&format!("hash_count\t{}\n", finalized.hash_count));
    info.push_str(&format!("bit_len\t{}\n", finalized.bit_len));
    info.push_str(&format!(
        "false_positive_rate\t{:.8}\n",
        finalized.false_positive_rate
    ));
    info.push_str(&format!(
        "expected_elements\t{}\n",
        finalized.expected_elements
    ));
    info.push_str(&format!("inserted_kmers\t{}\n", inserted));
    info.push_str("inputs\t");
    for (i, p) in cli.files.iter().enumerate() {
        if i > 0 {
            info.push(';');
        }
        info.push_str(&p.display().to_string());
    }
    info.push('\n');
    fs::write(&info_path, info)
        .with_context(|| format!("failed to write {}", info_path.display()))?;

    println!("wrote filter: {}", bf_path.display());
    println!("wrote info:   {}", info_path.display());
    println!(
        "k={} hashes={} bits={} inserted_kmers={}",
        finalized.kmer_size, finalized.hash_count, finalized.bit_len, inserted
    );
    println!("insert_kmer_ops\t{}", inserted);
    if inserted > 0 {
        println!(
            "insert_mean_ns_per_kmer\t{:.3}",
            (insert_elapsed.as_nanos() as f64) / (inserted as f64)
        );
    } else {
        println!("insert_mean_ns_per_kmer\t0.000");
    }

    Ok(())
}
