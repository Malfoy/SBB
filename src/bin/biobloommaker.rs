use anyhow::{Context, Result, bail};
use biobloom_rs::bloom::{
    BloomFilter, ConcurrentBloomFilter, DEFAULT_BLOCK_WORDS, bit_len_for_fpr_with_hash_count,
    optimal_bit_len, optimal_hash_count,
};
use biobloom_rs::fastx::{
    ProgressiveStats, count_kmers_in_file, insert_file_into_filter, progressive_insert_file,
};
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

    #[arg(long = "bit_len")]
    bit_len: Option<u64>,

    #[arg(long = "blocked", action = ArgAction::SetTrue)]
    blocked: bool,

    #[arg(long = "block_words", default_value_t = DEFAULT_BLOCK_WORDS)]
    block_words: u16,

    #[arg(short = 'r', long = "progressive")]
    progressive: Option<f64>,

    #[arg(short = 's', long = "subtract")]
    subtract: Option<PathBuf>,

    #[arg(short = 'e', long = "iterations", default_value_t = 1)]
    iterations: usize,

    #[arg(long = "seed_files", default_value_t = 1)]
    seed_files: usize,

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
    if let Some(bit_len) = cli.bit_len
        && bit_len == 0
    {
        bail!("bit_len must be > 0 when provided");
    }
    if cli.blocked && (cli.block_words == 0 || !cli.block_words.is_power_of_two()) {
        bail!("--block_words must be a non-zero power-of-two with --blocked");
    }
    if let Some(threshold) = cli.progressive {
        if !(0.0..=1.0).contains(&threshold) {
            bail!("--progressive threshold must be in [0, 1]");
        }
        if cli.seed_files == 0 || cli.seed_files >= cli.files.len() {
            bail!("--seed_files must be in [1, num_files-1] when --progressive is enabled");
        }
        if cli.iterations == 0 {
            bail!("--iterations must be > 0");
        }
    }
    if cli.iterations == 0 {
        bail!("--iterations must be > 0");
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
    let (hash_count, mut bit_len) = if let Some(k) = cli.hash_num {
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
    if let Some(bit_len_override) = cli.bit_len {
        bit_len = bit_len_override.max(64);
    }
    bit_len = bit_len
        .checked_next_power_of_two()
        .unwrap_or(bit_len)
        .max(64);
    if cli.blocked {
        let min_bits = u64::from(cli.block_words) * 64;
        if bit_len < min_bits {
            bit_len = min_bits;
        }
    }

    if cli.verbose {
        eprintln!(
            "building filter: k={} hashes={} bits={} expected_elements={} layout={} block_words={}",
            cli.kmer_size,
            hash_count,
            bit_len,
            expected_elements,
            if cli.blocked { "blocked" } else { "classic" },
            if cli.blocked { cli.block_words } else { 0 }
        );
    }

    let filter = ConcurrentBloomFilter::new(
        cli.kmer_size as u8,
        hash_count,
        bit_len,
        expected_elements,
        cli.false_positive_rate,
        cli.blocked,
        cli.block_words,
    )?;

    let insert_start = Instant::now();
    let mut inserted = 0_u64;
    let mut progressive_total_reads = 0_u64;
    let mut progressive_matched_reads = 0_u64;
    let mut progressive_queried_kmers = 0_u64;
    let mut progressive_passes = 0_u64;
    let mut subtract_loaded: Option<BloomFilter> = None;

    if let Some(progressive_threshold) = cli.progressive {
        let seed_files = &cli.files[..cli.seed_files];
        inserted += seed_files
            .par_iter()
            .map(|p| insert_file_into_filter(p, &filter))
            .try_reduce(|| 0_u64, |a, b| Ok(a + b))?;

        if let Some(path) = cli.subtract.as_ref() {
            subtract_loaded = Some(
                BloomFilter::load(path)
                    .with_context(|| format!("failed to load {}", path.display()))?,
            );
            if cli.verbose {
                eprintln!("loaded subtract filter: {}", path.display());
            }
        }

        let recruit_files = &cli.files[cli.seed_files..];
        for pass in 0..cli.iterations {
            let mut pass_stats = ProgressiveStats::default();
            for file in recruit_files {
                let stats = progressive_insert_file(
                    file,
                    &filter,
                    progressive_threshold,
                    subtract_loaded.as_ref(),
                )?;
                pass_stats.total_reads += stats.total_reads;
                pass_stats.matched_reads += stats.matched_reads;
                pass_stats.queried_kmers += stats.queried_kmers;
                pass_stats.inserted_kmers += stats.inserted_kmers;
            }
            progressive_passes += 1;
            progressive_total_reads += pass_stats.total_reads;
            progressive_matched_reads += pass_stats.matched_reads;
            progressive_queried_kmers += pass_stats.queried_kmers;
            inserted += pass_stats.inserted_kmers;
            if cli.verbose {
                eprintln!(
                    "progressive pass {}: total_reads={} matched_reads={} inserted_kmers={} queried_kmers={}",
                    pass + 1,
                    pass_stats.total_reads,
                    pass_stats.matched_reads,
                    pass_stats.inserted_kmers,
                    pass_stats.queried_kmers
                );
            }
            if pass_stats.matched_reads == 0 || pass_stats.inserted_kmers == 0 {
                break;
            }
        }
    } else {
        inserted += cli
            .files
            .par_iter()
            .map(|p| insert_file_into_filter(p, &filter))
            .try_reduce(|| 0_u64, |a, b| Ok(a + b))?;
    }
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
    info.push_str(&format!("bloom_layout\t{}\n", finalized.layout_name()));
    info.push_str(&format!("block_words\t{}\n", finalized.block_words));
    info.push_str(&format!(
        "false_positive_rate\t{:.8}\n",
        finalized.false_positive_rate
    ));
    info.push_str(&format!(
        "expected_elements\t{}\n",
        finalized.expected_elements
    ));
    if let Some(threshold) = cli.progressive {
        info.push_str(&format!("progressive_threshold\t{:.6}\n", threshold));
        info.push_str(&format!("progressive_seed_files\t{}\n", cli.seed_files));
        info.push_str(&format!("progressive_passes\t{}\n", progressive_passes));
        info.push_str(&format!(
            "progressive_total_reads\t{}\n",
            progressive_total_reads
        ));
        info.push_str(&format!(
            "progressive_matched_reads\t{}\n",
            progressive_matched_reads
        ));
        info.push_str(&format!(
            "progressive_queried_kmers\t{}\n",
            progressive_queried_kmers
        ));
        if let Some(path) = cli.subtract.as_ref() {
            info.push_str(&format!("progressive_subtract\t{}\n", path.display()));
        }
    }
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
    if let Some(threshold) = cli.progressive {
        println!("progressive_threshold\t{:.6}", threshold);
        println!("progressive_seed_files\t{}", cli.seed_files);
        println!("progressive_passes\t{}", progressive_passes);
        println!("progressive_total_reads\t{}", progressive_total_reads);
        println!("progressive_matched_reads\t{}", progressive_matched_reads);
        println!("progressive_queried_kmers\t{}", progressive_queried_kmers);
    }
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
