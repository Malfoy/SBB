use anyhow::{Context, Result, bail};
use biobloom_rs::bloom::{
    BloomFilter, ConcurrentBloomFilter, DEFAULT_BLOCK_WORDS, bit_len_for_fpr_with_hash_count,
    optimal_bit_len, optimal_hash_count,
};
use biobloom_rs::classify::{
    ClassifyConfig, FilterDb, apply_inclusive_pair, classify_seq, classify_seq_with_metrics,
    match_seq_against_filter_with_counts,
};
use biobloom_rs::fastx::{
    ProgressiveStats, ReadRecord, count_kmers_in_file, for_each_batch, for_each_paired_batch,
    for_each_seq_batch, insert_file_into_filter, progressive_insert_file,
};
use biobloom_rs::writer::{
    CategorizerWriters, OutputFormat, format_score_suffix, open_writer, write_record,
};
use clap::{ArgAction, Args, Parser, Subcommand};
use rayon::prelude::*;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

const AUTO_MAX_HASHES: u8 = 3;

#[derive(Parser, Debug)]
#[command(
    name = "sbb",
    about = "SBB: Bloom-filter maker/categorizer/recruit CLI"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    Maker(MakerArgs),
    Categorizer(CategorizerArgs),
    Recruit(RecruitArgs),
}

#[derive(Args, Debug)]
#[command(about = "Build Bloom filters from FASTA/FASTQ (optionally compressed)")]
struct MakerArgs {
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

    #[arg(long = "classic", action = ArgAction::SetTrue)]
    classic: bool,

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

#[derive(Args, Debug)]
#[command(about = "Categorize FASTA/FASTQ reads using one or more Bloom filters")]
struct CategorizerArgs {
    #[arg(short = 'p', long = "prefix")]
    prefix: Option<String>,

    #[arg(short = 'f', long = "filter_files", required = true)]
    filter_files: Vec<String>,

    #[arg(short = 'e', long = "paired_mode", action = ArgAction::SetTrue)]
    paired_mode: bool,

    #[arg(short = 'i', long = "inclusive", action = ArgAction::SetTrue)]
    inclusive: bool,

    #[arg(short = 's', long = "score", default_value_t = 0.15)]
    score: f64,

    #[arg(short = 'b', long = "best_hit", action = ArgAction::SetTrue)]
    best_hit: bool,

    #[arg(short = 'w', long = "with_score", action = ArgAction::SetTrue)]
    with_score: bool,

    #[arg(short = 't', long = "threads")]
    threads: Option<usize>,

    #[arg(short = 'g', long = "gz_output", action = ArgAction::SetTrue)]
    gz_output: bool,

    #[arg(long = "fa", action = ArgAction::SetTrue)]
    out_fa: bool,

    #[arg(long = "fq", action = ArgAction::SetTrue)]
    out_fq: bool,

    #[arg(long = "verbose", action = ArgAction::SetTrue)]
    verbose: bool,

    #[arg(required = true)]
    files: Vec<PathBuf>,
}

#[derive(Args, Debug)]
#[command(about = "Recruit reads matching a single Bloom filter")]
struct RecruitArgs {
    #[arg(short = 'p', long = "prefix")]
    prefix: Option<String>,

    #[arg(short = 'f', long = "filter_files", required = true)]
    filter_files: Vec<String>,

    #[arg(short = 's', long = "score", default_value_t = 0.15)]
    score: f64,

    #[arg(short = 'b', long = "best_hit", action = ArgAction::SetTrue)]
    best_hit: bool,

    #[arg(short = 't', long = "threads")]
    threads: Option<usize>,

    #[arg(short = 'g', long = "gz_output", action = ArgAction::SetTrue)]
    gz_output: bool,

    #[arg(long = "fa", action = ArgAction::SetTrue)]
    out_fa: bool,

    #[arg(long = "fq", action = ArgAction::SetTrue)]
    out_fq: bool,

    #[arg(required = true)]
    files: Vec<PathBuf>,
}

#[derive(Default, Debug)]
struct Stats {
    total_reads: u64,
    nomatch_reads: u64,
    multimatch_reads: u64,
    per_filter_reads: Vec<u64>,
}

#[derive(Default, Debug)]
struct QueryMetrics {
    query_kmers: u64,
    query_nanos: u128,
}

impl Stats {
    fn new(filter_count: usize) -> Self {
        Self {
            per_filter_reads: vec![0; filter_count],
            ..Self::default()
        }
    }

    fn add_matches(&mut self, matches: &[usize]) {
        self.total_reads += 1;
        if matches.is_empty() {
            self.nomatch_reads += 1;
        }
        if matches.len() > 1 {
            self.multimatch_reads += 1;
        }
        for &idx in matches {
            self.per_filter_reads[idx] += 1;
        }
    }
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Commands::Maker(args) => run_maker(args),
        Commands::Categorizer(args) => run_categorizer(args),
        Commands::Recruit(args) => run_recruit(args),
    }
}

fn configure_threads(threads: Option<usize>) -> Result<()> {
    if let Some(t) = threads {
        if t == 0 {
            bail!("threads must be > 0");
        }
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(t)
            .build_global();
    }
    Ok(())
}

fn run_maker(cli: MakerArgs) -> Result<()> {
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
    if cli.blocked && cli.classic {
        bail!("choose only one layout override: --blocked or --classic");
    }

    let blocked_layout = !cli.classic;
    if blocked_layout && (cli.block_words == 0 || !cli.block_words.is_power_of_two()) {
        bail!("--block_words must be a non-zero power-of-two in blocked mode");
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

    configure_threads(cli.threads)?;

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
    if blocked_layout {
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
            if blocked_layout { "blocked" } else { "classic" },
            if blocked_layout { cli.block_words } else { 0 }
        );
    }

    let filter = ConcurrentBloomFilter::new(
        cli.kmer_size as u8,
        hash_count,
        bit_len,
        expected_elements,
        cli.false_positive_rate,
        blocked_layout,
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

fn run_categorizer(cli: CategorizerArgs) -> Result<()> {
    if cli.out_fa && cli.out_fq {
        bail!("choose only one output format: --fa or --fq");
    }
    if !cli.best_hit && !(0.0..=1.0).contains(&cli.score) {
        bail!("--score must be in [0,1] unless --best_hit is used");
    }

    configure_threads(cli.threads)?;

    let db = FilterDb::load(&cli.filter_files)?;
    let filter_ids: Vec<String> = db.filters.iter().map(|f| f.id.clone()).collect();
    let prefix = cli.prefix.unwrap_or_else(|| default_prefix(&cli.files));

    if cli.verbose {
        eprintln!(
            "loaded {} filters: {}",
            db.filters.len(),
            filter_ids.join(", ")
        );
    }

    let mut writers = match (cli.out_fa, cli.out_fq) {
        (true, false) => Some(CategorizerWriters::new(
            &prefix,
            filter_ids.clone(),
            OutputFormat::Fasta,
            cli.gz_output,
        )?),
        (false, true) => Some(CategorizerWriters::new(
            &prefix,
            filter_ids.clone(),
            OutputFormat::Fastq,
            cli.gz_output,
        )?),
        _ => None,
    };

    let cfg = ClassifyConfig {
        threshold: cli.score,
        best_hit: cli.best_hit,
    };

    let mut stats = Stats::new(db.filters.len());
    let mut metrics = QueryMetrics::default();
    let stats_only_single_filter =
        !cli.with_score && writers.is_none() && !cli.paired_mode && db.filters.len() == 1;

    if cli.paired_mode {
        if cli.files.len() != 2 {
            bail!("paired mode currently requires exactly two input files");
        }
        process_paired(
            &cli.files[0],
            &cli.files[1],
            &db,
            cfg,
            cli.inclusive,
            cli.with_score,
            writers.as_mut(),
            &mut stats,
            &mut metrics,
        )?;
    } else {
        for file in &cli.files {
            if stats_only_single_filter {
                process_single_stats_only_single_filter(file, &db, cfg, &mut stats, &mut metrics)?;
            } else {
                process_single(
                    file,
                    &db,
                    cfg,
                    cli.with_score,
                    writers.as_mut(),
                    &mut stats,
                    &mut metrics,
                )?;
            }
        }
    }

    if let Some(w) = writers.as_mut() {
        w.flush_all()?;
    }

    println!("total_reads\t{}", stats.total_reads);
    println!("nomatch_reads\t{}", stats.nomatch_reads);
    println!("multimatch_reads\t{}", stats.multimatch_reads);
    for (id, count) in filter_ids.iter().zip(stats.per_filter_reads.iter()) {
        println!("filter_reads\t{}\t{}", id, count);
    }
    println!("query_kmer_queries\t{}", metrics.query_kmers);
    if metrics.query_kmers > 0 {
        println!(
            "query_mean_ns_per_kmer\t{:.3}",
            (metrics.query_nanos as f64) / (metrics.query_kmers as f64)
        );
    } else {
        println!("query_mean_ns_per_kmer\t0.000");
    }

    Ok(())
}

fn run_recruit(cli: RecruitArgs) -> Result<()> {
    if cli.out_fa && cli.out_fq {
        bail!("choose only one output format: --fa or --fq");
    }
    if !cli.best_hit && !(0.0..=1.0).contains(&cli.score) {
        bail!("--score must be in [0,1] unless --best_hit is used");
    }

    configure_threads(cli.threads)?;

    let db = FilterDb::load(&cli.filter_files)?;
    if db.filters.len() != 1 {
        bail!("sbb recruit requires exactly one filter");
    }

    let prefix = cli.prefix.unwrap_or_else(|| default_prefix(&cli.files));
    let format = if cli.out_fa {
        OutputFormat::Fasta
    } else {
        OutputFormat::Fastq
    };
    let ext = match format {
        OutputFormat::Fasta => "fa",
        OutputFormat::Fastq => "fq",
    };
    let gz_ext = if cli.gz_output { ".gz" } else { "" };

    let matched_path = format!("{prefix}_recruited.{ext}{gz_ext}");
    let unmatched_path = format!("{prefix}_unmatched.{ext}{gz_ext}");

    let mut matched = open_writer(Path::new(&matched_path), cli.gz_output)?;
    let mut unmatched = open_writer(Path::new(&unmatched_path), cli.gz_output)?;

    let cfg = ClassifyConfig {
        threshold: cli.score,
        best_hit: cli.best_hit,
    };

    let mut total = 0_u64;
    let mut recruited = 0_u64;

    for file in &cli.files {
        for_each_batch(file, 2048, |batch: Vec<ReadRecord>| {
            let classes: Vec<_> = batch
                .par_iter()
                .map(|r| classify_seq(&r.seq, &db, cfg))
                .collect();

            for (r, c) in batch.iter().zip(classes.iter()) {
                total += 1;
                if c.matches.is_empty() {
                    write_record(unmatched.as_mut(), r, format, None)?;
                } else {
                    recruited += 1;
                    write_record(matched.as_mut(), r, format, None)?;
                }
            }
            Ok(())
        })?;
    }

    matched.flush()?;
    unmatched.flush()?;

    println!("total_reads\t{}", total);
    println!("recruited_reads\t{}", recruited);
    println!("unmatched_reads\t{}", total - recruited);
    println!("output_recruited\t{}", matched_path);
    println!("output_unmatched\t{}", unmatched_path);

    Ok(())
}

fn process_single(
    input: &PathBuf,
    db: &FilterDb,
    cfg: ClassifyConfig,
    with_score: bool,
    mut writers: Option<&mut CategorizerWriters>,
    stats: &mut Stats,
    metrics: &mut QueryMetrics,
) -> Result<()> {
    const BATCH_SIZE: usize = 2048;

    for_each_batch(input, BATCH_SIZE, |batch: Vec<ReadRecord>| {
        let query_start = Instant::now();
        let results: Vec<_> = batch
            .par_iter()
            .map(|rec| classify_seq_with_metrics(&rec.seq, db, cfg))
            .collect();
        metrics.query_nanos += query_start.elapsed().as_nanos();

        for (rec, (classification, queried_kmers)) in batch.iter().zip(results.iter()) {
            metrics.query_kmers += *queried_kmers;
            stats.add_matches(&classification.matches);

            if let Some(w) = writers.as_deref_mut() {
                let suffix = if with_score {
                    Some(format_score_suffix(&w.filter_ids, &classification.scores))
                } else {
                    None
                };
                w.write_assignment(rec, &classification.matches, suffix.as_deref())?;
            }
        }

        Ok(())
    })
}

fn process_single_stats_only_single_filter(
    input: &PathBuf,
    db: &FilterDb,
    cfg: ClassifyConfig,
    stats: &mut Stats,
    metrics: &mut QueryMetrics,
) -> Result<()> {
    const BATCH_SIZE: usize = 2048;
    let bloom = &db.filters[0].bloom;

    for_each_seq_batch(input, BATCH_SIZE, |batch: Vec<Vec<u8>>| {
        let query_start = Instant::now();
        let results: Vec<_> = batch
            .par_iter()
            .map(|seq| {
                let (matched, queried_kmers) =
                    match_seq_against_filter_with_counts(seq, bloom, cfg);
                (u64::from(matched), queried_kmers)
            })
            .collect();
        metrics.query_nanos += query_start.elapsed().as_nanos();

        let mut matched_reads = 0_u64;
        for (matched, queried_kmers) in results {
            matched_reads += matched;
            metrics.query_kmers += queried_kmers;
        }

        let batch_total = batch.len() as u64;
        stats.total_reads += batch_total;
        stats.nomatch_reads += batch_total - matched_reads;
        stats.per_filter_reads[0] += matched_reads;
        Ok(())
    })
}

fn process_paired(
    r1: &PathBuf,
    r2: &PathBuf,
    db: &FilterDb,
    cfg: ClassifyConfig,
    inclusive: bool,
    with_score: bool,
    mut writers: Option<&mut CategorizerWriters>,
    stats: &mut Stats,
    metrics: &mut QueryMetrics,
) -> Result<()> {
    const BATCH_SIZE: usize = 1024;

    for_each_paired_batch(
        r1,
        r2,
        BATCH_SIZE,
        |pairs: Vec<(ReadRecord, ReadRecord)>| {
            let query_start = Instant::now();
            let results: Vec<_> = pairs
                .par_iter()
                .map(|(a, b)| {
                    let (ca, qa) = classify_seq_with_metrics(&a.seq, db, cfg);
                    let (cb, qb) = classify_seq_with_metrics(&b.seq, db, cfg);
                    (ca, cb, qa + qb)
                })
                .collect();
            metrics.query_nanos += query_start.elapsed().as_nanos();

            for ((a, b), (ca, cb, queried_kmers)) in pairs.iter().zip(results.iter()) {
                metrics.query_kmers += *queried_kmers;

                let (ma, mb) = apply_inclusive_pair(ca, cb, inclusive);
                stats.add_matches(&ma);
                stats.add_matches(&mb);

                if let Some(w) = writers.as_deref_mut() {
                    let sa = if with_score {
                        Some(format_score_suffix(&w.filter_ids, &ca.scores))
                    } else {
                        None
                    };
                    let sb = if with_score {
                        Some(format_score_suffix(&w.filter_ids, &cb.scores))
                    } else {
                        None
                    };
                    w.write_assignment(a, &ma, sa.as_deref())?;
                    w.write_assignment(b, &mb, sb.as_deref())?;
                }
            }

            Ok(())
        },
    )
}

fn default_prefix(files: &[PathBuf]) -> String {
    files
        .first()
        .and_then(|p| p.file_stem())
        .and_then(|s| s.to_str())
        .unwrap_or("sbb")
        .to_string()
}
