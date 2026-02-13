use anyhow::{Result, bail};
use biobloom_rs::classify::{
    ClassifyConfig, FilterDb, apply_inclusive_pair, classify_seq_with_metrics,
    match_seq_against_filter_with_counts,
};
use biobloom_rs::fastx::{ReadRecord, for_each_batch, for_each_paired_batch, for_each_seq_batch};
use biobloom_rs::writer::{CategorizerWriters, OutputFormat, format_score_suffix};
use clap::{ArgAction, Parser};
use rayon::prelude::*;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Parser, Debug)]
#[command(
    name = "biobloomcategorizer",
    about = "Categorize FASTA/FASTQ reads using one or more Bloom filters"
)]
struct Cli {
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

    if cli.out_fa && cli.out_fq {
        bail!("choose only one output format: --fa or --fq");
    }
    if !cli.best_hit && !(0.0..=1.0).contains(&cli.score) {
        bail!("--score must be in [0,1] unless --best_hit is used");
    }

    if let Some(t) = cli.threads {
        if t == 0 {
            bail!("threads must be > 0");
        }
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(t)
            .build_global();
    }

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
    input1: &PathBuf,
    input2: &PathBuf,
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
        input1,
        input2,
        BATCH_SIZE,
        |batch: Vec<(ReadRecord, ReadRecord)>| {
            let query_start = Instant::now();
            let results: Vec<_> = batch
                .par_iter()
                .map(|(r1, r2)| {
                    let (c1, q1) = classify_seq_with_metrics(&r1.seq, db, cfg);
                    let (c2, q2) = classify_seq_with_metrics(&r2.seq, db, cfg);
                    let (m1, m2) = apply_inclusive_pair(&c1, &c2, inclusive);
                    (c1, c2, m1, m2, q1, q2)
                })
                .collect();
            metrics.query_nanos += query_start.elapsed().as_nanos();

            for ((r1, r2), (c1, c2, m1, m2, q1, q2)) in batch.iter().zip(results.iter()) {
                metrics.query_kmers += *q1 + *q2;
                stats.add_matches(m1);
                stats.add_matches(m2);

                if let Some(w) = writers.as_deref_mut() {
                    let s1 = if with_score {
                        Some(format_score_suffix(&w.filter_ids, &c1.scores))
                    } else {
                        None
                    };
                    let s2 = if with_score {
                        Some(format_score_suffix(&w.filter_ids, &c2.scores))
                    } else {
                        None
                    };

                    w.write_assignment(r1, m1, s1.as_deref())?;
                    w.write_assignment(r2, m2, s2.as_deref())?;
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
        .unwrap_or("biobloom")
        .to_string()
}
