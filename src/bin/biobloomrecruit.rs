use anyhow::{Result, bail};
use biobloom_rs::classify::{ClassifyConfig, FilterDb, classify_seq};
use biobloom_rs::fastx::{ReadRecord, for_each_batch};
use biobloom_rs::writer::{OutputFormat, open_writer, write_record};
use clap::{ArgAction, Parser};
use rayon::prelude::*;
use std::path::{Path, PathBuf};

#[derive(Parser, Debug)]
#[command(
    name = "biobloomrecruit",
    about = "Recruit reads matching a single Bloom filter"
)]
struct Cli {
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
    if db.filters.len() != 1 {
        bail!("biobloomrecruit requires exactly one filter");
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

fn default_prefix(files: &[PathBuf]) -> String {
    files
        .first()
        .and_then(|p| p.file_stem())
        .and_then(|s| s.to_str())
        .unwrap_or("biobloom_recruit")
        .to_string()
}
