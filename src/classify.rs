use crate::bloom::{BloomFilter, for_each_canonical_kmer};
use anyhow::{Context, Result, bail};
use std::path::{Path, PathBuf};

#[derive(Debug)]
pub struct NamedFilter {
    pub id: String,
    pub path: PathBuf,
    pub bloom: BloomFilter,
}

#[derive(Debug)]
pub struct FilterDb {
    pub filters: Vec<NamedFilter>,
}

#[derive(Debug, Clone, Copy)]
pub struct ClassifyConfig {
    pub threshold: f64,
    pub best_hit: bool,
}

#[derive(Debug, Clone)]
pub struct Classification {
    pub scores: Vec<f64>,
    pub matches: Vec<usize>,
}

impl FilterDb {
    pub fn load(filter_specs: &[String]) -> Result<Self> {
        let mut paths = Vec::new();
        for spec in filter_specs {
            for token in spec.split_whitespace() {
                if !token.is_empty() {
                    paths.push(PathBuf::from(token));
                }
            }
        }

        if paths.is_empty() {
            bail!("no filter files were provided");
        }

        let mut filters = Vec::with_capacity(paths.len());
        for path in paths {
            let bloom = BloomFilter::load(&path)
                .with_context(|| format!("failed to load filter {}", path.display()))?;
            let id = filter_id_from_path(&path);
            filters.push(NamedFilter { id, path, bloom });
        }

        Ok(Self { filters })
    }
}

pub fn filter_id_from_path(path: &Path) -> String {
    path.file_stem()
        .and_then(|x| x.to_str())
        .unwrap_or("filter")
        .to_string()
}

pub fn classify_seq(seq: &[u8], db: &FilterDb, cfg: ClassifyConfig) -> Classification {
    let (classification, _) = classify_seq_with_metrics(seq, db, cfg);
    classification
}

pub fn classify_seq_with_metrics(
    seq: &[u8],
    db: &FilterDb,
    cfg: ClassifyConfig,
) -> (Classification, u64) {
    let mut scores = Vec::with_capacity(db.filters.len());
    let mut queried_kmers = 0_u64;

    for f in &db.filters {
        let (score, filter_kmers) = score_seq_against_filter_with_counts(seq, &f.bloom);
        scores.push(score);
        queried_kmers += filter_kmers;
    }

    let matches = if cfg.best_hit {
        let mut best_idx = None;
        let mut best_score = f64::MIN;
        for (i, score) in scores.iter().copied().enumerate() {
            if score > best_score {
                best_idx = Some(i);
                best_score = score;
            }
        }
        if let Some(i) = best_idx {
            if best_score > 0.0 {
                vec![i]
            } else {
                Vec::new()
            }
        } else {
            Vec::new()
        }
    } else {
        scores
            .iter()
            .enumerate()
            .filter_map(|(i, score)| {
                if *score >= cfg.threshold {
                    Some(i)
                } else {
                    None
                }
            })
            .collect()
    };

    (Classification { scores, matches }, queried_kmers)
}

pub fn score_seq_against_filter(seq: &[u8], bloom: &BloomFilter) -> f64 {
    let (score, _) = score_seq_against_filter_with_counts(seq, bloom);
    score
}

pub fn score_seq_against_filter_with_counts(seq: &[u8], bloom: &BloomFilter) -> (f64, u64) {
    let mut total = 0_u64;
    let mut hits = 0_u64;
    let k = bloom.kmer_size as usize;

    for_each_canonical_kmer(seq, k, |canonical| {
        total += 1;
        if bloom.contains_kmer(canonical) {
            hits += 1;
        }
    });

    if total == 0 {
        (0.0, 0)
    } else {
        (hits as f64 / total as f64, total)
    }
}

pub fn apply_inclusive_pair(
    a: &Classification,
    b: &Classification,
    inclusive: bool,
) -> (Vec<usize>, Vec<usize>) {
    if !inclusive {
        return (a.matches.clone(), b.matches.clone());
    }

    let mut union = a.matches.clone();
    for &idx in &b.matches {
        if !union.contains(&idx) {
            union.push(idx);
        }
    }

    (union.clone(), union)
}
