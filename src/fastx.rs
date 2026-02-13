use crate::bloom::{BloomFilter, ConcurrentBloomFilter, for_each_canonical_kmer};
use anyhow::{Context, Result, bail};
use needletail::parse_fastx_file;
use rayon::prelude::*;
use std::path::Path;

#[derive(Clone, Debug)]
pub struct ReadRecord {
    pub id: Vec<u8>,
    pub seq: Vec<u8>,
    pub qual: Option<Vec<u8>>,
}

#[derive(Clone, Copy, Debug, Default)]
pub struct ProgressiveStats {
    pub total_reads: u64,
    pub matched_reads: u64,
    pub queried_kmers: u64,
    pub inserted_kmers: u64,
}

pub fn count_kmers_in_file(path: &Path, kmer_size: usize) -> Result<u64> {
    let mut reader = parse_fastx_file(path)
        .with_context(|| format!("failed to open input file {}", path.display()))?;

    let mut total = 0_u64;
    while let Some(record) = reader.next() {
        let record = record.with_context(|| format!("failed to parse {}", path.display()))?;
        total += for_each_canonical_kmer(record.seq().as_ref(), kmer_size, |_| {}) as u64;
    }

    Ok(total)
}

pub fn insert_file_into_filter(path: &Path, filter: &ConcurrentBloomFilter) -> Result<u64> {
    let mut reader = parse_fastx_file(path)
        .with_context(|| format!("failed to open input file {}", path.display()))?;

    let k = filter.kmer_size as usize;
    let mut inserted_kmers = 0_u64;

    while let Some(record) = reader.next() {
        let record = record.with_context(|| format!("failed to parse {}", path.display()))?;
        inserted_kmers += insert_seq_into_filter(record.seq().as_ref(), filter, k);
    }

    Ok(inserted_kmers)
}

pub fn progressive_insert_file(
    path: &Path,
    filter: &ConcurrentBloomFilter,
    threshold: f64,
    subtract: Option<&BloomFilter>,
) -> Result<ProgressiveStats> {
    let mut reader = parse_fastx_file(path)
        .with_context(|| format!("failed to open input file {}", path.display()))?;
    let k = filter.kmer_size as usize;
    let mut stats = ProgressiveStats::default();

    while let Some(record) = reader.next() {
        let record = record.with_context(|| format!("failed to parse {}", path.display()))?;
        let seq = record.seq().as_ref().to_vec();
        stats.total_reads += 1;

        let mut total = 0_u64;
        let mut hits = 0_u64;
        for_each_canonical_kmer(&seq, k, |canonical| {
            total += 1;
            if filter.contains_kmer(canonical) {
                hits += 1;
            }
        });
        stats.queried_kmers += total;

        if total == 0 {
            continue;
        }

        let score = hits as f64 / total as f64;
        if score < threshold {
            continue;
        }

        stats.matched_reads += 1;
        for_each_canonical_kmer(&seq, k, |canonical| {
            if subtract.is_some_and(|sub| sub.contains_kmer(canonical)) {
                return;
            }
            filter.insert_kmer(canonical);
            stats.inserted_kmers += 1;
        });
    }

    Ok(stats)
}

#[inline]
fn nuc2bits(b: u8) -> Option<u8> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' | b'U' | b'u' => Some(3),
        _ => None,
    }
}

fn insert_seq_into_filter(seq: &[u8], filter: &ConcurrentBloomFilter, k: usize) -> u64 {
    const PAR_MIN_STARTS: usize = 1_000_000;
    const CHUNK_STARTS: usize = 1_000_000;

    if seq.len() < k {
        return 0;
    }
    let starts = seq.len() - k + 1;
    if starts < PAR_MIN_STARTS || rayon::current_num_threads() <= 1 {
        return for_each_canonical_kmer(seq, k, |canonical| filter.insert_kmer(canonical)) as u64;
    }

    let chunks = starts.div_ceil(CHUNK_STARTS);
    (0..chunks)
        .into_par_iter()
        .map(|chunk_idx| {
            let start_pos = chunk_idx * CHUNK_STARTS;
            let end_pos = (start_pos + CHUNK_STARTS).min(starts);
            let slice_end = end_pos + k - 1;
            let slice = &seq[start_pos..slice_end];

            let mut forward = 0_u64;
            let mut reverse = 0_u64;
            let mask = if k == 32 {
                u64::MAX
            } else {
                (1_u64 << (2 * k)) - 1
            };
            let rev_shift = 2 * (k - 1);
            let mut run = 0_usize;
            let mut count = 0_u64;

            for &b in slice {
                match nuc2bits(b) {
                    Some(code) => {
                        forward = ((forward << 2) | code as u64) & mask;
                        let comp = (3_u8 - code) as u64;
                        reverse = (reverse >> 2) | (comp << rev_shift);

                        run += 1;
                        if run >= k {
                            count += 1;
                            filter.insert_kmer(forward.min(reverse));
                        }
                    }
                    None => {
                        run = 0;
                        forward = 0;
                        reverse = 0;
                    }
                }
            }
            count
        })
        .sum()
}

pub fn for_each_batch<F>(path: &Path, batch_size: usize, mut f: F) -> Result<()>
where
    F: FnMut(Vec<ReadRecord>) -> Result<()>,
{
    let mut reader = parse_fastx_file(path)
        .with_context(|| format!("failed to open input file {}", path.display()))?;

    let mut batch = Vec::with_capacity(batch_size);

    while let Some(record) = reader.next() {
        let record = record.with_context(|| format!("failed to parse {}", path.display()))?;
        batch.push(ReadRecord {
            id: record.id().to_vec(),
            seq: record.seq().as_ref().to_vec(),
            qual: record.qual().map(|q| q.to_vec()),
        });

        if batch.len() >= batch_size {
            f(std::mem::take(&mut batch))?;
        }
    }

    if !batch.is_empty() {
        f(batch)?;
    }

    Ok(())
}

pub fn for_each_seq_batch<F>(path: &Path, batch_size: usize, mut f: F) -> Result<()>
where
    F: FnMut(Vec<Vec<u8>>) -> Result<()>,
{
    let mut reader = parse_fastx_file(path)
        .with_context(|| format!("failed to open input file {}", path.display()))?;

    let mut batch = Vec::with_capacity(batch_size);

    while let Some(record) = reader.next() {
        let record = record.with_context(|| format!("failed to parse {}", path.display()))?;
        batch.push(record.seq().as_ref().to_vec());

        if batch.len() >= batch_size {
            f(std::mem::take(&mut batch))?;
        }
    }

    if !batch.is_empty() {
        f(batch)?;
    }

    Ok(())
}

pub fn for_each_paired_batch<F>(
    path1: &Path,
    path2: &Path,
    batch_size: usize,
    mut f: F,
) -> Result<()>
where
    F: FnMut(Vec<(ReadRecord, ReadRecord)>) -> Result<()>,
{
    let mut r1 = parse_fastx_file(path1)
        .with_context(|| format!("failed to open input file {}", path1.display()))?;
    let mut r2 = parse_fastx_file(path2)
        .with_context(|| format!("failed to open input file {}", path2.display()))?;

    let mut batch = Vec::with_capacity(batch_size);

    loop {
        let a = r1.next();
        let b = r2.next();

        match (a, b) {
            (None, None) => break,
            (Some(_), None) | (None, Some(_)) => {
                bail!(
                    "paired files have different record counts: {} vs {}",
                    path1.display(),
                    path2.display()
                );
            }
            (Some(ra), Some(rb)) => {
                let ra = ra.with_context(|| format!("failed to parse {}", path1.display()))?;
                let rb = rb.with_context(|| format!("failed to parse {}", path2.display()))?;

                batch.push((
                    ReadRecord {
                        id: ra.id().to_vec(),
                        seq: ra.seq().as_ref().to_vec(),
                        qual: ra.qual().map(|q| q.to_vec()),
                    },
                    ReadRecord {
                        id: rb.id().to_vec(),
                        seq: rb.seq().as_ref().to_vec(),
                        qual: rb.qual().map(|q| q.to_vec()),
                    },
                ));

                if batch.len() >= batch_size {
                    f(std::mem::take(&mut batch))?;
                }
            }
        }
    }

    if !batch.is_empty() {
        f(batch)?;
    }

    Ok(())
}
