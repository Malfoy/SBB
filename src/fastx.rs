use crate::bloom::{BloomFilter, BloomLayout, ConcurrentBloomFilter, for_each_canonical_kmer};
use anyhow::{Context, Result, bail};
use helicase::input::FromFile;
use helicase::{Config, Event, FastxParser, Parser, ParserOptions};
use rayon::prelude::*;
use std::path::Path;

const FASTX_CONFIG: Config = ParserOptions::default().compute_quality().config();
type FastxReader = FastxParser<'static, FASTX_CONFIG>;

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

fn open_fastx_reader(path: &Path) -> Result<FastxReader> {
    FastxParser::<FASTX_CONFIG>::from_file(path)
        .with_context(|| format!("failed to open input file {}", path.display()))
}

fn next_record(reader: &mut FastxReader) -> Option<ReadRecord> {
    loop {
        match reader.next()? {
            Event::Record(_) => {
                return Some(ReadRecord {
                    id: reader.get_header().to_vec(),
                    seq: reader.get_dna_string().to_vec(),
                    qual: reader.get_quality().map(|q| q.to_vec()),
                });
            }
            Event::DnaChunk(_) => {}
        }
    }
}

pub fn count_kmers_in_file(path: &Path, kmer_size: usize) -> Result<u64> {
    let mut reader = open_fastx_reader(path)?;

    let mut total = 0_u64;
    while let Some(event) = reader.next() {
        if !matches!(event, Event::Record(_)) {
            continue;
        }
        total += count_valid_kmers(reader.get_dna_string(), kmer_size);
    }

    Ok(total)
}

pub fn insert_file_into_filter(path: &Path, filter: &ConcurrentBloomFilter) -> Result<u64> {
    let mut reader = open_fastx_reader(path)?;

    let k = filter.kmer_size as usize;
    let mut inserted_kmers = 0_u64;

    while let Some(event) = reader.next() {
        if !matches!(event, Event::Record(_)) {
            continue;
        }
        let seq = reader.get_dna_string();
        inserted_kmers += if filter.layout == BloomLayout::BloomyBloom {
            filter.insert_sequence(seq)
        } else {
            insert_seq_into_filter(seq, filter, k)
        };
    }

    Ok(inserted_kmers)
}

pub fn progressive_insert_file(
    path: &Path,
    filter: &ConcurrentBloomFilter,
    threshold: f64,
    subtract: Option<&BloomFilter>,
) -> Result<ProgressiveStats> {
    if filter.layout == BloomLayout::BloomyBloom && subtract.is_some() {
        bail!("progressive subtract is not supported for bloomybloom filters");
    }

    let mut reader = open_fastx_reader(path)?;
    let k = filter.kmer_size as usize;
    let mut stats = ProgressiveStats::default();

    while let Some(event) = reader.next() {
        if !matches!(event, Event::Record(_)) {
            continue;
        }
        let seq = reader.get_dna_string();
        stats.total_reads += 1;

        let (hits, total) = filter.score_sequence_with_counts(seq);
        stats.queried_kmers += total;

        if total == 0 {
            continue;
        }

        let score = hits as f64 / total as f64;
        if score < threshold {
            continue;
        }

        stats.matched_reads += 1;
        if filter.layout == BloomLayout::BloomyBloom {
            stats.inserted_kmers += filter.insert_sequence(seq);
            continue;
        }

        for_each_canonical_kmer(seq, k, |canonical| {
            if !filter.should_keep_kmer(canonical) {
                return;
            }
            if subtract.is_some_and(|sub| sub.contains_kmer(canonical)) {
                return;
            }
            if filter.insert_kmer(canonical) {
                stats.inserted_kmers += 1;
            }
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
        let mut inserted = 0_u64;
        for_each_canonical_kmer(seq, k, |canonical| {
            if filter.insert_kmer(canonical) {
                inserted += 1;
            }
        });
        return inserted;
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
                            if filter.insert_kmer(forward.min(reverse)) {
                                count += 1;
                            }
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

fn count_valid_kmers(seq: &[u8], k: usize) -> u64 {
    if k == 0 {
        return 0;
    }
    if k <= 32 {
        return for_each_canonical_kmer(seq, k, |_| {}) as u64;
    }

    let mut total = 0_u64;
    let mut run = 0_usize;
    for &b in seq {
        if nuc2bits(b).is_some() {
            run += 1;
            if run >= k {
                total += 1;
            }
        } else {
            run = 0;
        }
    }
    total
}

pub fn for_each_batch<F>(path: &Path, batch_size: usize, mut f: F) -> Result<()>
where
    F: FnMut(Vec<ReadRecord>) -> Result<()>,
{
    let mut reader = open_fastx_reader(path)?;

    let mut batch = Vec::with_capacity(batch_size);

    while let Some(record) = next_record(&mut reader) {
        batch.push(record);

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
    let mut reader = open_fastx_reader(path)?;

    let mut batch = Vec::with_capacity(batch_size);

    while let Some(event) = reader.next() {
        if !matches!(event, Event::Record(_)) {
            continue;
        }
        batch.push(reader.get_dna_string().to_vec());

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
    let mut r1 = open_fastx_reader(path1)?;
    let mut r2 = open_fastx_reader(path2)?;

    let mut batch = Vec::with_capacity(batch_size);

    loop {
        let a = next_record(&mut r1);
        let b = next_record(&mut r2);

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
                batch.push((ra, rb));

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
