use crate::bloom::{ConcurrentBloomFilter, for_each_canonical_kmer};
use anyhow::{Context, Result, bail};
use needletail::parse_fastx_file;
use std::path::Path;

#[derive(Clone, Debug)]
pub struct ReadRecord {
    pub id: Vec<u8>,
    pub seq: Vec<u8>,
    pub qual: Option<Vec<u8>>,
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
        inserted_kmers += for_each_canonical_kmer(record.seq().as_ref(), k, |canonical| {
            filter.insert_kmer(canonical)
        }) as u64;
    }

    Ok(inserted_kmers)
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
