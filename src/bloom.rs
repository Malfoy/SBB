use anyhow::{Context, Result, bail};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};

const MAGIC: [u8; 8] = *b"BBRFILT1";
const HEADER_SIZE: usize = 8 + 1 + 1 + 2 + 8 + 8 + 8 + 8;

#[derive(Clone, Debug)]
pub struct BloomFilter {
    pub kmer_size: u8,
    pub hash_count: u8,
    pub bit_len: u64,
    pub expected_elements: u64,
    pub false_positive_rate: f64,
    bit_mask: u64,
    words: Vec<u64>,
}

#[derive(Debug)]
pub struct ConcurrentBloomFilter {
    pub kmer_size: u8,
    pub hash_count: u8,
    pub bit_len: u64,
    pub expected_elements: u64,
    pub false_positive_rate: f64,
    bit_mask: u64,
    words: Vec<AtomicU64>,
}

impl ConcurrentBloomFilter {
    pub fn new(
        kmer_size: u8,
        hash_count: u8,
        bit_len: u64,
        expected_elements: u64,
        false_positive_rate: f64,
    ) -> Result<Self> {
        if bit_len == 0 {
            bail!("bit_len must be > 0");
        }
        if hash_count == 0 {
            bail!("hash_count must be > 0");
        }
        let word_len = bit_len.div_ceil(64) as usize;
        let mut words = Vec::with_capacity(word_len);
        words.resize_with(word_len, || AtomicU64::new(0));
        let bit_mask = if bit_len.is_power_of_two() {
            bit_len - 1
        } else {
            0
        };

        Ok(Self {
            kmer_size,
            hash_count,
            bit_len,
            expected_elements,
            false_positive_rate,
            bit_mask,
            words,
        })
    }

    #[inline]
    pub fn insert_kmer(&self, canonical_kmer: u64) {
        let (h1, h2) = hash_pair(canonical_kmer);
        let mut h = h1;
        for _ in 0..self.hash_count {
            let idx = if self.bit_mask != 0 {
                (h & self.bit_mask) as usize
            } else {
                (h % self.bit_len) as usize
            };
            let word_idx = idx >> 6;
            let bit = idx & 63;
            let mask = 1_u64 << bit;
            self.words[word_idx].fetch_or(mask, Ordering::Relaxed);
            h = h.wrapping_add(h2);
        }
    }

    pub fn finalize(self) -> BloomFilter {
        let words = self
            .words
            .into_iter()
            .map(|x| x.load(Ordering::Relaxed))
            .collect();

        BloomFilter {
            kmer_size: self.kmer_size,
            hash_count: self.hash_count,
            bit_len: self.bit_len,
            expected_elements: self.expected_elements,
            false_positive_rate: self.false_positive_rate,
            bit_mask: self.bit_mask,
            words,
        }
    }
}

impl BloomFilter {
    #[inline]
    pub fn contains_kmer(&self, canonical_kmer: u64) -> bool {
        let (h1, h2) = hash_pair(canonical_kmer);
        let mut h = h1;
        for _ in 0..self.hash_count {
            let idx = if self.bit_mask != 0 {
                (h & self.bit_mask) as usize
            } else {
                (h % self.bit_len) as usize
            };
            let word_idx = idx >> 6;
            let bit = idx & 63;
            let mask = 1_u64 << bit;
            if self.words[word_idx] & mask == 0 {
                return false;
            }
            h = h.wrapping_add(h2);
        }
        true
    }

    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let file = File::create(path.as_ref())
            .with_context(|| format!("failed to create {}", path.as_ref().display()))?;
        let mut w = BufWriter::new(file);

        w.write_all(&MAGIC)?;
        w.write_all(&[self.kmer_size])?;
        w.write_all(&[self.hash_count])?;
        w.write_all(&(HEADER_SIZE as u16).to_le_bytes())?;
        w.write_all(&self.bit_len.to_le_bytes())?;
        w.write_all(&self.expected_elements.to_le_bytes())?;
        w.write_all(&self.false_positive_rate.to_le_bytes())?;
        w.write_all(&(self.words.len() as u64).to_le_bytes())?;

        for word in &self.words {
            w.write_all(&word.to_le_bytes())?;
        }
        w.flush()?;
        Ok(())
    }

    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path.as_ref())
            .with_context(|| format!("failed to open {}", path.as_ref().display()))?;
        let mut r = BufReader::new(file);

        let mut magic = [0_u8; 8];
        r.read_exact(&mut magic)?;
        if magic != MAGIC {
            bail!("unsupported filter format for {}", path.as_ref().display());
        }

        let mut b1 = [0_u8; 1];
        let mut b2 = [0_u8; 2];
        let mut b8 = [0_u8; 8];

        r.read_exact(&mut b1)?;
        let kmer_size = b1[0];

        r.read_exact(&mut b1)?;
        let hash_count = b1[0];

        r.read_exact(&mut b2)?;
        let header_size = u16::from_le_bytes(b2) as usize;
        if header_size < HEADER_SIZE {
            bail!("invalid header size in {}", path.as_ref().display());
        }

        r.read_exact(&mut b8)?;
        let bit_len = u64::from_le_bytes(b8);

        r.read_exact(&mut b8)?;
        let expected_elements = u64::from_le_bytes(b8);

        r.read_exact(&mut b8)?;
        let false_positive_rate = f64::from_le_bytes(b8);

        r.read_exact(&mut b8)?;
        let word_len = u64::from_le_bytes(b8) as usize;

        if header_size > HEADER_SIZE {
            let mut skip = vec![0_u8; header_size - HEADER_SIZE];
            r.read_exact(&mut skip)?;
        }

        let mut words = Vec::with_capacity(word_len);
        for _ in 0..word_len {
            r.read_exact(&mut b8)?;
            words.push(u64::from_le_bytes(b8));
        }

        Ok(Self {
            kmer_size,
            hash_count,
            bit_len,
            expected_elements,
            false_positive_rate,
            bit_mask: if bit_len.is_power_of_two() {
                bit_len - 1
            } else {
                0
            },
            words,
        })
    }
}

pub fn optimal_bit_len(expected_elements: u64, false_positive_rate: f64) -> u64 {
    let n = expected_elements.max(1) as f64;
    let fpr = false_positive_rate.clamp(1e-12, 0.999_999_999_999);
    (-(n * fpr.ln()) / (std::f64::consts::LN_2 * std::f64::consts::LN_2)).ceil() as u64
}

pub fn optimal_hash_count(bit_len: u64, expected_elements: u64) -> u8 {
    let n = expected_elements.max(1) as f64;
    let m = bit_len.max(1) as f64;
    ((m / n) * std::f64::consts::LN_2).round().clamp(1.0, 255.0) as u8
}

pub fn bit_len_for_fpr_with_hash_count(
    expected_elements: u64,
    false_positive_rate: f64,
    hash_count: u8,
) -> u64 {
    let n = expected_elements.max(1) as f64;
    let p = false_positive_rate.clamp(1e-12, 0.999_999_999_999);
    let k = hash_count.max(1) as f64;
    let one_minus = 1.0 - p.powf(1.0 / k);
    if !(0.0..1.0).contains(&one_minus) {
        return optimal_bit_len(expected_elements, false_positive_rate);
    }
    let denom = one_minus.ln();
    if !denom.is_finite() || denom >= 0.0 {
        return optimal_bit_len(expected_elements, false_positive_rate);
    }
    ((-k * n) / denom).ceil() as u64
}

#[inline]
fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E3779B97F4A7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D049BB133111EB);
    x ^ (x >> 31)
}

#[inline]
fn hash_pair(key: u64) -> (u64, u64) {
    let h1 = splitmix64(key ^ 0xA24BAED4963EE407);
    let mut h2 = splitmix64(key ^ 0x9FB21C651E98DF25);
    if h2 == 0 {
        h2 = 0x9E3779B97F4A7C15;
    }
    (h1, h2 | 1)
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

pub fn for_each_canonical_kmer(seq: &[u8], k: usize, mut f: impl FnMut(u64)) -> usize {
    if !(1..=32).contains(&k) {
        return 0;
    }

    let mut forward = 0_u64;
    let mut reverse = 0_u64;
    let mask = if k == 32 {
        u64::MAX
    } else {
        (1_u64 << (2 * k)) - 1
    };
    let rev_shift = 2 * (k - 1);
    let mut run = 0_usize;
    let mut count = 0_usize;

    for &b in seq {
        match nuc2bits(b) {
            Some(code) => {
                forward = ((forward << 2) | code as u64) & mask;
                let comp = (3_u8 - code) as u64;
                reverse = (reverse >> 2) | (comp << rev_shift);

                run += 1;
                if run >= k {
                    count += 1;
                    let canonical = forward.min(reverse);
                    f(canonical);
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
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn canonical_kmers_count_and_values() {
        let mut out = Vec::new();
        let count = for_each_canonical_kmer(b"ACGTNACG", 3, |k| out.push(k));
        assert_eq!(count, 3);
        assert_eq!(out.len(), 3);
    }

    #[test]
    fn bloom_roundtrip() {
        let bf = ConcurrentBloomFilter::new(21, 4, 10_000, 1000, 0.01)
            .expect("concurrent bf init should succeed");
        bf.insert_kmer(123);
        bf.insert_kmer(456);

        let bf = bf.finalize();
        assert!(bf.contains_kmer(123));
        assert!(bf.contains_kmer(456));

        let d = tempdir().expect("tempdir should be creatable");
        let p = d.path().join("t.bf");
        bf.save(&p).expect("save should succeed");

        let loaded = BloomFilter::load(&p).expect("load should succeed");
        assert!(loaded.contains_kmer(123));
        assert!(loaded.contains_kmer(456));
        assert_eq!(loaded.kmer_size, 21);
        assert_eq!(loaded.hash_count, 4);
    }
}
