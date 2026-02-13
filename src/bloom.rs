use anyhow::{Context, Result, bail};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};

const MAGIC: [u8; 8] = *b"BBRFILT1";
const HEADER_V1_SIZE: usize = 8 + 1 + 1 + 2 + 8 + 8 + 8 + 8;
const HEADER_V2_SIZE: usize = HEADER_V1_SIZE + 1 + 2 + 1;
pub const DEFAULT_BLOCK_WORDS: u16 = 8;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum BloomLayout {
    Classic = 0,
    Blocked = 1,
}

impl BloomLayout {
    fn from_u8(v: u8) -> Result<Self> {
        match v {
            0 => Ok(Self::Classic),
            1 => Ok(Self::Blocked),
            _ => bail!("unsupported bloom layout {}", v),
        }
    }
}

#[derive(Clone, Debug)]
pub struct BloomFilter {
    pub kmer_size: u8,
    pub hash_count: u8,
    pub bit_len: u64,
    pub expected_elements: u64,
    pub false_positive_rate: f64,
    pub layout: BloomLayout,
    pub block_words: u16,
    bit_mask: u64,
    block_count: u64,
    block_count_mask: u64,
    block_bits: u64,
    block_bits_mask: u64,
    words: Vec<u64>,
}

#[derive(Debug)]
pub struct ConcurrentBloomFilter {
    pub kmer_size: u8,
    pub hash_count: u8,
    pub bit_len: u64,
    pub expected_elements: u64,
    pub false_positive_rate: f64,
    pub layout: BloomLayout,
    pub block_words: u16,
    bit_mask: u64,
    block_count: u64,
    block_count_mask: u64,
    block_bits: u64,
    block_bits_mask: u64,
    words: Vec<AtomicU64>,
}

impl ConcurrentBloomFilter {
    pub fn new(
        kmer_size: u8,
        hash_count: u8,
        bit_len: u64,
        expected_elements: u64,
        false_positive_rate: f64,
        blocked: bool,
        block_words: u16,
    ) -> Result<Self> {
        if bit_len == 0 {
            bail!("bit_len must be > 0");
        }
        if hash_count == 0 {
            bail!("hash_count must be > 0");
        }
        let word_len = bit_len.div_ceil(64) as usize;
        let layout = if blocked {
            BloomLayout::Blocked
        } else {
            BloomLayout::Classic
        };
        let (block_words, block_count, block_count_mask, block_bits, block_bits_mask) = match layout
        {
            BloomLayout::Classic => (0_u16, 0_u64, 0_u64, 0_u64, 0_u64),
            BloomLayout::Blocked => {
                if block_words == 0 || !block_words.is_power_of_two() {
                    bail!("block_words must be a non-zero power-of-two");
                }
                if word_len < block_words as usize {
                    bail!("bit_len is too small for blocked mode");
                }
                if word_len % block_words as usize != 0 {
                    bail!("bit_len/64 must be a multiple of block_words");
                }
                let block_count = (word_len / block_words as usize) as u64;
                let block_bits = u64::from(block_words) * 64;
                (
                    block_words,
                    block_count,
                    if block_count.is_power_of_two() {
                        block_count - 1
                    } else {
                        0
                    },
                    block_bits,
                    if block_bits.is_power_of_two() {
                        block_bits - 1
                    } else {
                        0
                    },
                )
            }
        };
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
            layout,
            block_words,
            bit_mask,
            block_count,
            block_count_mask,
            block_bits,
            block_bits_mask,
            words,
        })
    }

    #[inline]
    pub fn insert_kmer(&self, canonical_kmer: u64) {
        let (h1, h2) = hash_pair(canonical_kmer);
        if self.layout == BloomLayout::Blocked {
            let block_idx = if self.block_count_mask != 0 {
                (h1 & self.block_count_mask) as usize
            } else {
                (h1 % self.block_count) as usize
            };
            let base = block_idx * self.block_words as usize;
            let mut h = h1 ^ h2.rotate_left(17);
            for _ in 0..self.hash_count {
                let in_block = if self.block_bits_mask != 0 {
                    (h & self.block_bits_mask) as usize
                } else {
                    (h % self.block_bits) as usize
                };
                let word_idx = base + (in_block >> 6);
                let bit = in_block & 63;
                let mask = 1_u64 << bit;
                self.words[word_idx].fetch_or(mask, Ordering::Relaxed);
                h = h.wrapping_add(h2);
            }
        } else {
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
    }

    #[inline]
    pub fn contains_kmer(&self, canonical_kmer: u64) -> bool {
        let (h1, h2) = hash_pair(canonical_kmer);
        if self.layout == BloomLayout::Blocked {
            let block_idx = if self.block_count_mask != 0 {
                (h1 & self.block_count_mask) as usize
            } else {
                (h1 % self.block_count) as usize
            };
            let base = block_idx * self.block_words as usize;
            let mut h = h1 ^ h2.rotate_left(17);
            for _ in 0..self.hash_count {
                let in_block = if self.block_bits_mask != 0 {
                    (h & self.block_bits_mask) as usize
                } else {
                    (h % self.block_bits) as usize
                };
                let word_idx = base + (in_block >> 6);
                let bit = in_block & 63;
                let mask = 1_u64 << bit;
                if self.words[word_idx].load(Ordering::Relaxed) & mask == 0 {
                    return false;
                }
                h = h.wrapping_add(h2);
            }
        } else {
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
                if self.words[word_idx].load(Ordering::Relaxed) & mask == 0 {
                    return false;
                }
                h = h.wrapping_add(h2);
            }
        }
        true
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
            layout: self.layout,
            block_words: self.block_words,
            bit_mask: self.bit_mask,
            block_count: self.block_count,
            block_count_mask: self.block_count_mask,
            block_bits: self.block_bits,
            block_bits_mask: self.block_bits_mask,
            words,
        }
    }
}

impl BloomFilter {
    pub fn layout_name(&self) -> &'static str {
        match self.layout {
            BloomLayout::Classic => "classic",
            BloomLayout::Blocked => "blocked",
        }
    }

    #[inline]
    pub fn contains_kmer(&self, canonical_kmer: u64) -> bool {
        let (h1, h2) = hash_pair(canonical_kmer);
        if self.layout == BloomLayout::Blocked {
            let block_idx = if self.block_count_mask != 0 {
                (h1 & self.block_count_mask) as usize
            } else {
                (h1 % self.block_count) as usize
            };
            let base = block_idx * self.block_words as usize;
            let mut h = h1 ^ h2.rotate_left(17);
            for _ in 0..self.hash_count {
                let in_block = if self.block_bits_mask != 0 {
                    (h & self.block_bits_mask) as usize
                } else {
                    (h % self.block_bits) as usize
                };
                let word_idx = base + (in_block >> 6);
                let bit = in_block & 63;
                let mask = 1_u64 << bit;
                if self.words[word_idx] & mask == 0 {
                    return false;
                }
                h = h.wrapping_add(h2);
            }
        } else {
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
        w.write_all(&(HEADER_V2_SIZE as u16).to_le_bytes())?;
        w.write_all(&self.bit_len.to_le_bytes())?;
        w.write_all(&self.expected_elements.to_le_bytes())?;
        w.write_all(&self.false_positive_rate.to_le_bytes())?;
        w.write_all(&(self.words.len() as u64).to_le_bytes())?;
        w.write_all(&[self.layout as u8])?;
        w.write_all(&self.block_words.to_le_bytes())?;
        w.write_all(&[0_u8])?;

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
        if header_size < HEADER_V1_SIZE {
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

        let mut parsed = HEADER_V1_SIZE;
        let mut layout = BloomLayout::Classic;
        let mut block_words = 0_u16;
        if parsed < header_size {
            r.read_exact(&mut b1)?;
            layout = BloomLayout::from_u8(b1[0])?;
            parsed += 1;
        }
        if parsed + 1 < header_size {
            r.read_exact(&mut b2)?;
            block_words = u16::from_le_bytes(b2);
            parsed += 2;
        }
        if parsed < header_size {
            let mut skip = vec![0_u8; header_size - parsed];
            r.read_exact(&mut skip)?;
        }

        let (block_words, block_count, block_count_mask, block_bits, block_bits_mask) = match layout
        {
            BloomLayout::Classic => (0_u16, 0_u64, 0_u64, 0_u64, 0_u64),
            BloomLayout::Blocked => {
                if block_words == 0 || !block_words.is_power_of_two() {
                    bail!(
                        "invalid blocked layout header in {}",
                        path.as_ref().display()
                    );
                }
                if word_len < block_words as usize || word_len % block_words as usize != 0 {
                    bail!("blocked layout is inconsistent with word length");
                }
                let block_count = (word_len / block_words as usize) as u64;
                let block_bits = u64::from(block_words) * 64;
                (
                    block_words,
                    block_count,
                    if block_count.is_power_of_two() {
                        block_count - 1
                    } else {
                        0
                    },
                    block_bits,
                    if block_bits.is_power_of_two() {
                        block_bits - 1
                    } else {
                        0
                    },
                )
            }
        };

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
            layout,
            block_words,
            bit_mask: if bit_len.is_power_of_two() {
                bit_len - 1
            } else {
                0
            },
            block_count,
            block_count_mask,
            block_bits,
            block_bits_mask,
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
    for_each_canonical_kmer_until(seq, k, |canonical| {
        f(canonical);
        true
    })
}

pub fn for_each_canonical_kmer_until(
    seq: &[u8],
    k: usize,
    mut f: impl FnMut(u64) -> bool,
) -> usize {
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
                    if !f(canonical) {
                        return count;
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
        let bf = ConcurrentBloomFilter::new(21, 4, 10_000, 1000, 0.01, false, 0)
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
        assert_eq!(loaded.layout, BloomLayout::Classic);
    }

    #[test]
    fn blocked_bloom_roundtrip() {
        let bf = ConcurrentBloomFilter::new(21, 4, 32_768, 1000, 0.01, true, 8)
            .expect("blocked bf init should succeed");
        bf.insert_kmer(123);
        bf.insert_kmer(456);

        let bf = bf.finalize();
        assert!(bf.contains_kmer(123));
        assert!(bf.contains_kmer(456));
        assert_eq!(bf.layout, BloomLayout::Blocked);
        assert_eq!(bf.block_words, 8);

        let d = tempdir().expect("tempdir should be creatable");
        let p = d.path().join("t_blocked.bf");
        bf.save(&p).expect("save should succeed");

        let loaded = BloomFilter::load(&p).expect("load should succeed");
        assert!(loaded.contains_kmer(123));
        assert!(loaded.contains_kmer(456));
        assert_eq!(loaded.layout, BloomLayout::Blocked);
        assert_eq!(loaded.block_words, 8);
    }
}
