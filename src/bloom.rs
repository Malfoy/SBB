use crate::bloomybloom_adapter::{
    BloomyConcurrentState, BloomyFrozenState, BloomyParams, MIN_BLOOMY_BITS,
};
use anyhow::{Context, Result, bail};
use rayon::prelude::*;
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

const MAGIC: [u8; 8] = *b"BBRFILT1";
const SHARDED_MAGIC: [u8; 8] = *b"BBRSHRD1";
const HEADER_V1_SIZE: usize = 8 + 1 + 1 + 2 + 8 + 8 + 8 + 8;
const HEADER_V2_SIZE: usize = HEADER_V1_SIZE + 1 + 2 + 1;
const HEADER_V3_SIZE: usize = HEADER_V2_SIZE + 8;
const HEADER_V4_SIZE: usize = HEADER_V3_SIZE + 2 + 2 + 8 + 8 + 1 + 2;
const MIN_SHARD_BYTES: usize = 64 * 1024 * 1024;
pub const DEFAULT_BLOCK_WORDS: u16 = 8;

#[derive(Clone, Debug)]
struct BloomHeader {
    kmer_size: u8,
    hash_count: u8,
    bit_len: u64,
    expected_elements: u64,
    false_positive_rate: f64,
    word_len: usize,
    sample_rate: f64,
    layout: BloomLayout,
    block_words: u16,
    bloomy_params: Option<BloomyParams>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum BloomLayout {
    Classic = 0,
    Blocked = 1,
    BloomyBloom = 2,
}

impl BloomLayout {
    fn from_u8(v: u8) -> Result<Self> {
        match v {
            0 => Ok(Self::Classic),
            1 => Ok(Self::Blocked),
            2 => Ok(Self::BloomyBloom),
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
    pub sample_rate: f64,
    pub layout: BloomLayout,
    pub block_words: u16,
    sample_hash_bound: u64,
    bit_mask: u64,
    block_count: u64,
    block_count_mask: u64,
    block_bits: u64,
    block_bits_mask: u64,
    bloomy_params: Option<BloomyParams>,
    bloomy: Option<BloomyFrozenState>,
    words: Vec<u64>,
}

#[derive(Debug)]
pub struct ConcurrentBloomFilter {
    pub kmer_size: u8,
    pub hash_count: u8,
    pub bit_len: u64,
    pub expected_elements: u64,
    pub false_positive_rate: f64,
    pub sample_rate: f64,
    pub layout: BloomLayout,
    pub block_words: u16,
    sample_hash_bound: u64,
    bit_mask: u64,
    block_count: u64,
    block_count_mask: u64,
    block_bits: u64,
    block_bits_mask: u64,
    bloomy_params: Option<BloomyParams>,
    bloomy: Option<BloomyConcurrentState>,
    words: Vec<AtomicU64>,
}

#[inline]
fn blocked_word_masks_bw8(
    mut h: u64,
    h2: u64,
    hash_count: u8,
    block_bits_mask: u64,
    block_bits: u64,
) -> [u64; 8] {
    let mut masks = [0_u64; 8];
    for _ in 0..hash_count {
        let in_block = if block_bits_mask != 0 {
            (h & block_bits_mask) as usize
        } else {
            (h % block_bits) as usize
        };
        masks[in_block >> 6] |= 1_u64 << (in_block & 63);
        h = h.wrapping_add(h2);
    }
    masks
}

#[inline]
fn blocked_word_masks_bw8_h8(mut h: u64, h2: u64, block_bits_mask: u64) -> [u64; 8] {
    let mut masks = [0_u64; 8];
    macro_rules! step {
        () => {{
            let in_block = (h & block_bits_mask) as usize;
            masks[in_block >> 6] |= 1_u64 << (in_block & 63);
            h = h.wrapping_add(h2);
        }};
    }
    step!();
    step!();
    step!();
    step!();
    step!();
    step!();
    step!();
    step!();
    let _ = h;
    masks
}

impl ConcurrentBloomFilter {
    pub fn new(
        kmer_size: u8,
        hash_count: u8,
        bit_len: u64,
        expected_elements: u64,
        false_positive_rate: f64,
        sample_rate: f64,
        blocked: bool,
        block_words: u16,
    ) -> Result<Self> {
        if bit_len == 0 {
            bail!("bit_len must be > 0");
        }
        if hash_count == 0 {
            bail!("hash_count must be > 0");
        }
        if !(0.0..=1.0).contains(&sample_rate) {
            bail!("sample_rate must be in [0, 1]");
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
            BloomLayout::BloomyBloom => unreachable!("bloomybloom filters must use new_bloomybloom"),
        };
        let mut words = Vec::with_capacity(word_len);
        words.resize_with(word_len, || AtomicU64::new(0));
        let bit_mask = if bit_len.is_power_of_two() {
            bit_len - 1
        } else {
            0
        };

        let sample_hash_bound = sampling_hash_bound(sample_rate);
        Ok(Self {
            kmer_size,
            hash_count,
            bit_len,
            expected_elements,
            false_positive_rate,
            sample_rate,
            layout,
            block_words,
            sample_hash_bound,
            bit_mask,
            block_count,
            block_count_mask,
            block_bits,
            block_bits_mask,
            bloomy_params: None,
            bloomy: None,
            words,
        })
    }

    pub fn new_bloomybloom(
        kmer_size: u8,
        bit_len: u64,
        expected_elements: u64,
        false_positive_rate: f64,
        sample_rate: f64,
        requested_hash_count: Option<u8>,
    ) -> Result<Self> {
        if bit_len == 0 {
            bail!("bit_len must be > 0");
        }
        if !(0.0..=1.0).contains(&sample_rate) {
            bail!("sample_rate must be in [0, 1]");
        }

        let (bloomy, params) =
            BloomyConcurrentState::new(kmer_size, requested_hash_count, bit_len.max(MIN_BLOOMY_BITS))?;
        Ok(Self {
            kmer_size,
            hash_count: params.hash_count,
            bit_len: params.bit_len,
            expected_elements,
            false_positive_rate,
            sample_rate,
            layout: BloomLayout::BloomyBloom,
            block_words: 0,
            sample_hash_bound: sampling_hash_bound(sample_rate),
            bit_mask: 0,
            block_count: 0,
            block_count_mask: 0,
            block_bits: 0,
            block_bits_mask: 0,
            bloomy_params: Some(params),
            bloomy: Some(bloomy),
            words: Vec::new(),
        })
    }

    #[inline]
    pub fn should_keep_kmer(&self, canonical_kmer: u64) -> bool {
        if self.sample_rate >= 1.0 {
            return true;
        }
        if self.sample_rate <= 0.0 {
            return false;
        }
        sampling_hash(canonical_kmer) < self.sample_hash_bound
    }

    #[inline]
    pub fn insert_kmer(&self, canonical_kmer: u64) -> bool {
        if self.layout == BloomLayout::BloomyBloom {
            panic!("insert_kmer is not supported for bloomybloom filters; use insert_sequence");
        }
        if !self.should_keep_kmer(canonical_kmer) {
            return false;
        }
        let (h1, h2) = hash_pair(canonical_kmer);
        if self.layout == BloomLayout::Blocked {
            let block_idx = if self.block_count_mask != 0 {
                (h1 & self.block_count_mask) as usize
            } else {
                (h1 % self.block_count) as usize
            };
            let base = block_idx * self.block_words as usize;
            if self.block_words == 8 && self.hash_count == 8 && self.block_bits_mask != 0 {
                let masks = blocked_word_masks_bw8_h8(
                    h1 ^ h2.rotate_left(17),
                    h2,
                    self.block_bits_mask,
                );
                // SAFETY: block_words == 8 guarantees base..base+7 are valid indexes.
                let words = unsafe { self.words.get_unchecked(base..base + 8) };
                if masks[0] != 0 {
                    words[0].fetch_or(masks[0], Ordering::Relaxed);
                }
                if masks[1] != 0 {
                    words[1].fetch_or(masks[1], Ordering::Relaxed);
                }
                if masks[2] != 0 {
                    words[2].fetch_or(masks[2], Ordering::Relaxed);
                }
                if masks[3] != 0 {
                    words[3].fetch_or(masks[3], Ordering::Relaxed);
                }
                if masks[4] != 0 {
                    words[4].fetch_or(masks[4], Ordering::Relaxed);
                }
                if masks[5] != 0 {
                    words[5].fetch_or(masks[5], Ordering::Relaxed);
                }
                if masks[6] != 0 {
                    words[6].fetch_or(masks[6], Ordering::Relaxed);
                }
                if masks[7] != 0 {
                    words[7].fetch_or(masks[7], Ordering::Relaxed);
                }
            } else {
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
        true
    }

    #[inline]
    pub fn contains_kmer(&self, canonical_kmer: u64) -> bool {
        if self.layout == BloomLayout::BloomyBloom {
            panic!("contains_kmer is not supported for bloomybloom filters; use sequence queries");
        }
        if !self.should_keep_kmer(canonical_kmer) {
            return false;
        }
        self.contains_kmer_kept(canonical_kmer)
    }

    #[inline]
    pub(crate) fn contains_kmer_kept(&self, canonical_kmer: u64) -> bool {
        if self.layout == BloomLayout::BloomyBloom {
            panic!("contains_kmer is not supported for bloomybloom filters; use sequence queries");
        }
        let (h1, h2) = hash_pair(canonical_kmer);
        if self.layout == BloomLayout::Blocked {
            let block_idx = if self.block_count_mask != 0 {
                (h1 & self.block_count_mask) as usize
            } else {
                (h1 % self.block_count) as usize
            };
            let base = block_idx * self.block_words as usize;
            if self.block_words == 8 && self.hash_count == 8 && self.block_bits_mask != 0 {
                let mut h = h1 ^ h2.rotate_left(17);
                // SAFETY: block_words == 8 guarantees base..base+7 are valid indexes.
                let words = unsafe { self.words.get_unchecked(base..base + 8) };
                macro_rules! step {
                    () => {{
                        let in_block = (h & self.block_bits_mask) as usize;
                        let word =
                            unsafe { words.get_unchecked(in_block >> 6) }.load(Ordering::Relaxed);
                        if (word & (1_u64 << (in_block & 63))) == 0 {
                            return false;
                        }
                        h = h.wrapping_add(h2);
                    }};
                }
                step!();
                step!();
                step!();
                step!();
                step!();
                step!();
                step!();
                step!();
                let _ = h;
            } else if self.block_words == 8 {
                let masks = blocked_word_masks_bw8(
                    h1 ^ h2.rotate_left(17),
                    h2,
                    self.hash_count,
                    self.block_bits_mask,
                    self.block_bits,
                );
                for (i, mask) in masks.into_iter().enumerate() {
                    if mask != 0
                        && (self.words[base + i].load(Ordering::Relaxed) & mask) != mask
                    {
                        return false;
                    }
                }
            } else {
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
        if self.layout == BloomLayout::BloomyBloom {
            let params = self
                .bloomy_params
                .expect("bloomybloom params should be present");
            let bloomy = self
                .bloomy
                .expect("bloomybloom state should be present")
                .finalize();

            return BloomFilter {
                kmer_size: self.kmer_size,
                hash_count: self.hash_count,
                bit_len: self.bit_len,
                expected_elements: self.expected_elements,
                false_positive_rate: self.false_positive_rate,
                sample_rate: self.sample_rate,
                layout: self.layout,
                block_words: 0,
                sample_hash_bound: self.sample_hash_bound,
                bit_mask: 0,
                block_count: 0,
                block_count_mask: 0,
                block_bits: 0,
                block_bits_mask: 0,
                bloomy_params: Some(params),
                bloomy: Some(bloomy),
                words: Vec::new(),
            };
        }

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
            sample_rate: self.sample_rate,
            layout: self.layout,
            block_words: self.block_words,
            sample_hash_bound: self.sample_hash_bound,
            bit_mask: self.bit_mask,
            block_count: self.block_count,
            block_count_mask: self.block_count_mask,
            block_bits: self.block_bits,
            block_bits_mask: self.block_bits_mask,
            bloomy_params: None,
            bloomy: None,
            words,
        }
    }

    pub fn insert_sequence(&self, seq: &[u8]) -> u64 {
        if self.layout == BloomLayout::BloomyBloom {
            if self.sample_rate < 1.0 {
                panic!("sample_alpha < 1.0 is not supported for bloomybloom filters");
            }
            return self
                .bloomy
                .as_ref()
                .expect("bloomybloom state should be present")
                .insert_sequence(seq);
        }

        let mut inserted = 0_u64;
        let k = self.kmer_size as usize;
        for_each_canonical_kmer(seq, k, |canonical| {
            if self.insert_kmer(canonical) {
                inserted += 1;
            }
        });
        inserted
    }

    pub fn score_sequence_with_counts(&self, seq: &[u8]) -> (u64, u64) {
        if self.layout == BloomLayout::BloomyBloom {
            if self.sample_rate < 1.0 {
                panic!("sample_alpha < 1.0 is not supported for bloomybloom filters");
            }
            return self
                .bloomy
                .as_ref()
                .expect("bloomybloom state should be present")
                .score_sequence_with_counts(seq);
        }

        let mut total = 0_u64;
        let mut hits = 0_u64;
        let k = self.kmer_size as usize;
        for_each_canonical_kmer(seq, k, |canonical| {
            if !self.should_keep_kmer(canonical) {
                return;
            }
            total += 1;
            if self.contains_kmer_kept(canonical) {
                hits += 1;
            }
        });
        (hits, total)
    }

}

impl BloomFilter {
    fn word_len(&self) -> usize {
        if self.layout == BloomLayout::BloomyBloom {
            self.bloomy
                .as_ref()
                .expect("bloomybloom state should be present")
                .serialized_word_len()
        } else {
            self.words.len()
        }
    }

    pub fn layout_name(&self) -> &'static str {
        match self.layout {
            BloomLayout::Classic => "classic",
            BloomLayout::Blocked => "blocked",
            BloomLayout::BloomyBloom => "bloomybloom",
        }
    }

    #[inline]
    pub fn should_keep_kmer(&self, canonical_kmer: u64) -> bool {
        if self.sample_rate >= 1.0 {
            return true;
        }
        if self.sample_rate <= 0.0 {
            return false;
        }
        sampling_hash(canonical_kmer) < self.sample_hash_bound
    }

    #[inline]
    pub fn contains_kmer(&self, canonical_kmer: u64) -> bool {
        if self.layout == BloomLayout::BloomyBloom {
            panic!("contains_kmer is not supported for bloomybloom filters; use sequence queries");
        }
        if !self.should_keep_kmer(canonical_kmer) {
            return false;
        }
        self.contains_kmer_kept(canonical_kmer)
    }

    #[inline]
    pub(crate) fn contains_kmer_kept(&self, canonical_kmer: u64) -> bool {
        if self.layout == BloomLayout::BloomyBloom {
            panic!("contains_kmer is not supported for bloomybloom filters; use sequence queries");
        }
        let (h1, h2) = hash_pair(canonical_kmer);
        if self.layout == BloomLayout::Blocked {
            let block_idx = if self.block_count_mask != 0 {
                (h1 & self.block_count_mask) as usize
            } else {
                (h1 % self.block_count) as usize
            };
            let base = block_idx * self.block_words as usize;
            if self.block_words == 8 && self.hash_count == 8 && self.block_bits_mask != 0 {
                let mut h = h1 ^ h2.rotate_left(17);
                // SAFETY: block_words == 8 guarantees base..base+7 are valid indexes.
                let words = unsafe { self.words.get_unchecked(base..base + 8) };
                macro_rules! step {
                    () => {{
                        let in_block = (h & self.block_bits_mask) as usize;
                        let word = unsafe { *words.get_unchecked(in_block >> 6) };
                        if (word & (1_u64 << (in_block & 63))) == 0 {
                            return false;
                        }
                        h = h.wrapping_add(h2);
                    }};
                }
                step!();
                step!();
                step!();
                step!();
                step!();
                step!();
                step!();
                step!();
                let _ = h;
            } else if self.block_words == 8 {
                let masks = blocked_word_masks_bw8(
                    h1 ^ h2.rotate_left(17),
                    h2,
                    self.hash_count,
                    self.block_bits_mask,
                    self.block_bits,
                );
                for (i, mask) in masks.into_iter().enumerate() {
                    if mask != 0 && (self.words[base + i] & mask) != mask {
                        return false;
                    }
                }
            } else {
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

    pub fn score_sequence_with_counts(&self, seq: &[u8]) -> (u64, u64) {
        if self.layout == BloomLayout::BloomyBloom {
            if self.sample_rate < 1.0 {
                panic!("sample_alpha < 1.0 is not supported for bloomybloom filters");
            }
            return self
                .bloomy
                .as_ref()
                .expect("bloomybloom state should be present")
                .score_sequence_with_counts(seq);
        }

        let mut total = 0_u64;
        let mut hits = 0_u64;
        let k = self.kmer_size as usize;
        for_each_canonical_kmer(seq, k, |canonical| {
            if !self.should_keep_kmer(canonical) {
                return;
            }
            total += 1;
            if self.contains_kmer_kept(canonical) {
                hits += 1;
            }
        });
        (hits, total)
    }

    pub(crate) fn match_sequence_threshold_with_counts(
        &self,
        seq: &[u8],
        threshold: f64,
        best_hit: bool,
    ) -> (bool, u64) {
        let k = self.kmer_size as usize;
        if seq.len() < k {
            return (if best_hit { false } else { 0.0 >= threshold }, 0);
        }

        if self.layout == BloomLayout::BloomyBloom {
            if self.sample_rate < 1.0 {
                panic!("sample_alpha < 1.0 is not supported for bloomybloom filters");
            }
            return self
                .bloomy
                .as_ref()
                .expect("bloomybloom state should be present")
                .match_sequence_threshold_with_counts(seq, threshold, best_hit);
        }

        if best_hit {
            let mut queried = 0_u64;
            let mut matched = false;
            for_each_canonical_kmer_until(seq, k, |canonical| {
                if !self.should_keep_kmer(canonical) {
                    return true;
                }
                queried += 1;
                if self.contains_kmer_kept(canonical) {
                    matched = true;
                    false
                } else {
                    true
                }
            });
            return (matched, queried);
        }

        if threshold <= 0.0 {
            let queried = if self.sample_rate >= 1.0 {
                for_each_canonical_kmer(seq, k, |_| {}) as u64
            } else {
                let mut kept = 0_u64;
                for_each_canonical_kmer(seq, k, |canonical| {
                    if self.should_keep_kmer(canonical) {
                        kept += 1;
                    }
                });
                kept
            };
            return (true, queried);
        }

        if threshold >= 1.0 {
            let mut queried = 0_u64;
            let mut missed = false;
            for_each_canonical_kmer_until(seq, k, |canonical| {
                if !self.should_keep_kmer(canonical) {
                    return true;
                }
                queried += 1;
                if self.contains_kmer_kept(canonical) {
                    true
                } else {
                    missed = true;
                    false
                }
            });
            return (queried > 0 && !missed, queried);
        }

        let total_candidates = for_each_canonical_kmer(seq, k, |_| {}) as u64;
        if total_candidates == 0 {
            return (false, 0);
        }

        let mut queried = 0_u64;
        let mut hits = 0_u64;
        let mut processed = 0_u64;
        let mut decided: Option<bool> = None;
        for_each_canonical_kmer_until(seq, k, |canonical| {
            processed += 1;
            if self.should_keep_kmer(canonical) {
                queried += 1;
                if self.contains_kmer_kept(canonical) {
                    hits += 1;
                }
            }

            if queried == 0 {
                return true;
            }

            let remaining = total_candidates - processed;
            let denom_upper = queried + remaining;
            if (hits as f64 / denom_upper as f64) >= threshold {
                decided = Some(true);
                return false;
            }
            if ((hits + remaining) as f64 / denom_upper as f64) < threshold {
                decided = Some(false);
                return false;
            }
            true
        });

        if let Some(matched) = decided {
            return (matched, queried);
        }
        if queried == 0 {
            (false, 0)
        } else {
            ((hits as f64 / queried as f64) >= threshold, queried)
        }
    }

    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let path = path.as_ref();
        let shard_count = preferred_shard_count(self.word_len());
        if shard_count <= 1 {
            self.save_single_stream(path)?;
            cleanup_shard_files(path, 0)?;
            return Ok(());
        }

        self.save_sharded(path, shard_count)
    }

    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let file =
            File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
        let reader = BufReader::new(file);
        let mut reader = reader;
        let mut magic = [0_u8; 8];
        reader.read_exact(&mut magic)?;
        if magic == SHARDED_MAGIC {
            return Self::load_sharded(&mut reader, path);
        }

        let file =
            File::open(path).with_context(|| format!("failed to reopen {}", path.display()))?;
        let reader = BufReader::new(file);
        let mut decoder = zstd::stream::read::Decoder::new(reader).with_context(|| {
            format!(
                "failed to create zstd decoder for {} (expected zstd-compressed filter)",
                path.display()
            )
        })?;
        Self::load_payload(&mut decoder, path)
    }

    fn header_size(&self) -> usize {
        if self.layout == BloomLayout::BloomyBloom {
            HEADER_V4_SIZE
        } else {
            HEADER_V3_SIZE
        }
    }

    fn write_header<W: Write>(&self, w: &mut W) -> Result<()> {
        let header_size = self.header_size();
        w.write_all(&MAGIC)?;
        w.write_all(&[self.kmer_size])?;
        w.write_all(&[self.hash_count])?;
        w.write_all(&(header_size as u16).to_le_bytes())?;
        w.write_all(&self.bit_len.to_le_bytes())?;
        w.write_all(&self.expected_elements.to_le_bytes())?;
        w.write_all(&self.false_positive_rate.to_le_bytes())?;
        w.write_all(&(self.word_len() as u64).to_le_bytes())?;
        w.write_all(&[self.layout as u8])?;
        w.write_all(&self.block_words.to_le_bytes())?;
        w.write_all(&[0_u8])?;
        w.write_all(&self.sample_rate.to_le_bytes())?;
        if self.layout == BloomLayout::BloomyBloom {
            let params = self
                .bloomy_params
                .expect("bloomybloom params should be present");
            w.write_all(&params.m.to_le_bytes())?;
            w.write_all(&params.s.to_le_bytes())?;
            w.write_all(&params.block_size.to_le_bytes())?;
            w.write_all(&params.nb_blocks.to_le_bytes())?;
            w.write_all(&[params.mode_tag])?;
            w.write_all(&params.open_closed_t.to_le_bytes())?;
        }
        Ok(())
    }

    fn write_words<W: Write>(&self, w: &mut W) -> Result<()> {
        self.write_word_range(0..self.word_len(), w)
    }

    fn write_word_range<W: Write>(&self, range: std::ops::Range<usize>, w: &mut W) -> Result<()> {
        if self.layout == BloomLayout::BloomyBloom {
            self.bloomy
                .as_ref()
                .expect("bloomybloom state should be present")
                .write_word_range(range, w)?;
        } else {
            for word in &self.words[range] {
                w.write_all(&word.to_le_bytes())?;
            }
        }
        Ok(())
    }

    fn write_raw<W: Write>(&self, w: &mut W) -> Result<()> {
        self.write_header(w)?;
        self.write_words(w)?;
        Ok(())
    }

    fn save_single_stream(&self, path: &Path) -> Result<()> {
        let file =
            File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
        let writer = BufWriter::new(file);
        let mut encoder = zstd::stream::write::Encoder::new(writer, 0)
            .with_context(|| format!("failed to create zstd stream for {}", path.display()))?;
        let workers = rayon::current_num_threads();
        if workers > 1 {
            encoder.multithread(workers.try_into().unwrap_or(u32::MAX)).with_context(|| {
                format!(
                    "failed to enable zstd multithreaded compression for {}",
                    path.display()
                )
            })?;
        }
        self.write_raw(&mut encoder)?;
        let mut writer = encoder
            .finish()
            .with_context(|| format!("failed to finalize zstd stream for {}", path.display()))?;
        writer.flush()?;
        Ok(())
    }

    fn save_sharded(&self, path: &Path, shard_count: usize) -> Result<()> {
        (0..shard_count)
            .into_par_iter()
            .map(|shard_idx| -> Result<()> {
                let shard_path = shard_file_path(path, shard_idx)?;
                let file = File::create(&shard_path)
                    .with_context(|| format!("failed to create {}", shard_path.display()))?;
                let writer = BufWriter::new(file);
                let mut encoder =
                    zstd::stream::write::Encoder::new(writer, 0).with_context(|| {
                        format!(
                            "failed to create zstd stream for shard {}",
                            shard_path.display()
                        )
                    })?;
                let range = shard_word_range(self.word_len(), shard_count, shard_idx);
                self.write_word_range(range, &mut encoder)?;
                let mut writer = encoder.finish().with_context(|| {
                    format!(
                        "failed to finalize zstd stream for shard {}",
                        shard_path.display()
                    )
                })?;
                writer.flush()?;
                Ok(())
            })
            .collect::<Vec<_>>()
            .into_iter()
            .collect::<Result<Vec<_>>>()?;

        let manifest_path = manifest_temp_path(path)?;
        let file = File::create(&manifest_path)
            .with_context(|| format!("failed to create {}", manifest_path.display()))?;
        let mut writer = BufWriter::new(file);
        writer.write_all(&SHARDED_MAGIC)?;
        writer.write_all(&(shard_count as u32).to_le_bytes())?;
        self.write_header(&mut writer)?;
        writer.flush()?;
        drop(writer);
        fs::rename(&manifest_path, path).with_context(|| {
            format!(
                "failed to move manifest {} into place at {}",
                manifest_path.display(),
                path.display()
            )
        })?;
        cleanup_shard_files(path, shard_count)?;
        Ok(())
    }

    fn load_sharded<R: Read>(r: &mut R, path: &Path) -> Result<Self> {
        let mut b4 = [0_u8; 4];
        r.read_exact(&mut b4)?;
        let shard_count = u32::from_le_bytes(b4) as usize;
        if shard_count == 0 {
            bail!("invalid shard count in {}", path.display());
        }

        let header = Self::read_header(r, path)?;
        let word_len = header.word_len;
        let shard_words = (0..shard_count)
            .into_par_iter()
            .map(|shard_idx| {
                let shard_path = shard_file_path(path, shard_idx)?;
                let file = File::open(&shard_path)
                    .with_context(|| format!("failed to open {}", shard_path.display()))?;
                let reader = BufReader::new(file);
                let mut decoder = zstd::stream::read::Decoder::new(reader).with_context(|| {
                    format!(
                        "failed to create zstd decoder for shard {}",
                        shard_path.display()
                    )
                })?;
                let range = shard_word_range(word_len, shard_count, shard_idx);
                Self::read_words(&mut decoder, range.len(), &shard_path)
            })
            .collect::<Vec<_>>();

        if header.layout == BloomLayout::BloomyBloom {
            let params = header
                .bloomy_params
                .expect("bloomybloom params should be present");
            let mut words = Vec::with_capacity(word_len);
            for shard in shard_words {
                words.extend(shard?);
            }
            let bloomy = BloomyFrozenState::from_header(words, header.kmer_size, params)?;
            return Ok(Self {
                kmer_size: header.kmer_size,
                hash_count: header.hash_count,
                bit_len: header.bit_len,
                expected_elements: header.expected_elements,
                false_positive_rate: header.false_positive_rate,
                sample_rate: header.sample_rate,
                layout: header.layout,
                block_words: 0,
                sample_hash_bound: sampling_hash_bound(header.sample_rate),
                bit_mask: 0,
                block_count: 0,
                block_count_mask: 0,
                block_bits: 0,
                block_bits_mask: 0,
                bloomy_params: Some(params),
                bloomy: Some(bloomy),
                words: Vec::new(),
            });
        }

        let mut words = Vec::with_capacity(word_len);
        for shard in shard_words {
            words.extend(shard?);
        }
        if words.len() != word_len {
            bail!(
                "sharded filter {} reconstructed {} words, expected {}",
                path.display(),
                words.len(),
                word_len
            );
        }

        Self::from_parts(header, words, path)
    }

    fn read_header<R: Read>(r: &mut R, path: &Path) -> Result<BloomHeader> {
        let mut magic = [0_u8; 8];
        r.read_exact(&mut magic)?;
        if magic != MAGIC {
            bail!("unsupported filter format for {}", path.display());
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
            bail!("invalid header size in {}", path.display());
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
        let mut sample_rate = 1.0_f64;
        let mut bloomy_params = None;
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
            r.read_exact(&mut b1)?;
            parsed += 1;
        }
        if parsed + 7 < header_size {
            r.read_exact(&mut b8)?;
            sample_rate = f64::from_le_bytes(b8);
            parsed += 8;
            if !(0.0..=1.0).contains(&sample_rate) {
                bail!("invalid sampling rate in {}", path.display());
            }
        }
        if layout == BloomLayout::BloomyBloom {
            if header_size < HEADER_V4_SIZE {
                bail!("invalid bloomy layout header in {}", path.display());
            }

            r.read_exact(&mut b2)?;
            let m = u16::from_le_bytes(b2);
            parsed += 2;

            r.read_exact(&mut b2)?;
            let s = u16::from_le_bytes(b2);
            parsed += 2;

            r.read_exact(&mut b8)?;
            let block_size = u64::from_le_bytes(b8);
            parsed += 8;

            r.read_exact(&mut b8)?;
            let nb_blocks = u64::from_le_bytes(b8);
            parsed += 8;

            r.read_exact(&mut b1)?;
            let mode_tag = b1[0];
            parsed += 1;

            r.read_exact(&mut b2)?;
            let open_closed_t = u16::from_le_bytes(b2);
            parsed += 2;

            bloomy_params = Some(BloomyParams {
                hash_count,
                bit_len,
                m,
                s,
                block_size,
                nb_blocks,
                mode_tag,
                open_closed_t,
            });
        }
        if parsed < header_size {
            let mut skip = vec![0_u8; header_size - parsed];
            r.read_exact(&mut skip)?;
        }

        Ok(BloomHeader {
            kmer_size,
            hash_count,
            bit_len,
            expected_elements,
            false_positive_rate,
            word_len,
            sample_rate,
            layout,
            block_words,
            bloomy_params,
        })
    }

    fn read_words<R: Read>(r: &mut R, word_len: usize, path: &Path) -> Result<Vec<u64>> {
        let mut b8 = [0_u8; 8];
        let mut words = Vec::with_capacity(word_len);
        for _ in 0..word_len {
            r.read_exact(&mut b8)
                .with_context(|| format!("failed to read words from {}", path.display()))?;
            words.push(u64::from_le_bytes(b8));
        }
        Ok(words)
    }

    fn from_parts(header: BloomHeader, words: Vec<u64>, path: &Path) -> Result<Self> {
        let BloomHeader {
            kmer_size,
            hash_count,
            bit_len,
            expected_elements,
            false_positive_rate,
            word_len,
            sample_rate,
            layout,
            block_words,
            bloomy_params,
        } = header;
        if words.len() != word_len {
            bail!(
                "filter {} loaded {} words, expected {}",
                path.display(),
                words.len(),
                word_len
            );
        }

        let (block_words, block_count, block_count_mask, block_bits, block_bits_mask) = match layout
        {
            BloomLayout::Classic => (0_u16, 0_u64, 0_u64, 0_u64, 0_u64),
            BloomLayout::Blocked => {
                if block_words == 0 || !block_words.is_power_of_two() {
                    bail!("invalid blocked layout header in {}", path.display());
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
            BloomLayout::BloomyBloom => (0_u16, 0_u64, 0_u64, 0_u64, 0_u64),
        };

        if let Some(params) = bloomy_params {
            let bloomy = BloomyFrozenState::from_header(words, kmer_size, params)?;
            return Ok(Self {
                kmer_size,
                hash_count,
                bit_len,
                expected_elements,
                false_positive_rate,
                sample_rate,
                layout,
                block_words: 0,
                sample_hash_bound: sampling_hash_bound(sample_rate),
                bit_mask: 0,
                block_count: 0,
                block_count_mask: 0,
                block_bits: 0,
                block_bits_mask: 0,
                bloomy_params: Some(params),
                bloomy: Some(bloomy),
                words: Vec::new(),
            });
        }

        Ok(Self {
            kmer_size,
            hash_count,
            bit_len,
            expected_elements,
            false_positive_rate,
            sample_rate,
            layout,
            block_words,
            sample_hash_bound: sampling_hash_bound(sample_rate),
            bit_mask: if bit_len.is_power_of_two() {
                bit_len - 1
            } else {
                0
            },
            block_count,
            block_count_mask,
            block_bits,
            block_bits_mask,
            bloomy_params: None,
            bloomy: None,
            words,
        })
    }

    fn load_payload<R: Read>(r: &mut R, path: &Path) -> Result<Self> {
        let header = Self::read_header(r, path)?;
        let words = Self::read_words(r, header.word_len, path)?;
        Self::from_parts(header, words, path)
    }
}

fn preferred_shard_count(word_len: usize) -> usize {
    if word_len == 0 {
        return 1;
    }
    let workers = rayon::current_num_threads().max(1);
    let total_bytes = word_len.saturating_mul(std::mem::size_of::<u64>());
    let size_limited = total_bytes.div_ceil(MIN_SHARD_BYTES).max(1);
    workers.min(size_limited).min(word_len)
}

fn shard_word_range(word_len: usize, shard_count: usize, shard_idx: usize) -> std::ops::Range<usize> {
    let start = shard_idx * word_len / shard_count;
    let end = (shard_idx + 1) * word_len / shard_count;
    start..end
}

fn shard_file_path(path: &Path, shard_idx: usize) -> Result<PathBuf> {
    let file_name = path
        .file_name()
        .and_then(|x| x.to_str())
        .with_context(|| format!("filter path {} has no valid file name", path.display()))?;
    Ok(path.with_file_name(format!("{file_name}.shard{shard_idx:04}.zst")))
}

fn manifest_temp_path(path: &Path) -> Result<PathBuf> {
    let file_name = path
        .file_name()
        .and_then(|x| x.to_str())
        .with_context(|| format!("filter path {} has no valid file name", path.display()))?;
    Ok(path.with_file_name(format!("{file_name}.manifest.tmp")))
}

fn cleanup_shard_files(path: &Path, keep_count: usize) -> Result<()> {
    let parent = path.parent().unwrap_or_else(|| Path::new("."));
    let file_name = path
        .file_name()
        .and_then(|x| x.to_str())
        .with_context(|| format!("filter path {} has no valid file name", path.display()))?;
    let prefix = format!("{file_name}.shard");

    for entry in fs::read_dir(parent)
        .with_context(|| format!("failed to scan {}", parent.display()))?
    {
        let entry = entry?;
        let name = match entry.file_name().into_string() {
            Ok(name) => name,
            Err(_) => continue,
        };
        if !name.starts_with(&prefix) || !name.ends_with(".zst") {
            continue;
        }
        let index_str = &name[prefix.len()..name.len() - 4];
        let Ok(index) = index_str.parse::<usize>() else {
            continue;
        };
        if index >= keep_count {
            fs::remove_file(entry.path())
                .with_context(|| format!("failed to remove stale shard {}", name))?;
        }
    }

    Ok(())
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
fn sampling_hash(key: u64) -> u64 {
    splitmix64(key ^ 0xC6BC279692B5CC8B)
}

#[inline]
fn sampling_hash_bound(sample_rate: f64) -> u64 {
    if sample_rate <= 0.0 {
        0
    } else if sample_rate >= 1.0 {
        u64::MAX
    } else {
        (sample_rate * (u64::MAX as f64)).floor() as u64
    }
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
    use rayon::ThreadPoolBuilder;
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
        let bf = ConcurrentBloomFilter::new(21, 4, 10_000, 1000, 0.01, 1.0, false, 0)
            .expect("concurrent bf init should succeed");
        bf.insert_kmer(123);
        bf.insert_kmer(456);

        let bf = bf.finalize();
        assert!(bf.contains_kmer(123));
        assert!(bf.contains_kmer(456));

        let d = tempdir().expect("tempdir should be creatable");
        let p = d.path().join("t.bf.zst");
        bf.save(&p).expect("save should succeed");

        let loaded = BloomFilter::load(&p).expect("load should succeed");
        assert!(loaded.contains_kmer(123));
        assert!(loaded.contains_kmer(456));
        assert_eq!(loaded.kmer_size, 21);
        assert_eq!(loaded.hash_count, 4);
        assert_eq!(loaded.sample_rate, 1.0);
        assert_eq!(loaded.layout, BloomLayout::Classic);
    }

    #[test]
    fn blocked_bloom_roundtrip() {
        let bf = ConcurrentBloomFilter::new(21, 4, 32_768, 1000, 0.01, 1.0, true, 8)
            .expect("blocked bf init should succeed");
        bf.insert_kmer(123);
        bf.insert_kmer(456);

        let bf = bf.finalize();
        assert!(bf.contains_kmer(123));
        assert!(bf.contains_kmer(456));
        assert_eq!(bf.layout, BloomLayout::Blocked);
        assert_eq!(bf.block_words, 8);

        let d = tempdir().expect("tempdir should be creatable");
        let p = d.path().join("t_blocked.bf.zst");
        bf.save(&p).expect("save should succeed");

        let loaded = BloomFilter::load(&p).expect("load should succeed");
        assert!(loaded.contains_kmer(123));
        assert!(loaded.contains_kmer(456));
        assert_eq!(loaded.layout, BloomLayout::Blocked);
        assert_eq!(loaded.block_words, 8);
    }

    #[test]
    fn large_filters_roundtrip_through_sharded_format() {
        let pool = ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .expect("test pool should build");
        pool.install(|| {
            let bf = ConcurrentBloomFilter::new(21, 4, 1 << 30, 1000, 0.01, 1.0, false, 0)
                .expect("large concurrent bf init should succeed");
            bf.insert_kmer(123);
            bf.insert_kmer(456);

            let bf = bf.finalize();
            let expected_shards = preferred_shard_count(bf.words.len());
            assert!(expected_shards > 1, "test should force sharded persistence");

            let d = tempdir().expect("tempdir should be creatable");
            let p = d.path().join("t_large.bf.zst");
            bf.save(&p).expect("save should succeed");
            assert!(p.exists(), "manifest should exist");
            for shard_idx in 0..expected_shards {
                let shard_path = shard_file_path(&p, shard_idx).expect("shard path should build");
                assert!(shard_path.exists(), "missing shard {}", shard_path.display());
            }

            let loaded = BloomFilter::load(&p).expect("load should succeed");
            assert!(loaded.contains_kmer(123));
            assert!(loaded.contains_kmer(456));
            assert_eq!(loaded.layout, BloomLayout::Classic);
        });
    }
}
