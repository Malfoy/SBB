use crate::bloomybloom_vendor::{
    BlockShard, ConcurrentBloomFilter as RawConcurrentBloomFilter,
    FrozenBlockShard, FrozenBloomFilter as RawFrozenBloomFilter, export_filter_words,
    import_filter_from_shards, import_filter_from_words, into_frozen,
};
use crate::decyclers::Decycler;
use crate::minimizers::{MinimizerMode, selected_mins_x_pos};
use crate::utils::{hash_u128_to_u64, roll_u128_kmer, xorshift_u64};
use anyhow::{Result, bail};
use packed_seq::{PackedSeq, PackedSeqVec, Seq, SeqVec};
use rayon::prelude::*;

const SHARD_COUNT: usize = 1024;
const DEFAULT_BLOOMY_B: u16 = 4;
const DEFAULT_BLOOMY_BLOCK_SIZE: usize = 1 << 12;
pub(crate) const MIN_BLOOMY_BITS: u64 = (SHARD_COUNT * DEFAULT_BLOOMY_BLOCK_SIZE) as u64;
const BLOOMY_PAR_MIN_STARTS: usize = 1_000_000;
const BLOOMY_PAR_CHUNK_STARTS: usize = 1_000_000;

#[derive(Debug, Clone, Copy)]
pub(crate) struct BloomyParams {
    pub(crate) hash_count: u8,
    pub(crate) bit_len: u64,
    pub(crate) m: u16,
    pub(crate) s: u16,
    pub(crate) block_size: u64,
    pub(crate) nb_blocks: u64,
    pub(crate) mode_tag: u8,
    pub(crate) open_closed_t: u16,
}

pub(crate) struct BloomyConcurrentState {
    raw: RawConcurrentBloomFilter,
    k: u16,
    m: u16,
    s: u16,
    mode: MinimizerMode,
    decycler: Decycler,
}

pub(crate) struct BloomyFrozenState {
    raw: RawFrozenBloomFilter,
    k: u16,
    m: u16,
    s: u16,
    mode: MinimizerMode,
    decycler: Decycler,
}

impl std::fmt::Debug for BloomyConcurrentState {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BloomyConcurrentState")
            .field("k", &self.k)
            .field("m", &self.m)
            .field("s", &self.s)
            .field("mode", &self.mode)
            .field("block_size", &self.raw.block_size)
            .field("nb_blocks", &self.raw.nb_blocks)
            .field("n_hashes", &self.raw.n_hashes)
            .finish()
    }
}

impl std::fmt::Debug for BloomyFrozenState {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BloomyFrozenState")
            .field("k", &self.k)
            .field("m", &self.m)
            .field("s", &self.s)
            .field("mode", &self.mode)
            .field("block_size", &self.raw.block_size)
            .field("nb_blocks", &self.raw.nb_blocks)
            .field("n_hashes", &self.raw.n_hashes)
            .finish()
    }
}

impl Clone for BloomyFrozenState {
    fn clone(&self) -> Self {
        Self::from_parts(
            export_filter_words(&self.raw),
            self.k,
            self.m,
            self.s,
            self.raw.block_size,
            self.raw.nb_blocks,
            self.raw.n_hashes,
            self.mode,
        )
        .expect("cloning bloomy filter should succeed")
    }
}

fn default_smer_length(k: u16) -> u16 {
    k.saturating_sub(10).max(1)
}

fn default_minimizer_length(nb_blocks: usize) -> u16 {
    ((nb_blocks.ilog2() as u16) / 2) + DEFAULT_BLOOMY_B
}

fn resolve_smer_length(k: u16) -> Result<u16> {
    let Some(raw) = std::env::var_os("SBB_BLOOMY_SMER") else {
        return Ok(default_smer_length(k));
    };

    let raw = raw
        .to_str()
        .ok_or_else(|| anyhow::anyhow!("SBB_BLOOMY_SMER must be valid UTF-8"))?;
    let s: u16 = raw
        .parse()
        .map_err(|_| anyhow::anyhow!("SBB_BLOOMY_SMER must be an integer, got {}", raw))?;
    if s == 0 || s > k {
        bail!("SBB_BLOOMY_SMER must be in [1, {}], got {}", k, s);
    }
    Ok(s)
}

fn resolve_minimizer_length(k: u16, nb_blocks: usize) -> Result<u16> {
    let default_m = default_minimizer_length(nb_blocks);
    let Some(raw) = std::env::var_os("SBB_BLOOMY_M") else {
        return Ok(default_m);
    };

    let raw = raw
        .to_str()
        .ok_or_else(|| anyhow::anyhow!("SBB_BLOOMY_M must be valid UTF-8"))?;
    let m: u16 = raw
        .parse()
        .map_err(|_| anyhow::anyhow!("SBB_BLOOMY_M must be an integer, got {}", raw))?;
    if m == 0 || m >= 32 || m > k {
        bail!("SBB_BLOOMY_M must be in [1, min(31, {})], got {}", k, m);
    }
    Ok(m)
}

fn max_superkmer_kmers(k: u16, m: u16) -> usize {
    k as usize - m as usize + 1
}

fn superkmer_smer_count(k: u16, m: u16, s: u16) -> usize {
    max_superkmer_kmers(k, m) + k as usize - s as usize
}

fn optimal_hash_count(block_size: usize, k: u16, m: u16, s: u16) -> usize {
    let smer_count = superkmer_smer_count(k, m, s) as f64;
    ((block_size as f64 * std::f64::consts::LN_2) / smer_count)
        .round()
        .max(1.0) as usize
}

fn default_mode() -> MinimizerMode {
    MinimizerMode::Simd
}

fn mode_to_header(mode: MinimizerMode) -> (u8, u16) {
    match mode {
        MinimizerMode::Simd => (0, 0),
        MinimizerMode::Decycling => (1, 0),
        MinimizerMode::OpenClosed { t } => (2, t),
    }
}

fn mode_from_header(tag: u8, t: u16) -> Result<MinimizerMode> {
    match tag {
        0 => Ok(MinimizerMode::Simd),
        1 => Ok(MinimizerMode::Decycling),
        2 => Ok(MinimizerMode::OpenClosed { t }),
        _ => bail!("unsupported bloomy minimizer mode {}", tag),
    }
}

fn build_decycler(m: u16, mode: MinimizerMode) -> Decycler {
    if matches!(mode, MinimizerMode::Decycling) {
        let mut decycler = Decycler::new(m);
        decycler.compute_blocks();
        decycler
    } else {
        Decycler::new(1)
    }
}

fn is_valid_base(b: u8) -> bool {
    matches!(b, b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't' | b'U' | b'u')
}

fn for_each_valid_slice(seq: &[u8], mut f: impl FnMut(&[u8])) {
    let mut start = None;
    for (idx, &b) in seq.iter().enumerate() {
        if is_valid_base(b) {
            if start.is_none() {
                start = Some(idx);
            }
        } else if let Some(slice_start) = start.take() {
            f(&seq[slice_start..idx]);
        }
    }
    if let Some(slice_start) = start {
        f(&seq[slice_start..]);
    }
}

fn count_valid_kmers(seq: &[u8], k: usize) -> u64 {
    let mut total = 0_u64;
    for_each_valid_slice(seq, |slice| {
        if slice.len() >= k {
            total += (slice.len() + 1 - k) as u64;
        }
    });
    total
}

impl BloomyConcurrentState {
    pub(crate) fn new(k: u8, requested_hash_count: Option<u8>, bit_len: u64) -> Result<(Self, BloomyParams)> {
        if k % 2 == 0 {
            bail!("bloomybloom requires an odd k-mer size");
        }
        let bit_len = bit_len.max(MIN_BLOOMY_BITS).next_power_of_two();
        let block_size = DEFAULT_BLOOMY_BLOCK_SIZE;
        let nb_blocks = (bit_len as usize) / block_size;
        if nb_blocks < SHARD_COUNT {
            bail!("bloomybloom requires at least {} blocks", SHARD_COUNT);
        }

        let m = resolve_minimizer_length(u16::from(k), nb_blocks)?;
        if m == 0 || m >= 32 || m > u16::from(k) {
            bail!("invalid bloomy minimizer length {}", m);
        }
        let s = resolve_smer_length(u16::from(k))?;
        let hash_count =
            requested_hash_count.map(usize::from).unwrap_or_else(|| optimal_hash_count(block_size, u16::from(k), m, s));
        let mode = default_mode();
        let raw = RawConcurrentBloomFilter::new(bit_len as usize, hash_count, usize::from(k), block_size, nb_blocks);

        let state = Self {
            raw,
            k: u16::from(k),
            m,
            s,
            mode,
            decycler: build_decycler(m, mode),
        };
        let (mode_tag, open_closed_t) = mode_to_header(mode);
        let params = BloomyParams {
            hash_count: hash_count as u8,
            bit_len,
            m,
            s,
            block_size: block_size as u64,
            nb_blocks: nb_blocks as u64,
            mode_tag,
            open_closed_t,
        };
        Ok((state, params))
    }

    pub(crate) fn insert_sequence(&self, seq: &[u8]) -> u64 {
        let mut inserted = 0_u64;
        let k = self.k as usize;
        for_each_valid_slice(seq, |slice| {
            if slice.len() < k {
                return;
            }
            let starts = slice.len() + 1 - k;
            inserted += starts as u64;

            if starts < BLOOMY_PAR_MIN_STARTS || rayon::current_num_threads() <= 1 {
                handle_sequence(
                    &self.raw,
                    PackedSeqVec::from_ascii(slice),
                    self.k,
                    self.m,
                    self.raw.nb_blocks,
                    &self.decycler,
                    self.mode,
                    self.s,
                );
                return;
            }

            let chunk_count = starts.div_ceil(BLOOMY_PAR_CHUNK_STARTS);
            (0..chunk_count).into_par_iter().for_each(|chunk_idx| {
                // Partition by k-mer start positions so each k-mer is inserted exactly once.
                let start_kmer = chunk_idx * BLOOMY_PAR_CHUNK_STARTS;
                let end_kmer = (start_kmer + BLOOMY_PAR_CHUNK_STARTS).min(starts);
                let chunk = &slice[start_kmer..end_kmer + k - 1];
                handle_sequence(
                    &self.raw,
                    PackedSeqVec::from_ascii(chunk),
                    self.k,
                    self.m,
                    self.raw.nb_blocks,
                    &self.decycler,
                    self.mode,
                    self.s,
                );
            });
        });
        inserted
    }

    pub(crate) fn score_sequence_with_counts(&self, seq: &[u8]) -> (u64, u64) {
        let mut hits = 0_u64;
        let mut total = 0_u64;
        for_each_valid_slice(seq, |slice| {
            if slice.len() < self.k as usize {
                return;
            }
            let sequence = PackedSeqVec::from_ascii(slice);
            let presence = score_sequence_on_concurrent(
                &self.raw,
                sequence,
                self.k,
                self.m,
                self.s,
                &self.decycler,
                self.mode,
            );
            total += presence.len() as u64;
            hits += presence.iter().filter(|&&present| present).count() as u64;
        });
        (hits, total)
    }

    pub(crate) fn finalize(self) -> BloomyFrozenState {
        let raw = into_frozen(self.raw.into_frozen());
        BloomyFrozenState {
            raw,
            k: self.k,
            m: self.m,
            s: self.s,
            mode: self.mode,
            decycler: build_decycler(self.m, self.mode),
        }
    }
}

fn score_sequence_on_concurrent(
    bloom: &RawConcurrentBloomFilter,
    original_sequence: PackedSeqVec,
    k: u16,
    m: u16,
    s: u16,
    decycler_set: &Decycler,
    minimizer_mode: MinimizerMode,
) -> Vec<bool> {
    if original_sequence.len() < k as usize {
        return vec![];
    }

    let mut presence_vec = Vec::with_capacity(original_sequence.len() - k as usize + 1);
    let address_mask = (bloom.nb_blocks - 1) >> 10;
    let (super_kmers_positions, minimizers, sequence): (Vec<u32>, Vec<u64>, PackedSeqVec) =
        selected_mins_x_pos(original_sequence, k, m, decycler_set, minimizer_mode);

    for i in 0..super_kmers_positions.len() {
        let hashed_minimizer = xorshift_u64(minimizers[i]);
        let start_pos = super_kmers_positions[i] as usize;
        let end_pos = if i == super_kmers_positions.len() - 1 {
            sequence.len() + 1 - k as usize
        } else {
            super_kmers_positions[i + 1] as usize
        };
        let blocknum = (hashed_minimizer as usize) % SHARD_COUNT;
        let subblocknum = ((hashed_minimizer as usize) >> 10) & address_mask;
        let block = bloom.filter[blocknum].lock().unwrap();

        if s <= 31 {
            for j in start_pos..end_pos {
                let kmer: PackedSeq = sequence.slice(j..j + k as usize);
                presence_vec.push(check_packed_kmer(&block, subblocknum, kmer, s, bloom.block_size - 1, bloom.n_hashes));
            }
        } else {
            check_super_kmer_u128_concurrent(
                &block,
                subblocknum,
                &sequence,
                start_pos,
                end_pos,
                k,
                s,
                bloom.n_hashes,
                bloom.block_size - 1,
                &mut presence_vec,
            );
        }
    }

    presence_vec
}

fn check_packed_kmer(
    block: &BlockShard,
    subblock: usize,
    kmer: PackedSeq,
    s: u16,
    block_size_mask: usize,
    n_hashes: usize,
) -> bool {
    for i in 0..kmer.len() - s as usize + 1 {
        let smer = kmer.slice(i..i + s as usize);
        let mut hash = xorshift_u64(smer.as_u64());
        for _ in 0..n_hashes {
            let address = hash as usize & block_size_mask;
            if !block.get(subblock, address) {
                return false;
            }
            hash = xorshift_u64(hash);
        }
    }
    true
}

fn check_super_kmer_u128_concurrent(
    block: &BlockShard,
    subblock: usize,
    sequence: &PackedSeqVec,
    start_pos: usize,
    end_pos: usize,
    k: u16,
    s: u16,
    n_hashes: usize,
    block_size_mask: usize,
    presence_vec: &mut Vec<bool>,
) {
    let window_size = (k - s + 1) as usize;
    let smer_count = end_pos + window_size - 1 - start_pos;
    let mut current_smer = sequence.slice(start_pos..start_pos + s as usize).as_u128();
    let mut rolling_presence = vec![true; window_size];
    let mut missing_in_window = 0_usize;

    for offset in 0..smer_count {
        if offset >= window_size && !rolling_presence[offset % window_size] {
            missing_in_window -= 1;
        }

        let mut hash = hash_u128_to_u64(current_smer);
        let mut smer_present = true;
        for _ in 0..n_hashes {
            let address = hash as usize & block_size_mask;
            if !block.get(subblock, address) {
                smer_present = false;
                break;
            }
            hash = xorshift_u64(hash);
        }

        rolling_presence[offset % window_size] = smer_present;
        if !smer_present {
            missing_in_window += 1;
        }
        if offset + 1 >= window_size {
            presence_vec.push(missing_in_window == 0);
        }
        if offset + 1 < smer_count {
            let next_base = sequence.as_slice().get(start_pos + offset + s as usize);
            current_smer = roll_u128_kmer(current_smer, next_base, s);
        }
    }
}

impl BloomyFrozenState {
    pub(crate) fn from_parts(
        words: Vec<u64>,
        k: u16,
        m: u16,
        s: u16,
        block_size: usize,
        nb_blocks: usize,
        n_hashes: usize,
        mode: MinimizerMode,
    ) -> Result<Self> {
        let raw = import_filter_from_words(words, block_size, nb_blocks, n_hashes);
        Ok(Self {
            raw,
            k,
            m,
            s,
            mode,
            decycler: build_decycler(m, mode),
        })
    }

    pub(crate) fn from_shard_words(
        shard_words: Vec<Vec<u64>>,
        k: u8,
        params: BloomyParams,
    ) -> Result<Self> {
        let mode = mode_from_header(params.mode_tag, params.open_closed_t)?;
        let raw = import_filter_from_shards(
            shard_words,
            params.block_size as usize,
            params.nb_blocks as usize,
            params.hash_count as usize,
        );
        Ok(Self {
            raw,
            k: u16::from(k),
            m: params.m,
            s: params.s,
            mode,
            decycler: build_decycler(params.m, mode),
        })
    }

    pub(crate) fn from_header(words: Vec<u64>, k: u8, params: BloomyParams) -> Result<Self> {
        let mode = mode_from_header(params.mode_tag, params.open_closed_t)?;
        Self::from_parts(
            words,
            u16::from(k),
            params.m,
            params.s,
            params.block_size as usize,
            params.nb_blocks as usize,
            params.hash_count as usize,
            mode,
        )
    }

    pub(crate) fn score_sequence_with_counts(&self, seq: &[u8]) -> (u64, u64) {
        let mut hits = 0_u64;
        let mut total = 0_u64;
        for_each_valid_slice(seq, |slice| {
            if slice.len() < self.k as usize {
                return;
            }
            let sequence = PackedSeqVec::from_ascii(slice);
            let presence = score_sequence_on_frozen(
                &self.raw,
                sequence,
                self.k,
                self.m,
                self.s,
                &self.decycler,
                self.mode,
            );
            total += presence.len() as u64;
            hits += presence.iter().filter(|&&present| present).count() as u64;
        });
        (hits, total)
    }

    pub(crate) fn match_sequence_threshold_with_counts(
        &self,
        seq: &[u8],
        threshold: f64,
        best_hit: bool,
    ) -> (bool, u64) {
        let k = self.k as usize;
        if best_hit {
            let mut queried = 0_u64;
            let mut matched = false;
            for_each_valid_slice(seq, |slice| {
                if matched || slice.len() < k {
                    return;
                }
                let sequence = PackedSeqVec::from_ascii(slice);
                let _ = for_each_presence_on_frozen_until(
                    &self.raw,
                    sequence,
                    self.k,
                    self.m,
                    self.s,
                    &self.decycler,
                    self.mode,
                    |present| {
                        queried += 1;
                        if present {
                            matched = true;
                            false
                        } else {
                            true
                        }
                    },
                );
            });
            return (matched, queried);
        }

        if threshold <= 0.0 {
            return (true, count_valid_kmers(seq, k));
        }

        if threshold >= 1.0 {
            let mut queried = 0_u64;
            let mut missed = false;
            for_each_valid_slice(seq, |slice| {
                if missed || slice.len() < k {
                    return;
                }
                let sequence = PackedSeqVec::from_ascii(slice);
                let _ = for_each_presence_on_frozen_until(
                    &self.raw,
                    sequence,
                    self.k,
                    self.m,
                    self.s,
                    &self.decycler,
                    self.mode,
                    |present| {
                        queried += 1;
                        if present {
                            true
                        } else {
                            missed = true;
                            false
                        }
                    },
                );
            });
            return (queried > 0 && !missed, queried);
        }

        let total = count_valid_kmers(seq, k);
        if total == 0 {
            return (false, 0);
        }

        let mut queried = 0_u64;
        let mut hits = 0_u64;
        let mut decided: Option<bool> = None;
        for_each_valid_slice(seq, |slice| {
            if decided.is_some() || slice.len() < k {
                return;
            }
            let sequence = PackedSeqVec::from_ascii(slice);
            let _ = for_each_presence_on_frozen_until(
                &self.raw,
                sequence,
                self.k,
                self.m,
                self.s,
                &self.decycler,
                self.mode,
                |present| {
                    queried += 1;
                    if present {
                        hits += 1;
                    }
                    let remaining = total - queried;
                    if (hits as f64 / total as f64) >= threshold {
                        decided = Some(true);
                        return false;
                    }
                    if ((hits + remaining) as f64 / total as f64) < threshold {
                        decided = Some(false);
                        return false;
                    }
                    true
                },
            );
        });

        if let Some(matched) = decided {
            return (matched, queried);
        }
        ((hits as f64 / total as f64) >= threshold, queried)
    }

    pub(crate) fn serialized_word_len(&self) -> usize {
        self.raw.word_len()
    }

    pub(crate) fn write_word_range(
        &self,
        range: std::ops::Range<usize>,
        w: &mut impl std::io::Write,
    ) -> Result<()> {
        let mut cursor = 0_usize;
        for shard in &self.raw.filter {
            let shard_words = shard.as_slice();
            let shard_start = cursor;
            let shard_end = cursor + shard_words.len();
            let overlap_start = range.start.max(shard_start);
            let overlap_end = range.end.min(shard_end);
            if overlap_start < overlap_end {
                let local_start = overlap_start - shard_start;
                let local_end = overlap_end - shard_start;
                for word in &shard_words[local_start..local_end] {
                    w.write_all(&word.to_le_bytes())?;
                }
            }
            if shard_end >= range.end {
                break;
            }
            cursor = shard_end;
        }
        Ok(())
    }
}

fn for_each_presence_on_frozen_until(
    bloom: &RawFrozenBloomFilter,
    original_sequence: PackedSeqVec,
    k: u16,
    m: u16,
    s: u16,
    decycler_set: &Decycler,
    minimizer_mode: MinimizerMode,
    mut f: impl FnMut(bool) -> bool,
) -> usize {
    if original_sequence.len() < k as usize {
        return 0;
    }

    let mut count = 0_usize;
    let address_mask = (bloom.nb_blocks - 1) >> 10;
    let (super_kmers_positions, minimizers, sequence): (Vec<u32>, Vec<u64>, PackedSeqVec) =
        selected_mins_x_pos(original_sequence, k, m, decycler_set, minimizer_mode);

    for i in 0..super_kmers_positions.len() {
        let hashed_minimizer = xorshift_u64(minimizers[i]);
        let start_pos = super_kmers_positions[i] as usize;
        let end_pos = if i == super_kmers_positions.len() - 1 {
            sequence.len() + 1 - k as usize
        } else {
            super_kmers_positions[i + 1] as usize
        };
        let blocknum = (hashed_minimizer as usize) % SHARD_COUNT;
        let subblocknum = ((hashed_minimizer as usize) >> 10) & address_mask;
        let block = &bloom.filter[blocknum];

        if s <= 31 {
            for j in start_pos..end_pos {
                let kmer: PackedSeq = sequence.slice(j..j + k as usize);
                let present =
                    check_packed_kmer_frozen(block, subblocknum, kmer, s, bloom.block_size_mask, bloom.n_hashes);
                count += 1;
                if !f(present) {
                    return count;
                }
            }
        } else {
            let mut presence = Vec::with_capacity(end_pos.saturating_sub(start_pos));
            check_super_kmer_u128_frozen(
                block,
                subblocknum,
                &sequence,
                start_pos,
                end_pos,
                k,
                s,
                bloom.n_hashes,
                bloom.block_size_mask,
                &mut presence,
            );
            for present in presence {
                count += 1;
                if !f(present) {
                    return count;
                }
            }
        }
    }

    count
}

fn score_sequence_on_frozen(
    bloom: &RawFrozenBloomFilter,
    original_sequence: PackedSeqVec,
    k: u16,
    m: u16,
    s: u16,
    decycler_set: &Decycler,
    minimizer_mode: MinimizerMode,
) -> Vec<bool> {
    if original_sequence.len() < k as usize {
        return vec![];
    }

    let mut presence_vec = Vec::with_capacity(original_sequence.len() - k as usize + 1);
    let address_mask = (bloom.nb_blocks - 1) >> 10;
    let (super_kmers_positions, minimizers, sequence): (Vec<u32>, Vec<u64>, PackedSeqVec) =
        selected_mins_x_pos(original_sequence, k, m, decycler_set, minimizer_mode);

    for i in 0..super_kmers_positions.len() {
        let hashed_minimizer = xorshift_u64(minimizers[i]);
        let start_pos = super_kmers_positions[i] as usize;
        let end_pos = if i == super_kmers_positions.len() - 1 {
            sequence.len() + 1 - k as usize
        } else {
            super_kmers_positions[i + 1] as usize
        };
        let blocknum = (hashed_minimizer as usize) % SHARD_COUNT;
        let subblocknum = ((hashed_minimizer as usize) >> 10) & address_mask;
        let block = &bloom.filter[blocknum];

        if s <= 31 {
            for j in start_pos..end_pos {
                let kmer: PackedSeq = sequence.slice(j..j + k as usize);
                presence_vec.push(check_packed_kmer_frozen(
                    block,
                    subblocknum,
                    kmer,
                    s,
                    bloom.block_size_mask,
                    bloom.n_hashes,
                ));
            }
        } else {
            check_super_kmer_u128_frozen(
                block,
                subblocknum,
                &sequence,
                start_pos,
                end_pos,
                k,
                s,
                bloom.n_hashes,
                bloom.block_size_mask,
                &mut presence_vec,
            );
        }
    }

    presence_vec
}

fn check_packed_kmer_frozen(
    block: &FrozenBlockShard,
    subblock: usize,
    kmer: PackedSeq,
    s: u16,
    block_size_mask: usize,
    n_hashes: usize,
) -> bool {
    for i in 0..kmer.len() - s as usize + 1 {
        let smer = kmer.slice(i..i + s as usize);
        let mut hash = xorshift_u64(smer.as_u64());
        for _ in 0..n_hashes {
            let address = hash as usize & block_size_mask;
            if !block.get(subblock, address) {
                return false;
            }
            hash = xorshift_u64(hash);
        }
    }
    true
}

fn check_super_kmer_u128_frozen(
    block: &FrozenBlockShard,
    subblock: usize,
    sequence: &PackedSeqVec,
    start_pos: usize,
    end_pos: usize,
    k: u16,
    s: u16,
    n_hashes: usize,
    block_size_mask: usize,
    presence_vec: &mut Vec<bool>,
) {
    let window_size = (k - s + 1) as usize;
    let smer_count = end_pos + window_size - 1 - start_pos;
    let mut current_smer = sequence.slice(start_pos..start_pos + s as usize).as_u128();
    let mut rolling_presence = vec![true; window_size];
    let mut missing_in_window = 0_usize;

    for offset in 0..smer_count {
        if offset >= window_size && !rolling_presence[offset % window_size] {
            missing_in_window -= 1;
        }

        let mut hash = hash_u128_to_u64(current_smer);
        let mut smer_present = true;
        for _ in 0..n_hashes {
            let address = hash as usize & block_size_mask;
            if !block.get(subblock, address) {
                smer_present = false;
                break;
            }
            hash = xorshift_u64(hash);
        }

        rolling_presence[offset % window_size] = smer_present;
        if !smer_present {
            missing_in_window += 1;
        }
        if offset + 1 >= window_size {
            presence_vec.push(missing_in_window == 0);
        }
        if offset + 1 < smer_count {
            let next_base = sequence.as_slice().get(start_pos + offset + s as usize);
            current_smer = roll_u128_kmer(current_smer, next_base, s);
        }
    }
}

fn handle_sequence(
    bloom: &RawConcurrentBloomFilter,
    original_sequence: PackedSeqVec,
    k: u16,
    m: u16,
    nb_blocks: usize,
    decycler_set: &Decycler,
    minimizer_mode: MinimizerMode,
    s: u16,
) {
    if original_sequence.len() < k as usize {
        return;
    }

    let address_mask: usize = bloom.block_size - 1;
    let (super_kmers_positions, minimizer_values, sequence): (Vec<u32>, Vec<u64>, PackedSeqVec) =
        selected_mins_x_pos(original_sequence, k, m, decycler_set, minimizer_mode);

    let mut all_addresses = Vec::new();
    let mut kmer_number: usize = 0;
    for i in 0..super_kmers_positions.len().saturating_sub(1) {
        let hashed_minimizer: u64 = xorshift_u64(minimizer_values[i]) & (nb_blocks as u64 - 1);
        if s <= 31 {
            kmer_number = handle_super_kmer(
                super_kmers_positions[i],
                super_kmers_positions[i + 1],
                &sequence,
                bloom,
                k,
                hashed_minimizer,
                kmer_number,
                &mut all_addresses,
                address_mask,
                s,
            );
        } else {
            kmer_number = handle_super_kmer_u128(
                super_kmers_positions[i],
                super_kmers_positions[i + 1],
                &sequence,
                bloom,
                k,
                hashed_minimizer,
                kmer_number,
                &mut all_addresses,
                address_mask,
                s,
            );
        }
    }

    if !super_kmers_positions.is_empty() {
        let last = super_kmers_positions.len() - 1;
        let hashed_minimizer: u64 = xorshift_u64(minimizer_values[last]) & (nb_blocks as u64 - 1);
        if s <= 31 {
            let _ = handle_super_kmer(
                super_kmers_positions[last],
                (sequence.len() + 1 - k as usize) as u32,
                &sequence,
                bloom,
                k,
                hashed_minimizer,
                kmer_number,
                &mut all_addresses,
                address_mask,
                s,
            );
        } else {
            let _ = handle_super_kmer_u128(
                super_kmers_positions[last],
                (sequence.len() + 1 - k as usize) as u32,
                &sequence,
                bloom,
                k,
                hashed_minimizer,
                kmer_number,
                &mut all_addresses,
                address_mask,
                s,
            );
        }
    }
}

fn handle_super_kmer(
    start_pos: u32,
    end_pos: u32,
    sequence: &PackedSeqVec,
    bloom: &RawConcurrentBloomFilter,
    k: u16,
    hashed_minimizer: u64,
    mut kmer_number: usize,
    all_addresses: &mut Vec<usize>,
    address_mask: usize,
    s: u16,
) -> usize {
    let smer_count = end_pos as usize + (k - s) as usize - start_pos as usize;
    let required_addresses = smer_count * bloom.n_hashes;
    if all_addresses.len() < required_addresses {
        all_addresses.resize(required_addresses, 0);
    }
    let mut last_relevant_index: usize = 0;
    for j in (start_pos as usize)..(end_pos as usize) + (k - s) as usize {
        let smer: PackedSeq = sequence.slice(j..j + s as usize);
        let mut hash: u64 = xorshift_u64(smer.as_u64());
        for _ in 0..bloom.n_hashes {
            let address = hash as usize & address_mask;
            all_addresses[last_relevant_index] = address;
            last_relevant_index += 1;
            hash = xorshift_u64(hash);
        }
        kmer_number += 1;
    }

    let relevant_addresses = &mut all_addresses[..last_relevant_index];
    let blocknum: usize = (hashed_minimizer as usize) & 1023;
    let subblocknum: usize = ((hashed_minimizer as usize) >> 10) & ((bloom.nb_blocks >> 10) - 1);
    let mut block = bloom.filter[blocknum].lock().unwrap();
    for address in relevant_addresses {
        if !block.get(subblocknum, *address) {
            block.set(subblocknum, *address, true);
        }
    }
    drop(block);
    kmer_number
}

fn handle_super_kmer_u128(
    start_pos: u32,
    end_pos: u32,
    sequence: &PackedSeqVec,
    bloom: &RawConcurrentBloomFilter,
    k: u16,
    hashed_minimizer: u64,
    mut kmer_number: usize,
    _all_addresses: &mut Vec<usize>,
    address_mask: usize,
    s: u16,
) -> usize {
    let smer_count = end_pos as usize + (k - s) as usize - start_pos as usize;
    let blocknum: usize = (hashed_minimizer as usize) & 1023;
    let subblocknum: usize = ((hashed_minimizer as usize) >> 10) & ((bloom.nb_blocks >> 10) - 1);
    let mut block = bloom.filter[blocknum].lock().unwrap();

    let mut current_smer = sequence
        .slice(start_pos as usize..start_pos as usize + s as usize)
        .as_u128();
    for offset in 0..smer_count {
        let mut hash = hash_u128_to_u64(current_smer);
        for _ in 0..bloom.n_hashes {
            let address = hash as usize & address_mask;
            if !block.get(subblocknum, address) {
                block.set(subblocknum, address, true);
            }
            hash = xorshift_u64(hash);
        }

        kmer_number += 1;
        if offset + 1 < smer_count {
            let next_base = sequence.as_slice().get(start_pos as usize + offset + s as usize);
            current_smer = roll_u128_kmer(current_smer, next_base, s);
        }
    }
    kmer_number
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bloomybloom_vendor::export_filter_words;
    use rayon::ThreadPoolBuilder;

    fn patterned_sequence(len: usize) -> Vec<u8> {
        const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
        (0..len).map(|i| BASES[i & 3]).collect()
    }

    #[test]
    fn parallel_chunked_insert_matches_serial_bits() {
        let k = 25_u8;
        let bit_len = 1_u64 << 24;
        let hash_count = Some(3_u8);
        let seq = patterned_sequence(BLOOMY_PAR_MIN_STARTS + usize::from(k) + 65_536);

        let serial_pool = ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .expect("serial pool should build");
        let parallel_pool = ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .expect("parallel pool should build");

        let (serial_state, _) =
            BloomyConcurrentState::new(k, hash_count, bit_len).expect("serial state should build");
        let (parallel_state, _) =
            BloomyConcurrentState::new(k, hash_count, bit_len).expect("parallel state should build");

        let serial_inserted = serial_pool.install(|| serial_state.insert_sequence(&seq));
        let parallel_inserted = parallel_pool.install(|| parallel_state.insert_sequence(&seq));
        assert_eq!(
            serial_inserted, parallel_inserted,
            "parallel chunking should report the same inserted k-mer count"
        );

        let serial_words = export_filter_words(&serial_state.finalize().raw);
        let parallel_words = export_filter_words(&parallel_state.finalize().raw);
        assert_eq!(
            serial_words, parallel_words,
            "parallel chunking should produce the same bloom bits as serial insertion"
        );
    }
}
