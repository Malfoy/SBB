/// Utilities to compute super-kmer boundaries and minimizer values for different
/// grouping orders.

use crate::decyclers::Decycler;

use std::cell::RefCell;
use std::collections::VecDeque;

use packed_seq::{PackedSeq, PackedSeqVec, Seq, SeqVec};
use seq_hash::{KmerHasher, NtHasher};
use simd_minimizers::canonical_minimizers;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum MinimizerMode {
    Simd,
    Decycling,
    OpenClosed { t: u16 },
}

#[derive(Default)]
struct SimdMinimizerScratch {
    minimizer_positions: Vec<u32>,
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct OpenClosedKey {
    priority: u8,
    hash: u64,
}

thread_local! {
    static SIMD_MINIMIZER_SCRATCH: RefCell<SimdMinimizerScratch> =
        RefCell::new(SimdMinimizerScratch::default());
}

#[inline(always)]
fn open_closed_priority(m: u16, t: u16, min_t_offset: usize) -> u8 {
    let offset = (m - t) as usize / 2;
    let last = (m - t) as usize;
    if min_t_offset == offset {
        0
    } else if min_t_offset == 0 || min_t_offset == last {
        1
    } else {
        2
    }
}

#[inline(always)]
fn build_superkmers_from_positions(
    packed_seq: &PackedSeqVec,
    k: u16,
    m: u16,
    minimizer_positions: &[usize],
    minimizer_hashes: &[u64],
) -> (Vec<u32>, Vec<u64>) {
    debug_assert_eq!(minimizer_positions.len(), packed_seq.len() + 1 - k as usize);

    let mut super_kmers = Vec::new();
    let mut minimizer_vals = Vec::new();
    let mut first_superkmer_addr: u32 = 0;
    let mut last_minimizer_pos = minimizer_positions[0];

    for i in 1..minimizer_positions.len() {
        if minimizer_positions[i] != last_minimizer_pos {
            last_minimizer_pos = minimizer_positions[i];
            super_kmers.push(first_superkmer_addr);
            minimizer_vals.push(minimizer_hashes[minimizer_positions[first_superkmer_addr as usize]]);
            first_superkmer_addr = i as u32;
        }
    }

    super_kmers.push(first_superkmer_addr);
    minimizer_vals.push(minimizer_hashes[minimizer_positions[first_superkmer_addr as usize]]);

    debug_assert!(super_kmers.iter().all(|&p| p < (packed_seq.len() + 1 - k as usize) as u32));
    debug_assert_eq!(super_kmers.len(), minimizer_vals.len());
    debug_assert!(m <= k);
    (super_kmers, minimizer_vals)
}

/// Dispatch to the requested minimizer order while preserving the existing
/// `(positions, values, sequence)` interface used by the bloom code.
pub fn selected_mins_x_pos(
    packed_seq: PackedSeqVec,
    k: u16,
    m: u16,
    decycler_set: &Decycler,
    mode: MinimizerMode,
) -> (Vec<u32>, Vec<u64>, PackedSeqVec) {
    match mode {
        MinimizerMode::Simd => minimizers_x_positions(packed_seq, k, m),
        MinimizerMode::Decycling => decycling_mins_x_pos(packed_seq, k, m, decycler_set),
        MinimizerMode::OpenClosed { t } => open_closed_mins_x_pos(packed_seq, k, m, t),
    }
}

/// Standard SIMD minimizers.
pub fn minimizers_x_positions(
    packed_seq: PackedSeqVec,
    k: u16,
    m: u16,
) -> (Vec<u32>, Vec<u64>, PackedSeqVec) {
    let max_windows = packed_seq.len().saturating_sub(k as usize) + 1;
    let mut super_kmers = Vec::with_capacity(max_windows);
    let window_size = k - m + 1;
    let minimizer_length = m;

    let minimizer_vals = SIMD_MINIMIZER_SCRATCH.with(|scratch| {
        let mut scratch = scratch.borrow_mut();
        scratch.minimizer_positions.clear();
        let current_capacity = scratch.minimizer_positions.capacity();
        scratch
            .minimizer_positions
            .reserve(max_windows.saturating_sub(current_capacity));

        let seed: u32 = 42;
        let hasher = <NtHasher>::new_with_seed(minimizer_length.into(), seed);
        let output = canonical_minimizers(minimizer_length.into(), window_size.into())
            .hasher(&hasher)
            .super_kmers(&mut super_kmers)
            .run(packed_seq.as_slice(), &mut scratch.minimizer_positions);
        output.values_u64().collect()
    });

    (super_kmers, minimizer_vals, packed_seq)
}

/// Open-closed minimizers using the priority order from jermp/minimizers:
/// open syncmers first, then closed syncmers, then ordinary m-mers.
pub fn open_closed_mins_x_pos(
    packed_seq: PackedSeqVec,
    k: u16,
    m: u16,
    t: u16,
) -> (Vec<u32>, Vec<u64>, PackedSeqVec) {
    if packed_seq.len() < k as usize {
        return (Vec::new(), Vec::new(), packed_seq);
    }

    assert!(t > 0);
    assert!(t <= m);

    let num_kmers = packed_seq.len() + 1 - k as usize;
    let num_mmers = packed_seq.len() + 1 - m as usize;
    let outer_window = (k - m + 1) as usize;
    let inner_window = (m - t + 1) as usize;

    let m_hasher = <NtHasher>::new_with_seed(m as usize, 42);
    let t_hasher = <NtHasher>::new_with_seed(t as usize, 42);
    let m_hashes: Vec<u64> = m_hasher
        .hash_kmers_scalar(packed_seq.as_slice())
        .map(|hash| hash as u64)
        .collect();
    let t_hashes: Vec<u64> = t_hasher
        .hash_kmers_scalar(packed_seq.as_slice())
        .map(|hash| hash as u64)
        .collect();

    let mut keys = Vec::with_capacity(num_mmers);
    let mut inner_minima = VecDeque::<usize>::new();
    for i in 0..num_mmers {
        let new_idx = i + inner_window - 1;
        while let Some(&back) = inner_minima.back() {
            if t_hashes[new_idx] < t_hashes[back] {
                inner_minima.pop_back();
            } else {
                break;
            }
        }
        inner_minima.push_back(new_idx);

        while let Some(&front) = inner_minima.front() {
            if front < i {
                inner_minima.pop_front();
            } else {
                break;
            }
        }

        let min_t_pos = inner_minima.front().copied().expect("non-empty inner window") - i;
        keys.push(OpenClosedKey {
            priority: open_closed_priority(m, t, min_t_pos),
            hash: m_hashes[i],
        });
    }

    let mut minimizer_positions = Vec::with_capacity(num_kmers);
    let mut outer_minima = VecDeque::<usize>::new();
    for i in 0..num_kmers {
        let new_idx = i + outer_window - 1;
        while let Some(&back) = outer_minima.back() {
            if keys[new_idx] < keys[back] {
                outer_minima.pop_back();
            } else {
                break;
            }
        }
        outer_minima.push_back(new_idx);

        while let Some(&front) = outer_minima.front() {
            if front < i {
                outer_minima.pop_front();
            } else {
                break;
            }
        }

        minimizer_positions.push(*outer_minima.front().expect("non-empty outer window"));
    }

    let (super_kmers, minimizer_vals) =
        build_superkmers_from_positions(&packed_seq, k, m, &minimizer_positions, &m_hashes);
    (super_kmers, minimizer_vals, packed_seq)
}

/// Decycling-set minimizers.
pub fn decycling_mins_x_pos(
    packed_seq: PackedSeqVec,
    k: u16,
    m: u16,
    decycler_set: &Decycler,
) -> (Vec<u32>, Vec<u64>, PackedSeqVec) {
    if packed_seq.len() < k as usize {
        return (Vec::new(), Vec::new(), packed_seq);
    }

    let mut minimizer_vals = Vec::new();
    let mut super_kmers = Vec::new();

    let mut is_decycler = Vec::with_capacity(packed_seq.len() - m as usize + 1);
    for i in 0..packed_seq.len() - m as usize + 1 {
        is_decycler.push(decycler_set.lookup(packed_seq.slice(i..i + m as usize)));
    }

    let mut mini_addrs = Vec::with_capacity(packed_seq.len() - k as usize + 1);
    let (mut min_addr, mut _is_decyc, mut _min_lexic) =
        mins_from_kmer(packed_seq.as_slice(), &is_decycler, 0, m, k);
    mini_addrs.push(min_addr);

    for i in 1..packed_seq.len() + 1 - k as usize {
        if min_addr < i {
            (min_addr, _is_decyc, _min_lexic) =
                mins_from_kmer(packed_seq.as_slice(), &is_decycler, i, m, k);
        }

        mini_addrs.push(min_addr);
    }

    let mut first_superkmer_addr: u32 = 0;
    let mut last_minimizer_pos = mini_addrs[0];
    for i in 1..packed_seq.len() + 1 - k as usize {
        if mini_addrs[i] != last_minimizer_pos {
            last_minimizer_pos = mini_addrs[i];
            super_kmers.push(first_superkmer_addr);
            minimizer_vals.push(
                packed_seq
                    .slice(
                        mini_addrs[first_superkmer_addr as usize]
                            ..mini_addrs[first_superkmer_addr as usize] + m as usize,
                    )
                    .as_u64(),
            );
            first_superkmer_addr = i as u32;
        }
    }

    super_kmers.push(first_superkmer_addr);
    minimizer_vals.push(
        packed_seq
            .slice(
                mini_addrs[first_superkmer_addr as usize]
                    ..mini_addrs[first_superkmer_addr as usize] + m as usize,
            )
            .as_u64(),
    );

    (super_kmers, minimizer_vals, packed_seq)
}

fn mins_from_kmer<'a>(
    packed_seq: PackedSeq<'a>,
    is_decycler: &[bool],
    i: usize,
    m: u16,
    k: u16,
) -> (usize, bool, PackedSeq<'a>) {
    for j in i..i + k as usize - m as usize + 1 {
        if is_decycler[j] {
            return (j, true, packed_seq.slice(j..j + m as usize));
        }
    }
    (i, false, packed_seq.slice(i..i + m as usize))
}

#[cfg(test)]
mod tests {
    use super::{
        decycling_mins_x_pos, minimizers_x_positions, mins_from_kmer, open_closed_mins_x_pos,
        open_closed_priority, selected_mins_x_pos, MinimizerMode,
    };
    use crate::decyclers::Decycler;
    use packed_seq::{PackedSeqVec, Seq, SeqVec};

    #[test]
    fn simd_minimizers_return_aligned_positions_and_values() {
        let sequence = PackedSeqVec::from_ascii(b"ACGTACGTAC");
        let (super_kmers, minimizers, original) = minimizers_x_positions(sequence.clone(), 5, 3);

        assert_eq!(original.len(), sequence.len());
        assert!(!super_kmers.is_empty());
        assert_eq!(super_kmers.len(), minimizers.len());
        assert_eq!(super_kmers[0], 0);
        assert!(super_kmers.windows(2).all(|w| w[0] < w[1]));
    }

    #[test]
    fn open_closed_priority_prefers_open_then_closed_then_ordinary() {
        assert_eq!(open_closed_priority(7, 3, 2), 0);
        assert_eq!(open_closed_priority(7, 3, 0), 1);
        assert_eq!(open_closed_priority(7, 3, 4), 1);
        assert_eq!(open_closed_priority(7, 3, 1), 2);
    }

    #[test]
    fn open_closed_minimizers_return_monotonic_super_kmer_starts() {
        let sequence = PackedSeqVec::from_ascii(b"ACGTACGTACGTACGT");
        let (super_kmers, minimizers, original) = open_closed_mins_x_pos(sequence.clone(), 7, 5, 3);

        assert_eq!(original.len(), sequence.len());
        assert!(!super_kmers.is_empty());
        assert_eq!(super_kmers.len(), minimizers.len());
        assert_eq!(super_kmers[0], 0);
        assert!(super_kmers.windows(2).all(|w| w[0] < w[1]));
    }

    #[test]
    fn open_closed_minimizers_support_exact_kmer_length() {
        let sequence = PackedSeqVec::from_ascii(b"ACGTA");
        let (super_kmers, minimizers, _) = open_closed_mins_x_pos(sequence, 5, 3, 2);

        assert_eq!(super_kmers.len(), 1);
        assert_eq!(minimizers.len(), 1);
        assert_eq!(super_kmers[0], 0);
    }

    #[test]
    fn selected_mins_dispatches_open_closed() {
        let sequence = PackedSeqVec::from_ascii(b"ACGTACGTAC");
        let decycler = Decycler::new(1);
        let (super_kmers, minimizers, _) =
            selected_mins_x_pos(sequence, 5, 3, &decycler, MinimizerMode::OpenClosed { t: 2 });

        assert!(!super_kmers.is_empty());
        assert_eq!(super_kmers.len(), minimizers.len());
    }

    #[test]
    fn decycling_minimizers_return_empty_for_short_sequences() {
        let mut decycler = Decycler::new(3);
        decycler.compute_blocks();
        let sequence = PackedSeqVec::from_ascii(b"ACG");
        let (super_kmers, minimizers, original) =
            decycling_mins_x_pos(sequence.clone(), 4, 3, &decycler);

        assert!(super_kmers.is_empty());
        assert!(minimizers.is_empty());
        assert_eq!(original.len(), sequence.len());
    }

    #[test]
    fn decycling_minimizers_return_monotonic_super_kmer_starts() {
        let mut decycler = Decycler::new(3);
        decycler.compute_blocks();
        let sequence = PackedSeqVec::from_ascii(b"ACGTACGTACGT");
        let (super_kmers, minimizers, _) = decycling_mins_x_pos(sequence, 5, 3, &decycler);

        assert!(!super_kmers.is_empty());
        assert_eq!(super_kmers.len(), minimizers.len());
        assert_eq!(super_kmers[0], 0);
        assert!(super_kmers.windows(2).all(|w| w[0] < w[1]));
    }

    #[test]
    fn mins_from_kmer_prefers_first_decycling_candidate() {
        let sequence = PackedSeqVec::from_ascii(b"ACGTAC");
        let is_decycler = vec![false, true, true, false];
        let (addr, is_member, minimizer) =
            mins_from_kmer(sequence.as_slice(), &is_decycler, 0, 3, 5);

        assert_eq!(addr, 1);
        assert!(is_member);
        assert_eq!(minimizer.as_u64(), sequence.slice(1..4).as_u64());
    }

    #[test]
    fn mins_from_kmer_falls_back_to_window_start_when_needed() {
        let sequence = PackedSeqVec::from_ascii(b"ACGTAC");
        let is_decycler = vec![false, false, false, false];
        let (addr, is_member, minimizer) =
            mins_from_kmer(sequence.as_slice(), &is_decycler, 0, 3, 5);

        assert_eq!(addr, 0);
        assert!(!is_member);
        assert_eq!(minimizer.as_u64(), sequence.slice(0..3).as_u64());
    }

    #[test]
    fn mins_from_kmer_ignores_candidates_outside_window() {
        let sequence = PackedSeqVec::from_ascii(b"ACGTAC");
        let is_decycler = vec![false, false, false, true];
        let (addr, is_member, minimizer) =
            mins_from_kmer(sequence.as_slice(), &is_decycler, 0, 3, 5);

        assert_eq!(addr, 0);
        assert!(!is_member);
        assert_eq!(minimizer.as_u64(), sequence.slice(0..3).as_u64());
    }

    #[test]
    fn simd_minimizers_support_exact_kmer_length() {
        let sequence = PackedSeqVec::from_ascii(b"ACGTA");
        let (super_kmers, minimizers, _) = minimizers_x_positions(sequence, 5, 3);

        assert_eq!(super_kmers.len(), 1);
        assert_eq!(minimizers.len(), 1);
        assert_eq!(super_kmers[0], 0);
    }

    #[test]
    fn simd_minimizers_return_original_sequence_unchanged() {
        let sequence = PackedSeqVec::from_ascii(b"ACGTACGT");
        let expected = sequence.as_slice().as_u64();
        let (_, _, original) = minimizers_x_positions(sequence, 5, 3);

        assert_eq!(original.as_slice().as_u64(), expected);
    }

    #[test]
    fn decycling_minimizers_become_non_empty_at_k_plus_three() {
        let mut decycler = Decycler::new(3);
        decycler.compute_blocks();
        let sequence = PackedSeqVec::from_ascii(b"ACGTACG");
        let (super_kmers, minimizers, _) = decycling_mins_x_pos(sequence, 4, 3, &decycler);

        assert!(!super_kmers.is_empty());
        assert_eq!(super_kmers.len(), minimizers.len());
    }

    #[test]
    fn decycling_minimizer_values_fit_requested_mmer_width() {
        let mut decycler = Decycler::new(3);
        decycler.compute_blocks();
        let sequence = PackedSeqVec::from_ascii(b"ACGTACGTACGT");
        let (super_kmers, minimizers, original) =
            decycling_mins_x_pos(sequence, 5, 3, &decycler);

        assert_eq!(super_kmers.len(), minimizers.len());
        for minimizer in minimizers {
            assert!(minimizer < (1 << 6));
        }
        assert_eq!(original.len(), 12);
    }
}
