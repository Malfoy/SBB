use biobloom_rs::bloom::for_each_canonical_kmer;

#[test]
fn ambiguous_base_resets_kmer_run() {
    let count = for_each_canonical_kmer(b"AAAAANAAAAA", 5, |_| {});
    assert_eq!(count, 2);
}

#[test]
fn k1_counts_only_unambiguous_bases() {
    let mut out = Vec::new();
    let count = for_each_canonical_kmer(b"ACGTN", 1, |k| out.push(k));
    assert_eq!(count, 4);
    assert_eq!(out.len(), 4);
}

#[test]
fn k32_handles_max_kmer_size() {
    let seq = [b'A'; 32];
    let mut out = Vec::new();
    let count = for_each_canonical_kmer(&seq, 32, |k| out.push(k));
    assert_eq!(count, 1);
    assert_eq!(out, vec![0_u64]);
}

#[test]
fn short_sequence_has_zero_kmers() {
    let count = for_each_canonical_kmer(b"ACGT", 5, |_| {});
    assert_eq!(count, 0);
}
