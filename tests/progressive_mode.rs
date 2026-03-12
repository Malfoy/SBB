mod common;

use biobloom_rs::bloom::{BloomLayout, ConcurrentBloomFilter};
use biobloom_rs::fastx::{insert_file_into_filter, progressive_insert_file};
use tempfile::tempdir;

#[test]
fn progressive_insert_adds_kmers_and_subtract_can_suppress_them() {
    let d = tempdir().expect("tempdir should be creatable");
    let seed = d.path().join("seed.fa");
    let recruit = d.path().join("recruit.fa");

    common::write_fasta(&seed, &[("s1", "ACGTACGTACGTACGTACGTACGT")])
        .expect("seed fasta should be writable");
    common::write_fasta(
        &recruit,
        &[
            ("r1", "ACGTACGTACGTACGTACGTACGT"),
            ("r2", "TGCATGCATGCATGCATGCATGCA"),
        ],
    )
    .expect("recruit fasta should be writable");

    let filter = ConcurrentBloomFilter::new(5, 3, 1 << 16, 100, 0.01, 1.0, true, 8)
        .expect("filter init should succeed");
    let seeded = insert_file_into_filter(&seed, &filter).expect("seeding should succeed");
    assert!(seeded > 0);

    let stats = progressive_insert_file(&recruit, &filter, 1.0, None)
        .expect("progressive insert should succeed");
    assert_eq!(stats.total_reads, 2);
    assert!(stats.queried_kmers > 0);
    assert!(stats.matched_reads >= 1);
    assert!(stats.inserted_kmers > 0);

    let subtract_filter = ConcurrentBloomFilter::new(5, 3, 1 << 16, 100, 0.01, 1.0, false, 0)
        .expect("subtract filter init should succeed");
    let subtract_inserted =
        insert_file_into_filter(&recruit, &subtract_filter).expect("subtract load should succeed");
    assert!(subtract_inserted > 0);
    let subtract_filter = subtract_filter.finalize();

    let filter_with_subtract = ConcurrentBloomFilter::new(5, 3, 1 << 16, 100, 0.01, 1.0, true, 8)
        .expect("filter init should succeed");
    let _ = insert_file_into_filter(&seed, &filter_with_subtract).expect("seeding should succeed");
    let stats_sub =
        progressive_insert_file(&recruit, &filter_with_subtract, 1.0, Some(&subtract_filter))
            .expect("progressive insert with subtract should succeed");
    assert_eq!(stats_sub.inserted_kmers, 0);
}

#[test]
fn progressive_bloomy_subtract_is_rejected() {
    let d = tempdir().expect("tempdir should be creatable");
    let seed = d.path().join("seed.fa");
    let recruit = d.path().join("recruit.fa");

    common::write_fasta(
        &seed,
        &[(
            "s1",
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA",
        )],
    )
    .expect("seed fasta should be writable");
    common::write_fasta(
        &recruit,
        &[(
            "r1",
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA",
        )],
    )
    .expect("recruit fasta should be writable");

    let filter = ConcurrentBloomFilter::new_bloomybloom(41, 1 << 23, 100, 0.01, 1.0, None)
        .expect("bloomy filter init should succeed");
    assert_eq!(filter.layout, BloomLayout::BloomyBloom);
    let seeded = insert_file_into_filter(&seed, &filter).expect("seeding should succeed");
    assert!(seeded > 0);

    let subtract_filter = ConcurrentBloomFilter::new(5, 3, 1 << 16, 100, 0.01, 1.0, false, 0)
        .expect("subtract filter init should succeed")
        .finalize();

    let err = progressive_insert_file(&recruit, &filter, 1.0, Some(&subtract_filter))
        .expect_err("bloomy progressive subtract should be rejected");
    assert!(err.to_string().contains("not supported"));
}
