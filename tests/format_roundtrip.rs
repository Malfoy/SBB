use biobloom_rs::bloom::{BloomFilter, BloomLayout, ConcurrentBloomFilter};
use tempfile::tempdir;

fn assert_roundtrip(blocked: bool, block_words: u16, sample_rate: f64) {
    let d = tempdir().expect("tempdir should be creatable");
    let path = d.path().join(if blocked {
        "blocked.bf.zst"
    } else {
        "classic.bf.zst"
    });
    let expected_layout = if blocked {
        BloomLayout::Blocked
    } else {
        BloomLayout::Classic
    };

    let filter = ConcurrentBloomFilter::new(
        21,
        4,
        1 << 18,
        10_000,
        0.001,
        sample_rate,
        blocked,
        block_words,
    )
    .expect("filter init should succeed");
    for i in 0_u64..10_000 {
        let key = i.wrapping_mul(0xD6E8FEB86659FD93).rotate_left(11);
        filter.insert_kmer(key);
    }
    let filter = filter.finalize();
    filter.save(&path).expect("save should succeed");

    let loaded = BloomFilter::load(&path).expect("load should succeed");
    assert_eq!(loaded.layout, expected_layout);
    if blocked {
        assert_eq!(loaded.block_words, block_words);
    } else {
        assert_eq!(loaded.block_words, 0);
    }
    assert_eq!(loaded.hash_count, 4);
    assert_eq!(loaded.kmer_size, 21);
    assert_eq!(loaded.bit_len, 1 << 18);
    assert_eq!(loaded.sample_rate, sample_rate);

    for i in 0_u64..10_000 {
        let key = i.wrapping_mul(0xD6E8FEB86659FD93).rotate_left(11);
        if loaded.should_keep_kmer(key) {
            assert!(loaded.contains_kmer(key), "false negative for key={key}");
        } else {
            assert!(!loaded.contains_kmer(key), "unexpected hit for key={key}");
        }
    }
}

#[test]
fn classic_roundtrip_preserves_metadata_and_membership() {
    assert_roundtrip(false, 0, 1.0);
}

#[test]
fn blocked_roundtrip_preserves_metadata_and_membership() {
    assert_roundtrip(true, 8, 0.5);
}

#[test]
fn bloomy_roundtrip_preserves_metadata_and_sequence_hits() {
    let d = tempdir().expect("tempdir should be creatable");
    let path = d.path().join("bloomy.bf.zst");
    let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA";

    let filter = ConcurrentBloomFilter::new_bloomybloom(41, 1 << 23, 1_000, 0.001, 1.0, None)
        .expect("bloomy filter init should succeed");
    assert!(filter.insert_sequence(seq) > 0);

    let filter = filter.finalize();
    filter.save(&path).expect("save should succeed");

    let loaded = BloomFilter::load(&path).expect("load should succeed");
    assert_eq!(loaded.layout, BloomLayout::BloomyBloom);
    assert_eq!(loaded.kmer_size, 41);
    assert_eq!(loaded.sample_rate, 1.0);

    let (hits, total) = loaded.score_sequence_with_counts(seq);
    assert!(total > 0);
    assert_eq!(hits, total);
}
