use biobloom_rs::bloom::{BloomFilter, BloomLayout, ConcurrentBloomFilter};
use tempfile::tempdir;

fn assert_roundtrip(blocked: bool, block_words: u16) {
    let d = tempdir().expect("tempdir should be creatable");
    let path = d
        .path()
        .join(if blocked { "blocked.bf" } else { "classic.bf" });
    let expected_layout = if blocked {
        BloomLayout::Blocked
    } else {
        BloomLayout::Classic
    };

    let filter = ConcurrentBloomFilter::new(21, 4, 1 << 18, 10_000, 0.001, blocked, block_words)
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

    for i in 0_u64..10_000 {
        let key = i.wrapping_mul(0xD6E8FEB86659FD93).rotate_left(11);
        assert!(loaded.contains_kmer(key), "false negative for key={key}");
    }
}

#[test]
fn classic_roundtrip_preserves_metadata_and_membership() {
    assert_roundtrip(false, 0);
}

#[test]
fn blocked_roundtrip_preserves_metadata_and_membership() {
    assert_roundtrip(true, 8);
}
