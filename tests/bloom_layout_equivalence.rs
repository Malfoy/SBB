use biobloom_rs::bloom::{BloomLayout, ConcurrentBloomFilter};

#[test]
fn classic_and_blocked_have_no_false_negatives_on_inserted_set() {
    let bit_len = 1 << 20;
    let hash_count = 4;
    let kmer_size = 21;
    let mut inserted = Vec::with_capacity(20_000);

    let classic =
        ConcurrentBloomFilter::new(kmer_size, hash_count, bit_len, 20_000, 0.01, false, 0)
            .expect("classic init should succeed");
    let blocked = ConcurrentBloomFilter::new(kmer_size, hash_count, bit_len, 20_000, 0.01, true, 8)
        .expect("blocked init should succeed");

    for i in 0_u64..20_000 {
        let key = i.wrapping_mul(0x9E3779B97F4A7C15).rotate_left(17);
        inserted.push(key);
        classic.insert_kmer(key);
        blocked.insert_kmer(key);
    }

    let classic = classic.finalize();
    let blocked = blocked.finalize();
    assert_eq!(classic.layout, BloomLayout::Classic);
    assert_eq!(blocked.layout, BloomLayout::Blocked);

    for key in inserted {
        assert!(
            classic.contains_kmer(key),
            "classic false negative for key={key}"
        );
        assert!(
            blocked.contains_kmer(key),
            "blocked false negative for key={key}"
        );
    }
}
