use biobloom_rs::bloom::{bit_len_for_fpr_with_hash_count, optimal_bit_len, optimal_hash_count};

#[test]
fn lower_fpr_requires_more_bits() {
    let n = 1_000_000;
    let looser = optimal_bit_len(n, 1e-2);
    let tighter = optimal_bit_len(n, 1e-5);
    assert!(tighter > looser);
}

#[test]
fn more_elements_require_more_bits() {
    let fpr = 1e-3;
    let small = optimal_bit_len(100_000, fpr);
    let large = optimal_bit_len(1_000_000, fpr);
    assert!(large > small);
}

#[test]
fn optimal_hash_count_is_bounded_and_nonzero() {
    let k = optimal_hash_count(8_000_000, 1_000_000);
    assert!((1..=255).contains(&k));
}

#[test]
fn fixed_hash_bit_len_gets_tighter_with_lower_fpr() {
    let n = 1_000_000;
    let k = 4;
    let looser = bit_len_for_fpr_with_hash_count(n, 1e-2, k);
    let tighter = bit_len_for_fpr_with_hash_count(n, 1e-5, k);
    assert!(tighter > looser);
}
