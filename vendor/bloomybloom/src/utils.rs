///module with small utility fucntions

//use bitvec::prelude::*;
//use packed_seq::{PackedSeq, Seq};

///64-bit mixing function used by bloomybloom for bucket and address routing.
///The original xorshift form has weak low-bit avalanche; this splitmix64-style
///finalizer provides substantially better diffusion for masked-index use.
pub fn xorshift_u64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^ (x >> 31)
}

#[inline(always)]
pub fn fold_u128_to_u64(x: u128) -> u64 {
    (x as u64) ^ ((x >> 64) as u64)
}

#[inline(always)]
pub fn hash_u128_to_u64(x: u128) -> u64 {
    xorshift_u64(fold_u128_to_u64(x))
}

#[inline(always)]
pub fn roll_u128_kmer(window: u128, next_base: u8, len: u16) -> u128 {
    (window >> 2) | ((next_base as u128) << (2 * (len as usize - 1)))
}

pub fn _xorshift_u32(mut x: u32) -> u32 {
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    x
}

///since no implementation directly on u128 exists, I just used whatever numbers of shifts
pub fn xorshift_u128(mut x: u128) -> u128 {
    x ^= x << 17;
    x ^= x >> 23;
    x ^= x << 5;
    x
}

///returns the amount of positive bits in a boolean vector
pub fn sum_vec_bool(boolean_vector: &Vec<bool>) -> usize {
    let mut counter: usize = 0;
    for val in boolean_vector {
        if *val {counter += 1};
    }
    counter
}

#[cfg(test)]
mod tests {
    use super::{
        _xorshift_u32, fold_u128_to_u64, hash_u128_to_u64, roll_u128_kmer, sum_vec_bool,
        xorshift_u64, xorshift_u128,
    };
    use packed_seq::{Seq, SeqVec};

    #[test]
    fn xorshift_u64_is_deterministic() {
        assert_eq!(xorshift_u64(42), xorshift_u64(42));
        assert_ne!(xorshift_u64(42), 42);
    }

    #[test]
    fn xorshift_u128_is_deterministic() {
        assert_eq!(xorshift_u128(42), xorshift_u128(42));
        assert_ne!(xorshift_u128(42), 42);
    }

    #[test]
    fn sum_vec_bool_counts_true_values() {
        let values = vec![true, false, true, true, false];
        assert_eq!(sum_vec_bool(&values), 3);
    }

    #[test]
    fn xorshift_u64_zero_is_well_mixed() {
        assert_ne!(xorshift_u64(0), 0);
    }

    #[test]
    fn xorshift_u128_zero_stays_zero() {
        assert_eq!(xorshift_u128(0), 0);
    }

    #[test]
    fn fold_u128_to_u64_xors_high_and_low_halves() {
        let value = ((0x0123_4567_89ab_cdef_u128) << 64) | 0xfedc_ba98_7654_3210_u128;
        assert_eq!(fold_u128_to_u64(value), 0xffff_ffff_ffff_ffff_u64);
    }

    #[test]
    fn hash_u128_to_u64_is_deterministic() {
        assert_eq!(hash_u128_to_u64(123456789), hash_u128_to_u64(123456789));
    }

    #[test]
    fn roll_u128_kmer_drops_first_base_and_appends_new_one() {
        let start = packed_seq::PackedSeqVec::from_ascii(b"ACGT").as_slice().as_u128();
        let rolled = roll_u128_kmer(start, 0, 4);
        let expected = packed_seq::PackedSeqVec::from_ascii(b"CGTA").as_slice().as_u128();
        assert_eq!(rolled, expected);
    }

    #[test]
    fn xorshift_u32_is_deterministic() {
        assert_eq!(_xorshift_u32(12345), _xorshift_u32(12345));
        assert_ne!(_xorshift_u32(12345), 12345);
    }

    #[test]
    fn sum_vec_bool_empty_vector_is_zero() {
        let values = vec![];
        assert_eq!(sum_vec_bool(&values), 0);
    }

    #[test]
    fn sum_vec_bool_all_true_counts_every_entry() {
        let values = vec![true, true, true, true];
        assert_eq!(sum_vec_bool(&values), 4);
    }
}
