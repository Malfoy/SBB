// module to implement decycling minimizers
// because we will always use (relatively) small minimizers
// the strategy will always be to first compute the complete list of all decycled minimizers, and
// then to query it in read only accross the different threads

use rayon::prelude::*;
//use packed_seq::{PackedSeqVec, SeqVec, PackedSeq, Seq};
use packed_seq::{PackedSeq, Seq};

const CYCLER_BLOCK_SIZE: usize = 512; //just always keep it a power of 2
const CYCLER_BLOCK_MASK_MULT_64: usize = 64*CYCLER_BLOCK_SIZE-1;
const CYCLER_BLOCK_SHIFTER: usize = 15;

//const CYCLER_BLOCK_SIZE: usize = 1; //just always keep it a power of 2

pub struct Decycler {
    pub m: u16, //minimizer length
    pub direct_list: Vec<Vec<u64>>, //stores in "booleans" if a string is a direct decycling minimizer
                                //based on its address
    //indirect_lit, to add
}

impl Decycler {

    ///initialization, creating enough space for all the minimizers
    pub fn new(m: u16) -> Self {
        let block_count = usize::max(1, (1usize << (2 * m as usize)) >> CYCLER_BLOCK_SHIFTER);
        let direct_list: Vec<Vec<u64>> = vec![vec![0; CYCLER_BLOCK_SIZE]; block_count];
        init_vec_ci(m);
        Self {
            m,
            direct_list,
        }
    }

    ///computes the belonging (or not) of all the kmers
    pub fn compute_blocks(&mut self) {
        let vec_ci: Vec<f64> = init_vec_ci(self.m);
        self.direct_list.par_iter_mut().enumerate().for_each(|(i, mut block)| {
            //println!("ca compute un block en léééégende");
            compute_block(i, &mut block, self.m, &vec_ci);
        })
    }

    pub fn lookup(&self, minimizer: PackedSeq) -> bool {
        //start by converting the kmer to an address we can use
        let address = minimizer.as_u64() as usize;

        //finnding the block
        //let block_adress: usize = address/(CYCLER_BLOCK_SIZE*64);
        //REMOVED MODULO
        let block_adress: usize = address >> CYCLER_BLOCK_SHIFTER;
        let block = &self.direct_list[block_adress];

        //lookup the corresponding u64
        //let integer_adress: usize = (address%(CYCLER_BLOCK_SIZE*64))/64;
        //REMOVED MODULO
        let integer_adress: usize = (address&(CYCLER_BLOCK_MASK_MULT_64))>>6;
        let integer: u64 = block[integer_adress];

        //reading the correct bit using bitshifting
        //let boolean: bool = if (integer>>(63-address%64))%2 == 1 {true} else {false};
        //REMOVED MODULO
        let boolean: bool = if (integer>>(63-(address&63)))&1 == 1 {true} else {false};

        //return
        boolean
    }
}


fn compute_block(i: usize, block: &mut Vec<u64>, m: u16, vec_ci: &Vec<f64>) {
    let mut kmer: u64 = (i*CYCLER_BLOCK_SIZE*64) as u64;

    for j in 0..CYCLER_BLOCK_SIZE {
        let mut to_insert: u64 = 0;
        for k in 0..64 {
            //let is_decycler: bool = compute_membership(4_u64.pow(m as u32)-kmer, m, vec_ci);
            let is_decycler: bool = compute_membership(kmer, m, vec_ci);
            if is_decycler {
                to_insert += 1<<(63-k);
            }

            kmer +=1;
        }
        block[j] = to_insert;
    }
}

///algorithm to check if a kmer is member of a minimum decysling set
///see "Efficient minimizer orders for large values of k using minimum decycling sets" by
///David Pellow, Lianrong Pu, Baris Ekim, Kior Kotlar, Bonnie Berger, Ron Shamir and Yaron
///Orenstein; page 3
pub fn compute_membership(kmer: u64, m: u16, vec_ci: &Vec<f64>) -> bool {
    let epsilon: f64 = 0.00001;
    let mut imaginary_x: f64 = 0.0; 
    let mut imaginary_x_prime: f64 = 0.0;
    for i in 0..m {
        //let x_i: u64 = (kmer>>(2*i))%4;
        //REMOVED MODULO
        let x_i: u64 = (kmer>>(2*i))&3;
        imaginary_x += vec_ci[i as usize]*x_i as f64;
        let i_prime: usize = if i<m-1 {i as usize+1} else {0};
        imaginary_x_prime += vec_ci[i_prime]*x_i as f64;
    }
    //println!("partie imaginaire : {imaginary_x}, du précédent : {imaginary_x_prime}");

    if imaginary_x > epsilon {
        if imaginary_x_prime <= epsilon {
            return true
        }
    }else if imaginary_x >= -epsilon && imaginary_x <= epsilon { //testing equality actually
        if imaginary_x_prime >= -epsilon && imaginary_x_prime <= epsilon {
            let mut k: u16 = 0;
            for l in 0..2*m {
                //REMOVED MODULO
                let x_l_mod_m: u64 = (kmer>>(2*(m-(l%m)-1)))&3;
                let x_k: u64 = (kmer>>(2*(m-k-1)))&3;
                if x_l_mod_m < x_k {return false};
                if x_l_mod_m > x_k {k = 0} else {k += 1};
                if (l >= m-1) && (k%m == 0) {return true};
            }
        }
    }
    return false
}

pub fn init_vec_ci(m: u16) -> Vec<f64> {
    let mut vec_ci: Vec<f64> = Vec::with_capacity(m as usize);
    for i in 0..m {
        vec_ci.push((2.0*std::f64::consts::PI*i as f64/m as f64).sin());
    }
    vec_ci
}

#[cfg(test)]
mod tests {
    use super::{compute_membership, init_vec_ci, Decycler};
    use packed_seq::{PackedSeqVec, Seq, SeqVec};

    fn decode_kmer(mut encoded: u64, m: u16) -> String {
        let mut chars = vec!['A'; m as usize];
        for idx in (0..m as usize).rev() {
            chars[idx] = match encoded & 3 {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                3 => 'T',
                _ => unreachable!(),
            };
            encoded >>= 2;
        }
        chars.into_iter().collect()
    }

    #[test]
    fn init_vec_ci_has_expected_shape() {
        let values = init_vec_ci(4);
        assert_eq!(values.len(), 4);
        assert!((values[0] - 0.0).abs() < 1e-12);
        assert!((values[1] - 1.0).abs() < 1e-12);
        assert!(values[2].abs() < 1e-12);
        assert!((values[3] + 1.0).abs() < 1e-12);
    }

    #[test]
    fn compute_blocks_matches_direct_membership_for_all_mers_of_size_three() {
        let m = 3;
        let vec_ci = init_vec_ci(m);
        let mut decycler = Decycler::new(m);
        decycler.compute_blocks();

        for encoded in 0..(1_u64 << (2 * m)) {
            let kmer = decode_kmer(encoded, m);
            let packed = PackedSeqVec::from_ascii(kmer.as_bytes());
            assert_eq!(
                decycler.lookup(packed.as_slice()),
                compute_membership(packed.as_slice().as_u64(), m, &vec_ci),
                "mismatch for {kmer}"
            );
        }
    }

    #[test]
    fn known_block_contents_are_stable() {
        let mut decycler = Decycler::new(3);
        decycler.compute_blocks();

        assert_eq!(
            decycler.direct_list[0][0],
            0b1000110011101111000001001110111100000000001011110000000000000001
        );
    }

    #[test]
    fn new_small_m_allocates_at_least_one_block() {
        let decycler = Decycler::new(1);
        assert_eq!(decycler.direct_list.len(), 1);
        assert_eq!(decycler.direct_list[0].len(), 512);
    }

    #[test]
    fn init_vec_ci_for_m_one_is_zero() {
        let values = init_vec_ci(1);
        assert_eq!(values.len(), 1);
        assert!(values[0].abs() < 1e-12);
    }

    #[test]
    fn compute_membership_is_deterministic() {
        let vec_ci = init_vec_ci(4);
        let kmer = PackedSeqVec::from_ascii(b"ACGT").as_slice().as_u64();
        assert_eq!(
            compute_membership(kmer, 4, &vec_ci),
            compute_membership(kmer, 4, &vec_ci)
        );
    }

    #[test]
    fn lookup_matches_direct_membership_for_all_mers_of_size_one() {
        let m = 1;
        let vec_ci = init_vec_ci(m);
        let mut decycler = Decycler::new(m);
        decycler.compute_blocks();

        for kmer in [b"A", b"C", b"G", b"T"] {
            let seq = PackedSeqVec::from_ascii(kmer);
            let packed = seq.as_slice();
            assert_eq!(
                decycler.lookup(packed),
                compute_membership(packed.as_u64(), m, &vec_ci)
            );
        }
    }

    #[test]
    fn lookup_matches_direct_membership_for_all_mers_of_size_two() {
        let m = 2;
        let vec_ci = init_vec_ci(m);
        let mut decycler = Decycler::new(m);
        decycler.compute_blocks();

        for left in [b'A', b'C', b'G', b'T'] {
            for right in [b'A', b'C', b'G', b'T'] {
                let kmer = [left, right];
                let seq = PackedSeqVec::from_ascii(&kmer);
                let packed = seq.as_slice();
                assert_eq!(
                    decycler.lookup(packed),
                    compute_membership(packed.as_u64(), m, &vec_ci)
                );
            }
        }
    }

    #[test]
    fn lookup_matches_direct_membership_for_specific_size_four_examples() {
        let m = 4;
        let vec_ci = init_vec_ci(m);
        let mut decycler = Decycler::new(m);
        decycler.compute_blocks();

        for kmer in [b"AAAA", b"ACGT", b"TTTT", b"CGTA"] {
            let seq = PackedSeqVec::from_ascii(kmer);
            let packed = seq.as_slice();
            assert_eq!(
                decycler.lookup(packed),
                compute_membership(packed.as_u64(), m, &vec_ci)
            );
        }
    }

    #[test]
    fn computed_decycler_contains_at_least_one_member_for_m_three() {
        let mut decycler = Decycler::new(3);
        decycler.compute_blocks();
        assert!(decycler.direct_list[0].iter().any(|entry| *entry != 0));
    }
}
