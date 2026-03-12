///module taht implements my super u64 bitvecs
use std::fmt::Error;
use std::fmt::Display;
use std::fmt;

pub struct SuperBitVec {
    vector: Vec<u64>,
    size: usize,
}

impl SuperBitVec {

    ///creates a new superbitvec, filled with 0's, "makes the actual complete structure"
    pub fn new(size : usize) -> Self {
        let vec_size: usize = if size%64==0 {size/64} else {size/64+1};
        let vector: Vec<u64> = vec![0; vec_size];
        
        Self {
            vector,
            size,
        }
    }

    //sets the bit at "index" address to value
    pub fn set(&mut self, address: usize, value: bool) {
        //first check the address is legal
        if address >= self.size {
            panic!("Index out of SuperBitVec range");
        }

        //compute which bit to insert
        let block_num: usize = address>>6;
        let mut to_insert: u64 = 1<<(63-(address&63)) as u64; //trust the calculation
        if value == false {
            to_insert = u64::MAX-to_insert;
        }

        //now performing the insertion with an atomic or
        if value == true {
            self.vector[block_num] = self.vector[block_num] | to_insert;
        } else {
            self.vector[block_num] = self.vector[block_num] & to_insert;
        }
    }

    ///getter for a certain bit
    pub fn get(&self, address: usize) -> bool {
        let block = self.vector[address>>6];
        let boolean: bool = if (block>>(63-(address&63)))&1 == 1 {true} else {false};
        boolean
    }

    ///very standard len method
    pub fn len(&self) -> usize {
        self.size
    }
}


impl Display for SuperBitVec {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), Error> {

        //is meant for some debugging on smaller examples, so we have a max threshhold
        let write_threshhold: usize = 8192;
        if self.size > write_threshhold {
            return write!(f, "SuperBitVec of length over {write_threshhold} wasn't written.");
        }

        //special case for the first block depends on the actual length of the SuperBitVec
        let mut to_write = String::new();

        let to_push: usize = if self.size&63==0 {0} else {64-(self.size&63)};
        let first_block = self.vector[0] >> to_push;
        to_write += &format!("{:0width$b}", first_block, width = 64-to_push); 

        let last_block_number: usize = if self.size&63==0 {self.size>>6} else {(self.size>>6)+1};
        for i in 1..(last_block_number) {
            to_write += &format!("{:064b}", self.vector[i]);
        }
        write!(f, "{to_write}")
    }
}

impl Clone for SuperBitVec {
    fn clone(&self) -> Self {
        let vector = self.vector.clone();
        let size = self.size;
        Self {
            vector,
            size,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::SuperBitVec;

    #[test]
    fn set_and_get_bits_across_blocks() {
        let mut bitvec = SuperBitVec::new(130);
        bitvec.set(0, true);
        bitvec.set(64, true);
        bitvec.set(129, true);

        assert!(bitvec.get(0));
        assert!(bitvec.get(64));
        assert!(bitvec.get(129));
        assert!(!bitvec.get(1));
        assert!(!bitvec.get(128));
    }

    #[test]
    fn clearing_bit_restores_false_value() {
        let mut bitvec = SuperBitVec::new(10);
        bitvec.set(3, true);
        assert!(bitvec.get(3));
        bitvec.set(3, false);
        assert!(!bitvec.get(3));
    }

    #[test]
    #[should_panic(expected = "Index out of SuperBitVec range")]
    fn setting_out_of_range_panics() {
        let mut bitvec = SuperBitVec::new(8);
        bitvec.set(8, true);
    }

    #[test]
    fn display_matches_small_bit_pattern() {
        let mut bitvec = SuperBitVec::new(5);
        bitvec.set(0, true);
        bitvec.set(2, true);
        bitvec.set(4, true);

        assert_eq!(format!("{bitvec}"), "10101");
    }

    #[test]
    fn clone_preserves_content() {
        let mut bitvec = SuperBitVec::new(66);
        bitvec.set(1, true);
        bitvec.set(65, true);

        let cloned = bitvec.clone();
        assert!(cloned.get(1));
        assert!(cloned.get(65));
        assert_eq!(cloned.len(), 66);
    }

    #[test]
    fn fresh_bitvec_is_all_false() {
        let bitvec = SuperBitVec::new(70);
        for idx in 0..70 {
            assert!(!bitvec.get(idx));
        }
    }

    #[test]
    fn len_matches_requested_size() {
        let bitvec = SuperBitVec::new(127);
        assert_eq!(bitvec.len(), 127);
    }

    #[test]
    fn setting_last_bit_of_partial_block_works() {
        let mut bitvec = SuperBitVec::new(65);
        bitvec.set(64, true);
        assert!(bitvec.get(64));
    }

    #[test]
    fn clearing_unset_bit_keeps_it_false() {
        let mut bitvec = SuperBitVec::new(12);
        bitvec.set(7, false);
        assert!(!bitvec.get(7));
    }

    #[test]
    fn display_large_vectors_uses_summary_message() {
        let bitvec = SuperBitVec::new(9000);
        assert_eq!(
            format!("{bitvec}"),
            "SuperBitVec of length over 8192 wasn't written."
        );
    }
}
