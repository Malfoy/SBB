///module to read an inputed fasta file, and maybe later a .txt file of file
///output will always be some compressed sequences (packed_seq) from imartayan

use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
use packed_seq::{PackedSeqVec, SeqVec};
use needletail::FastxReader;

//pub fn read_fasta(fasta_file: String) -> packed_seq::packed_seq::PackedSeqVecBase<2> {
pub fn _read_fasta(fasta_file: String) -> PackedSeqVec {
    //var that will contain the concatenation of all lines before conversion to packed_seq
    let mut full_ascii: String = String::new();
    if let Ok(lines) = _read_lines(fasta_file) {
        for line in lines {
            let unwrapped_line = line.expect("Problem reading a FASTA");
            let line_bytes = unwrapped_line.as_bytes();
            //filter out all the comments
            if line_bytes.len() >0 && line_bytes[0] != b'>' && line_bytes[0] != b';' {
                full_ascii += &unwrapped_line;
            }
        }
    }
    //once we have everything in a single String, its time to turn it into a more efficient
    //packed_seq, to be used quickly later
    let packed_seq = PackedSeqVec::from_ascii(full_ascii.as_bytes());

    packed_seq
}

///function for reading a file of file to handle lots of fasta at once
///the fof should have the path to a single fasta on each line
pub fn _read_fof(fof_file: String) -> Vec<String> {
    let mut iter_files: Vec<String> = Vec::new();
    if let Ok(lines) = _read_lines(fof_file) {
        for line in lines {
            iter_files.push(line.unwrap());
        }
    }
    iter_files
}

///classic function to simply read any file line by line efficiently
pub fn _read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

/// I never hated a struct more than an Iterator with lifetimes, worst invention in Mankind history
pub struct Hell {
    pub fxreader: Box<dyn FastxReader>,
    pub chunk_size: usize,
}

impl Iterator for Hell {
    type Item = Vec<Vec<u8>>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut chunk = Vec::new();
        for _ in 0..self.chunk_size {
            let result = match self.fxreader.next() {
                Some(res) => res,
                None => break,
            };
            let seq_red = result.unwrap().seq().to_mut().clone();
            chunk.push(seq_red);
        }
        if chunk.is_empty() {
            None
        } else {
            return Some(chunk);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{_read_fasta, _read_fof, Hell};
    use needletail::parse_fastx_file;
    use packed_seq::{Seq, SeqVec};
    use std::fs;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temp_path(name: &str) -> PathBuf {
        let mut path = std::env::temp_dir();
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        path.push(format!("bloomybloom_{name}_{unique}_{}", std::process::id()));
        path
    }

    #[test]
    fn read_fasta_skips_headers_and_comments() {
        let path = temp_path("input.fasta");
        fs::write(&path, b">seq1\nACGT\n;comment\nTGCA\n").unwrap();

        let packed = _read_fasta(path.to_string_lossy().into_owned());

        assert_eq!(packed.len(), 8);
        assert_eq!(packed.as_slice().as_u64(), packed_seq::PackedSeqVec::from_ascii(b"ACGTTGCA").as_slice().as_u64());

        fs::remove_file(path).unwrap();
    }

    #[test]
    fn read_fof_returns_all_lines() {
        let fasta_a = temp_path("a.fasta");
        let fasta_b = temp_path("b.fasta");
        let fof = temp_path("files.txt");
        fs::write(&fof, format!("{}\n{}\n", fasta_a.display(), fasta_b.display())).unwrap();

        let files = _read_fof(fof.to_string_lossy().into_owned());
        assert_eq!(files.len(), 2);
        assert_eq!(files[0], fasta_a.to_string_lossy());
        assert_eq!(files[1], fasta_b.to_string_lossy());

        fs::remove_file(fof).unwrap();
    }

    #[test]
    fn hell_iterator_chunks_fastx_records() {
        let path = temp_path("reads.fasta");
        fs::write(&path, b">r1\nAAAA\n>r2\nCCCC\n>r3\nGGGG\n").unwrap();

        let reader = parse_fastx_file(&path).unwrap();
        let mut hell = Hell {
            fxreader: reader,
            chunk_size: 2,
        };

        let first = hell.next().unwrap();
        let second = hell.next().unwrap();
        let third = hell.next();

        assert_eq!(first, vec![b"AAAA".to_vec(), b"CCCC".to_vec()]);
        assert_eq!(second, vec![b"GGGG".to_vec()]);
        assert!(third.is_none());

        fs::remove_file(path).unwrap();
    }

    #[test]
    fn read_lines_returns_exact_lines() {
        let path = temp_path("lines.txt");
        fs::write(&path, b"alpha\nbeta\ngamma\n").unwrap();

        let lines: Vec<String> = super::_read_lines(&path)
            .unwrap()
            .map(|line| line.unwrap())
            .collect();

        assert_eq!(lines, vec!["alpha", "beta", "gamma"]);
        fs::remove_file(path).unwrap();
    }

    #[test]
    fn read_fasta_empty_file_returns_empty_sequence() {
        let path = temp_path("empty.fasta");
        fs::write(&path, b"").unwrap();

        let packed = _read_fasta(path.to_string_lossy().into_owned());
        assert_eq!(packed.len(), 0);

        fs::remove_file(path).unwrap();
    }

    #[test]
    fn read_fasta_comment_only_returns_empty_sequence() {
        let path = temp_path("comment_only.fasta");
        fs::write(&path, b">seq1\n;comment\n").unwrap();

        let packed = _read_fasta(path.to_string_lossy().into_owned());
        assert_eq!(packed.len(), 0);

        fs::remove_file(path).unwrap();
    }

    #[test]
    fn hell_iterator_with_chunk_size_one_yields_single_records() {
        let path = temp_path("chunk1.fasta");
        fs::write(&path, b">r1\nAAAA\n>r2\nCCCC\n").unwrap();

        let reader = parse_fastx_file(&path).unwrap();
        let mut hell = Hell {
            fxreader: reader,
            chunk_size: 1,
        };

        assert_eq!(hell.next().unwrap(), vec![b"AAAA".to_vec()]);
        assert_eq!(hell.next().unwrap(), vec![b"CCCC".to_vec()]);
        assert!(hell.next().is_none());

        fs::remove_file(path).unwrap();
    }

    #[test]
    fn parse_fastx_file_rejects_empty_input() {
        let path = temp_path("empty_reads.fasta");
        fs::write(&path, b"").unwrap();

        assert!(parse_fastx_file(&path).is_err());

        fs::remove_file(path).unwrap();
    }
}
