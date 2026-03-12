use bzip2::Compression as BzCompression;
use bzip2::write::BzEncoder;
use flate2::Compression as GzCompression;
use flate2::write::GzEncoder;
use liblzma::write::XzEncoder;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;
use tempfile::tempdir;

fn write_bytes(path: &Path, data: &[u8]) {
    fs::write(path, data).expect("test input should be writable");
}

fn gzip_bytes(data: &[u8]) -> Vec<u8> {
    let mut enc = GzEncoder::new(Vec::new(), GzCompression::default());
    enc.write_all(data).expect("gzip write should succeed");
    enc.finish().expect("gzip finish should succeed")
}

fn zstd_bytes(data: &[u8]) -> Vec<u8> {
    let mut enc =
        zstd::stream::write::Encoder::new(Vec::new(), 0).expect("zstd encoder should be creatable");
    enc.write_all(data).expect("zstd write should succeed");
    enc.finish().expect("zstd finish should succeed")
}

fn xz_bytes(data: &[u8]) -> Vec<u8> {
    let mut enc = XzEncoder::new(Vec::new(), 6);
    enc.write_all(data).expect("xz write should succeed");
    enc.finish().expect("xz finish should succeed")
}

fn bz2_bytes(data: &[u8]) -> Vec<u8> {
    let mut enc = BzEncoder::new(Vec::new(), BzCompression::default());
    enc.write_all(data).expect("bz2 write should succeed");
    enc.finish().expect("bz2 finish should succeed")
}

fn run_maker(input: &Path, out_dir: &Path, prefix: &str) {
    let out = Command::new(env!("CARGO_BIN_EXE_sbb"))
        .args([
            "maker",
            "-p",
            prefix,
            "-o",
            out_dir.to_str().expect("temp path should be utf-8"),
            input.to_str().expect("input path should be utf-8"),
        ])
        .output()
        .expect("maker command should run");

    if !out.status.success() {
        let stderr = String::from_utf8_lossy(&out.stderr);
        panic!(
            "maker failed for {} with stderr:\n{}",
            input.display(),
            stderr
        );
    }
    let filter_path = out_dir.join(format!("{prefix}.bf.zst"));
    assert!(
        filter_path.exists(),
        "expected output filter {}",
        filter_path.display()
    );
}

#[test]
fn maker_accepts_all_supported_compressions_for_fasta() {
    let d = tempdir().expect("tempdir should be creatable");
    let fasta = b">r1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";

    let inputs: Vec<(&str, Vec<u8>)> = vec![
        ("fa", fasta.to_vec()),
        ("fa.gz", gzip_bytes(fasta)),
        ("fa.zst", zstd_bytes(fasta)),
        ("fa.xz", xz_bytes(fasta)),
        ("fa.bz2", bz2_bytes(fasta)),
    ];

    for (i, (ext, bytes)) in inputs.into_iter().enumerate() {
        let input_path: PathBuf = d.path().join(format!("in_{i}.{ext}"));
        write_bytes(&input_path, &bytes);
        run_maker(&input_path, d.path(), &format!("fasta_{i}"));
    }
}

#[test]
fn maker_accepts_all_supported_compressions_for_fastq() {
    let d = tempdir().expect("tempdir should be creatable");
    let fastq = b"@r1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";

    let inputs: Vec<(&str, Vec<u8>)> = vec![
        ("fq", fastq.to_vec()),
        ("fq.gz", gzip_bytes(fastq)),
        ("fq.zst", zstd_bytes(fastq)),
        ("fq.xz", xz_bytes(fastq)),
        ("fq.bz2", bz2_bytes(fastq)),
    ];

    for (i, (ext, bytes)) in inputs.into_iter().enumerate() {
        let input_path: PathBuf = d.path().join(format!("in_{i}.{ext}"));
        write_bytes(&input_path, &bytes);
        run_maker(&input_path, d.path(), &format!("fastq_{i}"));
    }
}
