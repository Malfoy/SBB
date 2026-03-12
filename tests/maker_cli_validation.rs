mod common;

use biobloom_rs::bloom::BloomFilter;
use std::process::Command;
use tempfile::tempdir;

#[test]
fn rejects_conflicting_layout_flags() {
    let d = tempdir().expect("tempdir should be creatable");
    let input = d.path().join("in.fa");
    common::write_fasta(&input, &[("r1", "ACGTACGTACGT")]).expect("input fasta should be writable");

    let out = Command::new(env!("CARGO_BIN_EXE_sbb"))
        .args([
            "maker",
            "-p",
            "t",
            "-o",
            d.path().to_str().expect("temp path should be utf-8"),
            "--blocked",
            "--classic",
            input.to_str().expect("input path should be utf-8"),
        ])
        .output()
        .expect("maker command should run");

    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("choose only one layout override"));
}

#[test]
fn rejects_non_power_of_two_block_words_in_blocked_mode() {
    let d = tempdir().expect("tempdir should be creatable");
    let input = d.path().join("in.fa");
    common::write_fasta(&input, &[("r1", "ACGTACGTACGT")]).expect("input fasta should be writable");

    let out = Command::new(env!("CARGO_BIN_EXE_sbb"))
        .args([
            "maker",
            "-p",
            "t",
            "-o",
            d.path().to_str().expect("temp path should be utf-8"),
            "--block_words",
            "3",
            input.to_str().expect("input path should be utf-8"),
        ])
        .output()
        .expect("maker command should run");

    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("--block_words must be a non-zero power-of-two"));
}

#[test]
fn rejects_invalid_progressive_threshold() {
    let d = tempdir().expect("tempdir should be creatable");
    let a = d.path().join("a.fa");
    let b = d.path().join("b.fa");
    common::write_fasta(&a, &[("r1", "ACGTACGTACGT")]).expect("seed fasta should be writable");
    common::write_fasta(&b, &[("r2", "ACGTACGTACGT")]).expect("recruit fasta should be writable");

    let out = Command::new(env!("CARGO_BIN_EXE_sbb"))
        .args([
            "maker",
            "-p",
            "t",
            "-o",
            d.path().to_str().expect("temp path should be utf-8"),
            "--progressive",
            "1.2",
            "--seed_files",
            "1",
            a.to_str().expect("seed path should be utf-8"),
            b.to_str().expect("recruit path should be utf-8"),
        ])
        .output()
        .expect("maker command should run");

    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("--progressive threshold must be in [0, 1]"));
}

#[test]
fn rejects_invalid_sample_alpha() {
    let d = tempdir().expect("tempdir should be creatable");
    let input = d.path().join("in.fa");
    common::write_fasta(&input, &[("r1", "ACGTACGTACGT")]).expect("input fasta should be writable");

    let out = Command::new(env!("CARGO_BIN_EXE_sbb"))
        .args([
            "maker",
            "-p",
            "t",
            "-o",
            d.path().to_str().expect("temp path should be utf-8"),
            "--sample_alpha",
            "1.1",
            input.to_str().expect("input path should be utf-8"),
        ])
        .output()
        .expect("maker command should run");

    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("--sample_alpha must be in [0, 1]"));
}

#[test]
fn persists_sample_alpha_in_filter_header() {
    let d = tempdir().expect("tempdir should be creatable");
    let input = d.path().join("in.fa");
    common::write_fasta(&input, &[("r1", "ACGTACGTACGTACGTACGT")])
        .expect("input fasta should be writable");

    let out = Command::new(env!("CARGO_BIN_EXE_sbb"))
        .args([
            "maker",
            "-p",
            "t",
            "-o",
            d.path().to_str().expect("temp path should be utf-8"),
            "--sample_alpha",
            "0.5",
            input.to_str().expect("input path should be utf-8"),
        ])
        .output()
        .expect("maker command should run");

    assert!(out.status.success());
    let filter_path = d.path().join("t.bf.zst");
    assert!(filter_path.exists());
    assert!(!d.path().join("t.txt").exists());
    let filter = BloomFilter::load(&filter_path).expect("written filter should load");
    assert_eq!(filter.sample_rate, 0.5);
}

#[test]
fn maker_builds_bloomybloom_filters() {
    let d = tempdir().expect("tempdir should be creatable");
    let input = d.path().join("in.fa");
    common::write_fasta(
        &input,
        &[(
            "r1",
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA",
        )],
    )
    .expect("input fasta should be writable");

    let out = Command::new(env!("CARGO_BIN_EXE_sbb"))
        .args([
            "maker",
            "-p",
            "t",
            "-o",
            d.path().to_str().expect("temp path should be utf-8"),
            "--bloomybloom",
            "-k",
            "41",
            input.to_str().expect("input path should be utf-8"),
        ])
        .output()
        .expect("maker command should run");

    assert!(out.status.success());
    let filter_path = d.path().join("t.bf.zst");
    let filter = BloomFilter::load(&filter_path).expect("written filter should load");
    assert_eq!(filter.layout_name(), "bloomybloom");
    assert_eq!(filter.kmer_size, 41);
}

#[test]
fn rejects_even_k_for_bloomybloom() {
    let d = tempdir().expect("tempdir should be creatable");
    let input = d.path().join("in.fa");
    common::write_fasta(
        &input,
        &[(
            "r1",
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA",
        )],
    )
    .expect("input fasta should be writable");

    let out = Command::new(env!("CARGO_BIN_EXE_sbb"))
        .args([
            "maker",
            "-p",
            "t",
            "-o",
            d.path().to_str().expect("temp path should be utf-8"),
            "--bloomybloom",
            "-k",
            "40",
            input.to_str().expect("input path should be utf-8"),
        ])
        .output()
        .expect("maker command should run");

    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("odd k-mer size"));
}
