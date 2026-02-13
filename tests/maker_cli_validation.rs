mod common;

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
