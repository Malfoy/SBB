use anyhow::{Context, Result};
use std::fs;
use std::path::Path;

pub fn write_fasta(path: &Path, records: &[(&str, &str)]) -> Result<()> {
    let mut out = String::new();
    for (id, seq) in records {
        out.push('>');
        out.push_str(id);
        out.push('\n');
        out.push_str(seq);
        out.push('\n');
    }
    fs::write(path, out).with_context(|| format!("failed to write {}", path.display()))?;
    Ok(())
}
