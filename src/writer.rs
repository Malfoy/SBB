use crate::fastx::ReadRecord;
use anyhow::{Context, Result};
use flate2::Compression;
use flate2::write::GzEncoder;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

#[derive(Clone, Copy, Debug)]
pub enum OutputFormat {
    Fasta,
    Fastq,
}

pub fn open_writer(path: &Path, gzip: bool) -> Result<Box<dyn Write>> {
    let file =
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let buffered = BufWriter::new(file);

    if gzip {
        Ok(Box::new(GzEncoder::new(buffered, Compression::default())))
    } else {
        Ok(Box::new(buffered))
    }
}

pub fn format_score_suffix(filter_ids: &[String], scores: &[f64]) -> String {
    let mut out = String::from(" bbscores=");
    for (i, (id, score)) in filter_ids.iter().zip(scores).enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push_str(id);
        out.push(':');
        out.push_str(&format!("{score:.4}"));
    }
    out
}

pub fn write_record(
    writer: &mut dyn Write,
    record: &ReadRecord,
    format: OutputFormat,
    suffix: Option<&str>,
) -> Result<()> {
    match format {
        OutputFormat::Fasta => {
            writer.write_all(b">")?;
            writer.write_all(&record.id)?;
            if let Some(s) = suffix {
                writer.write_all(s.as_bytes())?;
            }
            writer.write_all(b"\n")?;
            writer.write_all(&record.seq)?;
            writer.write_all(b"\n")?;
        }
        OutputFormat::Fastq => {
            writer.write_all(b"@")?;
            writer.write_all(&record.id)?;
            if let Some(s) = suffix {
                writer.write_all(s.as_bytes())?;
            }
            writer.write_all(b"\n")?;
            writer.write_all(&record.seq)?;
            writer.write_all(b"\n+\n")?;

            if let Some(qual) = &record.qual {
                writer.write_all(qual)?;
            } else {
                let fallback = vec![b'I'; record.seq.len()];
                writer.write_all(&fallback)?;
            }
            writer.write_all(b"\n")?;
        }
    }

    Ok(())
}

pub struct CategorizerWriters {
    pub format: OutputFormat,
    pub filter_ids: Vec<String>,
    pub filter_writers: Vec<Box<dyn Write>>,
    pub nomatch: Box<dyn Write>,
    pub multimatch: Box<dyn Write>,
}

impl CategorizerWriters {
    pub fn new(
        prefix: &str,
        filter_ids: Vec<String>,
        format: OutputFormat,
        gzip: bool,
    ) -> Result<Self> {
        let ext = match format {
            OutputFormat::Fasta => "fa",
            OutputFormat::Fastq => "fq",
        };
        let gz_ext = if gzip { ".gz" } else { "" };

        let mut filter_writers = Vec::with_capacity(filter_ids.len());
        for id in &filter_ids {
            let path = format!("{prefix}_{id}.{ext}{gz_ext}");
            filter_writers.push(open_writer(Path::new(&path), gzip)?);
        }

        let nomatch_path = format!("{prefix}_nomatch.{ext}{gz_ext}");
        let multimatch_path = format!("{prefix}_multimatch.{ext}{gz_ext}");

        Ok(Self {
            format,
            filter_ids,
            filter_writers,
            nomatch: open_writer(Path::new(&nomatch_path), gzip)?,
            multimatch: open_writer(Path::new(&multimatch_path), gzip)?,
        })
    }

    pub fn write_assignment(
        &mut self,
        record: &ReadRecord,
        matches: &[usize],
        suffix: Option<&str>,
    ) -> Result<()> {
        if matches.is_empty() {
            return write_record(self.nomatch.as_mut(), record, self.format, suffix);
        }

        if matches.len() > 1 {
            write_record(self.multimatch.as_mut(), record, self.format, suffix)?;
        }

        for &idx in matches {
            write_record(
                self.filter_writers[idx].as_mut(),
                record,
                self.format,
                suffix,
            )?;
        }

        Ok(())
    }

    pub fn flush_all(&mut self) -> Result<()> {
        for writer in &mut self.filter_writers {
            writer.flush()?;
        }
        self.nomatch.flush()?;
        self.multimatch.flush()?;
        Ok(())
    }
}
