use anyhow::{bail, Context, Result};
use bzip2::read::BzDecoder;
use flate2::read::GzDecoder;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

use super::colorspace;
use super::Sequence;

pub struct FastqReader {
    reader: Box<dyn BufRead>,
    line_buf: String,
    colorspace_detected: Option<bool>,
}

impl FastqReader {
    /// Open a FASTQ file (optionally compressed with gzip or bzip2).
    pub fn open(path: &Path) -> Result<Self> {
        let name = path
            .file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .to_lowercase();

        let file =
            File::open(path).with_context(|| format!("Cannot open file: {}", path.display()))?;

        let reader: Box<dyn Read> = if name.ends_with(".gz") {
            Box::new(GzDecoder::new(file))
        } else if name.ends_with(".bz2") {
            Box::new(BzDecoder::new(file))
        } else {
            Box::new(file)
        };

        Ok(FastqReader {
            reader: Box::new(BufReader::with_capacity(1024 * 1024, reader)),
            line_buf: String::with_capacity(512),
            colorspace_detected: None,
        })
    }

    /// Create a reader from stdin.
    pub fn from_stdin() -> Self {
        FastqReader {
            reader: Box::new(BufReader::with_capacity(1024 * 1024, io::stdin())),
            line_buf: String::with_capacity(512),
            colorspace_detected: None,
        }
    }

    pub fn next_sequence(&mut self) -> Result<Option<Sequence>> {
        // Read header line (starts with @)
        self.line_buf.clear();
        if self.reader.read_line(&mut self.line_buf)? == 0 {
            return Ok(None);
        }

        // Skip blank lines
        while self.line_buf.trim().is_empty() {
            self.line_buf.clear();
            if self.reader.read_line(&mut self.line_buf)? == 0 {
                return Ok(None);
            }
        }

        let header = self.line_buf.trim().to_string();
        if !header.starts_with('@') {
            bail!("Expected FASTQ header starting with '@', got: {}", header);
        }

        // Check for Casava filtered flag
        let filtered = header.contains(":Y:");

        // Read sequence line
        self.line_buf.clear();
        if self.reader.read_line(&mut self.line_buf)? == 0 {
            bail!("Unexpected end of file reading sequence");
        }
        let sequence = self.line_buf.trim().as_bytes().to_vec();

        // Read separator line (+)
        self.line_buf.clear();
        if self.reader.read_line(&mut self.line_buf)? == 0 {
            bail!("Unexpected end of file reading separator");
        }
        if !self.line_buf.trim().starts_with('+') {
            bail!(
                "Expected '+' separator, got: {}",
                self.line_buf.trim()
            );
        }

        // Read quality line
        self.line_buf.clear();
        if self.reader.read_line(&mut self.line_buf)? == 0 {
            bail!("Unexpected end of file reading quality");
        }
        let quality = self.line_buf.trim().as_bytes().to_vec();

        // Auto-detect colorspace on first sequence
        if self.colorspace_detected.is_none() {
            self.colorspace_detected = Some(colorspace::is_colorspace(&sequence));
        }

        // Decode colorspace if detected
        let sequence = if self.colorspace_detected == Some(true) {
            colorspace::decode_colorspace(&sequence).unwrap_or(sequence)
        } else {
            sequence
        };

        Ok(Some(Sequence {
            header,
            sequence,
            quality,
            filtered,
        }))
    }
}
