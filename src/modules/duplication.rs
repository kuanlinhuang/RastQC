use crate::config::FastQCConfig;
use crate::io::Sequence;
use super::{QCModule, QCResult};
use std::collections::HashMap;

const OBSERVATION_CUTOFF: usize = 100_000;

/// Duplication level bins: 1-9 individual, then >10, >50, >100, >500, >1k, >5k, >10k
const BIN_LABELS: &[&str] = &[
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
    ">10", ">50", ">100", ">500", ">1k", ">5k", ">10k",
];

pub struct DuplicationLevel {
    sequences: HashMap<u64, u64>, // hash -> count
    total_sequences: u64,
    unique_count_at_limit: u64,
    count_at_limit: u64,
    dup_length: usize,
    reached_limit: bool,
    // Results
    total_percentages: [f64; 17],
    dedup_percentages: [f64; 17],
    percent_different: f64,
    qc_result: QCResult,
}

impl DuplicationLevel {
    pub fn new(dup_length: usize) -> Self {
        DuplicationLevel {
            sequences: HashMap::new(),
            total_sequences: 0,
            unique_count_at_limit: 0,
            count_at_limit: 0,
            dup_length,
            reached_limit: false,
            total_percentages: [0.0; 17],
            dedup_percentages: [0.0; 17],
            percent_different: 100.0,
            qc_result: QCResult::NotRun,
        }
    }

    fn hash_sequence(seq: &[u8]) -> u64 {
        // Simple FNV-1a hash
        let mut hash: u64 = 0xcbf29ce484222325;
        for &b in seq {
            hash ^= b.to_ascii_uppercase() as u64;
            hash = hash.wrapping_mul(0x100000001b3);
        }
        hash
    }

    fn dup_bin(count: u64) -> usize {
        match count {
            1 => 0,
            2 => 1,
            3 => 2,
            4 => 3,
            5 => 4,
            6 => 5,
            7 => 6,
            8 => 7,
            9 => 8,
            10 => 9,
            11..=50 => 10,
            51..=100 => 11,
            101..=500 => 12,
            501..=1000 => 13,
            1001..=5000 => 14,
            5001..=10000 => 15,
            _ => 16,
        }
    }
}

impl QCModule for DuplicationLevel {
    fn name(&self) -> &str {
        "Sequence Duplication Levels"
    }

    fn key(&self) -> &str {
        "duplication"
    }

    fn process_sequence(&mut self, seq: &Sequence) {
        self.total_sequences += 1;

        // Truncate for duplication detection
        let end = seq.sequence.len().min(self.dup_length);
        let trunc = &seq.sequence[..end];

        let hash = Self::hash_sequence(trunc);

        if self.reached_limit {
            // Only increment existing entries
            if let Some(count) = self.sequences.get_mut(&hash) {
                *count += 1;
            }
        } else {
            *self.sequences.entry(hash).or_insert(0) += 1;

            if self.sequences.len() >= OBSERVATION_CUTOFF && !self.reached_limit {
                self.reached_limit = true;
                self.unique_count_at_limit = self.sequences.len() as u64;
                self.count_at_limit = self.total_sequences;
            }
        }
    }

    fn calculate_results(&mut self, config: &FastQCConfig) {
        if self.sequences.is_empty() {
            return;
        }

        // Collate duplication counts
        let mut dedup_counts = [0.0_f64; 17]; // unique sequences per bin
        let mut total_counts = [0.0_f64; 17]; // total sequences per bin

        let correction_factor = if self.reached_limit && self.unique_count_at_limit > 0 {
            self.total_sequences as f64 / self.count_at_limit as f64
        } else {
            1.0
        };

        for (&_hash, &count) in &self.sequences {
            let corrected_count = if self.reached_limit {
                // Apply binomial correction for unseen sequences
                let p_seen = self.count_at_limit as f64 / self.total_sequences as f64;
                let corrected = (count as f64) / (1.0 - (1.0 - p_seen).powi(count as i32).max(0.0001));
                corrected.round() as u64
            } else {
                count
            };

            let bin = Self::dup_bin(corrected_count);
            dedup_counts[bin] += 1.0;
            total_counts[bin] += corrected_count as f64;
        }

        // Scale if we hit the limit
        if self.reached_limit {
            for i in 0..17 {
                dedup_counts[i] *= correction_factor;
                total_counts[i] *= correction_factor;
            }
        }

        let total_seqs: f64 = total_counts.iter().sum();
        let total_unique: f64 = dedup_counts.iter().sum();

        if total_seqs > 0.0 {
            for i in 0..17 {
                self.total_percentages[i] = total_counts[i] / total_seqs * 100.0;
                self.dedup_percentages[i] = dedup_counts[i] / total_unique.max(1.0) * 100.0;
            }
        }

        self.percent_different = if total_seqs > 0.0 {
            total_unique / total_seqs * 100.0
        } else {
            100.0
        };

        let warn = config.get_limit("duplication").map(|l| l.warn).unwrap_or(70.0);
        let error = config.get_limit("duplication").map(|l| l.error).unwrap_or(50.0);

        self.qc_result = if self.percent_different < error {
            QCResult::Fail
        } else if self.percent_different < warn {
            QCResult::Warn
        } else {
            QCResult::Pass
        };
    }

    fn result(&self) -> QCResult {
        self.qc_result
    }

    fn text_data(&self) -> String {
        let mut out = format!(
            ">>Sequence Duplication Levels\t{}\n\
             #Total Deduplicated Percentage\t{:.2}\n\
             #Duplication Level\tPercentage of deduplicated\tPercentage of total\n",
            self.result().label(),
            self.percent_different
        );

        for i in 0..17 {
            out.push_str(&format!(
                "{}\t{:.2}\t{:.2}\n",
                BIN_LABELS[i], self.dedup_percentages[i], self.total_percentages[i]
            ));
        }
        out.push_str(">>END_MODULE\n");
        out
    }

    fn svg_chart(&self) -> String {
        let width = 800.0_f64;
        let height = 400.0_f64;
        let ml = 60.0;
        let mr = 20.0;
        let mt = 30.0;
        let mb = 80.0;
        let pw = width - ml - mr;
        let ph = height - mt - mb;
        let n = 17;

        let max_pct = self
            .total_percentages
            .iter()
            .chain(self.dedup_percentages.iter())
            .copied()
            .fold(0.0_f64, f64::max)
            .max(1.0);

        let mut svg = format!(
            r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {} {}" font-family="Arial, sans-serif" font-size="11">"##,
            width, height
        );

        // Total line (blue)
        svg.push_str(r##"<polyline points=""##);
        for i in 0..n {
            let x = ml + (i as f64 + 0.5) / n as f64 * pw;
            let y = mt + ph * (1.0 - self.total_percentages[i] / max_pct);
            if i > 0 { svg.push(' '); }
            svg.push_str(&format!("{:.1},{:.1}", x, y));
        }
        svg.push_str(r##"" fill="none" stroke="#0000ff" stroke-width="2" />"##);

        // Deduplicated line (red)
        svg.push_str(r##"<polyline points=""##);
        for i in 0..n {
            let x = ml + (i as f64 + 0.5) / n as f64 * pw;
            let y = mt + ph * (1.0 - self.dedup_percentages[i] / max_pct);
            if i > 0 { svg.push(' '); }
            svg.push_str(&format!("{:.1},{:.1}", x, y));
        }
        svg.push_str(r##"" fill="none" stroke="#ff0000" stroke-width="2" />"##);

        // Axes
        svg.push_str(&format!(
            r##"<line x1="{ml}" y1="{mt}" x2="{ml}" y2="{}" stroke="black" />"##,
            mt + ph
        ));
        svg.push_str(&format!(
            r##"<line x1="{ml}" y1="{}" x2="{}" y2="{}" stroke="black" />"##,
            mt + ph, ml + pw, mt + ph
        ));

        // X labels
        for i in 0..n {
            let x = ml + (i as f64 + 0.5) / n as f64 * pw;
            let y = mt + ph + 15.0;
            svg.push_str(&format!(
                r##"<text x="{x}" y="{y}" text-anchor="end" transform="rotate(-45 {x} {y})" font-size="9">{}</text>"##,
                BIN_LABELS[i]
            ));
        }

        // Title
        svg.push_str(&format!(
            r##"<text x="{}" y="18" text-anchor="middle" font-size="13" font-weight="bold">Sequence Duplication Levels [{:.2}% remaining after dedup]</text>"##,
            width / 2.0, self.percent_different
        ));

        svg.push_str("</svg>");
        svg
    }
}
