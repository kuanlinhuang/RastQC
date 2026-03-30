use crate::config::FastQCConfig;
use crate::io::Sequence;
use super::{QCModule, QCResult};
use std::collections::HashMap;

pub struct SequenceLengthDist {
    length_counts: HashMap<usize, u64>,
    min_length: usize,
    max_length: usize,
    has_zero_length: bool,
    // Results
    labels: Vec<String>,
    counts: Vec<f64>,
    multiple_lengths: bool,
    qc_result: QCResult,
}

impl SequenceLengthDist {
    pub fn new() -> Self {
        SequenceLengthDist {
            length_counts: HashMap::new(),
            min_length: usize::MAX,
            max_length: 0,
            has_zero_length: false,
            labels: Vec::new(),
            counts: Vec::new(),
            multiple_lengths: false,
            qc_result: QCResult::NotRun,
        }
    }

    /// Find a nice bin interval to get ~50 categories
    fn find_interval(min: usize, max: usize) -> usize {
        if min == max {
            return 1;
        }
        let range = max - min;
        let divisors = [1, 2, 5];
        let mut power = 1usize;
        loop {
            for &d in &divisors {
                let interval = d * power;
                if range / interval <= 50 {
                    return interval;
                }
            }
            power *= 10;
            if power > 1_000_000 {
                return power;
            }
        }
    }
}

impl QCModule for SequenceLengthDist {
    fn name(&self) -> &str {
        "Sequence Length Distribution"
    }

    fn key(&self) -> &str {
        "sequence_length"
    }

    fn process_sequence(&mut self, seq: &Sequence) {
        let len = seq.len();
        *self.length_counts.entry(len).or_insert(0) += 1;

        if len == 0 {
            self.has_zero_length = true;
        }
        if len < self.min_length {
            self.min_length = len;
        }
        if len > self.max_length {
            self.max_length = len;
        }
    }

    fn calculate_results(&mut self, config: &FastQCConfig) {
        if self.length_counts.is_empty() {
            return;
        }

        let min = if self.min_length == usize::MAX { 0 } else { self.min_length };
        let max = self.max_length;

        self.multiple_lengths = self.length_counts.len() > 1;

        let interval = Self::find_interval(min, max);

        // Build binned distribution
        let bin_start = (min / interval) * interval;
        let mut pos = bin_start;
        while pos <= max {
            let bin_end = (pos + interval - 1).min(max);
            let mut count = 0u64;
            for len in pos..=bin_end {
                count += self.length_counts.get(&len).copied().unwrap_or(0);
            }

            let label = if interval == 1 || pos == bin_end {
                format!("{}", pos)
            } else {
                format!("{}-{}", pos, bin_end)
            };

            self.labels.push(label);
            self.counts.push(count as f64);
            pos += interval;
        }

        let warn_on = config.get_limit("sequence_length").map(|l| l.warn != 0.0).unwrap_or(true);
        let error_on = config.get_limit("sequence_length").map(|l| l.error != 0.0).unwrap_or(true);

        self.qc_result = if self.has_zero_length && error_on {
            QCResult::Fail
        } else if self.multiple_lengths && warn_on {
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
            ">>Sequence Length Distribution\t{}\n\
             #Length\tCount\n",
            self.result().label()
        );

        for (i, label) in self.labels.iter().enumerate() {
            out.push_str(&format!("{}\t{:.1}\n", label, self.counts[i]));
        }
        out.push_str(">>END_MODULE\n");
        out
    }

    fn svg_chart(&self) -> String {
        if self.labels.is_empty() {
            return String::new();
        }

        let width = 800.0_f64;
        let height = 400.0_f64;
        let ml = 60.0;
        let mr = 20.0;
        let mt = 30.0;
        let mb = 80.0;
        let pw = width - ml - mr;
        let ph = height - mt - mb;
        let n = self.labels.len();

        let max_count = self.counts.iter().copied().fold(0.0_f64, f64::max).max(1.0);

        let mut svg = format!(
            r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {} {}" font-family="Arial, sans-serif" font-size="11">"##,
            width, height
        );

        svg.push_str(r##"<polyline points=""##);
        for (i, &count) in self.counts.iter().enumerate() {
            let x = ml + (i as f64 + 0.5) / n as f64 * pw;
            let y = mt + ph * (1.0 - count / max_count);
            if i > 0 { svg.push(' '); }
            svg.push_str(&format!("{:.1},{:.1}", x, y));
        }
        svg.push_str(r##"" fill="none" stroke="#ff0000" stroke-width="1.5" />"##);

        svg.push_str(&format!(
            r##"<line x1="{ml}" y1="{mt}" x2="{ml}" y2="{}" stroke="black" />"##,
            mt + ph
        ));
        svg.push_str(&format!(
            r##"<line x1="{ml}" y1="{}" x2="{}" y2="{}" stroke="black" />"##,
            mt + ph, ml + pw, mt + ph
        ));

        // X labels
        let step = (n / 15).max(1);
        for i in (0..n).step_by(step) {
            let x = ml + (i as f64 + 0.5) / n as f64 * pw;
            let y = mt + ph + 15.0;
            svg.push_str(&format!(
                r##"<text x="{x}" y="{y}" text-anchor="end" transform="rotate(-45 {x} {y})" font-size="9">{}</text>"##,
                self.labels[i]
            ));
        }

        svg.push_str(&format!(
            r##"<text x="{}" y="18" text-anchor="middle" font-size="13" font-weight="bold">Distribution of sequence lengths over all sequences</text>"##,
            width / 2.0
        ));

        svg.push_str("</svg>");
        svg
    }
}
