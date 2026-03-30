use crate::config::FastQCConfig;
use crate::io::Sequence;
use super::{QCModule, QCResult};

pub struct PerSequenceGC {
    /// Histogram of GC percentages (0..=100)
    gc_distribution: [f64; 101],
    total_count: u64,
    // Results
    theoretical: [f64; 101],
    deviation_percent: f64,
    qc_result: QCResult,
}

impl PerSequenceGC {
    pub fn new() -> Self {
        PerSequenceGC {
            gc_distribution: [0.0; 101],
            total_count: 0,
            theoretical: [0.0; 101],
            deviation_percent: 0.0,
            qc_result: QCResult::NotRun,
        }
    }
}

impl QCModule for PerSequenceGC {
    fn name(&self) -> &str {
        "Per sequence GC content"
    }

    fn key(&self) -> &str {
        "gc_sequence"
    }

    fn process_sequence(&mut self, seq: &Sequence) {
        if seq.sequence.is_empty() {
            return;
        }

        let mut gc = 0u64;
        let mut total = 0u64;
        for &b in &seq.sequence {
            match b {
                b'G' | b'C' | b'g' | b'c' => {
                    gc += 1;
                    total += 1;
                }
                b'A' | b'T' | b'a' | b't' => {
                    total += 1;
                }
                _ => {} // Skip N
            }
        }

        if total > 0 {
            let pct = ((gc as f64 / total as f64) * 100.0).round() as usize;
            let pct = pct.min(100);
            self.gc_distribution[pct] += 1.0;
            self.total_count += 1;
        }
    }

    fn calculate_results(&mut self, config: &FastQCConfig) {
        if self.total_count == 0 {
            return;
        }

        // Find mode (peak of distribution)
        let max_count = self
            .gc_distribution
            .iter()
            .copied()
            .fold(0.0_f64, f64::max);

        if max_count == 0.0 {
            return;
        }

        // Average values near the peak to find the mode
        let threshold = max_count * 0.9;
        let mut mode_sum = 0.0_f64;
        let mut mode_weight = 0.0_f64;
        for i in 0..=100 {
            if self.gc_distribution[i] >= threshold {
                mode_sum += i as f64 * self.gc_distribution[i];
                mode_weight += self.gc_distribution[i];
            }
        }
        let mode = if mode_weight > 0.0 {
            mode_sum / mode_weight
        } else {
            50.0
        };

        // Calculate standard deviation
        let total = self.total_count as f64;
        let mut variance_sum = 0.0_f64;
        for i in 0..=100 {
            let diff = i as f64 - mode;
            variance_sum += diff * diff * self.gc_distribution[i];
        }
        let stdev = (variance_sum / (total - 1.0).max(1.0)).sqrt();

        // Generate theoretical normal distribution
        if stdev > 0.0 {
            let two_var = 2.0 * stdev * stdev;
            let norm_factor = 1.0 / (stdev * (2.0 * std::f64::consts::PI).sqrt());

            for i in 0..=100 {
                let diff = i as f64 - mode;
                self.theoretical[i] = norm_factor * (-diff * diff / two_var).exp() * total;
            }
        }

        // Calculate deviation
        let mut deviation_sum = 0.0_f64;
        for i in 0..=100 {
            deviation_sum += (self.gc_distribution[i] - self.theoretical[i]).abs();
        }
        self.deviation_percent = deviation_sum / total * 100.0;

        let warn = config.get_limit("gc_sequence").map(|l| l.warn).unwrap_or(15.0);
        let error = config.get_limit("gc_sequence").map(|l| l.error).unwrap_or(30.0);

        self.qc_result = if self.deviation_percent > error {
            QCResult::Fail
        } else if self.deviation_percent > warn {
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
            ">>Per sequence GC content\t{}\n\
             #GC Content\tCount\n",
            self.result().label()
        );

        for i in 0..=100 {
            out.push_str(&format!("{}\t{:.1}\n", i, self.gc_distribution[i]));
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
        let mb = 50.0;
        let pw = width - ml - mr;
        let ph = height - mt - mb;

        let max_count = self
            .gc_distribution
            .iter()
            .chain(self.theoretical.iter())
            .copied()
            .fold(0.0_f64, f64::max)
            .max(1.0);

        let mut svg = format!(
            r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {} {}" font-family="Arial, sans-serif" font-size="11">"##,
            width, height
        );

        // Observed distribution
        svg.push_str(r##"<polyline points=""##);
        for i in 0..=100 {
            let x = ml + i as f64 / 100.0 * pw;
            let y = mt + ph * (1.0 - self.gc_distribution[i] / max_count);
            if i > 0 { svg.push(' '); }
            svg.push_str(&format!("{:.1},{:.1}", x, y));
        }
        svg.push_str(r##"" fill="none" stroke="#ff0000" stroke-width="2" />"##);

        // Theoretical distribution
        svg.push_str(r##"<polyline points=""##);
        for i in 0..=100 {
            let x = ml + i as f64 / 100.0 * pw;
            let y = mt + ph * (1.0 - self.theoretical[i] / max_count);
            if i > 0 { svg.push(' '); }
            svg.push_str(&format!("{:.1},{:.1}", x, y));
        }
        svg.push_str(r##"" fill="none" stroke="#0000ff" stroke-width="1.5" stroke-dasharray="5,3" />"##);

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
        for v in (0..=100).step_by(10) {
            let x = ml + v as f64 / 100.0 * pw;
            svg.push_str(&format!(
                r##"<text x="{}" y="{}" text-anchor="middle" font-size="10">{}</text>"##,
                x, mt + ph + 15.0, v
            ));
        }

        // Legend
        svg.push_str(&format!(
            r##"<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="#ff0000" stroke-width="2" />"##,
            ml + pw - 150.0, mt + 10.0, ml + pw - 135.0, mt + 10.0
        ));
        svg.push_str(&format!(
            r##"<text x="{}" y="{}" font-size="10">GC count per read</text>"##,
            ml + pw - 130.0, mt + 14.0
        ));
        svg.push_str(&format!(
            r##"<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="#0000ff" stroke-width="1.5" stroke-dasharray="5,3" />"##,
            ml + pw - 150.0, mt + 25.0, ml + pw - 135.0, mt + 25.0
        ));
        svg.push_str(&format!(
            r##"<text x="{}" y="{}" font-size="10">Theoretical distribution</text>"##,
            ml + pw - 130.0, mt + 29.0
        ));

        svg.push_str(&format!(
            r##"<text x="{}" y="18" text-anchor="middle" font-size="13" font-weight="bold">GC distribution over all sequences</text>"##,
            width / 2.0
        ));

        svg.push_str("</svg>");
        svg
    }
}
