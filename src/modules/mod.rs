pub mod basic_stats;
pub mod per_base_quality;
pub mod per_tile_quality;
pub mod per_sequence_quality;
pub mod per_base_content;
pub mod per_sequence_gc;
pub mod n_content;
pub mod sequence_length;
pub mod duplication;
pub mod overrepresented;
pub mod adapter_content;
pub mod kmer_content;

use crate::config::FastQCConfig;
use crate::io::Sequence;

/// QC result status
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum QCResult {
    Pass,
    Warn,
    Fail,
    NotRun,
}

impl QCResult {
    pub fn icon(&self) -> &'static str {
        match self {
            QCResult::Pass => "pass",
            QCResult::Warn => "warn",
            QCResult::Fail => "fail",
            QCResult::NotRun => "na",
        }
    }

    pub fn label(&self) -> &'static str {
        match self {
            QCResult::Pass => "PASS",
            QCResult::Warn => "WARN",
            QCResult::Fail => "FAIL",
            QCResult::NotRun => "N/A",
        }
    }
}

/// Trait for all QC modules
pub trait QCModule: Send {
    /// Module name for display
    fn name(&self) -> &str;

    /// Short key for config lookup
    #[allow(dead_code)]
    fn key(&self) -> &str;

    /// Process one sequence
    fn process_sequence(&mut self, seq: &Sequence);

    /// Calculate final results after all sequences processed
    fn calculate_results(&mut self, config: &FastQCConfig);

    /// Get pass/warn/fail status
    fn result(&self) -> QCResult;

    /// Generate tab-separated text data for this module
    fn text_data(&self) -> String;

    /// Generate SVG chart content
    fn svg_chart(&self) -> String;

    /// Whether this module has a chart
    fn has_chart(&self) -> bool {
        true
    }
}

/// Base group for aggregating positions
#[derive(Debug, Clone)]
pub struct BaseGroup {
    pub start: usize,
    pub end: usize,
}

impl BaseGroup {
    pub fn label(&self) -> String {
        if self.start == self.end {
            format!("{}", self.start + 1)
        } else {
            format!("{}-{}", self.start + 1, self.end + 1)
        }
    }

    /// Create base groups for a given max length, using adaptive binning
    pub fn make_groups(max_length: usize) -> Vec<BaseGroup> {
        if max_length == 0 {
            return vec![];
        }

        let mut groups = Vec::new();

        if max_length <= 75 {
            // Individual positions
            for i in 0..max_length {
                groups.push(BaseGroup { start: i, end: i });
            }
        } else {
            // First 9 individually
            let individual = 9.min(max_length);
            for i in 0..individual {
                groups.push(BaseGroup { start: i, end: i });
            }

            // Then groups of increasing size
            let mut pos = individual;
            let mut group_size = if max_length <= 200 {
                2
            } else if max_length <= 500 {
                5
            } else {
                10
            };

            while pos < max_length {
                let end = (pos + group_size - 1).min(max_length - 1);
                groups.push(BaseGroup {
                    start: pos,
                    end,
                });
                pos = end + 1;

                // Increase group size at thresholds
                if pos >= 50 && group_size < 5 {
                    group_size = 5;
                }
                if pos >= 100 && group_size < 10 {
                    group_size = 10;
                }
                if pos >= 500 && group_size < 50 {
                    group_size = 50;
                }
            }
        }

        groups
    }
}

/// Phred encoding detection
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PhredEncoding {
    Sanger,      // Phred+33
    Illumina1_3, // Phred+64
    Illumina1_5, // Phred+64
}

impl PhredEncoding {
    pub fn detect(lowest_char: u8) -> Self {
        if lowest_char < 64 {
            PhredEncoding::Sanger
        } else if lowest_char == 65 {
            PhredEncoding::Illumina1_3
        } else {
            PhredEncoding::Illumina1_5
        }
    }

    pub fn offset(&self) -> u8 {
        match self {
            PhredEncoding::Sanger => 33,
            PhredEncoding::Illumina1_3 | PhredEncoding::Illumina1_5 => 64,
        }
    }

    pub fn name(&self) -> &str {
        match self {
            PhredEncoding::Sanger => "Sanger / Illumina 1.9",
            PhredEncoding::Illumina1_3 => "Illumina 1.3",
            PhredEncoding::Illumina1_5 => "Illumina 1.5",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_base_group_single_position() {
        let groups = BaseGroup::make_groups(5);
        assert_eq!(groups.len(), 5);
        assert_eq!(groups[0].label(), "1");
        assert_eq!(groups[4].label(), "5");
    }

    #[test]
    fn test_base_group_short_reads() {
        let groups = BaseGroup::make_groups(50);
        assert_eq!(groups.len(), 50); // all individual
    }

    #[test]
    fn test_base_group_long_reads() {
        let groups = BaseGroup::make_groups(150);
        assert!(groups.len() < 150); // should bin
        assert_eq!(groups[0].start, 0);
        assert_eq!(groups[0].end, 0); // first 9 individual
        assert_eq!(groups.last().unwrap().end, 149);
    }

    #[test]
    fn test_base_group_empty() {
        let groups = BaseGroup::make_groups(0);
        assert!(groups.is_empty());
    }

    #[test]
    fn test_phred_sanger() {
        assert_eq!(PhredEncoding::detect(33), PhredEncoding::Sanger);
        assert_eq!(PhredEncoding::detect(50), PhredEncoding::Sanger);
        assert_eq!(PhredEncoding::Sanger.offset(), 33);
    }

    #[test]
    fn test_phred_illumina() {
        assert_eq!(PhredEncoding::detect(65), PhredEncoding::Illumina1_3);
        assert_eq!(PhredEncoding::detect(66), PhredEncoding::Illumina1_5);
        assert_eq!(PhredEncoding::Illumina1_3.offset(), 64);
    }

    #[test]
    fn test_qc_result_labels() {
        assert_eq!(QCResult::Pass.label(), "PASS");
        assert_eq!(QCResult::Warn.label(), "WARN");
        assert_eq!(QCResult::Fail.label(), "FAIL");
    }
}

/// Factory to create all modules
pub struct ModuleFactory;

impl ModuleFactory {
    pub fn create_modules(config: &FastQCConfig) -> Vec<Box<dyn QCModule>> {
        let mut modules: Vec<Box<dyn QCModule>> = Vec::new();

        modules.push(Box::new(basic_stats::BasicStats::new()));

        if !config.is_ignored("quality_base") {
            modules.push(Box::new(per_base_quality::PerBaseQuality::new()));
        }
        if !config.is_ignored("tile") {
            modules.push(Box::new(per_tile_quality::PerTileQuality::new()));
        }
        if !config.is_ignored("quality_sequence") {
            modules.push(Box::new(per_sequence_quality::PerSequenceQuality::new()));
        }
        if !config.is_ignored("sequence") {
            modules.push(Box::new(per_base_content::PerBaseContent::new()));
        }
        if !config.is_ignored("gc_sequence") {
            modules.push(Box::new(per_sequence_gc::PerSequenceGC::new()));
        }
        if !config.is_ignored("n_content") {
            modules.push(Box::new(n_content::NContent::new()));
        }
        if !config.is_ignored("sequence_length") {
            modules.push(Box::new(sequence_length::SequenceLengthDist::new()));
        }
        if !config.is_ignored("duplication") {
            modules.push(Box::new(duplication::DuplicationLevel::new(config.dup_length)));
        }
        if !config.is_ignored("overrepresented") {
            modules.push(Box::new(overrepresented::OverrepresentedSeqs::new(config.dup_length)));
        }
        if !config.is_ignored("adapter") {
            modules.push(Box::new(adapter_content::AdapterContent::new(config)));
        }
        if !config.is_ignored("kmer") {
            modules.push(Box::new(kmer_content::KmerContent::new(config.kmer_size)));
        }

        modules
    }
}
