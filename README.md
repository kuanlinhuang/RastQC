# RastQC

A fast quality control tool for high-throughput sequencing data, written in Rust. Drop-in replacement for [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) with identical QC modules, matching algorithms, and compatible output formats.

## Features

- **12 QC modules** matching FastQC's analysis pipeline
- **Fast**: multi-threaded parallel file processing via rayon
- **Portable**: single static binary, no Java runtime needed
- **Compatible output**: HTML reports, tab-separated data files, ZIP archives
- **Multi-file summary**: generates an overview report when processing many files
- **Input formats**: FASTQ (`.fastq`, `.fq`), gzip (`.gz`), bzip2 (`.bz2`), BAM, SAM

## Installation

### From source

```bash
# Requires Rust 1.70+
cargo install --path .
```

### Build manually

```bash
git clone https://github.com/huangk06/RastQC.git
cd RastQC
cargo build --release
# Binary at ./target/release/rastqc
```

## Quick start

```bash
# Single file
rastqc sample.fastq.gz

# Multiple files (processed in parallel)
rastqc *.fastq.gz

# Specify output directory
rastqc -o results/ sample_R1.fastq.gz sample_R2.fastq.gz

# HTML only (no ZIP)
rastqc --nozip -o results/ sample.fastq.gz

# Use 8 threads
rastqc -t 8 -o results/ *.fastq.gz
```

## Usage

```
rastqc [OPTIONS] <FILES>...

Arguments:
  <FILES>...  Input files (FASTQ, BAM, SAM)

Options:
  -o, --outdir <DIR>            Output directory [default: current directory]
  -t, --threads <N>             Number of threads [default: all CPUs]
  -c, --contaminants <FILE>     Custom contaminant list (tab-separated: name\tsequence)
  -a, --adapters <FILE>         Custom adapter list (tab-separated: name\tsequence)
  -l, --limits <FILE>           Custom pass/warn/fail thresholds
  -k, --kmer-size <N>           Kmer size for enrichment analysis [default: 7]
      --nofilter                Include all reads (don't skip QC-failed reads)
      --extract                 Extract ZIP contents after creation
      --nozip                   Write HTML report only, skip ZIP archive
      --summary                 Write multi-file summary report (summary.html + summary.tsv)
  -q, --quiet                   Suppress progress output
      --dup-length <N>          Truncation length for duplication detection [default: 50]
  -h, --help                    Print help
  -V, --version                 Print version
```

## Output files

For each input file `sample.fastq.gz`, RastQC produces:

| File | Description |
|------|-------------|
| `sample_fastqc.zip` | ZIP archive containing all outputs below |
| `sample_fastqc/fastqc_report.html` | Self-contained HTML report with SVG charts |
| `sample_fastqc/fastqc_data.txt` | Tab-separated data for each module |
| `sample_fastqc/summary.txt` | One-line PASS/WARN/FAIL per module |

When processing multiple files with `--summary`:

| File | Description |
|------|-------------|
| `summary.tsv` | Tab-separated matrix: rows = files, columns = modules, cells = PASS/WARN/FAIL |
| `summary.html` | Overview dashboard linking to all individual reports |

## QC modules

| # | Module | What it checks | Pass/Warn/Fail criteria |
|---|--------|---------------|------------------------|
| 1 | **Basic Statistics** | Sequence count, length, %GC, encoding | Informational only |
| 2 | **Per Base Sequence Quality** | Quality score distribution at each position | Median < 25 (warn) / < 20 (fail); LQ < 10 (warn) / < 5 (fail) |
| 3 | **Per Tile Sequence Quality** | Quality variation between flowcell tiles | Max deviation > 5 (warn) / > 10 (fail) |
| 4 | **Per Sequence Quality Scores** | Distribution of mean quality per read | Mode <= 27 (warn) / <= 20 (fail) |
| 5 | **Per Base Sequence Content** | A/T/G/C proportions at each position | |A-T| or |G-C| > 10% (warn) / > 20% (fail) |
| 6 | **Per Sequence GC Content** | GC% distribution vs theoretical normal | Deviation > 15% (warn) / > 30% (fail) |
| 7 | **Per Base N Content** | Unknown base (N) frequency per position | N% > 5 (warn) / > 20 (fail) |
| 8 | **Sequence Length Distribution** | Read length variability | Variable lengths (warn); zero-length reads (fail) |
| 9 | **Sequence Duplication Levels** | Library complexity estimate | < 70% unique (warn) / < 50% unique (fail) |
| 10 | **Overrepresented Sequences** | Frequently occurring sequences + contaminant matching | Any seq > 0.1% (warn) / > 1% (fail) |
| 11 | **Adapter Content** | Known adapter sequence contamination | > 5% (warn) / > 10% (fail) |
| 12 | **Kmer Content** | Positionally biased k-mers | -log10(p) > 2 (warn) / > 5 (fail) |

## Working with many files

### Batch processing

Process hundreds of files in parallel. RastQC automatically uses all available CPU cores:

```bash
# Process all FASTQ files in a directory
rastqc -o qc_results/ data/*.fastq.gz

# Process with summary dashboard
rastqc -o qc_results/ --summary data/*.fastq.gz

# Use find for recursive discovery
find data/ -name "*.fastq.gz" | xargs rastqc -o qc_results/ --summary
```

### Summary report

The `--summary` flag generates two files for multi-file review:

**`summary.tsv`** -- machine-readable matrix for scripting:
```
Sample	Basic Statistics	Per Base Quality	...	Adapter Content
sample_A	PASS	PASS	...	WARN
sample_B	PASS	FAIL	...	PASS
```

**`summary.html`** -- browser-friendly dashboard:
- Color-coded PASS/WARN/FAIL table for all samples and modules
- Click any cell to jump to the full report for that sample
- Quickly spot problematic samples across a run

### Filtering results

```bash
# Find all failing samples
grep "FAIL" qc_results/summary.tsv

# Count warnings per sample
awk -F'\t' '{n=0; for(i=2;i<=NF;i++) if($i=="WARN") n++; print $1, n}' qc_results/summary.tsv
```

## Custom configuration

### Adapter list

Tab-separated file with adapter name and 12bp sequence:

```
My Custom Adapter	AGATCGGAAGAG
Another Adapter		CTGTCTCTTATA
```

### Contaminant list

Tab-separated file with contaminant name and full sequence:

```
PhiX Control	GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACT
Custom Primer	AATGATACGGCGACCACCGA
```

### Limits file

Controls pass/warn/fail thresholds and which modules run:

```
# Disable a module
kmer    ignore  1

# Adjust thresholds
quality_base_lower  warn    10
quality_base_lower  error   5
quality_base_median warn    25
quality_base_median error   20
adapter             warn    5
adapter             error   10
```

## Compatibility with FastQC

RastQC produces output compatible with tools that consume FastQC results:

- **MultiQC**: `fastqc_data.txt` files are compatible with MultiQC's FastQC module
- **summary.txt**: same PASS/WARN/FAIL format per module
- **Identical module names and data headers** in text output

## Performance

| Metric | RastQC | FastQC (Java) |
|--------|--------|---------------|
| Binary size | ~2 MB | ~300 MB (with JRE) |
| Startup time | instant | ~2s JVM warmup |
| Threading | per-file parallel | per-file parallel |
| Memory | streaming (low) | streaming (low) |

## License

MIT
