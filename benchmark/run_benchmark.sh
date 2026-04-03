#!/usr/bin/env bash
set -euo pipefail

# Benchmark: FastQC vs RastQC on real human exome FASTQ files
# ============================================================

RASTQC="/Users/kuan-lin.huang/Projects/RastQC/target/release/rastqc"
DATADIR="/Users/kuan-lin.huang/Projects/RastQC/benchmark/data"
RESULTSDIR="/Users/kuan-lin.huang/Projects/RastQC/benchmark/results"
THREADS=4

mkdir -p "$RESULTSDIR/fastqc" "$RESULTSDIR/rastqc"

# Collect all .fastq.gz files
FILES=($(ls "$DATADIR"/*.fastq.gz 2>/dev/null))
if [ ${#FILES[@]} -eq 0 ]; then
    echo "ERROR: No FASTQ files found in $DATADIR"
    exit 1
fi

echo "=============================================="
echo "  FastQC vs RastQC Benchmark"
echo "  Date: $(date)"
echo "  Threads: $THREADS"
echo "  Files: ${#FILES[@]}"
echo "=============================================="
echo ""

# Print file info
echo "--- File Sizes ---"
for f in "${FILES[@]}"; do
    fname=$(basename "$f")
    size_bytes=$(stat -f%z "$f")
    size_mb=$((size_bytes / 1024 / 1024))
    echo "  $fname: ${size_mb} MB"
done
echo ""

# CSV output
CSV="$RESULTSDIR/benchmark_results.csv"
echo "tool,file,reads,wall_sec,user_sec,sys_sec,max_rss_kb" > "$CSV"

# Function to count reads in a gzipped FASTQ
count_reads() {
    local f="$1"
    # Count lines / 4 for FASTQ
    gzip -dc "$f" | wc -l | awk '{print int($1/4)}'
}

# Benchmark function using /usr/bin/time
bench() {
    local tool="$1"
    local label="$2"
    shift 2
    local cmd=("$@")

    echo "  Running $tool on $label..."

    # Use GNU time for memory stats (macOS /usr/bin/time -l)
    local timefile=$(mktemp)
    local start=$(date +%s.%N 2>/dev/null || python3 -c "import time; print(f'{time.time():.3f}')")

    /usr/bin/time -l "${cmd[@]}" > /dev/null 2> "$timefile"
    local exit_code=$?

    local end=$(date +%s.%N 2>/dev/null || python3 -c "import time; print(f'{time.time():.3f}')")
    local wall=$(python3 -c "print(f'{$end - $start:.2f}')")

    # Parse /usr/bin/time output (macOS format)
    local user_sec=$(grep "user" "$timefile" | head -1 | awk '{print $1}')
    local sys_sec=$(grep "sys" "$timefile" | head -1 | awk '{print $1}')
    local max_rss=$(grep "maximum resident set size" "$timefile" | awk '{print $1}')

    echo "    Wall: ${wall}s | User: ${user_sec}s | Sys: ${sys_sec}s | MaxRSS: $((max_rss / 1024))MB"

    # Return results
    echo "$tool,$label,,$wall,$user_sec,$sys_sec,$max_rss" >> "$CSV"
    rm -f "$timefile"
    return $exit_code
}

# ---- Run benchmarks per file ----
for f in "${FILES[@]}"; do
    fname=$(basename "$f")
    echo ""
    echo "=== $fname ==="

    # Count reads
    echo "  Counting reads..."
    nreads=$(count_reads "$f")
    echo "  Reads: $nreads"

    # FastQC
    bench "fastqc" "$fname" fastqc -t "$THREADS" -o "$RESULTSDIR/fastqc" --quiet "$f"

    # RastQC
    bench "rastqc" "$fname" "$RASTQC" -t "$THREADS" -o "$RESULTSDIR/rastqc" -q "$f"

    echo ""
done

# ---- Run multi-file benchmark (all files at once) ----
echo ""
echo "=== ALL FILES TOGETHER ==="
echo ""

bench "fastqc" "ALL_FILES" fastqc -t "$THREADS" -o "$RESULTSDIR/fastqc" --quiet "${FILES[@]}"
bench "rastqc" "ALL_FILES" "$RASTQC" -t "$THREADS" -o "$RESULTSDIR/rastqc" -q "${FILES[@]}"

echo ""
echo "=============================================="
echo "  Results saved to: $CSV"
echo "=============================================="
echo ""
cat "$CSV" | column -t -s ','

echo ""
echo "--- Speedup Summary ---"
python3 << 'PYEOF'
import csv

results = {}
with open("/Users/kuan-lin.huang/Projects/RastQC/benchmark/results/benchmark_results.csv") as f:
    reader = csv.DictReader(f)
    for row in reader:
        key = row["file"]
        tool = row["tool"]
        if key not in results:
            results[key] = {}
        results[key][tool] = float(row["wall_sec"])

print(f"{'File':<30} {'FastQC (s)':>12} {'RastQC (s)':>12} {'Speedup':>10}")
print("-" * 68)
for fname, tools in results.items():
    if "fastqc" in tools and "rastqc" in tools:
        speedup = tools["fastqc"] / tools["rastqc"] if tools["rastqc"] > 0 else float('inf')
        print(f"{fname:<30} {tools['fastqc']:>12.2f} {tools['rastqc']:>12.2f} {speedup:>9.1f}x")
PYEOF
