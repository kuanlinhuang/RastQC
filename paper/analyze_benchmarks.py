#!/usr/bin/env python3
"""
Analyze RastQC benchmark results and generate:
1. Markdown tables for manuscript
2. SVG figures for key performance metrics
3. Concordance summary
"""

import os
import sys

PAPER_DIR = os.path.dirname(os.path.abspath(__file__))
BENCH_DIR = os.path.join(PAPER_DIR, "benchmarks")

# ─── Load Data ────────────────────────────────────────────────────────────────

def load_tsv(path):
    rows = []
    with open(path) as f:
        header = f.readline().strip().split("\t")
        for line in f:
            vals = line.strip().split("\t")
            if len(vals) >= len(header):
                rows.append(dict(zip(header, vals)))
    return rows

# ─── Synthetic Benchmark Table ────────────────────────────────────────────────

def generate_synthetic_table():
    path = os.path.join(BENCH_DIR, "benchmark_results_v2.tsv")
    rows = load_tsv(path)

    datasets = [
        ("startup", "Startup (1 read)", 1),
        ("100K", "100K synthetic", 100_000),
        ("1M", "1M synthetic", 1_000_000),
        ("1M_gz", "1M gzipped", 1_000_000),
    ]

    print("## Table 3: Performance comparison on synthetic datasets (single-threaded)")
    print()
    print("| Dataset | Reads | RastQC Time | FastQC Time | Speedup | RastQC RSS | FastQC RSS |")
    print("|---------|-------|-------------|-------------|---------|------------|------------|")

    for key, label, reads in datasets:
        r = [x for x in rows if x["Tool"] == "RastQC" and x["Dataset"] == key]
        f = [x for x in rows if x["Tool"] == "FastQC" and x["Dataset"] == key]
        if not r or not f:
            continue
        r, f = r[0], f[0]
        rt = float(r["Time_s"])
        ft = float(f["Time_s"])
        speedup = ft / rt if rt > 0 else float("inf")
        r_mb = int(r["RSS_bytes"]) / 1_048_576
        f_mb = int(f["RSS_bytes"]) / 1_048_576

        if rt < 0.01:
            rt_str = f"{rt*1000:.0f} ms"
        else:
            rt_str = f"{rt:.2f} s"
        ft_str = f"{ft:.2f} s"
        sp_str = f"{speedup:.1f}x" if speedup < 1000 else f"{speedup:.0f}x"

        print(f"| {label} | {reads:,} | {rt_str} | {ft_str} | {sp_str} | {r_mb:.0f} MB | {f_mb:.0f} MB |")

    print()

# ─── Real Genome Benchmark Table ─────────────────────────────────────────────

def generate_real_genome_table():
    path = os.path.join(BENCH_DIR, "real_genomes", "real_genome_benchmarks.tsv")
    rows = load_tsv(path)

    organisms = [
        ("fly", "*D. melanogaster*"),
        ("yeast", "*S. cerevisiae*"),
        ("ecoli", "*E. coli* K-12"),
        ("mouse", "*M. musculus*"),
        ("human", "*H. sapiens*"),
    ]

    print("## Table 4: Performance comparison on real genome data")
    print()
    print("| Organism | Reads | Read Len | RastQC Time | FastQC Time | Speedup | RastQC RSS | FastQC RSS | Mem. Ratio |")
    print("|----------|-------|----------|-------------|-------------|---------|------------|------------|------------|")

    for key, label in organisms:
        r = [x for x in rows if x["Tool"] == "RastQC" and x["Organism"] == key]
        f = [x for x in rows if x["Tool"] == "FastQC" and x["Organism"] == key]
        if not r or not f:
            continue
        r, f = r[0], f[0]
        rt = float(r["Time_s"])
        ft = float(f["Time_s"])
        speedup = ft / rt if rt > 0 else 0
        r_mb = int(r["RSS_bytes"]) / 1_048_576
        f_mb = int(f["RSS_bytes"]) / 1_048_576
        mem_ratio = f_mb / r_mb if r_mb > 0 else 0
        reads = int(r["Reads"])
        readlen = r["ReadLength"]

        print(f"| {label} | {reads:,} | {readlen} bp | {rt:.2f} s | {ft:.2f} s | {speedup:.1f}x | {r_mb:.0f} MB | {f_mb:.0f} MB | {mem_ratio:.0f}x |")

    print()

# ─── Concordance Table ────────────────────────────────────────────────────────

def generate_concordance_table():
    path = os.path.join(BENCH_DIR, "concordance_v2", "concordance_results.tsv")
    rows = load_tsv(path)

    print("## Table 6: Module-level concordance across five model organisms")
    print()
    print("| Module | Concordant | Discordant |")
    print("|--------|-----------|------------|")

    for row in rows:
        print(f"| {row['Module']} | {row['Concordant']} | {row['Discordant']} |")

    print()

# ─── SVG Figure Generation ───────────────────────────────────────────────────

def generate_speedup_figure():
    """Generate SVG bar chart comparing RastQC vs FastQC speed."""
    path = os.path.join(BENCH_DIR, "real_genomes", "real_genome_benchmarks.tsv")
    rows = load_tsv(path)

    organisms = ["fly", "yeast", "ecoli", "mouse", "human"]
    labels = ["Fly", "Yeast", "E. coli", "Mouse", "Human"]

    rastqc_times = []
    fastqc_times = []
    for org in organisms:
        r = [x for x in rows if x["Tool"] == "RastQC" and x["Organism"] == org][0]
        f = [x for x in rows if x["Tool"] == "FastQC" and x["Organism"] == org][0]
        rastqc_times.append(float(r["Time_s"]))
        fastqc_times.append(float(f["Time_s"]))

    max_time = max(max(rastqc_times), max(fastqc_times))

    w, h = 700, 400
    margin_left, margin_bottom, margin_top, margin_right = 80, 60, 40, 120
    plot_w = w - margin_left - margin_right
    plot_h = h - margin_top - margin_bottom
    n = len(organisms)
    group_w = plot_w / n
    bar_w = group_w * 0.35

    svg = []
    svg.append(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" font-family="Arial, sans-serif">')
    svg.append(f'<rect width="{w}" height="{h}" fill="white"/>')

    # Title
    svg.append(f'<text x="{w/2}" y="22" text-anchor="middle" font-size="14" font-weight="bold">Wall-Clock Time: RastQC vs FastQC (Real Genomes)</text>')

    # Y axis
    y_ticks = [0, 5, 10, 15, 20, 25, 30]
    for t in y_ticks:
        y = margin_top + plot_h - (t / max(y_ticks) * plot_h)
        svg.append(f'<line x1="{margin_left}" y1="{y}" x2="{margin_left + plot_w}" y2="{y}" stroke="#e0e0e0" stroke-width="0.5"/>')
        svg.append(f'<text x="{margin_left - 8}" y="{y + 4}" text-anchor="end" font-size="11" fill="#555">{t}s</text>')

    # Bars
    colors_r, colors_f = "#2980b9", "#e74c3c"
    for i in range(n):
        x_center = margin_left + group_w * i + group_w / 2
        x_r = x_center - bar_w - 2
        x_f = x_center + 2

        h_r = (rastqc_times[i] / max(y_ticks)) * plot_h
        h_f = (fastqc_times[i] / max(y_ticks)) * plot_h
        y_base = margin_top + plot_h

        svg.append(f'<rect x="{x_r}" y="{y_base - h_r}" width="{bar_w}" height="{h_r}" fill="{colors_r}" rx="2"/>')
        svg.append(f'<rect x="{x_f}" y="{y_base - h_f}" width="{bar_w}" height="{h_f}" fill="{colors_f}" rx="2"/>')

        # Value labels
        svg.append(f'<text x="{x_r + bar_w/2}" y="{y_base - h_r - 4}" text-anchor="middle" font-size="9" fill="{colors_r}">{rastqc_times[i]:.1f}s</text>')
        svg.append(f'<text x="{x_f + bar_w/2}" y="{y_base - h_f - 4}" text-anchor="middle" font-size="9" fill="{colors_f}">{fastqc_times[i]:.1f}s</text>')

        # X label
        svg.append(f'<text x="{x_center}" y="{y_base + 18}" text-anchor="middle" font-size="11" fill="#333">{labels[i]}</text>')

    # Axes
    svg.append(f'<line x1="{margin_left}" y1="{margin_top}" x2="{margin_left}" y2="{margin_top + plot_h}" stroke="#333" stroke-width="1.5"/>')
    svg.append(f'<line x1="{margin_left}" y1="{margin_top + plot_h}" x2="{margin_left + plot_w}" y2="{margin_top + plot_h}" stroke="#333" stroke-width="1.5"/>')

    # Y axis label
    svg.append(f'<text x="18" y="{margin_top + plot_h/2}" text-anchor="middle" font-size="12" fill="#333" transform="rotate(-90 18 {margin_top + plot_h/2})">Time (seconds)</text>')

    # Legend
    lx = w - margin_right + 10
    ly = margin_top + 20
    svg.append(f'<rect x="{lx}" y="{ly}" width="14" height="14" fill="{colors_r}" rx="2"/>')
    svg.append(f'<text x="{lx + 20}" y="{ly + 12}" font-size="11" fill="#333">RastQC</text>')
    svg.append(f'<rect x="{lx}" y="{ly + 22}" width="14" height="14" fill="{colors_f}" rx="2"/>')
    svg.append(f'<text x="{lx + 20}" y="{ly + 34}" font-size="11" fill="#333">FastQC</text>')

    svg.append('</svg>')

    out_path = os.path.join(PAPER_DIR, "figures", "fig1_speed_comparison.svg")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        f.write("\n".join(svg))
    print(f"Figure 1 saved to {out_path}")


def generate_memory_figure():
    """Generate SVG bar chart comparing memory usage."""
    path = os.path.join(BENCH_DIR, "real_genomes", "real_genome_benchmarks.tsv")
    rows = load_tsv(path)

    organisms = ["fly", "yeast", "ecoli", "mouse", "human"]
    labels = ["Fly", "Yeast", "E. coli", "Mouse", "Human"]

    rastqc_mem = []
    fastqc_mem = []
    for org in organisms:
        r = [x for x in rows if x["Tool"] == "RastQC" and x["Organism"] == org][0]
        f = [x for x in rows if x["Tool"] == "FastQC" and x["Organism"] == org][0]
        rastqc_mem.append(int(r["RSS_bytes"]) / 1_048_576)
        fastqc_mem.append(int(f["RSS_bytes"]) / 1_048_576)

    max_mem = max(max(fastqc_mem), 700)

    w, h = 700, 400
    margin_left, margin_bottom, margin_top, margin_right = 80, 60, 40, 120
    plot_w = w - margin_left - margin_right
    plot_h = h - margin_top - margin_bottom
    n = len(organisms)
    group_w = plot_w / n
    bar_w = group_w * 0.35

    svg = []
    svg.append(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" font-family="Arial, sans-serif">')
    svg.append(f'<rect width="{w}" height="{h}" fill="white"/>')

    svg.append(f'<text x="{w/2}" y="22" text-anchor="middle" font-size="14" font-weight="bold">Peak Memory (RSS): RastQC vs FastQC</text>')

    y_ticks = [0, 100, 200, 300, 400, 500, 600, 700]
    for t in y_ticks:
        y = margin_top + plot_h - (t / max(y_ticks) * plot_h)
        svg.append(f'<line x1="{margin_left}" y1="{y}" x2="{margin_left + plot_w}" y2="{y}" stroke="#e0e0e0" stroke-width="0.5"/>')
        svg.append(f'<text x="{margin_left - 8}" y="{y + 4}" text-anchor="end" font-size="11" fill="#555">{t}</text>')

    colors_r, colors_f = "#2980b9", "#e74c3c"
    for i in range(n):
        x_center = margin_left + group_w * i + group_w / 2
        x_r = x_center - bar_w - 2
        x_f = x_center + 2
        y_base = margin_top + plot_h

        h_r = (rastqc_mem[i] / max(y_ticks)) * plot_h
        h_f = (fastqc_mem[i] / max(y_ticks)) * plot_h

        svg.append(f'<rect x="{x_r}" y="{y_base - h_r}" width="{bar_w}" height="{h_r}" fill="{colors_r}" rx="2"/>')
        svg.append(f'<rect x="{x_f}" y="{y_base - h_f}" width="{bar_w}" height="{h_f}" fill="{colors_f}" rx="2"/>')

        svg.append(f'<text x="{x_r + bar_w/2}" y="{y_base - h_r - 4}" text-anchor="middle" font-size="9" fill="{colors_r}">{rastqc_mem[i]:.0f}</text>')
        svg.append(f'<text x="{x_f + bar_w/2}" y="{y_base - h_f - 4}" text-anchor="middle" font-size="9" fill="{colors_f}">{fastqc_mem[i]:.0f}</text>')

        ratio = fastqc_mem[i] / rastqc_mem[i] if rastqc_mem[i] > 0 else 0
        svg.append(f'<text x="{x_center}" y="{y_base + 18}" text-anchor="middle" font-size="11" fill="#333">{labels[i]}</text>')
        svg.append(f'<text x="{x_center}" y="{y_base + 32}" text-anchor="middle" font-size="9" fill="#888">{ratio:.0f}x less</text>')

    svg.append(f'<line x1="{margin_left}" y1="{margin_top}" x2="{margin_left}" y2="{margin_top + plot_h}" stroke="#333" stroke-width="1.5"/>')
    svg.append(f'<line x1="{margin_left}" y1="{margin_top + plot_h}" x2="{margin_left + plot_w}" y2="{margin_top + plot_h}" stroke="#333" stroke-width="1.5"/>')

    svg.append(f'<text x="18" y="{margin_top + plot_h/2}" text-anchor="middle" font-size="12" fill="#333" transform="rotate(-90 18 {margin_top + plot_h/2})">Memory (MB)</text>')

    lx = w - margin_right + 10
    ly = margin_top + 20
    svg.append(f'<rect x="{lx}" y="{ly}" width="14" height="14" fill="{colors_r}" rx="2"/>')
    svg.append(f'<text x="{lx + 20}" y="{ly + 12}" font-size="11" fill="#333">RastQC</text>')
    svg.append(f'<rect x="{lx}" y="{ly + 22}" width="14" height="14" fill="{colors_f}" rx="2"/>')
    svg.append(f'<text x="{lx + 20}" y="{ly + 34}" font-size="11" fill="#333">FastQC</text>')

    svg.append('</svg>')

    out_path = os.path.join(PAPER_DIR, "figures", "fig2_memory_comparison.svg")
    with open(out_path, "w") as f:
        f.write("\n".join(svg))
    print(f"Figure 2 saved to {out_path}")


def generate_startup_figure():
    """Generate SVG showing startup time comparison."""
    path = os.path.join(BENCH_DIR, "benchmark_results_v2.tsv")
    rows = load_tsv(path)

    datasets = ["startup", "100K", "1M", "1M_gz"]
    labels = ["Startup\n(1 read)", "100K\nreads", "1M\nreads", "1M\ngzipped"]

    rastqc_times = []
    fastqc_times = []
    for ds in datasets:
        r = [x for x in rows if x["Tool"] == "RastQC" and x["Dataset"] == ds]
        f = [x for x in rows if x["Tool"] == "FastQC" and x["Dataset"] == ds]
        if r and f:
            rastqc_times.append(float(r[0]["Time_s"]))
            fastqc_times.append(float(f[0]["Time_s"]))

    max_time = max(max(fastqc_times), 10)

    w, h = 600, 380
    margin_left, margin_bottom, margin_top, margin_right = 70, 70, 40, 120
    plot_w = w - margin_left - margin_right
    plot_h = h - margin_top - margin_bottom
    n = len(datasets)
    group_w = plot_w / n
    bar_w = group_w * 0.35

    svg = []
    svg.append(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" font-family="Arial, sans-serif">')
    svg.append(f'<rect width="{w}" height="{h}" fill="white"/>')
    svg.append(f'<text x="{w/2}" y="22" text-anchor="middle" font-size="14" font-weight="bold">Synthetic Benchmarks: RastQC vs FastQC</text>')

    y_ticks = [0, 2, 4, 6, 8, 10]
    for t in y_ticks:
        y = margin_top + plot_h - (t / max(y_ticks) * plot_h)
        svg.append(f'<line x1="{margin_left}" y1="{y}" x2="{margin_left + plot_w}" y2="{y}" stroke="#e0e0e0" stroke-width="0.5"/>')
        svg.append(f'<text x="{margin_left - 8}" y="{y + 4}" text-anchor="end" font-size="11" fill="#555">{t}s</text>')

    colors_r, colors_f = "#2980b9", "#e74c3c"
    display_labels = ["Startup", "100K", "1M", "1M gz"]
    for i in range(n):
        x_center = margin_left + group_w * i + group_w / 2
        x_r = x_center - bar_w - 2
        x_f = x_center + 2
        y_base = margin_top + plot_h

        rt = min(rastqc_times[i], max(y_ticks))
        ft = min(fastqc_times[i], max(y_ticks))
        h_r = (rt / max(y_ticks)) * plot_h
        h_f = (ft / max(y_ticks)) * plot_h

        svg.append(f'<rect x="{x_r}" y="{y_base - h_r}" width="{bar_w}" height="{max(h_r, 2)}" fill="{colors_r}" rx="2"/>')
        svg.append(f'<rect x="{x_f}" y="{y_base - h_f}" width="{bar_w}" height="{max(h_f, 2)}" fill="{colors_f}" rx="2"/>')

        rt_label = f"{rastqc_times[i]*1000:.0f}ms" if rastqc_times[i] < 0.1 else f"{rastqc_times[i]:.1f}s"
        ft_label = f"{fastqc_times[i]:.1f}s"
        svg.append(f'<text x="{x_r + bar_w/2}" y="{y_base - h_r - 4}" text-anchor="middle" font-size="8" fill="{colors_r}">{rt_label}</text>')
        svg.append(f'<text x="{x_f + bar_w/2}" y="{y_base - h_f - 4}" text-anchor="middle" font-size="8" fill="{colors_f}">{ft_label}</text>')

        svg.append(f'<text x="{x_center}" y="{y_base + 18}" text-anchor="middle" font-size="11" fill="#333">{display_labels[i]}</text>')

    svg.append(f'<line x1="{margin_left}" y1="{margin_top}" x2="{margin_left}" y2="{margin_top + plot_h}" stroke="#333" stroke-width="1.5"/>')
    svg.append(f'<line x1="{margin_left}" y1="{margin_top + plot_h}" x2="{margin_left + plot_w}" y2="{margin_top + plot_h}" stroke="#333" stroke-width="1.5"/>')
    svg.append(f'<text x="18" y="{margin_top + plot_h/2}" text-anchor="middle" font-size="12" fill="#333" transform="rotate(-90 18 {margin_top + plot_h/2})">Time (seconds)</text>')

    lx = w - margin_right + 10
    ly = margin_top + 20
    svg.append(f'<rect x="{lx}" y="{ly}" width="14" height="14" fill="{colors_r}" rx="2"/>')
    svg.append(f'<text x="{lx + 20}" y="{ly + 12}" font-size="11" fill="#333">RastQC</text>')
    svg.append(f'<rect x="{lx}" y="{ly + 22}" width="14" height="14" fill="{colors_f}" rx="2"/>')
    svg.append(f'<text x="{lx + 20}" y="{ly + 34}" font-size="11" fill="#333">FastQC</text>')

    svg.append('</svg>')

    out_path = os.path.join(PAPER_DIR, "figures", "fig3_synthetic_benchmarks.svg")
    with open(out_path, "w") as f:
        f.write("\n".join(svg))
    print(f"Figure 3 saved to {out_path}")


# ─── Main ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 70)
    print("RastQC Benchmark Analysis")
    print("=" * 70)
    print()

    generate_synthetic_table()
    generate_real_genome_table()
    generate_concordance_table()

    print("=" * 70)
    print("Generating figures...")
    print("=" * 70)
    generate_speedup_figure()
    generate_memory_figure()
    generate_startup_figure()

    print()
    print("Done. Tables printed above; figures saved to paper/figures/")
