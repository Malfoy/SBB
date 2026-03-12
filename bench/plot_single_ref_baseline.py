#!/usr/bin/env python3
import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt


def parse_float(v: str):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def mean_std(values):
    if not values:
        return (None, None)
    m = sum(values) / len(values)
    if len(values) == 1:
        return (m, 0.0)
    var = sum((x - m) ** 2 for x in values) / (len(values) - 1)
    return (m, math.sqrt(var))


def load_rows(tsv_path: Path):
    with tsv_path.open("r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader)


def aggregate(rows, metric):
    grouped = defaultdict(list)
    for r in rows:
        arm = f'{r["tool"]}:{r["layout"]}'
        v = parse_float(r.get(metric, ""))
        if v is not None:
            grouped[arm].append(v)
    out = []
    for arm in sorted(grouped.keys()):
        m, s = mean_std(grouped[arm])
        if m is not None:
            out.append((arm, m, s))
    return out


def bar_with_error(rows, metric, ylabel, out_png):
    stats = aggregate(rows, metric)
    if not stats:
        return
    labels = [x[0] for x in stats]
    means = [x[1] for x in stats]
    stds = [x[2] for x in stats]

    plt.figure(figsize=(8, 4.5))
    plt.bar(labels, means, yerr=stds, capsize=4)
    plt.ylabel(ylabel)
    plt.title(f"Single-Reference Baseline: {ylabel}")
    plt.tight_layout()
    plt.savefig(out_png, dpi=180)
    plt.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", required=True, type=Path)
    args = ap.parse_args()

    rows = load_rows(args.tsv)
    out_dir = args.tsv.parent
    stem = args.tsv.stem

    bar_with_error(rows, "build_sec", "Build Time (s)", out_dir / f"{stem}_build_time.png")
    bar_with_error(rows, "query_sec", "Query Time (s)", out_dir / f"{stem}_query_time.png")
    bar_with_error(rows, "index_bytes", "Index Size (bytes)", out_dir / f"{stem}_index_size.png")
    bar_with_error(rows, "build_rss_kb", "Build Max RSS (KB)", out_dir / f"{stem}_build_rss.png")
    bar_with_error(rows, "query_rss_kb", "Query Max RSS (KB)", out_dir / f"{stem}_query_rss.png")
    bar_with_error(rows, "reads_per_sec", "Reads/sec", out_dir / f"{stem}_reads_per_sec.png")

    print(f"plots_dir={out_dir}")


if __name__ == "__main__":
    main()
