#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt


def parse_float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def load_rows(path: Path):
    with path.open("r", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def mean_by_threads(rows, metric):
    table = defaultdict(list)
    for r in rows:
        arm = f'{r["tool"]}:{r["layout"]}'
        t = parse_float(r["threads"])
        y = parse_float(r.get(metric, ""))
        if t is None or y is None:
            continue
        table[(arm, int(t))].append(y)

    out = defaultdict(list)
    for (arm, t), vals in table.items():
        out[arm].append((t, sum(vals) / len(vals)))
    for arm in out:
        out[arm].sort(key=lambda p: p[0])
    return out


def plot_lines(data, ylabel, title, out_png):
    if not data:
        return
    plt.figure(figsize=(7.5, 4.8))
    for arm in sorted(data.keys()):
        xs = [p[0] for p in data[arm]]
        ys = [p[1] for p in data[arm]]
        plt.plot(xs, ys, marker="o", label=arm)
    plt.xlabel("Threads")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=180)
    plt.close()


def build_efficiency(data):
    out = defaultdict(list)
    for arm, pts in data.items():
        if not pts:
            continue
        base_threads = min(x for x, _ in pts)
        base_time = next(y for x, y in pts if x == base_threads)
        if base_time <= 0:
            continue
        for t, tm in pts:
            if tm <= 0:
                continue
            speedup = base_time / tm
            eff = speedup / (t / base_threads)
            out[arm].append((t, eff))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", required=True, type=Path)
    args = ap.parse_args()

    rows = load_rows(args.tsv)
    out_dir = args.tsv.parent
    stem = args.tsv.stem

    build_sec = mean_by_threads(rows, "build_sec")
    query_sec = mean_by_threads(rows, "query_sec")
    reads_per_sec = mean_by_threads(rows, "reads_per_sec")

    plot_lines(
        build_sec,
        "Build Time (s)",
        "Single-Reference Scaling: Build Time vs Threads",
        out_dir / f"{stem}_build_time.png",
    )
    plot_lines(
        query_sec,
        "Query Time (s)",
        "Single-Reference Scaling: Query Time vs Threads",
        out_dir / f"{stem}_query_time.png",
    )
    plot_lines(
        reads_per_sec,
        "Reads/sec",
        "Single-Reference Scaling: Throughput vs Threads",
        out_dir / f"{stem}_reads_per_sec.png",
    )

    eff = build_efficiency(query_sec)
    plot_lines(
        eff,
        "Parallel Efficiency",
        "Single-Reference Scaling: Query Efficiency vs Threads",
        out_dir / f"{stem}_query_efficiency.png",
    )

    print(f"plots_dir={out_dir}")


if __name__ == "__main__":
    main()
