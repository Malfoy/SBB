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


def grouped_means(rows, sweep, metric):
    table = defaultdict(list)
    for r in rows:
        if r["sweep"] != sweep:
            continue
        arm = f'{r["tool"]}:{r["layout"]}'
        x = parse_float(r["param"])
        y = parse_float(r.get(metric, ""))
        if x is None or y is None:
            continue
        table[(arm, x)].append(y)
    out = defaultdict(list)
    for (arm, x), vals in table.items():
        out[arm].append((x, sum(vals) / len(vals)))
    for arm in out:
        out[arm].sort(key=lambda p: p[0])
    return out


def plot_metric(rows, sweep, metric, ylabel, out_png):
    data = grouped_means(rows, sweep, metric)
    if not data:
        return
    plt.figure(figsize=(7.5, 4.8))
    for arm in sorted(data.keys()):
        xs = [p[0] for p in data[arm]]
        ys = [p[1] for p in data[arm]]
        plt.plot(xs, ys, marker="o", label=arm)
    plt.xlabel("Parameter")
    plt.ylabel(ylabel)
    plt.title(f"Single-Reference Sweep ({sweep}): {ylabel}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=180)
    plt.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", required=True, type=Path)
    args = ap.parse_args()

    rows = load_rows(args.tsv)
    sweeps = sorted({r["sweep"] for r in rows})
    out_dir = args.tsv.parent
    stem = args.tsv.stem

    for sweep in sweeps:
        plot_metric(
            rows,
            sweep,
            "query_sec",
            "Query Time (s)",
            out_dir / f"{stem}_{sweep}_query_time.png",
        )
        plot_metric(
            rows,
            sweep,
            "index_bytes",
            "Index Size (bytes)",
            out_dir / f"{stem}_{sweep}_index_size.png",
        )
        plot_metric(
            rows,
            sweep,
            "reads_per_sec",
            "Reads/sec",
            out_dir / f"{stem}_{sweep}_reads_per_sec.png",
        )
        plot_metric(
            rows,
            sweep,
            "build_sec",
            "Build Time (s)",
            out_dir / f"{stem}_{sweep}_build_time.png",
        )

    print(f"plots_dir={out_dir}")


if __name__ == "__main__":
    main()
