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
        return list(csv.DictReader(fh, delimiter="\t"))


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


def bar_metric(rows, metric, ylabel, out_png):
    stats = aggregate(rows, metric)
    if not stats:
        return
    labels = [x[0] for x in stats]
    means = [x[1] for x in stats]
    stds = [x[2] for x in stats]

    plt.figure(figsize=(8, 4.5))
    plt.bar(labels, means, yerr=stds, capsize=4)
    plt.ylabel(ylabel)
    plt.title(f"Single-Reference Accuracy: {ylabel}")
    plt.ylim(0, 1.05 if metric in {"precision", "recall", "f1", "specificity", "near_recall"} else None)
    plt.tight_layout()
    plt.savefig(out_png, dpi=180)
    plt.close()


def scatter_tradeoff(rows, out_png):
    points = defaultdict(list)
    for r in rows:
        arm = f'{r["tool"]}:{r["layout"]}'
        x = parse_float(r.get("reads_per_sec_total", ""))
        y = parse_float(r.get("f1", ""))
        if x is not None and y is not None:
            points[arm].append((x, y))

    if not points:
        return

    plt.figure(figsize=(7.5, 5))
    for arm in sorted(points.keys()):
        xs = [p[0] for p in points[arm]]
        ys = [p[1] for p in points[arm]]
        plt.scatter(xs, ys, label=arm)
    plt.xlabel("Reads/sec")
    plt.ylabel("F1")
    plt.title("Single-Reference Accuracy/Speed Tradeoff")
    plt.legend()
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

    for metric, label in [
        ("precision", "Precision"),
        ("recall", "Recall"),
        ("f1", "F1"),
        ("specificity", "Specificity"),
        ("near_recall", "Near-Neighbor Recall"),
        ("query_sec_total", "Total Query Time (s)"),
        ("reads_per_sec_total", "Total Reads/sec"),
    ]:
        bar_metric(rows, metric, label, out_dir / f"{stem}_{metric}.png")

    scatter_tradeoff(rows, out_dir / f"{stem}_tradeoff.png")
    print(f"plots_dir={out_dir}")


if __name__ == "__main__":
    main()
