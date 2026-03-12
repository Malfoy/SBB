#!/usr/bin/env python3

import argparse
import csv
import re
import statistics
import subprocess
from datetime import datetime
from pathlib import Path

try:
    import matplotlib.pyplot as plt
except ImportError as exc:
    raise SystemExit("matplotlib is required: pip install matplotlib") from exc


BENCHMARK_NAME = "bench_s"
S_VALUES = [17, 21, 25, 28]
K_VALUE = 31
M_VALUE = 13
N_HASHES = 3
SIZE_EXPONENT = 35
BLOCK_SIZE_EXPONENT = 13
THREADS = 32
REPEATS = 1
BUILD_FIRST = True
USE_INDEXED_FILE_FLAG = False
EXTRA_ARGS: list[str] = []

METRIC_PATTERNS = {
    "index_wall_s": r"total indexing time \(s\) : ([0-9eE+.\-]+)",
    "index_cpu_s": r"index cpu time \(s\) : ([0-9eE+.\-]+)",
    "query_wall_s": r"total query time \(s\) : ([0-9eE+.\-]+)",
    "query_cpu_s": r"query cpu time \(s\) : ([0-9eE+.\-]+)",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("index_file")
    parser.add_argument("query_file")
    return parser.parse_args()


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def results_dir() -> Path:
    path = Path(__file__).resolve().parent / "results"
    path.mkdir(parents=True, exist_ok=True)
    return path


def build_release(root: Path) -> None:
    if BUILD_FIRST:
        subprocess.run(["cargo", "build", "-r"], cwd=root, check=True)


def run_bloom(root: Path, command: list[str]) -> dict[str, float]:
    completed = subprocess.run(
        command,
        cwd=root,
        check=True,
        text=True,
        capture_output=True,
    )
    metrics: dict[str, float] = {}
    for key, pattern in METRIC_PATTERNS.items():
        match = re.search(pattern, completed.stdout)
        if match is None:
            raise RuntimeError(
                f"Missing metric {key} in output.\nSTDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}"
            )
        metrics[key] = float(match.group(1))
    return metrics


def build_command(index_file: str, query_file: str, s_value: int) -> list[str]:
    command = [
        "./target/release/bloomybloom",
        "--query-file",
        str(Path(query_file).expanduser()),
        "--threads",
        str(THREADS),
        "-k",
        str(K_VALUE),
        "-m",
        str(M_VALUE),
        "-s",
        str(s_value),
        "--size",
        str(SIZE_EXPONENT),
        "--block-size",
        str(BLOCK_SIZE_EXPONENT),
        "--n-hashes",
        str(N_HASHES),
    ]
    if USE_INDEXED_FILE_FLAG:
        command.extend(["--indexed-file", str(Path(index_file).expanduser())])
    else:
        command.append(str(Path(index_file).expanduser()))
    command.extend(EXTRA_ARGS)
    return command


def aggregate(metrics_list: list[dict[str, float]]) -> dict[str, float]:
    return {
        key: statistics.fmean(run[key] for run in metrics_list)
        for key in METRIC_PATTERNS
    }


def write_tsv(rows: list[dict[str, object]], output_path: Path) -> None:
    fieldnames = [
        "benchmark",
        "index_file",
        "query_file",
        "k",
        "m",
        "s",
        "n_hashes",
        "size_exponent",
        "block_size_exponent",
        "threads",
        "repeats",
        "index_wall_s",
        "index_cpu_s",
        "query_wall_s",
        "query_cpu_s",
    ]
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def plot_rows(rows: list[dict[str, object]], output_path: Path) -> None:
    x_values = [int(row["s"]) for row in rows]
    index_wall = [float(row["index_wall_s"]) for row in rows]
    index_cpu = [float(row["index_cpu_s"]) for row in rows]
    query_wall = [float(row["query_wall_s"]) for row in rows]
    query_cpu = [float(row["query_cpu_s"]) for row in rows]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharex=True)

    axes[0].plot(x_values, index_wall, marker="o", label="index wall")
    axes[0].plot(x_values, index_cpu, marker="o", label="index cpu")
    axes[0].set_title("Index Phase")
    axes[0].set_xlabel("s")
    axes[0].set_ylabel("seconds")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()

    axes[1].plot(x_values, query_wall, marker="o", label="query wall")
    axes[1].plot(x_values, query_cpu, marker="o", label="query cpu")
    axes[1].set_title("Query Phase")
    axes[1].set_xlabel("s")
    axes[1].set_ylabel("seconds")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()

    fig.suptitle("bloomybloom benchmark: s sweep")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    root = repo_root()
    build_release(root)

    rows: list[dict[str, object]] = []
    for s_value in S_VALUES:
        metrics_list = [
            run_bloom(root, build_command(args.index_file, args.query_file, s_value))
            for _ in range(REPEATS)
        ]
        metrics = aggregate(metrics_list)
        rows.append({
            "benchmark": BENCHMARK_NAME,
            "index_file": str(Path(args.index_file).expanduser()),
            "query_file": str(Path(args.query_file).expanduser()),
            "k": K_VALUE,
            "m": M_VALUE,
            "s": s_value,
            "n_hashes": N_HASHES,
            "size_exponent": SIZE_EXPONENT,
            "block_size_exponent": BLOCK_SIZE_EXPONENT,
            "threads": THREADS,
            "repeats": REPEATS,
            **metrics,
        })

    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    out_dir = results_dir()
    tsv_path = out_dir / f"{BENCHMARK_NAME}-{timestamp}.tsv"
    png_path = out_dir / f"{BENCHMARK_NAME}-{timestamp}.png"
    write_tsv(rows, tsv_path)
    plot_rows(rows, png_path)
    print(tsv_path)
    print(png_path)


if __name__ == "__main__":
    main()
