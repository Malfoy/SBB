#!/usr/bin/env python3
import argparse
import random
from pathlib import Path


def read_fasta(path: Path):
    seqs = []
    parts = []
    with path.open("r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if parts:
                    seqs.append("".join(parts))
                    parts = []
            else:
                parts.append(line.upper())
    if parts:
        seqs.append("".join(parts))
    return [s for s in seqs if s]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", required=True, type=Path)
    ap.add_argument("--out-reads", required=True, type=Path)
    ap.add_argument("--read-len", required=True, type=int)
    ap.add_argument("--read-count", required=True, type=int)
    ap.add_argument("--seed", required=True, type=int)
    ap.add_argument("--out-meta", type=Path)
    args = ap.parse_args()

    if args.read_len < 1 or args.read_count < 1:
        raise SystemExit("read-len and read-count must be >= 1")

    seqs = read_fasta(args.reference)
    if not seqs:
        raise SystemExit(f"no sequence found in {args.reference}")

    spans = []
    for s in seqs:
        span = len(s) - args.read_len + 1
        if span > 0:
            spans.append(span)
        else:
            spans.append(0)
    total_span = sum(spans)
    if total_span <= 0:
        raise SystemExit("reference is shorter than read-len across all contigs")

    rng = random.Random(args.seed)
    cumulative = []
    acc = 0
    for sp in spans:
        acc += sp
        cumulative.append(acc)

    args.out_reads.parent.mkdir(parents=True, exist_ok=True)
    with args.out_reads.open("w") as out:
        for i in range(args.read_count):
            pick = rng.randrange(total_span)
            cidx = 0
            while pick >= cumulative[cidx]:
                cidx += 1
            seq = seqs[cidx]
            start = rng.randrange(len(seq) - args.read_len + 1)
            read = seq[start : start + args.read_len]
            out.write(f">r{i}\n{read}\n")

    if args.out_meta:
        args.out_meta.parent.mkdir(parents=True, exist_ok=True)
        with args.out_meta.open("w") as meta:
            meta.write(f"reference={args.reference}\n")
            meta.write(f"read_len={args.read_len}\n")
            meta.write(f"read_count={args.read_count}\n")
            meta.write(f"seed={args.seed}\n")
            meta.write(f"total_reference_bases={sum(len(s) for s in seqs)}\n")


if __name__ == "__main__":
    main()
