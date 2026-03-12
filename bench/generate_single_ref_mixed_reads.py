#!/usr/bin/env python3
import argparse
import random
from pathlib import Path


BASES = ("A", "C", "G", "T")


def read_fasta_records(path: Path):
    records = []
    cur_id = None
    seq_parts = []
    with path.open("r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    records.append((cur_id, "".join(seq_parts)))
                cur_id = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
    if cur_id is not None:
        records.append((cur_id, "".join(seq_parts)))
    return records


def mutate_seq(seq: str, rate: float, rng: random.Random):
    out = []
    changed = 0
    for b in seq:
        if rng.random() < rate:
            choices = [x for x in BASES if x != b]
            out.append(rng.choice(choices))
            changed += 1
        else:
            out.append(b)
    if changed == 0 and seq:
        idx = rng.randrange(len(seq))
        b = out[idx]
        out[idx] = rng.choice([x for x in BASES if x != b])
    return "".join(out)


def random_seq(length: int, rng: random.Random):
    return "".join(rng.choice(BASES) for _ in range(length))


def write_fasta(path: Path, records):
    with path.open("w") as fh:
        for rid, seq in records:
            fh.write(f">{rid}\n{seq}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-reads", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--pos-count", required=True, type=int)
    ap.add_argument("--near-count", required=True, type=int)
    ap.add_argument("--neg-count", required=True, type=int)
    ap.add_argument("--near-mutation-rate", required=True, type=float)
    ap.add_argument("--seed", required=True, type=int)
    args = ap.parse_args()

    if args.pos_count < 1 or args.near_count < 0 or args.neg_count < 0:
        raise SystemExit("counts must be non-negative and pos-count >= 1")
    if not (0.0 < args.near_mutation_rate < 1.0):
        raise SystemExit("near-mutation-rate must be in (0,1)")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)

    src = read_fasta_records(args.source_reads)
    if not src:
        raise SystemExit(f"no source reads in {args.source_reads}")
    if len(src) < args.pos_count:
        raise SystemExit(
            f"source has {len(src)} reads but pos-count is {args.pos_count}; increase benchgen read count"
        )

    pos_src = src[: args.pos_count]
    read_len = len(pos_src[0][1])
    if read_len == 0:
        raise SystemExit("source reads have zero length")

    positives = [(f"pos_{i}", seq) for i, (_, seq) in enumerate(pos_src)]
    near = []
    for i in range(args.near_count):
        seq = pos_src[i % len(pos_src)][1]
        near.append((f"near_{i}", mutate_seq(seq, args.near_mutation_rate, rng)))
    negatives = [(f"neg_{i}", random_seq(read_len, rng)) for i in range(args.neg_count)]

    mixed = positives + near + negatives
    rng.shuffle(mixed)

    write_fasta(args.out_dir / "reads_pos.fa", positives)
    write_fasta(args.out_dir / "reads_near.fa", near)
    write_fasta(args.out_dir / "reads_neg.fa", negatives)
    write_fasta(args.out_dir / "reads_mixed.fa", mixed)

    with (args.out_dir / "truth.tsv").open("w") as fh:
        fh.write("read_id\tclass\n")
        for rid, _ in positives:
            fh.write(f"{rid}\tpos\n")
        for rid, _ in near:
            fh.write(f"{rid}\tnear\n")
        for rid, _ in negatives:
            fh.write(f"{rid}\tneg\n")

    with (args.out_dir / "dataset.meta").open("w") as fh:
        fh.write(f"pos_count={len(positives)}\n")
        fh.write(f"near_count={len(near)}\n")
        fh.write(f"neg_count={len(negatives)}\n")
        fh.write(f"read_len={read_len}\n")
        fh.write(f"near_mutation_rate={args.near_mutation_rate}\n")
        fh.write(f"seed={args.seed}\n")


if __name__ == "__main__":
    main()
