#!/usr/bin/env python3
import argparse
import os
import random
import sys
from pathlib import Path


def write_fasta(path: Path, seq: bytes, wrap: int = 80) -> None:
    with path.open("wb", buffering=16 * 1024 * 1024) as f:
        f.write(b">ref1\n")
        for i in range(0, len(seq), wrap):
            f.write(seq[i : i + wrap])
            f.write(b"\n")


def main() -> int:
    p = argparse.ArgumentParser(description="Generate deterministic benchmark dataset")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--ref-len", type=int, default=60_000_000)
    p.add_argument("--read-len", type=int, default=20_000)
    p.add_argument("--read-count", type=int, default=300_000)
    p.add_argument("--seed", type=int, default=1337)
    args = p.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    ref_path = out / "ref.fa"
    reads_path = out / "reads.fa"
    meta_path = out / "dataset.meta"

    rng = random.Random(args.seed)

    # Fast deterministic DNA generation: random bytes -> translated A/C/G/T.
    trans = bytes((b"ACGT"[i & 3] for i in range(256)))
    ref = rng.randbytes(args.ref_len).translate(trans)

    print(f"[gen] writing reference: {ref_path}", file=sys.stderr)
    write_fasta(ref_path, ref)

    print(f"[gen] writing reads (FASTA): {reads_path}", file=sys.stderr)
    progress_every = max(1, args.read_count // 20)

    with reads_path.open("wb", buffering=16 * 1024 * 1024) as f:
        for i in range(args.read_count):
            start = rng.randrange(0, args.ref_len - args.read_len + 1)
            seq = ref[start : start + args.read_len]
            f.write(b">r")
            f.write(str(i).encode("ascii"))
            f.write(b"\n")
            f.write(seq)
            f.write(b"\n")
            if i and i % progress_every == 0:
                pct = (100.0 * i) / args.read_count
                print(f"[gen] reads: {i}/{args.read_count} ({pct:.1f}%)", file=sys.stderr)

    print("[gen] done", file=sys.stderr)
    meta_path.write_text(
        "\n".join(
            [
                f"ref_len={args.ref_len}",
                f"read_len={args.read_len}",
                f"read_count={args.read_count}",
                f"seed={args.seed}",
                "reads_format=fasta",
            ]
        )
        + "\n",
        encoding="ascii",
    )
    print(f"ref_path={ref_path}")
    print(f"reads_path={reads_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
