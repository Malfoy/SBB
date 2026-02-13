#!/usr/bin/env bash
set -euo pipefail

BENCH_ROOT="${BENCH_ROOT:-/tmp/sbb-bench-big}"
DATA_DIR="$BENCH_ROOT/data"
RS_DIR="$BENCH_ROOT/rs"
THREADS="${THREADS:-$(nproc)}"
TAG="${TAG:-run}"
REF_LEN="${REF_LEN:-60000000}"
READ_LEN="${READ_LEN:-20000}"
READ_COUNT="${READ_COUNT:-300000}"
SEED="${SEED:-1337}"
MAKER_FPR="${MAKER_FPR:-0.00001}"

RS_MAKER="${RS_MAKER:-/home/nadine/Code/SBB/target/release/biobloommaker}"
RS_CAT="${RS_CAT:-/home/nadine/Code/SBB/target/release/biobloomcategorizer}"

mkdir -p "$DATA_DIR" "$RS_DIR"

DATA_META="$DATA_DIR/dataset.meta"
EXPECTED_META=$'ref_len='"$REF_LEN"$'\nread_len='"$READ_LEN"$'\nread_count='"$READ_COUNT"$'\nseed='"$SEED"$'\nreads_format=fasta\n'
CURRENT_META=""
if [[ -f "$DATA_META" ]]; then
  CURRENT_META="$(cat "$DATA_META")"
fi

if [[ ! -s "$DATA_DIR/ref.fa" || ! -s "$DATA_DIR/reads.fa" || "$CURRENT_META"$'\n' != "$EXPECTED_META" ]]; then
  /home/nadine/Code/SBB/bench/gen_dataset.py \
    --out-dir "$DATA_DIR" \
    --ref-len "$REF_LEN" \
    --read-len "$READ_LEN" \
    --read-count "$READ_COUNT" \
    --seed "$SEED"
fi

READS_INPUT="$DATA_DIR/reads.fa"
if "$RS_CAT" --help 2>&1 | rg -qi 'gz'; then
  if [[ ! -s "$DATA_DIR/reads.fa.gz" || "$DATA_DIR/reads.fa" -nt "$DATA_DIR/reads.fa.gz" ]]; then
    gzip -1 -c "$DATA_DIR/reads.fa" > "$DATA_DIR/reads.fa.gz"
  fi
  READS_INPUT="$DATA_DIR/reads.fa.gz"
fi

SUMMARY="$BENCH_ROOT/sbb_${TAG}.tsv"
echo -e "tool\tstep\tseconds\tmax_rss_kb" > "$SUMMARY"

/usr/bin/time -f "sbb\tmaker\t%e\t%M" -o "$SUMMARY" -a \
  "$RS_MAKER" -p rs_ref -o "$RS_DIR" -k 25 -f "$MAKER_FPR" -t "$THREADS" "$DATA_DIR/ref.fa" \
  >/tmp/sbb_maker_${TAG}.log 2>/tmp/sbb_maker_${TAG}.err

/usr/bin/time -f "sbb\tcategorizer\t%e\t%M" -o "$SUMMARY" -a \
  "$RS_CAT" -f "$RS_DIR/rs_ref.bf" -t "$THREADS" "$READS_INPUT" \
  >/tmp/sbb_cat_${TAG}.log 2>/tmp/sbb_cat_${TAG}.err

echo "threads\t$THREADS"
echo "maker_fpr\t$MAKER_FPR"
echo "reads_input\t$READS_INPUT"
wc -l "$DATA_DIR/ref.fa" "$DATA_DIR/reads.fa"
if [[ "$READS_INPUT" == *.gz ]]; then
  ls -lh "$DATA_DIR/reads.fa.gz"
fi
cat "$SUMMARY"
