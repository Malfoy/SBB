#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BENCH_ROOT="${BENCH_ROOT:-/tmp/sbb-bench-rust}"
DATA_DIR="$BENCH_ROOT/data"
RESULT_DIR="$BENCH_ROOT/results"
RUN_DIR="$BENCH_ROOT/runs"

THREADS="${THREADS:-$(nproc)}"
REF_LEN="${REF_LEN:-100000000}"
READ_LEN="${READ_LEN:-20000}"
READ_COUNT="${READ_COUNT:-150000}"
SEED="${SEED:-1337}"
KMER_SIZE="${KMER_SIZE:-25}"
FPR="${FPR:-0.000001}"
CAT_SCORE="${CAT_SCORE:-1.0}"
BLOCK_WORDS="${BLOCK_WORDS:-1}"
HASH_NUM="${HASH_NUM:-4}"
BITS_PER_ELEMENT_LIST="${BITS_PER_ELEMENT_LIST:-8,12,16,24,32}"

BENCHGEN="$ROOT_DIR/target/release/sbbbenchgen"
RS_MAKER="$ROOT_DIR/target/release/biobloommaker"
RS_CAT="$ROOT_DIR/target/release/biobloomcategorizer"

mkdir -p "$DATA_DIR" "$RESULT_DIR" "$RUN_DIR"

next_pow2() {
  local v="$1"
  local p=1
  while (( p < v )); do
    p=$((p << 1))
  done
  echo "$p"
}

echo "[bench] building release binaries"
cargo build --release --bin sbbbenchgen --bin biobloommaker --bin biobloomcategorizer >/tmp/sbb_bench_build.log

META_PATH="$DATA_DIR/dataset.meta"
EXPECTED_META=$'ref_len='"$REF_LEN"$'\nread_len='"$READ_LEN"$'\nread_count='"$READ_COUNT"$'\nseed='"$SEED"$'\nreads_format=fasta\n'
CURRENT_META=""
if [[ -f "$META_PATH" ]]; then
  CURRENT_META="$(cat "$META_PATH")"
fi

if [[ ! -s "$DATA_DIR/ref.fa" || ! -s "$DATA_DIR/reads.fa" || "$CURRENT_META"$'\n' != "$EXPECTED_META" ]]; then
  echo "[bench] generating dataset with sbbbenchgen"
  "$BENCHGEN" \
    --out_dir "$DATA_DIR" \
    --ref_len "$REF_LEN" \
    --read_len "$READ_LEN" \
    --read_count "$READ_COUNT" \
    --seed "$SEED" \
    --threads "$THREADS" \
    >/tmp/sbb_benchgen.log
fi

READS_INPUT="$DATA_DIR/reads.fa.gz"
if [[ ! -s "$READS_INPUT" || "$DATA_DIR/reads.fa" -nt "$READS_INPUT" ]]; then
  echo "[bench] compressing reads with gzip -1"
  gzip -1 -c "$DATA_DIR/reads.fa" > "$READS_INPUT"
fi

NUM_ELE=$((REF_LEN - KMER_SIZE + 1))
if (( NUM_ELE <= 0 )); then
  echo "invalid dimensions: REF_LEN ($REF_LEN) must be >= KMER_SIZE ($KMER_SIZE)" >&2
  exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
RESULTS="$RESULT_DIR/bloom_size_${STAMP}.tsv"
echo -e "suite\tparam\tlayout\tbit_len\thash_num\tmaker_sec\tmaker_rss_kb\tinsert_kmer_ops\tinsert_mean_ns_per_kmer\tcategorizer_sec\tcategorizer_rss_kb\tquery_kmer_queries\tquery_mean_ns_per_kmer" > "$RESULTS"

run_case() {
  local bpe="$1"
  local layout="$2"
  local bit_len="$3"
  local run_tag="bpe${bpe}_${layout}_b${bit_len}_h${HASH_NUM}"
  local out_dir="$RUN_DIR/$run_tag"
  local maker_log="$out_dir/maker.out"
  local cat_log="$out_dir/cat.out"
  local maker_time="$out_dir/maker.time"
  local cat_time="$out_dir/cat.time"
  mkdir -p "$out_dir"

  local layout_args=()
  if [[ "$layout" == "blocked" ]]; then
    layout_args+=(--block_words "$BLOCK_WORDS")
  else
    layout_args+=(--classic)
  fi

  /usr/bin/time -f "%e\t%M" -o "$maker_time" \
    "$RS_MAKER" \
      -p "ref_${run_tag}" \
      -o "$out_dir" \
      -k "$KMER_SIZE" \
      -n "$NUM_ELE" \
      -f "$FPR" \
      -t "$THREADS" \
      --hash_num "$HASH_NUM" \
      --bit_len "$bit_len" \
      "${layout_args[@]}" \
      "$DATA_DIR/ref.fa" \
      >"$maker_log" 2>"$out_dir/maker.err"

  /usr/bin/time -f "%e\t%M" -o "$cat_time" \
    "$RS_CAT" \
      -f "$out_dir/ref_${run_tag}.bf" \
      -t "$THREADS" \
      -s "$CAT_SCORE" \
      "$READS_INPUT" \
      >"$cat_log" 2>"$out_dir/cat.err"

  local maker_sec maker_rss cat_sec cat_rss insert_ops insert_ns query_ops query_ns
  maker_sec="$(cut -f1 "$maker_time")"
  maker_rss="$(cut -f2 "$maker_time")"
  cat_sec="$(cut -f1 "$cat_time")"
  cat_rss="$(cut -f2 "$cat_time")"
  insert_ops="$(awk -F'\t' '$1=="insert_kmer_ops"{print $2}' "$maker_log" | tail -n1)"
  insert_ns="$(awk -F'\t' '$1=="insert_mean_ns_per_kmer"{print $2}' "$maker_log" | tail -n1)"
  query_ops="$(awk -F'\t' '$1=="query_kmer_queries"{print $2}' "$cat_log" | tail -n1)"
  query_ns="$(awk -F'\t' '$1=="query_mean_ns_per_kmer"{print $2}' "$cat_log" | tail -n1)"

  echo -e "bloom_size\t$bpe\t$layout\t$bit_len\t$HASH_NUM\t$maker_sec\t$maker_rss\t$insert_ops\t$insert_ns\t$cat_sec\t$cat_rss\t$query_ops\t$query_ns" >> "$RESULTS"
}

IFS=',' read -r -a BPE_VALUES <<< "$BITS_PER_ELEMENT_LIST"
for bpe in "${BPE_VALUES[@]}"; do
  bpe="${bpe// /}"
  [[ -n "$bpe" ]] || continue
  bit_len="$(next_pow2 $((NUM_ELE * bpe)))"
  run_case "$bpe" classic "$bit_len"
  run_case "$bpe" blocked "$bit_len"
done

echo
echo "results_file=$RESULTS"
column -ts $'\t' "$RESULTS"
