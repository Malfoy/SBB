#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$ROOT_DIR/bench/common.sh"
bench_init_cpp_env

BENCH_ROOT="${BENCH_ROOT:-/tmp/sbb-single-ref-scaling}"
DATA_DIR="$BENCH_ROOT/data"
RUN_DIR="$BENCH_ROOT/runs"
RESULT_DIR="$BENCH_ROOT/results"

REPEATS="${REPEATS:-3}"
THREAD_LIST="${THREAD_LIST:-1,2,4,8,16,32}"
REF_LEN="${REF_LEN:-60000000}"
READ_LEN="${READ_LEN:-20000}"
READ_COUNT="${READ_COUNT:-300000}"
SEED="${SEED:-1337}"
KMER_SIZE="${KMER_SIZE:-25}"
FPR="${FPR:-0.00001}"
CAT_SCORE="${CAT_SCORE:-1.0}"
BLOCK_WORDS="${BLOCK_WORDS:-8}"
REFERENCE_FASTA="${REFERENCE_FASTA:-}"
PLOT="${PLOT:-1}"

SBB_BIN="${SBB_BIN:-$ROOT_DIR/target/release/sbb}"
BENCHGEN="${BENCHGEN:-$ROOT_DIR/target/release/sbbbenchgen}"
READ_GEN="${READ_GEN:-$ROOT_DIR/bench/generate_reads_from_reference.py}"
PLOT_SCRIPT="${PLOT_SCRIPT:-$ROOT_DIR/bench/plot_single_ref_scaling.py}"

if (( REPEATS < 1 )); then
  echo "REPEATS must be >= 1" >&2
  exit 1
fi

mkdir -p "$DATA_DIR" "$RUN_DIR" "$RESULT_DIR"

for bin in "$CPP_MAKER" "$CPP_CAT"; do
  if [[ ! -x "$bin" ]]; then
    echo "required executable not found: $bin" >&2
    exit 1
  fi
done
if [[ -z "$REFERENCE_FASTA" && -s "$ROOT_DIR/bench/refs/ecoli_k12_mg1655.fa" ]]; then
  REFERENCE_FASTA="$ROOT_DIR/bench/refs/ecoli_k12_mg1655.fa"
fi

echo "[scaling] building binaries"
cargo build --release --bin sbb --bin sbbbenchgen >/tmp/sbb_single_ref_scaling_build.log

for bin in "$SBB_BIN" "$BENCHGEN"; do
  if [[ ! -x "$bin" ]]; then
    echo "required executable not found after build: $bin" >&2
    exit 1
  fi
done

DATA_META="$DATA_DIR/dataset.meta"
EXPECTED_META=""
if [[ -n "$REFERENCE_FASTA" ]]; then
  EXPECTED_META=$'source=real_reference\nreference_fasta='"$(realpath "$REFERENCE_FASTA")"$'\nread_len='"$READ_LEN"$'\nread_count='"$READ_COUNT"$'\nseed='"$SEED"$'\nreads_format=fasta\n'
else
  EXPECTED_META=$'source=synthetic\nref_len='"$REF_LEN"$'\nread_len='"$READ_LEN"$'\nread_count='"$READ_COUNT"$'\nseed='"$SEED"$'\nreads_format=fasta\n'
fi
CURRENT_META=""
if [[ -f "$DATA_META" ]]; then
  CURRENT_META="$(cat "$DATA_META")"
fi

if [[ ! -s "$DATA_DIR/ref.fa" || ! -s "$DATA_DIR/reads.fa" || "$CURRENT_META"$'\n' != "$EXPECTED_META" ]]; then
  if [[ -n "$REFERENCE_FASTA" ]]; then
    echo "[scaling] preparing dataset from reference: $REFERENCE_FASTA"
    cp "$REFERENCE_FASTA" "$DATA_DIR/ref.fa"
    python3 "$READ_GEN" \
      --reference "$DATA_DIR/ref.fa" \
      --out-reads "$DATA_DIR/reads.fa" \
      --read-len "$READ_LEN" \
      --read-count "$READ_COUNT" \
      --seed "$SEED" \
      --out-meta "$DATA_DIR/reads.meta"
  else
    echo "[scaling] generating synthetic dataset"
    "$BENCHGEN" \
      --out_dir "$DATA_DIR" \
      --ref_len "$REF_LEN" \
      --read_len "$READ_LEN" \
      --read_count "$READ_COUNT" \
      --seed "$SEED" \
      --threads 1 \
      >/tmp/sbb_single_ref_scaling_benchgen.log
  fi
  printf "%s" "$EXPECTED_META" > "$DATA_META"
fi

READS_INPUT="$DATA_DIR/reads.fa"
if "$CPP_CAT" --help 2>&1 | rg -qi 'gz'; then
  if [[ ! -s "$DATA_DIR/reads.fa.gz" || "$DATA_DIR/reads.fa" -nt "$DATA_DIR/reads.fa.gz" ]]; then
    gzip -1 -c "$DATA_DIR/reads.fa" > "$DATA_DIR/reads.fa.gz"
  fi
  READS_INPUT="$DATA_DIR/reads.fa.gz"
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
RESULTS="$RESULT_DIR/single_ref_scaling_${STAMP}.tsv"
echo -e "benchmark\ttool\tlayout\trepeat\tthreads\tkmer_size\tfpr\tscore\tread_count\tbuild_sec\tbuild_rss_kb\tquery_sec\tquery_rss_kb\tindex_bytes\tinsert_kmer_ops\tinsert_mean_ns_per_kmer\tquery_kmer_queries\tquery_mean_ns_per_kmer\treads_per_sec\tquery_kmers_per_sec\tfilter_path\treads_path" > "$RESULTS"

to_num_or_na() {
  local v="${1:-}"
  if [[ -z "$v" ]]; then
    echo "na"
  elif [[ "$v" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    echo "$v"
  else
    echo "na"
  fi
}

tsv_value() {
  local file="$1"
  local key="$2"
  awk -F'\t' -v k="$key" '$1==k{print $2}' "$file" | tail -n1
}

run_case() {
  local tool="$1"
  local layout="$2"
  local rep="$3"
  local threads="$4"
  local run_tag="${tool}_${layout}_t${threads}_r${rep}"
  local out_dir="$RUN_DIR/$run_tag"
  mkdir -p "$out_dir"

  local maker_log="$out_dir/maker.out"
  local query_log="$out_dir/query.out"
  local maker_time="$out_dir/maker.time"
  local query_time="$out_dir/query.time"

  local maker_cmd=()
  local query_cmd=()
  local filter_path=""

  case "$tool/$layout" in
    cpp/classic)
      maker_cmd=("$CPP_MAKER" -p ref -o "$out_dir" -k "$KMER_SIZE" -f "$FPR" -t "$threads" "$DATA_DIR/ref.fa")
      filter_path="$out_dir/ref.bf"
      query_cmd=("$CPP_CAT" -f "$filter_path" -t "$threads" -s "$CAT_SCORE" "$READS_INPUT")
      ;;
    rust/classic)
      maker_cmd=("$SBB_BIN" maker -p ref -o "$out_dir" -k "$KMER_SIZE" -f "$FPR" -t "$threads" --classic "$DATA_DIR/ref.fa")
      filter_path="$out_dir/ref.bf.zst"
      query_cmd=("$SBB_BIN" categorizer -f "$filter_path" -t "$threads" -s "$CAT_SCORE" "$READS_INPUT")
      ;;
    rust/blocked)
      maker_cmd=("$SBB_BIN" maker -p ref -o "$out_dir" -k "$KMER_SIZE" -f "$FPR" -t "$threads" --block_words "$BLOCK_WORDS" "$DATA_DIR/ref.fa")
      filter_path="$out_dir/ref.bf.zst"
      query_cmd=("$SBB_BIN" categorizer -f "$filter_path" -t "$threads" -s "$CAT_SCORE" "$READS_INPUT")
      ;;
    *)
      echo "unknown case: $tool/$layout" >&2
      exit 1
      ;;
  esac

  /usr/bin/time -f "%e\t%M" -o "$maker_time" "${maker_cmd[@]}" >"$maker_log" 2>"$out_dir/maker.err"
  /usr/bin/time -f "%e\t%M" -o "$query_time" "${query_cmd[@]}" >"$query_log" 2>"$out_dir/query.err"

  local build_sec build_rss query_sec query_rss index_bytes
  build_sec="$(cut -f1 "$maker_time")"
  build_rss="$(cut -f2 "$maker_time")"
  query_sec="$(cut -f1 "$query_time")"
  query_rss="$(cut -f2 "$query_time")"
  index_bytes="$(stat -c%s "$filter_path" 2>/dev/null || echo na)"

  local insert_ops insert_ns query_ops query_ns reads_per_sec query_kmers_per_sec
  insert_ops="$(to_num_or_na "$(tsv_value "$maker_log" "insert_kmer_ops")")"
  insert_ns="$(to_num_or_na "$(tsv_value "$maker_log" "insert_mean_ns_per_kmer")")"
  query_ops="$(to_num_or_na "$(tsv_value "$query_log" "query_kmer_queries")")"
  query_ns="$(to_num_or_na "$(tsv_value "$query_log" "query_mean_ns_per_kmer")")"
  reads_per_sec="$(awk -v n="$READ_COUNT" -v t="$query_sec" 'BEGIN{if(t>0){printf "%.6f", n/t}else{print "na"}}')"
  if [[ "$query_ops" == "na" ]]; then
    query_kmers_per_sec="na"
  else
    query_kmers_per_sec="$(awk -v n="$query_ops" -v t="$query_sec" 'BEGIN{if(t>0){printf "%.6f", n/t}else{print "na"}}')"
  fi

  echo -e "single_ref_scaling\t$tool\t$layout\t$rep\t$threads\t$KMER_SIZE\t$FPR\t$CAT_SCORE\t$READ_COUNT\t$build_sec\t$build_rss\t$query_sec\t$query_rss\t$index_bytes\t$insert_ops\t$insert_ns\t$query_ops\t$query_ns\t$reads_per_sec\t$query_kmers_per_sec\t$filter_path\t$READS_INPUT" >> "$RESULTS"
}

IFS=',' read -r -a THREADS_ARR <<< "$THREAD_LIST"
for threads in "${THREADS_ARR[@]}"; do
  threads="${threads// /}"
  [[ -n "$threads" ]] || continue
  if (( threads < 1 )); then
    echo "invalid thread count: $threads" >&2
    exit 1
  fi
  for rep in $(seq 1 "$REPEATS"); do
    run_case cpp classic "$rep" "$threads"
    run_case rust classic "$rep" "$threads"
    run_case rust blocked "$rep" "$threads"
  done
done

echo "results_tsv=$RESULTS"

if [[ "$PLOT" == "1" ]]; then
  if command -v python3 >/dev/null 2>&1; then
    python3 "$PLOT_SCRIPT" --tsv "$RESULTS"
  else
    echo "[scaling] python3 not found; skipping plots" >&2
  fi
fi
