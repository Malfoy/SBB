#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$ROOT_DIR/bench/common.sh"
bench_init_cpp_env

BENCH_ROOT="${BENCH_ROOT:-/tmp/sbb-single-ref-accuracy}"
DATA_DIR="$BENCH_ROOT/data"
RUN_DIR="$BENCH_ROOT/runs"
RESULT_DIR="$BENCH_ROOT/results"

REPEATS="${REPEATS:-3}"
THREADS="${THREADS:-$(nproc)}"
REF_LEN="${REF_LEN:-60000000}"
READ_LEN="${READ_LEN:-20000}"
SEED="${SEED:-1337}"
KMER_SIZE="${KMER_SIZE:-25}"
FPR="${FPR:-0.00001}"
CAT_SCORE="${CAT_SCORE:-1.0}"
BLOCK_WORDS="${BLOCK_WORDS:-8}"
REFERENCE_FASTA="${REFERENCE_FASTA:-}"
POS_COUNT="${POS_COUNT:-100000}"
NEAR_COUNT="${NEAR_COUNT:-100000}"
NEG_COUNT="${NEG_COUNT:-100000}"
NEAR_MUTATION_RATE="${NEAR_MUTATION_RATE:-0.03}"
PLOT="${PLOT:-1}"

SBB_BIN="${SBB_BIN:-$ROOT_DIR/target/release/sbb}"
BENCHGEN="${BENCHGEN:-$ROOT_DIR/target/release/sbbbenchgen}"
MIX_GEN="${MIX_GEN:-$ROOT_DIR/bench/generate_single_ref_mixed_reads.py}"
READ_GEN="${READ_GEN:-$ROOT_DIR/bench/generate_reads_from_reference.py}"
PLOT_SCRIPT="${PLOT_SCRIPT:-$ROOT_DIR/bench/plot_single_ref_accuracy.py}"

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
if [[ ! -f "$MIX_GEN" ]]; then
  echo "required generator script not found: $MIX_GEN" >&2
  exit 1
fi
if [[ -z "$REFERENCE_FASTA" && -s "$ROOT_DIR/bench/refs/ecoli_k12_mg1655.fa" ]]; then
  REFERENCE_FASTA="$ROOT_DIR/bench/refs/ecoli_k12_mg1655.fa"
fi

echo "[accuracy] building binaries"
cargo build --release --bin sbb --bin sbbbenchgen >/tmp/sbb_single_ref_accuracy_build.log

for bin in "$SBB_BIN" "$BENCHGEN"; do
  if [[ ! -x "$bin" ]]; then
    echo "required executable not found after build: $bin" >&2
    exit 1
  fi
done

BASE_READ_COUNT="$POS_COUNT"
if (( NEAR_COUNT > BASE_READ_COUNT )); then
  BASE_READ_COUNT="$NEAR_COUNT"
fi

DATA_META="$DATA_DIR/source.meta"
EXPECTED_META=""
if [[ -n "$REFERENCE_FASTA" ]]; then
  EXPECTED_META=$'source=real_reference\nreference_fasta='"$(realpath "$REFERENCE_FASTA")"$'\nread_len='"$READ_LEN"$'\nread_count='"$BASE_READ_COUNT"$'\nseed='"$SEED"$'\nreads_format=fasta\n'
else
  EXPECTED_META=$'source=synthetic\nref_len='"$REF_LEN"$'\nread_len='"$READ_LEN"$'\nread_count='"$BASE_READ_COUNT"$'\nseed='"$SEED"$'\nreads_format=fasta\n'
fi
CURRENT_META=""
if [[ -f "$DATA_META" ]]; then
  CURRENT_META="$(cat "$DATA_META")"
fi

if [[ ! -s "$DATA_DIR/ref.fa" || ! -s "$DATA_DIR/source_reads.fa" || "$CURRENT_META"$'\n' != "$EXPECTED_META" ]]; then
  if [[ -n "$REFERENCE_FASTA" ]]; then
    echo "[accuracy] preparing source reads from reference: $REFERENCE_FASTA"
    cp "$REFERENCE_FASTA" "$DATA_DIR/ref.fa"
    python3 "$READ_GEN" \
      --reference "$DATA_DIR/ref.fa" \
      --out-reads "$DATA_DIR/source_reads.fa" \
      --read-len "$READ_LEN" \
      --read-count "$BASE_READ_COUNT" \
      --seed "$SEED" \
      --out-meta "$DATA_DIR/source_reads.meta"
  else
    echo "[accuracy] generating synthetic source reference/reads"
    "$BENCHGEN" \
      --out_dir "$DATA_DIR/source_tmp" \
      --ref_len "$REF_LEN" \
      --read_len "$READ_LEN" \
      --read_count "$BASE_READ_COUNT" \
      --seed "$SEED" \
      --threads "$THREADS" \
      >/tmp/sbb_single_ref_accuracy_benchgen.log
    mv "$DATA_DIR/source_tmp/ref.fa" "$DATA_DIR/ref.fa"
    mv "$DATA_DIR/source_tmp/reads.fa" "$DATA_DIR/source_reads.fa"
    rm -rf "$DATA_DIR/source_tmp"
  fi
  printf "%s" "$EXPECTED_META" > "$DATA_META"
fi

echo "[accuracy] generating labeled mixed read sets"
python3 "$MIX_GEN" \
  --source-reads "$DATA_DIR/source_reads.fa" \
  --out-dir "$DATA_DIR" \
  --pos-count "$POS_COUNT" \
  --near-count "$NEAR_COUNT" \
  --neg-count "$NEG_COUNT" \
  --near-mutation-rate "$NEAR_MUTATION_RATE" \
  --seed "$SEED"

POS_READS="$DATA_DIR/reads_pos.fa"
NEAR_READS="$DATA_DIR/reads_near.fa"
NEG_READS="$DATA_DIR/reads_neg.fa"
TOTAL_QUERY_READS=$((POS_COUNT + NEAR_COUNT + NEG_COUNT))

STAMP="$(date +%Y%m%d_%H%M%S)"
RESULTS="$RESULT_DIR/single_ref_accuracy_${STAMP}.tsv"
echo -e "benchmark\ttool\tlayout\trepeat\tthreads\tkmer_size\tfpr\tscore\tpos_total\tnear_total\tneg_total\tpos_hits\tnear_hits\tneg_hits\ttp\tfp\ttn\tfn\tprecision\trecall\tf1\tspecificity\tfalse_positive_rate\tnear_recall\tbuild_sec\tbuild_rss_kb\tquery_sec_total\tquery_rss_kb_max\treads_per_sec_total\tindex_bytes\tfilter_path" > "$RESULTS"

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

extract_positive_hits() {
  local log="$1"
  local v=""

  v="$(awk -F'\t' '$1=="filter_reads"{sum+=$3;seen=1} END{if(seen) print sum}' "$log")"
  if [[ -n "$v" ]]; then
    echo "$v"
    return
  fi
  v="$(awk -F'\t' '$1=="recruited_reads"{print $2}' "$log" | tail -n1)"
  if [[ -n "$v" ]]; then
    echo "$v"
    return
  fi
  v="$(grep -Eio '(recruited|matched|hits?)[^0-9]{0,20}[0-9]+' "$log" | grep -Eo '[0-9]+' | tail -n1 || true)"
  if [[ -n "$v" ]]; then
    echo "$v"
    return
  fi
  echo "na"
}

query_one_set() {
  local tool="$1"
  local filter_path="$2"
  local reads="$3"
  local out_prefix="$4"
  local log="$5"
  local tm="$6"

  local cmd=()
  if [[ "$tool" == "cpp" ]]; then
    cmd=("$CPP_CAT" -f "$filter_path" -t "$THREADS" -s "$CAT_SCORE" "$reads")
  else
    cmd=("$SBB_BIN" categorizer -f "$filter_path" -t "$THREADS" -s "$CAT_SCORE" "$reads")
  fi

  /usr/bin/time -f "%e\t%M" -o "$tm" "${cmd[@]}" >"$log" 2>"${out_prefix}.err"

  local sec rss hits
  sec="$(cut -f1 "$tm")"
  rss="$(cut -f2 "$tm")"
  hits="$(to_num_or_na "$(extract_positive_hits "$log")")"
  echo "$sec|$rss|$hits"
}

run_case() {
  local tool="$1"
  local layout="$2"
  local rep="$3"
  local run_tag="${tool}_${layout}_r${rep}"
  local out_dir="$RUN_DIR/$run_tag"
  mkdir -p "$out_dir"

  local maker_cmd=()
  local filter_path=""
  case "$tool/$layout" in
    cpp/classic)
      maker_cmd=("$CPP_MAKER" -p ref -o "$out_dir" -k "$KMER_SIZE" -f "$FPR" -t "$THREADS" "$DATA_DIR/ref.fa")
      filter_path="$out_dir/ref.bf"
      ;;
    rust/classic)
      maker_cmd=("$SBB_BIN" maker -p ref -o "$out_dir" -k "$KMER_SIZE" -f "$FPR" -t "$THREADS" --classic "$DATA_DIR/ref.fa")
      filter_path="$out_dir/ref.bf.zst"
      ;;
    rust/blocked)
      maker_cmd=("$SBB_BIN" maker -p ref -o "$out_dir" -k "$KMER_SIZE" -f "$FPR" -t "$THREADS" --block_words "$BLOCK_WORDS" "$DATA_DIR/ref.fa")
      filter_path="$out_dir/ref.bf.zst"
      ;;
    *)
      echo "unknown case: $tool/$layout" >&2
      exit 1
      ;;
  esac

  local maker_log="$out_dir/maker.out"
  local maker_time="$out_dir/maker.time"
  /usr/bin/time -f "%e\t%M" -o "$maker_time" "${maker_cmd[@]}" >"$maker_log" 2>"$out_dir/maker.err"

  local build_sec build_rss index_bytes
  build_sec="$(cut -f1 "$maker_time")"
  build_rss="$(cut -f2 "$maker_time")"
  index_bytes="$(stat -c%s "$filter_path" 2>/dev/null || echo na)"

  local qpos qnear qneg
  qpos="$(query_one_set "$tool" "$filter_path" "$POS_READS" "$out_dir/pos" "$out_dir/pos.out" "$out_dir/pos.time")"
  qnear="$(query_one_set "$tool" "$filter_path" "$NEAR_READS" "$out_dir/near" "$out_dir/near.out" "$out_dir/near.time")"
  qneg="$(query_one_set "$tool" "$filter_path" "$NEG_READS" "$out_dir/neg" "$out_dir/neg.out" "$out_dir/neg.time")"

  local pos_sec pos_rss pos_hits
  local near_sec near_rss near_hits
  local neg_sec neg_rss neg_hits
  IFS='|' read -r pos_sec pos_rss pos_hits <<< "$qpos"
  IFS='|' read -r near_sec near_rss near_hits <<< "$qnear"
  IFS='|' read -r neg_sec neg_rss neg_hits <<< "$qneg"

  pos_hits="$(to_num_or_na "$pos_hits")"
  near_hits="$(to_num_or_na "$near_hits")"
  neg_hits="$(to_num_or_na "$neg_hits")"

  local tp fp tn fn precision recall f1 specificity fprate near_recall
  tp="na"; fp="na"; tn="na"; fn="na"
  precision="na"; recall="na"; f1="na"; specificity="na"; fprate="na"; near_recall="na"

  if [[ "$pos_hits" != "na" && "$neg_hits" != "na" && "$near_hits" != "na" ]]; then
    tp="$pos_hits"
    fp="$neg_hits"
    fn="$((POS_COUNT - pos_hits))"
    tn="$((NEG_COUNT - neg_hits))"
    precision="$(awk -v tp="$tp" -v fp="$fp" 'BEGIN{d=tp+fp; if(d>0){printf "%.6f", tp/d}else{print "na"}}')"
    recall="$(awk -v tp="$tp" -v fn="$fn" 'BEGIN{d=tp+fn; if(d>0){printf "%.6f", tp/d}else{print "na"}}')"
    f1="$(awk -v p="$precision" -v r="$recall" 'BEGIN{if(p=="na" || r=="na"){print "na"}else{d=p+r; if(d>0){printf "%.6f", 2*p*r/d}else{print "na"}}}')"
    specificity="$(awk -v tn="$tn" -v fp="$fp" 'BEGIN{d=tn+fp; if(d>0){printf "%.6f", tn/d}else{print "na"}}')"
    fprate="$(awk -v fp="$fp" -v n="$NEG_COUNT" 'BEGIN{if(n>0){printf "%.6f", fp/n}else{print "na"}}')"
    near_recall="$(awk -v h="$near_hits" -v n="$NEAR_COUNT" 'BEGIN{if(n>0){printf "%.6f", h/n}else{print "na"}}')"
  fi

  local query_sec_total query_rss_kb_max reads_per_sec_total
  query_sec_total="$(awk -v a="$pos_sec" -v b="$near_sec" -v c="$neg_sec" 'BEGIN{printf "%.6f", a+b+c}')"
  query_rss_kb_max="$(awk -v a="$pos_rss" -v b="$near_rss" -v c="$neg_rss" 'BEGIN{m=a; if(b>m)m=b; if(c>m)m=c; print m}')"
  reads_per_sec_total="$(awk -v n="$TOTAL_QUERY_READS" -v t="$query_sec_total" 'BEGIN{if(t>0){printf "%.6f", n/t}else{print "na"}}')"

  echo -e "single_ref_accuracy\t$tool\t$layout\t$rep\t$THREADS\t$KMER_SIZE\t$FPR\t$CAT_SCORE\t$POS_COUNT\t$NEAR_COUNT\t$NEG_COUNT\t$pos_hits\t$near_hits\t$neg_hits\t$tp\t$fp\t$tn\t$fn\t$precision\t$recall\t$f1\t$specificity\t$fprate\t$near_recall\t$build_sec\t$build_rss\t$query_sec_total\t$query_rss_kb_max\t$reads_per_sec_total\t$index_bytes\t$filter_path" >> "$RESULTS"
}

for rep in $(seq 1 "$REPEATS"); do
  run_case cpp classic "$rep"
  run_case rust classic "$rep"
  run_case rust blocked "$rep"
done

echo "results_tsv=$RESULTS"

if [[ "$PLOT" == "1" ]]; then
  if command -v python3 >/dev/null 2>&1; then
    python3 "$PLOT_SCRIPT" --tsv "$RESULTS"
  else
    echo "[accuracy] python3 not found; skipping plots" >&2
  fi
fi
