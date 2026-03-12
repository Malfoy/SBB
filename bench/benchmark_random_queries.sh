#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$ROOT_DIR/bench/common.sh"
bench_init_cpp_env

BENCH_ROOT="${BENCH_ROOT:-/tmp/sbb-random-queries}"
DATA_DIR="$BENCH_ROOT/data"
RUN_DIR="$BENCH_ROOT/runs"
RESULT_DIR="$BENCH_ROOT/results"

INDEX_FASTA="${INDEX_FASTA:-$ROOT_DIR/bench/refs/ecoli_k12_mg1655.fa}"
QUERY_FASTA_INPUT="${QUERY_FASTA:-/home/nadine/Code/SBB/ecoli1d2.fa}"
QUERY_COUNT="${QUERY_COUNT:-1000000000}"
QUERY_LEN="${QUERY_LEN:-31}"
SEED="${SEED:-1337}"

THREADS="${THREADS:-$(nproc)}"
REPEATS="${REPEATS:-1}"
KMER_SIZE="${KMER_SIZE:-31}"
HASH_NUM="${HASH_NUM:-3}"
BIT_EXP="${BIT_EXP:-27}"
BIT_EXP_LIST="${BIT_EXP_LIST:-25,26,27,28,29,30,31,32,33,34,35}"
CPP_NUM_ELE="${CPP_NUM_ELE:-5000000}"
ALLOW_NON_COMPARABLE="${ALLOW_NON_COMPARABLE:-0}"
CAT_SCORE="${CAT_SCORE:-0.01}"
BLOCK_WORDS="${BLOCK_WORDS:-8}"
BLOOMY_S="${BLOOMY_S:-$KMER_SIZE}"
BLOOMY_S_LIST="${BLOOMY_S_LIST:-31,30,28,24}"
BLOOMY_M="${BLOOMY_M:-21}"
FIXED_MAKER_FPR="0.00001"

RUN_CPP="${RUN_CPP:-1}"
RUN_RUST_CLASSIC="${RUN_RUST_CLASSIC:-1}"
RUN_RUST_BLOCKED="${RUN_RUST_BLOCKED:-1}"
RUN_RUST_BLOOMY="${RUN_RUST_BLOOMY:-1}"

SBB_BIN="${SBB_BIN:-$ROOT_DIR/target/release/sbb}"
BENCHRAND="${BENCHRAND:-$ROOT_DIR/target/release/sbbbenchrand}"

mkdir -p "$DATA_DIR" "$RUN_DIR" "$RESULT_DIR"

require_file() {
  local path="$1"
  if [[ ! -f "$path" ]]; then
    echo "missing file: $path" >&2
    exit 1
  fi
}

require_bin() {
  local path="$1"
  local name="$2"
  if [[ ! -x "$path" ]]; then
    echo "missing executable for $name: $path" >&2
    exit 1
  fi
}

usage() {
  cat <<'EOF'
usage: benchmark_random_queries.sh

Environment variables:
  INDEX_FASTA        FASTA indexed by all tools
  QUERY_FASTA        optional FASTA query dataset; when set, QUERY_COUNT/QUERY_LEN/SEED are ignored
  QUERY_COUNT        number of random query sequences
  QUERY_LEN          query sequence length (default 31)
  CAT_SCORE          categorizer threshold (default 0.001)
  THREADS            worker threads
  REPEATS            repeats for each tool/layout
  KMER_SIZE          k-mer size for index creation
  HASH_NUM           optional fixed hash count passed to maker
  BIT_EXP            bloom size exponent (for example 10 => 2^10 bits)
  BIT_EXP_LIST       comma-separated exponents (for example 24,25,26); overrides BIT_EXP
  BLOOMY_M           bloomybloom minimizer length m (default 21)
  BLOOMY_S           bloomybloom s-mer length (single value)
  BLOOMY_S_LIST      comma-separated bloomybloom s-mer values (overrides BLOOMY_S)
  CPP_NUM_ELE        required with BIT_EXP+RUN_CPP for C++ FPR derivation
  ALLOW_NON_COMPARABLE  set to 1 to bypass strict comparable-memory checks
  RUN_CPP            1/0 enable C++ benchmark case
  RUN_RUST_CLASSIC   1/0 enable Rust classic case
  RUN_RUST_BLOCKED   1/0 enable Rust blocked case
  RUN_RUST_BLOOMY    1/0 enable Rust bloomybloom case
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

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

format_int_commas() {
  local v="${1:-}"
  if [[ ! "$v" =~ ^-?[0-9]+$ ]]; then
    echo "$v"
    return
  fi
  local sign=""
  if [[ "$v" == -* ]]; then
    sign="-"
    v="${v#-}"
  fi
  local out=""
  while [[ ${#v} -gt 3 ]]; do
    out=",${v: -3}${out}"
    v="${v:0:${#v}-3}"
  done
  echo "${sign}${v}${out}"
}

read_fasta_stats() {
  local fasta="$1"
  local count_var="$2"
  local len_var="$3"

  local stats
  stats="$(
    awk '
      BEGIN {count=0; seq_len=0; fixed_len=-1; mixed=0}
      /^>/ {
        if (count > 0) {
          if (fixed_len < 0) {
            fixed_len = seq_len
          } else if (seq_len != fixed_len) {
            mixed = 1
          }
        }
        count++
        seq_len = 0
        next
      }
      {
        gsub(/[[:space:]]/, "", $0)
        seq_len += length($0)
      }
      END {
        if (count > 0) {
          if (fixed_len < 0) {
            fixed_len = seq_len
          } else if (seq_len != fixed_len) {
            mixed = 1
          }
        }
        if (count == 0) {
          print "0\tna"
        } else if (mixed == 1) {
          print count "\tna"
        } else {
          print count "\t" fixed_len
        }
      }
    ' "$fasta"
  )"

  local count_value len_value
  IFS=$'\t' read -r count_value len_value <<< "$stats"
  if [[ "$count_value" == "0" ]]; then
    echo "QUERY_FASTA has no FASTA records: $fasta" >&2
    exit 1
  fi
  printf -v "$count_var" '%s' "$count_value"
  printf -v "$len_var" '%s' "$len_value"
}

extract_positive_hits() {
  local log="$1"
  local v=""
  local run_dir
  run_dir="$(dirname "$log")"
  local cpp_summary="$run_dir/_summary.tsv"

  if [[ -s "$cpp_summary" ]]; then
    v="$(awk -F'\t' 'NR>1 && $1!="multiMatch" && $1!="noMatch"{sum+=$2;seen=1} END{if(seen) print sum}' "$cpp_summary")"
    if [[ -n "$v" ]]; then
      echo "$v"
      return
    fi
  fi

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

compute_cpp_fpr() {
  local bit_len="$1"
  local num_ele="$2"
  local hash_num="$3"
  awk -v m="$bit_len" -v n="$num_ele" -v k="$hash_num" 'BEGIN {
    p = (1 - exp((-k * n) / m)) ^ k
    if (p <= 0) {
      p = 1e-300
    }
    printf "%.12g", p
  }'
}

if [[ ! "$ALLOW_NON_COMPARABLE" =~ ^[01]$ ]]; then
  echo "ALLOW_NON_COMPARABLE must be 0 or 1, got: $ALLOW_NON_COMPARABLE" >&2
  exit 1
fi

if [[ -z "$QUERY_FASTA_INPUT" ]]; then
  if [[ "$QUERY_COUNT" =~ ^[0-9]+$ ]] && (( QUERY_COUNT > 0 )); then
    :
  else
    echo "QUERY_COUNT must be a positive integer, got: $QUERY_COUNT" >&2
    exit 1
  fi

  if [[ "$QUERY_LEN" =~ ^[0-9]+$ ]] && (( QUERY_LEN > 0 )); then
    :
  else
    echo "QUERY_LEN must be a positive integer, got: $QUERY_LEN" >&2
    exit 1
  fi
fi

if [[ "$REPEATS" =~ ^[0-9]+$ ]] && (( REPEATS > 0 )); then
  :
else
  echo "REPEATS must be a positive integer, got: $REPEATS" >&2
  exit 1
fi

if [[ "$KMER_SIZE" =~ ^[0-9]+$ ]] && (( KMER_SIZE > 0 )); then
  :
else
  echo "KMER_SIZE must be a positive integer, got: $KMER_SIZE" >&2
  exit 1
fi

BIT_EXP_VALUES=()
if [[ -n "$BIT_EXP_LIST" ]]; then
  IFS=',' read -r -a raw_exp_values <<< "$BIT_EXP_LIST"
  for value in "${raw_exp_values[@]}"; do
    value="${value//[[:space:]]/}"
    [[ -n "$value" ]] || continue
    if [[ ! "$value" =~ ^[0-9]+$ ]] || (( value < 1 || value > 62 )); then
      echo "BIT_EXP_LIST contains invalid exponent (must be in [1,62]): $value" >&2
      exit 1
    fi
    BIT_EXP_VALUES+=("$value")
  done
  if (( ${#BIT_EXP_VALUES[@]} == 0 )); then
    echo "BIT_EXP_LIST did not contain any valid values" >&2
    exit 1
  fi
else
  if [[ ! "$BIT_EXP" =~ ^[0-9]+$ ]] || (( BIT_EXP < 1 || BIT_EXP > 62 )); then
    echo "BIT_EXP must be an integer in [1, 62], got: $BIT_EXP" >&2
    exit 1
  fi
  BIT_EXP_VALUES=("$BIT_EXP")
fi

if (( RUN_RUST_BLOOMY == 1 && (KMER_SIZE % 2) == 0 )); then
  echo "RUN_RUST_BLOOMY=1 requires an odd KMER_SIZE, got $KMER_SIZE" >&2
  exit 1
fi

BLOOMY_S_VALUES=()
if (( RUN_RUST_BLOOMY == 1 )); then
  if [[ -n "$BLOOMY_S_LIST" ]]; then
    IFS=',' read -r -a raw_bloomy_s_values <<< "$BLOOMY_S_LIST"
    for value in "${raw_bloomy_s_values[@]}"; do
      value="${value//[[:space:]]/}"
      [[ -n "$value" ]] || continue
      if [[ ! "$value" =~ ^[0-9]+$ ]] || (( value < 1 || value > KMER_SIZE )); then
        echo "BLOOMY_S_LIST contains invalid s (must be in [1, KMER_SIZE]): $value (KMER_SIZE=$KMER_SIZE)" >&2
        exit 1
      fi
      BLOOMY_S_VALUES+=("$value")
    done
    if (( ${#BLOOMY_S_VALUES[@]} == 0 )); then
      echo "BLOOMY_S_LIST did not contain any valid values" >&2
      exit 1
    fi
  else
    if [[ ! "$BLOOMY_S" =~ ^[0-9]+$ ]] || (( BLOOMY_S < 1 || BLOOMY_S > KMER_SIZE )); then
      echo "BLOOMY_S must be an integer in [1, KMER_SIZE]; got BLOOMY_S=$BLOOMY_S KMER_SIZE=$KMER_SIZE" >&2
      exit 1
    fi
    BLOOMY_S_VALUES=("$BLOOMY_S")
  fi
  if [[ ! "$BLOOMY_M" =~ ^[0-9]+$ ]] || (( BLOOMY_M < 1 || BLOOMY_M > 31 || BLOOMY_M > KMER_SIZE )); then
    echo "BLOOMY_M must be an integer in [1, min(31, KMER_SIZE)]; got BLOOMY_M=$BLOOMY_M KMER_SIZE=$KMER_SIZE" >&2
    exit 1
  fi
fi

if (( (RUN_CPP == 1 || RUN_RUST_CLASSIC == 1 || RUN_RUST_BLOCKED == 1) && KMER_SIZE > 32 )); then
  echo "KMER_SIZE > 32 is only supported by bloomybloom in this repo; disable classic/blocked/cpp or lower KMER_SIZE" >&2
  exit 1
fi

if [[ -n "$HASH_NUM" ]]; then
  if [[ ! "$HASH_NUM" =~ ^[0-9]+$ ]] || (( HASH_NUM <= 0 )); then
    echo "HASH_NUM must be a positive integer, got: $HASH_NUM" >&2
    exit 1
  fi
fi

if [[ -z "$HASH_NUM" ]]; then
  echo "HASH_NUM must be set when BIT_EXP is set" >&2
  exit 1
fi
if (( RUN_CPP == 1 )); then
  if [[ ! "$CPP_NUM_ELE" =~ ^[0-9]+$ ]] || (( CPP_NUM_ELE <= 0 )); then
    echo "CPP_NUM_ELE must be a positive integer when RUN_CPP=1" >&2
    exit 1
  fi
fi

if (( ALLOW_NON_COMPARABLE == 0 )); then
  if [[ -z "$HASH_NUM" ]]; then
    echo "HASH_NUM must be set for comparable-memory FPR benchmarking (or set ALLOW_NON_COMPARABLE=1)" >&2
    exit 1
  fi
  if (( RUN_CPP == 1 )) && [[ -z "$CPP_NUM_ELE" ]]; then
    echo "CPP_NUM_ELE must be set when RUN_CPP=1 and comparable-memory mode is enabled" >&2
    exit 1
  fi
fi

require_file "$INDEX_FASTA"
if [[ -n "$QUERY_FASTA_INPUT" ]]; then
  require_file "$QUERY_FASTA_INPUT"
fi

echo "[random] building benchmark binaries"
cargo build --release --bin sbb --bin sbbbenchrand >/tmp/sbb_random_queries_build.log

require_bin "$SBB_BIN" "sbb"
require_bin "$BENCHRAND" "sbbbenchrand"
if (( RUN_CPP == 1 )); then
  require_bin "$CPP_MAKER" "cpp maker"
  require_bin "$CPP_CAT" "cpp categorizer"
fi

GENERATED_QUERY_FASTA="$DATA_DIR/random_queries.fa"
QUERY_FASTA_PATH="$GENERATED_QUERY_FASTA"
QUERY_COUNT_VALUE="$QUERY_COUNT"
QUERY_LEN_VALUE="$QUERY_LEN"
if [[ -n "$QUERY_FASTA_INPUT" ]]; then
  QUERY_FASTA_PATH="$QUERY_FASTA_INPUT"
  read_fasta_stats "$QUERY_FASTA_PATH" QUERY_COUNT_VALUE QUERY_LEN_VALUE
  echo "[random] using provided query FASTA: $QUERY_FASTA_PATH (count=$QUERY_COUNT_VALUE len=$QUERY_LEN_VALUE)"
else
  QUERY_META="$DATA_DIR/random_queries.meta"
  EXPECTED_META=$'query_count='"$QUERY_COUNT"$'\nquery_len='"$QUERY_LEN"$'\nseed='"$SEED"$'\n'
  CURRENT_META=""
  if [[ -f "$QUERY_META" ]]; then
    CURRENT_META="$(cat "$QUERY_META")"
  fi
  if [[ ! -s "$GENERATED_QUERY_FASTA" || "$CURRENT_META"$'\n' != "$EXPECTED_META" ]]; then
    echo "[random] generating queries with Rust generator"
    "$BENCHRAND" \
      --out_fa "$GENERATED_QUERY_FASTA" \
      --count "$QUERY_COUNT" \
      --len "$QUERY_LEN" \
      --seed "$SEED" \
      --threads "$THREADS" \
      >/tmp/sbb_random_queries_gen.log
    printf "%s" "$EXPECTED_META" > "$QUERY_META"
  fi
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
RESULTS="$RESULT_DIR/random_queries_${STAMP}.tsv"
echo -e "benchmark\ttool\tlayout\trepeat\tthreads\tkmer_size\tfpr\thash_num\tbit_exp\ttarget_bit_len\tactual_bit_len\tscore\tbloomy_m\tbloomy_s\tquery_count\tquery_len\tpositive_hits\tbuild_sec\tbuild_rss_kb\tquery_sec\tquery_rss_kb\tindex_bytes\tserialized_index_bits\tfilter_path" > "$RESULTS"
CPP_BASE_ACTUAL_BITS=""

run_case() {
  local tool="$1"
  local layout="$2"
  local rep="$3"
  local bit_exp="$4"
  local bloomy_s_override="${5:-}"
  local bit_len="$((1 << bit_exp))"
  local out_dir_suffix="_e${bit_exp}_r${rep}"
  if [[ "$tool/$layout" == "rust/bloomybloom" ]]; then
    if [[ -z "$bloomy_s_override" ]]; then
      echo "internal error: bloomybloom run missing s-mer value" >&2
      exit 1
    fi
    out_dir_suffix="${out_dir_suffix}_s${bloomy_s_override}"
  fi
  local out_dir="$RUN_DIR/${tool}_${layout}${out_dir_suffix}"
  mkdir -p "$out_dir"

  local maker_cmd=()
  local query_cmd=()
  local filter_path=""
  local hash_value="${HASH_NUM:-na}"
  local bit_len_value="$bit_len"
  local fpr_value="$FIXED_MAKER_FPR"
  local bloomy_m_value="na"
  local bloomy_s_value="na"

  case "$tool/$layout" in
    cpp/classic)
      fpr_value="$(compute_cpp_fpr "$bit_len" "$CPP_NUM_ELE" "$HASH_NUM")"
      maker_cmd=(
        "$CPP_MAKER"
        -p ref
        -o "$out_dir"
        -k "$KMER_SIZE"
        -f "$fpr_value"
        -t "$THREADS"
      )
      if [[ -n "$CPP_NUM_ELE" ]]; then
        maker_cmd+=(-n "$CPP_NUM_ELE")
      fi
      if [[ -n "$HASH_NUM" ]]; then
        maker_cmd+=(-g "$HASH_NUM")
      fi
      maker_cmd+=("$INDEX_FASTA")
      filter_path="$out_dir/ref.bf"
      query_cmd=(
        "$CPP_CAT"
        -f "$filter_path"
        -t "$THREADS"
        -S binomial
        -s "$CAT_SCORE"
        "$QUERY_FASTA_PATH"
      )
      ;;
    rust/classic)
      maker_cmd=(
        "$SBB_BIN"
        maker
        -p ref
        -o "$out_dir"
        -k "$KMER_SIZE"
        -f "$FIXED_MAKER_FPR"
        -t "$THREADS"
        --classic
      )
      if [[ -n "$HASH_NUM" ]]; then
        maker_cmd+=(-g "$HASH_NUM")
      fi
      maker_cmd+=(--bit_len "$bit_len")
      maker_cmd+=("$INDEX_FASTA")
      filter_path="$out_dir/ref.bf.zst"
      query_cmd=(
        "$SBB_BIN"
        categorizer
        -f "$filter_path"
        -t "$THREADS"
        -s "$CAT_SCORE"
        "$QUERY_FASTA_PATH"
      )
      ;;
    rust/blocked)
      maker_cmd=(
        "$SBB_BIN"
        maker
        -p ref
        -o "$out_dir"
        -k "$KMER_SIZE"
        -f "$FIXED_MAKER_FPR"
        -t "$THREADS"
        --block_words "$BLOCK_WORDS"
      )
      if [[ -n "$HASH_NUM" ]]; then
        maker_cmd+=(-g "$HASH_NUM")
      fi
      maker_cmd+=(--bit_len "$bit_len")
      maker_cmd+=("$INDEX_FASTA")
      filter_path="$out_dir/ref.bf.zst"
      query_cmd=(
        "$SBB_BIN"
        categorizer
        -f "$filter_path"
        -t "$THREADS"
        -s "$CAT_SCORE"
        "$QUERY_FASTA_PATH"
      )
      ;;
    rust/bloomybloom)
      bloomy_m_value="$BLOOMY_M"
      bloomy_s_value="$bloomy_s_override"
      maker_cmd=(
        env
        "SBB_BLOOMY_M=$BLOOMY_M"
        "SBB_BLOOMY_SMER=$bloomy_s_override"
        "$SBB_BIN"
        maker
        -p ref
        -o "$out_dir"
        -k "$KMER_SIZE"
        -f "$FIXED_MAKER_FPR"
        -t "$THREADS"
        --bloomybloom
      )
      if [[ -n "$HASH_NUM" ]]; then
        maker_cmd+=(-g "$HASH_NUM")
      fi
      maker_cmd+=(--bit_len "$bit_len")
      maker_cmd+=("$INDEX_FASTA")
      filter_path="$out_dir/ref.bf.zst"
      query_cmd=(
        "$SBB_BIN"
        categorizer
        -f "$filter_path"
        -t "$THREADS"
        -s "$CAT_SCORE"
        "$QUERY_FASTA_PATH"
      )
      ;;
    *)
      echo "unknown benchmark case: $tool/$layout" >&2
      exit 1
      ;;
  esac

  local maker_time="$out_dir/maker.time"
  local query_time="$out_dir/query.time"
  local maker_log="$out_dir/maker.out"
  local query_log="$out_dir/query.out"
  local maker_err="$out_dir/maker.err"
  local query_err="$out_dir/query.err"

  local maker_cmd_q query_cmd_q
  printf -v maker_cmd_q '%q ' "${maker_cmd[@]}"
  printf -v query_cmd_q '%q ' "${query_cmd[@]}"

  /usr/bin/time -f "%e\t%M" -o "$maker_time" bash -lc "set -e; $maker_cmd_q > $(printf '%q' "$maker_log") 2> $(printf '%q' "$maker_err")"

  if [[ "$tool" == "cpp" ]]; then
    /usr/bin/time -f "%e\t%M" -o "$query_time" bash -lc "set -e; cd $(printf '%q' "$out_dir"); $query_cmd_q > $(printf '%q' "$query_log") 2> $(printf '%q' "$query_err")"
  else
    /usr/bin/time -f "%e\t%M" -o "$query_time" bash -lc "set -e; $query_cmd_q > $(printf '%q' "$query_log") 2> $(printf '%q' "$query_err")"
  fi

  local build_sec build_rss query_sec query_rss positive_hits index_bytes
  local serialized_index_bits="na"
  local actual_bit_len="na"
  build_sec="$(cut -f1 "$maker_time")"
  build_rss="$(cut -f2 "$maker_time")"
  query_sec="$(cut -f1 "$query_time")"
  query_rss="$(cut -f2 "$query_time")"
  positive_hits="$(to_num_or_na "$(extract_positive_hits "$query_log")")"
  index_bytes="$(stat -c%s "$filter_path" 2>/dev/null || echo na)"
  if [[ "$index_bytes" =~ ^[0-9]+$ ]]; then
    serialized_index_bits="$(( index_bytes * 8 ))"
  fi

  if [[ "$tool" == "cpp" ]]; then
    actual_bit_len="$(sed -n 's/^Allocating \([0-9][0-9]*\) bits.*/\1/p' "$maker_err" | tail -n1)"
  else
    actual_bit_len="$(sed -n 's/.* bits=\([0-9][0-9]*\).*/\1/p' "$maker_log" | tail -n1)"
  fi
  actual_bit_len="$(to_num_or_na "$actual_bit_len")"

  if [[ "$tool" == "cpp" && "$actual_bit_len" =~ ^[0-9]+$ ]]; then
    CPP_BASE_ACTUAL_BITS="$actual_bit_len"
  fi

  local overhead_pct="na"
  if [[ -n "$CPP_BASE_ACTUAL_BITS" && "$CPP_BASE_ACTUAL_BITS" =~ ^[0-9]+$ && "$actual_bit_len" =~ ^[0-9]+$ ]]; then
    overhead_pct="$(awk -v x="$actual_bit_len" -v c="$CPP_BASE_ACTUAL_BITS" 'BEGIN{printf "%.6f", ((x-c)*100.0)/c}')"
  fi

  echo -e "random_queries\t$tool\t$layout\t$rep\t$THREADS\t$KMER_SIZE\t$fpr_value\t$hash_value\t$bit_exp\t$bit_len_value\t$actual_bit_len\t$CAT_SCORE\t$bloomy_m_value\t$bloomy_s_value\t$QUERY_COUNT_VALUE\t$QUERY_LEN_VALUE\t$positive_hits\t$build_sec\t$build_rss\t$query_sec\t$query_rss\t$index_bytes\t$serialized_index_bits\t$filter_path" >> "$RESULTS"
  if [[ "$tool" == "cpp" ]]; then
    echo "[random] (${overhead_pct}% overhead vs cpp actual bits)"
  fi
  local positive_hits_pretty
  positive_hits_pretty="$(format_int_commas "$positive_hits")"
  if [[ "$tool/$layout" == "rust/bloomybloom" ]]; then
    echo "[random] result tool=$tool layout=$layout bit_exp=$bit_exp s=$bloomy_s_value repeat=$rep positive_hits=$positive_hits_pretty"
  else
    echo "[random] result tool=$tool layout=$layout bit_exp=$bit_exp repeat=$rep positive_hits=$positive_hits_pretty"
  fi
}

for bit_exp in "${BIT_EXP_VALUES[@]}"; do
  for rep in $(seq 1 "$REPEATS"); do
    CPP_BASE_ACTUAL_BITS=""
    if (( RUN_CPP == 1 )); then
      run_case cpp classic "$rep" "$bit_exp"
    fi
    if (( RUN_RUST_CLASSIC == 1 )); then
      run_case rust classic "$rep" "$bit_exp"
    fi
    if (( RUN_RUST_BLOCKED == 1 )); then
      run_case rust blocked "$rep" "$bit_exp"
    fi
    if (( RUN_RUST_BLOOMY == 1 )); then
      for bloomy_s in "${BLOOMY_S_VALUES[@]}"; do
        run_case rust bloomybloom "$rep" "$bit_exp" "$bloomy_s"
      done
    fi
  done
done

echo
echo "query_fasta=$QUERY_FASTA_PATH"
echo "results_tsv=$RESULTS"
