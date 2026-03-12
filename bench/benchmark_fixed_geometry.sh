#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$ROOT_DIR/bench/common.sh"
bench_init_cpp_env

BENCH_ROOT="${BENCH_ROOT:-/tmp/sbb-bench-fixed-geometry}"
RUN_DIR="$BENCH_ROOT/runs"
RESULT_DIR="$BENCH_ROOT/results"

INDEX_DATASET="${INDEX_DATASET:-/home/nadine/Downloads/HG.fa}"
QUERY_DATASET="${QUERY_DATASET:-/home/nadine/Code/SBB/p.HUMANHIFI10x.fa}"
BIT_LEN="${BIT_LEN:-137438953472}"
HASH_NUM="${HASH_NUM:-8}"
HASH_LIST="${HASH_LIST:-}"

KMER_SIZE="${KMER_SIZE:-31}"
THREADS="${THREADS:-$(nproc)}"
CAT_SCORE="${CAT_SCORE:-1.0}"
BLOCK_WORDS="${BLOCK_WORDS:-8}"
BLOOMY_S="${BLOOMY_S:-31}"
REPEATS="${REPEATS:-1}"
CPP_NUM_ELE="${CPP_NUM_ELE:-3000000000}"

RUN_CPP="${RUN_CPP:-1}"
RUN_RUST_CLASSIC="${RUN_RUST_CLASSIC:-1}"
RUN_RUST_BLOCKED="${RUN_RUST_BLOCKED:-1}"
RUN_RUST_BLOOMY="${RUN_RUST_BLOOMY:-1}"

SBB_BIN="${SBB_BIN:-$ROOT_DIR/target/release/sbb}"

mkdir -p "$RUN_DIR" "$RESULT_DIR"

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
  echo "usage: $0 [--cpp-num-ele N]"
  echo
  echo "options:"
  echo "  -n, --cpp-num-ele N  pass expected element count to BioBloomMaker only"
  echo "  -h, --help           show this help"
}

is_power_of_two() {
  local value="$1"
  (( value > 0 && (value & (value - 1)) == 0 ))
}

tsv_value() {
  local file="$1"
  local key="$2"
  awk -F'\t' -v k="$key" '$1==k{print $2}' "$file" | tail -n1
}

to_num_or_na() {
  local value="${1:-}"
  if [[ -n "$value" ]]; then
    echo "$value"
  else
    echo "na"
  fi
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

mean_ns_per_op() {
  local seconds="$1"
  local ops="$2"
  awk -v t="$seconds" -v n="$ops" 'BEGIN {
    if (n > 0 && t >= 0) {
      printf "%.3f", (t * 1000000000) / n
    } else {
      print "na"
    }
  }'
}

rate_per_sec() {
  local count="$1"
  local seconds="$2"
  awk -v n="$count" -v t="$seconds" 'BEGIN {
    if (n >= 0 && t > 0) {
      printf "%.6f", n / t
    } else {
      print "na"
    }
  }'
}

sum_cpu_seconds() {
  local user_sec="$1"
  local sys_sec="$2"
  awk -v u="$user_sec" -v s="$sys_sec" 'BEGIN { printf "%.2f", u + s }'
}

rss_kb_to_mb() {
  local rss_kb="$1"
  awk -v kb="$rss_kb" 'BEGIN { printf "%.2f", kb / 1024 }'
}

cleanup_serialized_filters() {
  local out_dir="$1"
  rm -f \
    "$out_dir/ref.bf" \
    "$out_dir/ref.bf.zst" \
    "$out_dir/ref.bf.zst.shard"*.zst \
    "$out_dir/ref.txt"
}

cpp_maker_insert_ops() {
  local file="$1"
  sed -n 's/^Approximated (due to false positives) total unique k-mers in reference files \([0-9][0-9]*\)$/\1/p' "$file" | tail -n1
}

cpp_query_total_reads() {
  local file="$1"
  sed -n 's/^Total Reads: \([0-9][0-9]*\)$/\1/p' "$file" | tail -n1
}

while (( $# > 0 )); do
  case "$1" in
    -n|--cpp-num-ele)
      if (( $# < 2 )); then
        echo "missing value for $1" >&2
        usage >&2
        exit 1
      fi
      CPP_NUM_ELE="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --)
      shift
      break
      ;;
    *)
      echo "unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

if (( $# > 0 )); then
  echo "unexpected positional arguments: $*" >&2
  usage >&2
  exit 1
fi

run_case() {
  local tool="$1"
  local layout="$2"
  local rep="$3"
  local impl="${tool}_${layout}"
  local out_dir="$RUN_DIR/${tool}_${layout}_h${HASH_NUM}_r${rep}"
  local maker_log="$out_dir/maker.out"
  local query_log="$out_dir/query.out"
  local maker_time="$out_dir/maker.time"
  local query_time="$out_dir/query.time"
  mkdir -p "$out_dir"
  cleanup_serialized_filters "$out_dir"

  local maker_cmd=()
  local query_cmd=()
  local cpp_fpr_for_run="na"
  local maker_err="$out_dir/maker.err"
  local query_err="$out_dir/query.err"
  local filter_path=""

  echo "[bench] running impl=$impl repeat=$rep/$REPEATS"

  case "$tool/$layout" in
    cpp/classic)
      cpp_fpr_for_run="$(compute_cpp_fpr "$BIT_LEN" "$CPP_NUM_ELE" "$HASH_NUM")"
      maker_cmd=(
        "$CPP_MAKER"
        -p ref
        -o "$out_dir"
        -k "$KMER_SIZE"
      )
      if [[ -n "$CPP_NUM_ELE" ]]; then
        maker_cmd+=(-n "$CPP_NUM_ELE")
      fi
      maker_cmd+=(
        -f "$cpp_fpr_for_run"
        -g "$HASH_NUM"
        -t "$THREADS"
        "$INDEX_DATASET"
      )
      filter_path="$out_dir/ref.bf"
      query_cmd=(
        "$CPP_CAT"
        -f "$filter_path"
        -t "$THREADS"
        -s "$CAT_SCORE"
        "$QUERY_DATASET"
      )
      ;;
    rust/classic)
      maker_cmd=(
        "$SBB_BIN"
        maker
        -p ref
        -o "$out_dir"
        -k "$KMER_SIZE"
        -f 0.5
        -g "$HASH_NUM"
        -t "$THREADS"
        --bit_len "$BIT_LEN"
        --classic
        "$INDEX_DATASET"
      )
      filter_path="$out_dir/ref.bf.zst"
      query_cmd=(
        "$SBB_BIN"
        categorizer
        -f "$filter_path"
        -t "$THREADS"
        -s "$CAT_SCORE"
        "$QUERY_DATASET"
      )
      ;;
    rust/blocked)
      maker_cmd=(
        "$SBB_BIN"
        maker
        -p ref
        -o "$out_dir"
        -k "$KMER_SIZE"
        -f 0.5
        -g "$HASH_NUM"
        -t "$THREADS"
        --bit_len "$BIT_LEN"
        --block_words "$BLOCK_WORDS"
        "$INDEX_DATASET"
      )
      filter_path="$out_dir/ref.bf.zst"
      query_cmd=(
        "$SBB_BIN"
        categorizer
        -f "$filter_path"
        -t "$THREADS"
        -s "$CAT_SCORE"
        "$QUERY_DATASET"
      )
      ;;
    rust/bloomybloom)
      maker_cmd=(
        env
        "SBB_BLOOMY_SMER=$BLOOMY_S"
        "$SBB_BIN"
        maker
        -p ref
        -o "$out_dir"
        -k "$KMER_SIZE"
        -f 0.5
        -g "$HASH_NUM"
        -t "$THREADS"
        --bit_len "$BIT_LEN"
        --bloomybloom
        "$INDEX_DATASET"
      )
      filter_path="$out_dir/ref.bf.zst"
      query_cmd=(
        "$SBB_BIN"
        categorizer
        -f "$filter_path"
        -t "$THREADS"
        -s "$CAT_SCORE"
        "$QUERY_DATASET"
      )
      ;;
    *)
      echo "unknown benchmark case: $tool/$layout" >&2
      exit 1
      ;;
  esac

  local maker_cmd_q query_cmd_q
  printf -v maker_cmd_q '%q ' "${maker_cmd[@]}"
  printf -v query_cmd_q '%q ' "${query_cmd[@]}"
  echo "[bench] phase=index impl=$impl repeat=$rep/$REPEATS hash_num=$HASH_NUM"
  /usr/bin/time -f "%e\t%U\t%S\t%M" -o "$maker_time" bash -lc "set -e; $maker_cmd_q > $(printf '%q' "$maker_log") 2> $(printf '%q' "$maker_err")"
  echo "[bench] phase=query impl=$impl repeat=$rep/$REPEATS hash_num=$HASH_NUM"
  if [[ "$tool" == "cpp" ]]; then
    /usr/bin/time -f "%e\t%U\t%S\t%M" -o "$query_time" bash -lc "set -e; cd $(printf '%q' "$out_dir"); $query_cmd_q > $(printf '%q' "$query_log") 2> $(printf '%q' "$query_err")"
  else
    /usr/bin/time -f "%e\t%U\t%S\t%M" -o "$query_time" bash -lc "set -e; $query_cmd_q > $(printf '%q' "$query_log") 2> $(printf '%q' "$query_err")"
  fi

  local maker_wall maker_user maker_sys maker_rss_kb maker_cpu maker_rss_mb
  local query_wall query_user query_sys query_rss_kb query_cpu query_rss_mb
  maker_wall="$(cut -f1 "$maker_time")"
  maker_user="$(cut -f2 "$maker_time")"
  maker_sys="$(cut -f3 "$maker_time")"
  maker_rss_kb="$(cut -f4 "$maker_time")"
  query_wall="$(cut -f1 "$query_time")"
  query_user="$(cut -f2 "$query_time")"
  query_sys="$(cut -f3 "$query_time")"
  query_rss_kb="$(cut -f4 "$query_time")"
  maker_cpu="$(sum_cpu_seconds "$maker_user" "$maker_sys")"
  query_cpu="$(sum_cpu_seconds "$query_user" "$query_sys")"
  maker_rss_mb="$(rss_kb_to_mb "$maker_rss_kb")"
  query_rss_mb="$(rss_kb_to_mb "$query_rss_kb")"

  local result_line
  result_line="$(printf "%s,%s,%s,%s,%s,%s,%s,%s,%s" \
    "$impl" \
    "$rep" \
    "$HASH_NUM" \
    "$maker_wall" \
    "$maker_cpu" \
    "$maker_rss_mb" \
    "$query_wall" \
    "$query_cpu" \
    "$query_rss_mb")"
  echo "$result_line" | tee -a "$RESULTS"

  cleanup_serialized_filters "$out_dir"
}

require_file "$INDEX_DATASET"
require_file "$QUERY_DATASET"

if ! is_power_of_two "$BIT_LEN"; then
  echo "BIT_LEN must be a power of two, got $BIT_LEN" >&2
  exit 1
fi
if [[ -n "$CPP_NUM_ELE" ]]; then
  if [[ ! "$CPP_NUM_ELE" =~ ^[0-9]+$ ]]; then
    echo "CPP_NUM_ELE must be a positive integer, got $CPP_NUM_ELE" >&2
    exit 1
  fi
  if (( CPP_NUM_ELE <= 0 )); then
    echo "CPP_NUM_ELE must be > 0, got $CPP_NUM_ELE" >&2
    exit 1
  fi
fi
if (( KMER_SIZE <= 0 )); then
  echo "KMER_SIZE must be > 0, got $KMER_SIZE" >&2
  exit 1
fi
if (( RUN_RUST_BLOOMY == 1 && (KMER_SIZE % 2) == 0 )); then
  echo "RUN_RUST_BLOOMY=1 requires an odd KMER_SIZE, got $KMER_SIZE" >&2
  exit 1
fi
if (( (RUN_CPP == 1 || RUN_RUST_CLASSIC == 1 || RUN_RUST_BLOCKED == 1) && KMER_SIZE > 32 )); then
  echo "KMER_SIZE > 32 is only supported by bloomybloom in this repo; disable classic/blocked/cpp or lower KMER_SIZE" >&2
  exit 1
fi

require_bin "$SBB_BIN" "sbb"
if (( RUN_CPP == 1 )); then
  require_bin "$CPP_MAKER" "cpp maker"
  require_bin "$CPP_CAT" "cpp categorizer"
  if [[ -z "$CPP_NUM_ELE" ]]; then
    echo "CPP_NUM_ELE must be set when RUN_CPP=1 so the script can derive BioBloomMaker FPR from BIT_LEN/HASH_NUM/CPP_NUM_ELE" >&2
    exit 1
  fi
fi

echo "[bench] building release sbb binary"
cargo build --release --bin sbb >/tmp/sbb_bench_fixed_build.log

STAMP="$(date +%Y%m%d_%H%M%S)"
RESULTS="$RESULT_DIR/fixed_geometry_${STAMP}.csv"
echo "impl,repeat,hash_num,maker_wall_sec,maker_cpu_sec,maker_max_mem_mb,query_wall_sec,query_cpu_sec,query_max_mem_mb" | tee "$RESULTS"

HASH_VALUES=()
if [[ -n "$HASH_LIST" ]]; then
  IFS=',' read -r -a raw_hash_values <<< "$HASH_LIST"
  for value in "${raw_hash_values[@]}"; do
    value="${value//[[:space:]]/}"
    [[ -n "$value" ]] || continue
    if [[ ! "$value" =~ ^[0-9]+$ ]]; then
      echo "HASH_LIST must contain positive integers, got: $value" >&2
      exit 1
    fi
    if (( value <= 0 )); then
      echo "HASH_LIST values must be > 0, got: $value" >&2
      exit 1
    fi
    HASH_VALUES+=("$value")
  done
  if (( ${#HASH_VALUES[@]} == 0 )); then
    echo "HASH_LIST did not contain any valid values" >&2
    exit 1
  fi
else
  if [[ ! "$HASH_NUM" =~ ^[0-9]+$ ]] || (( HASH_NUM <= 0 )); then
    echo "HASH_NUM must be a positive integer, got $HASH_NUM" >&2
    exit 1
  fi
  HASH_VALUES=("$HASH_NUM")
fi

for hash_num in "${HASH_VALUES[@]}"; do
  HASH_NUM="$hash_num"
  echo "[bench] hash_num=$HASH_NUM"
  for rep in $(seq 1 "$REPEATS"); do
    if (( RUN_CPP == 1 )); then
      run_case cpp classic "$rep"
    fi
    if (( RUN_RUST_CLASSIC == 1 )); then
      run_case rust classic "$rep"
    fi
    if (( RUN_RUST_BLOCKED == 1 )); then
      run_case rust blocked "$rep"
    fi
    if (( RUN_RUST_BLOOMY == 1 )); then
      run_case rust bloomybloom "$rep"
    fi
  done
done

echo
echo "results_file=$RESULTS"
