#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="${OUT_DIR:-$ROOT_DIR/bench/refs}"
FORCE="${FORCE:-0}"
ONLY="${ONLY:-all}"

mkdir -p "$OUT_DIR"

download_one() {
  local name="$1"
  local acc="$2"
  local out="$OUT_DIR/${name}.fa"
  local url="https://www.ebi.ac.uk/ena/browser/api/fasta/${acc}?download=true"

  if [[ "$ONLY" != "all" && "$ONLY" != "$name" ]]; then
    return
  fi
  if [[ -s "$out" && "$FORCE" != "1" ]]; then
    echo "skip (exists): $out"
    return
  fi

  echo "download: $name ($acc)"
  curl -L --fail --retry 5 --retry-delay 2 -o "$out" "$url"
  if [[ ! -s "$out" ]]; then
    echo "download failed (empty): $out" >&2
    exit 1
  fi
  if ! head -n1 "$out" | grep -q '^>'; then
    echo "download failed (not FASTA): $out" >&2
    exit 1
  fi
}

download_one "ecoli_k12_mg1655" "U00096.3"
download_one "salmonella_typhimurium_lt2" "AE006468.2"
download_one "bacillus_subtilis_168" "AL009126.3"
download_one "staphylococcus_aureus_n315" "BA000018.3"

MANIFEST="$OUT_DIR/manifest.tsv"
{
  echo -e "name\taccession\tpath\tbytes"
  for f in "$OUT_DIR"/*.fa; do
    [[ -f "$f" ]] || continue
    base="$(basename "$f" .fa)"
    acc="na"
    case "$base" in
      ecoli_k12_mg1655) acc="U00096.3" ;;
      salmonella_typhimurium_lt2) acc="AE006468.2" ;;
      bacillus_subtilis_168) acc="AL009126.3" ;;
      staphylococcus_aureus_n315) acc="BA000018.3" ;;
    esac
    echo -e "${base}\t${acc}\t${f}\t$(stat -c%s "$f")"
  done
} > "$MANIFEST"

echo "refs_dir=$OUT_DIR"
echo "manifest=$MANIFEST"
