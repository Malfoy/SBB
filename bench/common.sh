#!/usr/bin/env bash

bench_init_cpp_env() {
  local cpp_bin_dir="${CPP_BIN_DIR:-$ROOT_DIR/bench/bin}"
  CPP_MAKER="${CPP_MAKER:-$cpp_bin_dir/biobloommaker}"
  CPP_CAT="${CPP_CAT:-$cpp_bin_dir/biobloomcategorizer}"

  local cpp_lib_dir="${CPP_LIB_DIR:-$ROOT_DIR/bench/lib}"
  if [[ -d "$cpp_lib_dir" ]]; then
    if [[ -n "${LD_LIBRARY_PATH:-}" ]]; then
      export LD_LIBRARY_PATH="$cpp_lib_dir:$LD_LIBRARY_PATH"
    else
      export LD_LIBRARY_PATH="$cpp_lib_dir"
    fi
  fi
}
