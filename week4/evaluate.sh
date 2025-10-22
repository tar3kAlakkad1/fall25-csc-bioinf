#!/usr/bin/env bash
set -euo pipefail

# Resolve base directory (week4)
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PY_MAIN="$BASE_DIR/python/main.py"
CODON_MAIN="$BASE_DIR/codon/main.py"
CODON_ALIGN_MOD="$BASE_DIR/codon/alignment.py"
CODON_EXE="$BASE_DIR/codon/align"
DATA_DIR="$BASE_DIR/python/data"

# Ensure Python deps (numpy) are available
if ! python3 -c 'import numpy' >/dev/null 2>&1; then
  python3 -m pip install --user -q numpy >/dev/null 2>&1 || true
fi

# Ensure Codon is installed and compile the Codon CLI if needed
if ! command -v codon >/dev/null 2>&1; then
  echo "Codon compiler not found on PATH" >&2
  exit 1
fi

if [ ! -x "$CODON_EXE" ] || [ "$CODON_MAIN" -nt "$CODON_EXE" ] || [ "$CODON_ALIGN_MOD" -nt "$CODON_EXE" ]; then
  codon build -release -o "$CODON_EXE" "$CODON_MAIN"
fi

# On macOS, Codon/Numpy may require libomp.dylib at runtime. Add common Homebrew paths.
if [ "$(uname)" = "Darwin" ]; then
  # Ensure libomp is installed (Homebrew) if available
  if command -v brew >/dev/null 2>&1; then
    if ! brew list libomp >/dev/null 2>&1; then
      brew install -q libomp || true
    fi
  fi

  # Try to locate libomp from common paths and Homebrew prefix
  OMP_CANDIDATES=()
  if command -v brew >/dev/null 2>&1; then
    brew_prefix="$(brew --prefix libomp 2>/dev/null || true)"
    if [ -n "$brew_prefix" ]; then
      OMP_CANDIDATES+=("$brew_prefix/lib")
    fi
  fi
  OMP_CANDIDATES+=(
    "/opt/homebrew/opt/libomp/lib"
    "/usr/local/opt/libomp/lib"
  )

  for omp_dir in "${OMP_CANDIDATES[@]}"; do
    if [ -f "$omp_dir/libomp.dylib" ]; then
      if [ -z "${DYLD_LIBRARY_PATH+x}" ] || [ -z "${DYLD_LIBRARY_PATH}" ]; then
        export DYLD_LIBRARY_PATH="$omp_dir"
      else
        export DYLD_LIBRARY_PATH="$DYLD_LIBRARY_PATH:$omp_dir"
      fi
      # Also copy next to the executable so @loader_path/libomp.dylib resolves
      codon_dir="$(dirname "$CODON_EXE")"
      if [ -w "$codon_dir" ]; then
        cp -f "$omp_dir/libomp.dylib" "$codon_dir/libomp.dylib" || true
      fi
      break
    fi
  done
fi

# Helper: measure command wall time in milliseconds using Python's perf_counter
measure_ms() {
  python3 "$BASE_DIR/test.py" run "$@"
}

print_header() {
  printf "Method          Language    Runtime\n"
  printf '%s\n' "--------------------------------------"
}

run_case() {
  local algo="$1" dataset="$2" f1="$3" id1="$4" f2="$5" id2="$6"
  # Python
  local ms
  ms=$(measure_ms python "$algo" "$f1" "$id1" "$f2" "$id2")
  printf "%-16s %-10s %sms\n" "${algo}-${dataset}" "python" "$ms"
  # Codon
  ms=$(measure_ms codon "$algo" "$f1" "$id1" "$f2" "$id2")
  printf "%-16s %-10s %sms\n" "${algo}-${dataset}" "codon" "$ms"
}

main() {
  print_header

  local F_HUMAN="$DATA_DIR/MT-human.fa" ID_HUMAN="MT_human"
  local F_ORANG="$DATA_DIR/MT-orang.fa" ID_ORANG="MT_orang"
  local F_Q="$DATA_DIR/q1.fa"
  local F_T="$DATA_DIR/t1.fa"

  for algo in global local semi-global affine; do
    # MT_human vs MT_orang
    run_case "$algo" mt_human "$F_HUMAN" "$ID_HUMAN" "$F_ORANG" "$ID_ORANG"
    # q1..q5 vs t1..t5
    for n in 1 2 3 4 5; do
      run_case "$algo" "q$n" "$F_Q" "q$n" "$F_T" "t$n"
    done
  done
}

main "$@"


