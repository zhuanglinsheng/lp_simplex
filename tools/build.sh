#!/bin/sh
# Full Release build by default; forward extra CMake options after BUILD_DIR.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=${1:-"$root/build"}
if [ "$#" -gt 0 ]; then shift; fi
if [ "$(uname -s)" = Darwin ] && command -v brew >/dev/null 2>&1; then
    brew_prefix=$(brew --prefix)
    CMAKE_PREFIX_PATH="$brew_prefix${CMAKE_PREFIX_PATH:+:$CMAKE_PREFIX_PATH}"
    export CMAKE_PREFIX_PATH
    set -- -DBLA_VENDOR=Apple "$@"
fi
cmake -S "$root" -B "$build_dir" -DCMAKE_BUILD_TYPE=Release \
    -DLP_SIMPLEX_USE_KLU=ON -DLP_SIMPLEX_USE_SPQR=ON -DLP_SIMPLEX_USE_GMP=ON \
    -DLP_SIMPLEX_REQUIRE_OPTIONAL_BACKENDS=ON "$@"
cmake --build "$build_dir" --parallel "${CMAKE_BUILD_PARALLEL_LEVEL:-3}"
ctest --test-dir "$build_dir" --output-on-failure
