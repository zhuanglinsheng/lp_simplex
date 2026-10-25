#!/bin/sh
# Install build dependencies explicitly; CMake never changes the host system.
set -eu
case "$(uname -s)" in
    Darwin)
        command -v brew >/dev/null 2>&1 || {
            echo 'Install Homebrew first: https://brew.sh' >&2
            exit 1
        }
        brew install cmake ninja suite-sparse gmp
        ;;
    Linux)
        command -v apt-get >/dev/null 2>&1 || {
            echo 'Install C/C++, CMake, BLAS, LAPACK, SuiteSparse and GMP using your package manager.' >&2
            exit 1
        }
        if [ "$(id -u)" -eq 0 ]; then
            apt-get update
            apt-get install -y build-essential cmake ninja-build libopenblas-dev liblapack-dev libsuitesparse-dev libgmp-dev
        else
            sudo apt-get update
            sudo apt-get install -y build-essential cmake ninja-build libopenblas-dev liblapack-dev libsuitesparse-dev libgmp-dev
        fi
        ;;
    *) echo 'Unsupported dependency installer platform.' >&2; exit 1 ;;
esac
