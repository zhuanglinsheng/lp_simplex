# Release procedure

Keep reusable version, dependency, packaging and verification tools on `main`.
The dependency installer supports apt-based Linux, Homebrew macOS and pinned
vcpkg on Windows. Normal CMake configuration never installs system packages.
Linux/macOS packages require KLU, SPQR and GMP and use host-provided numerical
shared libraries. Windows validation defaults to BLAS/LAPACK only, with all
optional backends OFF; numerical DLLs required by the library are bundled.
Windows Release/Debug tests and relocated C/C++ ZIP consumers must pass before
Windows support is accepted. Optional Windows backends require separate validation.

Prepare the permanent `1.0.0` branch with CMake version `1.0.0`, matching release
notes in `CHANGELOG.md`, and installation documentation. Record acceptance
results before creating the annotated `v1.0.0` tag on the final branch commit.
Never move a published tag or replace its assets. Release fixes use a new
version branch and tag; `main` continues development independently.

Run source gates and the existing Release/Debug and backend-on/off CI matrix.
For a local native Release rehearsal, with Python 3.12+:

```sh
python3 data/netlib/prepare.py --all
sh tools/build.sh build-release
python3 tools/package_release.py --build-dir build-release --platform macos-arm64
python3 test/release/verify_package.py build-dist/lp_simplex-1.0.0-macos-arm64.tar.gz --version 1.0.0 --require-full
```

Select the actual native platform. Build directories default to three parallel
jobs; set `CMAKE_BUILD_PARALLEL_LEVEL` to override. Packaging refuses Debug
builds and an existing staging directory. Remove only the previous generated
staging directory before repeating a package build. BUILD-INFO.json records
whether the checkout was dirty; official archives must come from the clean tag.

The release workflow supports manual rehearsals without publishing. Tag pushes
validate the version, run source gates, build/test Linux x64, macOS x64/arm64 and native Windows x64,
verify the actual archives after extraction and relocation, then publish the
native TAR/ZIP archives, an explicit frozen source archive and SHA256SUMS.
Every platform must pass before publication. The tag gate rejects Unreleased
notes and mismatched CMake/vcpkg versions. Versioned
CMake discovery, loaded-library/header version agreement, AUTO/Dual/Tableau/BDA
and progress callbacks are exercised by the installed C/C++ consumers.

The existing installer relocation regression also checks exported CMake paths.
Archive checks run on the build runner with installed numerical libraries; they
do not prove compatibility with arbitrary distributions or dependency versions.
System dependency versions are not locked by the package installer, and builds
are not claimed to be bit-for-bit reproducible. Each version's released assets
and checksums are the frozen binary reference.

Do not infer a new performance result from release regression checks. Existing
benchmark tables keep their original version and measurement scope; no new
full performance sweep is required for 1.0.0.
