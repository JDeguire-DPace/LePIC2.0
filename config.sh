#!/usr/bin/env bash
set -euo pipefail

# Expect env files to set these; provide safe fallbacks
: "${FC:=mpiifx}"
: "${TYPE:=BASE}"
: "${HDF5_ROOT:=/opt/hdf5-ifx}"
: "${HDF5_DIR:=/opt/hdf5-ifx/lib/cmake/hdf5}"

echo "==> Cleaning build and stale .mod"
find Src -type f \( -name '*.mod' -o -name '*.smod' \) -delete || true
rm -rf build

echo "==> Configuring (FC=$FC, TYPE=$TYPE)"
cmake -S . -B build \
  -DCMAKE_Fortran_COMPILER="${FC}" \
  -DCMAKE_BUILD_TYPE="${TYPE}" \
  -DHDF5_ROOT="${HDF5_ROOT}" \
  -DHDF5_DIR="${HDF5_DIR}"

echo "==> Building"
cmake --build build -j

echo "==> Linked libs (sanity)"
ldd build/bin/run_pic3D | grep -E 'libmpi|hdf5' || true

echo "==> Done -> build/bin/run_pic3D"
