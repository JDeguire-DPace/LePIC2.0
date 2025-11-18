#!/usr/bin/env bash
rm -f Outputs/state.h5
set -euo pipefail
# Path to the executable built by CMake
BIN="build/bin/run_pic3D"

# Allow overriding mpirun (for Intel MPI etc.)
MPIRUN_BIN="${MPIRUN:-mpirun}"
NP="${NP:-4}"

echo "==> mpirun: $(command -v "${MPIRUN_BIN}" || echo "${MPIRUN_BIN} (not found)")"
echo "==> Linked (ldd):"
ldd "${BIN}" | grep -E 'hdf5|mpi' || true

# Disable HDF5 file locking in a *portable* way
export HDF5_USE_FILE_LOCKING=FALSE

echo "==> Running: ${MPIRUN_BIN} -np ${NP} ${BIN}"
exec "${MPIRUN_BIN}" -np 4 "${BIN}"