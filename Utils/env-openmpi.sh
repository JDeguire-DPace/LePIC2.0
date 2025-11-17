#!/usr/bin/env bash
# Use system OpenMPI + distro HDF5(openmpi)

# 1) Remove Intel MPI from this shell, if present
export PATH="$(printf "%s" "$PATH" | tr ':' '\n' | grep -v '/oneapi/mpi/' | paste -sd: -)"
export LD_LIBRARY_PATH="$(printf "%s" "${LD_LIBRARY_PATH:-}" | tr ':' '\n' | grep -v '/oneapi/mpi/' | paste -sd: -)"

# 2) Ensure OpenMPI is first
export PATH="/usr/bin:${PATH}"
export LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/openmpi/lib:/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH:-}"

# 3) Tell CMake what to use (this matches your current config)
export FC=mpifort
export TYPE="${TYPE:-Debug}"
export HDF5_ROOT="/usr/lib/x86_64-linux-gnu/hdf5/openmpi"
export HDF5_DIR="/usr/lib/x86_64-linux-gnu/cmake/hdf5"

echo "[env-openmpi] FC=$FC  TYPE=$TYPE"
command -v mpifort >/dev/null && mpifort --version | head -1 || true
