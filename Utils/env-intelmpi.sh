#!/usr/bin/env bash
# Intel oneAPI MPI + Intel-parallel HDF5 env for LePIC2.0

# Be defensive but allow setvars.sh to read unset env vars
set -eo pipefail

# ---- scrub Open MPI paths to avoid vendor mixing ----
export PATH="$(printf "%s" "${PATH}" | tr ':' '\n' | grep -v '/openmpi/' | paste -sd: -)"
export LD_LIBRARY_PATH="$(printf "%s" "${LD_LIBRARY_PATH:-}" | tr ':' '\n' | grep -v '/openmpi/' | paste -sd: -)"

SETVARS="/opt/intel/oneapi/setvars.sh"
if [[ ! -r "$SETVARS" ]]; then
  echo "[env-intelmpi] ERROR: $SETVARS not readable"; return 1 2>/dev/null || exit 1
fi

echo "[env-intelmpi] Sourcing $SETVARS with --force ..."
# Temporarily disable nounset while sourcing setvars.sh
# (setvars may reference unset vars like OCL_ICD_FILENAMES)
set +u
export SETVARS_ARGS="--force"
# shellcheck source=/opt/intel/oneapi/setvars.sh
source "$SETVARS"
set -u || true  # re-enable nounset if your shell supports it

# Ensure Intel MPI libs are early in the runtime path
I_MPI_LIB_CANDIDATES=(
  "${I_MPI_ROOT:-/opt/intel/oneapi/mpi/latest}/lib/release"
  "${I_MPI_ROOT:-/opt/intel/oneapi/mpi/latest}/lib"
)
for d in "${I_MPI_LIB_CANDIDATES[@]}"; do
  [[ -d "$d" ]] && export LD_LIBRARY_PATH="$d:${LD_LIBRARY_PATH:-}"
done

# CMake hints (adjust HDF5 path if different on your box)
export FC=mpiifx
export TYPE="${TYPE:-BASE}"
export HDF5_ROOT="${HDF5_ROOT:-/opt/hdf5-ifx}"
export HDF5_DIR="${HDF5_DIR:-$HDF5_ROOT/lib/cmake/hdf5}"

echo "[env-intelmpi] FC=$FC TYPE=$TYPE"
echo "[env-intelmpi] mpirun: $(command -v mpirun || echo 'not found')"
command -v mpiifx >/dev/null && mpiifx --version | head -1 || echo "[env-intelmpi] WARN: mpiifx not on PATH?"
