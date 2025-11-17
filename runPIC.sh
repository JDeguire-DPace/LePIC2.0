#!/usr/bin/env bash
rm -f Outputs/state.h5
set -euo pipefail
NP="${1:-4}"

python3 Utils/SetupGUI.py

BIN="build/bin/run_pic3D"

if [[ ! -x "$BIN" ]]; then
  echo "Error: $BIN not found. Build first: ./config.sh" >&2
  exit 1
fi

echo "==> mpirun: $(command -v mpirun || echo 'not found')"
echo "==> Linked (ldd):"
ldd "$BIN" | grep -E 'libmpi|hdf5' || true

exec mpirun -np 4 -genv HDF5_USE_FILE_LOCKING=FALSE "$BIN"
