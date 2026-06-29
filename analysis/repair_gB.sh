#!/bin/bash
# Finish the HCMV gB (5C6T chain A) repair once an AlphaFold3 model is available.
# gB is absent from AlphaFold DB (viral), so its model is generated with AF3 on Gemini:
#
#   1. On Gemini (repo root):  sbatch md/5C6T/apo/af3.sbatch
#      -> writes md/5C6T/apo/af3/out/gb_5c6t/..._model.cif
#   2. Pull md/5C6T/apo/af3/out back to this machine (the repair env is local).
#   3. Run this script from the repo root (it auto-finds the AF3 model cif, or pass it):
#         bash analysis/repair_gB.sh [path/to/model.cif]
#
# It grafts the unresolved internal loops (115-120, 219-220, 237-239, and the 32-residue
# 437-468) from the AF3 model onto the crystal, completes side chains, closes the
# junctions, and writes the frozen input structures/HCMVgB_5C6T_fixed.pdb that
# md/5C6T/apo/prep.sh reads. Inspect the printed global CA-RMSD: if AF3 disagrees with the
# crystal (>3 A) the graft is not trustworthy and we should instead simulate the AF3
# ectodomain directly (decide then). See manuscript sec:repair.
set -euo pipefail
MM="$HOME/.local/bin/micromamba run -n bcell-repair"
AF="${1:-}"
if [[ -z "$AF" ]]; then
  AF="$(find md/5C6T/apo/af3/out structures/repair -name '*_model.cif' -o -name 'AF3-gB_5C6T.cif' 2>/dev/null | head -1)"
fi
[[ -n "$AF" && -f "$AF" ]] || { echo "ERROR: AF3 model cif not found -- run md/5C6T/apo/af3.sbatch on Gemini and pull it back (see header)"; exit 1; }
echo "using AF3 model: $AF"

$MM python analysis/repair_structure.py --crystal structures/HCMVgB_5C6T.pdb --chain A \
    --af "$AF" --out structures/repair/gB_5C6T_grafted.pdb --max-gap 40
$MM python analysis/finalize_structure.py structures/repair/gB_5C6T_grafted.pdb \
    structures/repair/gB_5C6T_pre.pdb
# free a window around each grafted gap for closure
$MM python analysis/close_loops.py structures/repair/gB_5C6T_pre.pdb \
    structures/repair/gB_5C6T_closed.pdb 114-121,218-221,236-240,436-469
$MM python - structures/repair/gB_5C6T_closed.pdb structures/HCMVgB_5C6T_fixed.pdb <<'PY'
import sys
from Bio.PDB import PDBParser, PDBIO, Select
class NoH(Select):
    def accept_atom(self, a): return a.element != 'H'
io = PDBIO(); io.set_structure(PDBParser(QUIET=True).get_structure('x', sys.argv[1]))
io.save(sys.argv[2], NoH()); print("wrote", sys.argv[2])
PY
echo "Done -> structures/HCMVgB_5C6T_fixed.pdb. Validate, then sbatch md/5C6T/apo/prep.sh on Gemini."
