#!/bin/bash
# =============================================================================
# postprocess.sh — PBC-correct + fit the apo-HEL (1AKI) production trajectory
# for visualization/analysis. Submit FROM THE REPO ROOT:
#   sbatch md/1AKI/apo/postprocess.sh
#
# Pipeline (canonical GROMACS PBC recipe; -fit cannot be combined with -pbc,
# so it is a separate final pass):
#   1. -pbc whole    make molecules whole across periodic boundaries
#   2. -pbc nojump    stop atoms jumping box edges between frames
#   3. -pbc mol -center -ur compact    center protein, compact box
#   4. -fit rot+trans    remove overall rotation + translation (fit on Backbone)
# Final outputs (in md/1AKI/apo/out/):
#   md_fit.xtc   protein-only, PBC-clean, rotation/translation removed  (VMD trajectory)
#   md_ref.pdb   first frame of md_fit.xtc                              (VMD topology)
#   md_center.xtc  full-system PBC-clean, NOT fitted (kept for analysis; rm if unwanted)
# =============================================================================
#SBATCH -J 1AKI_apo_post
#SBATCH -p gpu-a100              # gmx requires a GPU node on Gemini (per env notes); trjconv
#SBATCH --gres=gpu:1            # itself is CPU-only — gpu-v100 works too if you prefer.
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH -o md/1AKI/apo/logs/slurm_post_%j.out
#SBATCH -e md/1AKI/apo/logs/slurm_post_%j.err

set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR:-$PWD}"     # repo root — submit from here
OUT="$ROOT/md/1AKI/apo/out"
cd "$OUT"

# --- group selections (default tpr groups; no .ndx needed) -------------------
CENTER=Protein      # group to center the box on
FITGRP=Backbone     # group to least-squares fit on (removes tumbling)
OUTGRP=Protein      # atoms written to the final fitted trajectory (no water)

TPR=md.tpr
TRAJ=md.xtc
[[ -f "$TPR"  ]] || { echo "ERROR: $OUT/$TPR not found (run md.sh first)"; exit 1; }
[[ -f "$TRAJ" ]] || { echo "ERROR: $OUT/$TRAJ not found (run md.sh first)"; exit 1; }

echo "=== 1AKI apo postprocess | $(date) | node $(hostname) | job ${SLURM_JOB_ID} ==="
module load Gromacs              # Gemini module = "Gromacs" (capital G); GROMACS 2023.2-dev
gmx --version | grep -i version || true

# 1. make molecules whole (output full System)
printf "System\n"        | gmx trjconv -s "$TPR" -f "$TRAJ"          -o md_whole.xtc  -pbc whole -nobackup
# 2. remove jumps across box edges (full System)
printf "System\n"        | gmx trjconv -s "$TPR" -f md_whole.xtc     -o md_nojump.xtc -pbc nojump -nobackup
# 3. center on protein, compact unit cell (center group, then output group)
printf "%s\nSystem\n" "$CENTER" | gmx trjconv -s "$TPR" -f md_nojump.xtc -o md_center.xtc -pbc mol -center -ur compact -nobackup
# 4. fit rot+trans (fit group, then output group = protein only)
printf "%s\n%s\n" "$FITGRP" "$OUTGRP" | gmx trjconv -s "$TPR" -f md_center.xtc -o md_fit.xtc -fit rot+trans -nobackup

# reference structure: first frame of the fitted protein-only trajectory
printf "%s\n" "$OUTGRP"  | gmx trjconv -s "$TPR" -f md_fit.xtc -o md_ref.pdb -dump 0 -nobackup

# tidy the bulky redundant intermediates (regenerable); keep md_center + md_fit + md_ref
rm -f md_whole.xtc md_nojump.xtc

echo "=== postprocess complete $(date) ==="
echo "VMD:  load md_ref.pdb as topology, then md_fit.xtc as trajectory."
echo "Outputs: $OUT/md_fit.xtc  $OUT/md_ref.pdb  (full-system PBC-clean: $OUT/md_center.xtc)"
