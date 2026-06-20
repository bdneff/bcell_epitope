#!/bin/bash
# =============================================================================
# md.sh — apo-HEL (1AKI) production MD (100 ns). Submit FROM THE REPO ROOT:
#   sbatch md/1AKI/apo/md.sh
# Reads prep outputs from md/1AKI/apo/out/ (npt.gro, npt.cpt, topol.top).
# Outputs -> md/1AKI/apo/out/ (gitignored); slurm logs -> md/1AKI/apo/logs/.
# =============================================================================
#SBATCH -J 1AKI_apo_prod
#SBATCH -p gpu-a100              # confirmed on Gemini 2026-06-18 (4-day limit, A100 fastest)
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00          # 100 ns; generous ceiling (expect ~10-20 h on one A100)
#SBATCH -o md/1AKI/apo/logs/slurm_prod_%j.out
#SBATCH -e md/1AKI/apo/logs/slurm_prod_%j.err

set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR:-$PWD}"     # repo root — submit from here
RUNDIR="$ROOT/md/1AKI/apo"
CFG="$RUNDIR/configs"
OUT="$RUNDIR/out"

[[ -d "$CFG" ]]         || { echo "ERROR: configs not found at $CFG (submit from repo root)"; exit 1; }
[[ -f "$OUT/npt.gro" ]] || { echo "ERROR: $OUT/npt.gro not found — run prep.sh first"; exit 1; }
cd "$OUT"

echo "=== apo-HEL (1AKI) production | $(date) | node $(hostname) | job ${SLURM_JOB_ID} ==="
module load Gromacs              # Gemini module = "Gromacs" (capital G); GROMACS 2023.2-dev
gmx --version | grep -i version || true

gmx grompp -f "$CFG/md.mdp" -c npt.gro -t npt.cpt -p topol.top -o md.tpr -nobackup

gmx mdrun -v -deffnm md -ntmpi 1 -ntomp 8 -gpu_id 0 -pin on -nobackup

echo "=== production finished $(date) ==="
echo "Record ns/day from md.log in the run log. Pull back md.xtc for analysis."
