#!/bin/bash
# =============================================================================
# md.sh — apo production MD (100 ns) for VEGF. Submit FROM THE REPO ROOT:
#   sbatch md/1BJ1/apo/md.sh
# Reads prep outputs from md/1BJ1/apo/out/ (npt.gro, npt.cpt, topol.top).
# =============================================================================
#SBATCH -J VEGF_apo_prod
#SBATCH -p gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH -o md/1BJ1/apo/logs/slurm_prod_%j.out
#SBATCH -e md/1BJ1/apo/logs/slurm_prod_%j.err

set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
RUNDIR="$ROOT/md/1BJ1/apo"; CFG="$RUNDIR/configs"; OUT="$RUNDIR/out"
[[ -d "$CFG" ]]         || { echo "ERROR: configs not found at $CFG (submit from repo root)"; exit 1; }
[[ -f "$OUT/npt.gro" ]] || { echo "ERROR: $OUT/npt.gro not found — run prep.sh first"; exit 1; }
cd "$OUT"

echo "=== VEGF apo production | $(date) | node $(hostname) | job ${SLURM_JOB_ID} ==="
module load Gromacs
gmx --version | grep -i version || true

gmx grompp -f "$CFG/md.mdp" -c npt.gro -t npt.cpt -p topol.top -o md.tpr -nobackup
gmx mdrun -v -deffnm md -ntmpi 1 -ntomp 8 -gpu_id 0 -pin on -nobackup

echo "=== VEGF production finished $(date) ==="
echo "Record ns/day from md.log in the run log. Pull back md.xtc for analysis."
