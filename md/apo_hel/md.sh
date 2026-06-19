#!/bin/bash
# =============================================================================
# md.sh — apo-HEL production MD (100 ns). Submit with: sbatch md.sh
# Run from the same directory as prep.sh outputs (npt.gro, npt.cpt, topol.top).
# =============================================================================
#SBATCH -J hel_prod
#SBATCH -p gpu-a100              # confirmed on Gemini 2026-06-18 (4-day limit, A100 fastest)
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00          # 100 ns; generous ceiling (expect ~10-20 h on one A100)
#SBATCH -o slurm_md_%j.out
#SBATCH -e slurm_md_%j.err

set -euo pipefail
# Slurm copies the batch script to a private spool dir, so $0/dirname won't find configs/.
# Anchor to the submit dir (repo root when you run `sbatch md/apo_hel/md.sh`).
CFG="${SLURM_SUBMIT_DIR:-$PWD}/md/apo_hel/configs"
[[ -d "$CFG" ]] || { echo "ERROR: configs dir not found at $CFG (submit from repo root)"; exit 1; }

echo "=== apo-HEL production | $(date) | node $(hostname) | job ${SLURM_JOB_ID} ==="
module load Gromacs              # Gemini module = "Gromacs" (capital G); GROMACS 2023.2-dev
gmx --version | grep -i version || true

gmx grompp -f "$CFG/md.mdp" -c npt.gro -t npt.cpt -p topol.top -o md.tpr -nobackup

gmx mdrun -v -deffnm md -ntmpi 1 -ntomp 8 -gpu_id 0 -pin on -nobackup

echo "=== production finished $(date) ==="
echo "Record ns/day from md.log in the manuscript run-log. Pull back md.xtc for analysis."
