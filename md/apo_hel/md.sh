#!/bin/bash
# =============================================================================
# md.sh — apo-HEL production MD (100 ns). Submit with: sbatch md.sh
# Run from the same directory as prep.sh outputs (npt.gro, npt.cpt, topol.top).
# =============================================================================
#SBATCH -J hel_prod
#SBATCH -p gpu-a100              # CONFIRM partition on Gemini
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH -o slurm_md_%j.out
#SBATCH -e slurm_md_%j.err

set -euo pipefail
CFG="$(dirname "$0")/configs"

echo "=== apo-HEL production | $(date) | node $(hostname) | job ${SLURM_JOB_ID} ==="
module load gromacs              # CONFIRM module name/version on Gemini
gmx --version | grep -i version || true

gmx grompp -f "$CFG/md.mdp" -c npt.gro -t npt.cpt -p topol.top -o md.tpr -nobackup

gmx mdrun -v -deffnm md -ntmpi 1 -ntomp 8 -gpu_id 0 -pin on -nobackup

echo "=== production finished $(date) ==="
echo "Record ns/day from md.log in the manuscript run-log. Pull back md.xtc for analysis."
