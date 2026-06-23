#!/bin/bash
# =============================================================================
# postprocess.sh — PBC-correct + rot/trans-fit the apo human growth hormone (1HGU chain A) trajectory.
# Submit FROM THE REPO ROOT:  sbatch md/1HGU/apo/postprocess.sh
# Pipeline: -pbc whole -> nojump -> mol/center/compact -> fit rot+trans.
# Outputs (md/1HGU/apo/out/): md_fit.xtc (protein-only, VMD traj), md_ref.pdb, md_center.xtc.
# =============================================================================
#SBATCH -J hGH_apo_post
#SBATCH -p gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH -o md/1HGU/apo/logs/slurm_post_%j.out
#SBATCH -e md/1HGU/apo/logs/slurm_post_%j.err

set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
OUT="$ROOT/md/1HGU/apo/out"
cd "$OUT"

CENTER=Protein; FITGRP=Backbone; OUTGRP=Protein
TPR=md.tpr; TRAJ=md.xtc
[[ -f "$TPR"  ]] || { echo "ERROR: $OUT/$TPR not found (run md.sh first)"; exit 1; }
[[ -f "$TRAJ" ]] || { echo "ERROR: $OUT/$TRAJ not found (run md.sh first)"; exit 1; }

echo "=== apo human growth hormone (1HGU chain A) postprocess | $(date) | node $(hostname) | job ${SLURM_JOB_ID} ==="
module load Gromacs
gmx --version | grep -i version || true

printf "System\n"               | gmx trjconv -s "$TPR" -f "$TRAJ"        -o md_whole.xtc  -pbc whole -nobackup
printf "System\n"               | gmx trjconv -s "$TPR" -f md_whole.xtc    -o md_nojump.xtc -pbc nojump -nobackup
printf "%s\nSystem\n" "$CENTER" | gmx trjconv -s "$TPR" -f md_nojump.xtc   -o md_center.xtc -pbc mol -center -ur compact -nobackup
printf "%s\n%s\n" "$FITGRP" "$OUTGRP" | gmx trjconv -s "$TPR" -f md_center.xtc -o md_fit.xtc -fit rot+trans -nobackup
printf "%s\n" "$OUTGRP"         | gmx trjconv -s "$TPR" -f md_fit.xtc      -o md_ref.pdb -dump 0 -nobackup
rm -f md_whole.xtc md_nojump.xtc

echo "=== apo human growth hormone (1HGU chain A) postprocess complete $(date) ==="
echo "VMD: load md_ref.pdb then md_fit.xtc."
