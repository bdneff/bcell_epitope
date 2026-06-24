#!/bin/bash
# postprocess.sh — PBC + rot/trans fit for apo Dengue-2 envelope ectodomain (sE; 1OKE chain A). Submit FROM REPO ROOT: sbatch md/1OKE/apo/postprocess.sh
#SBATCH -J DengueE_apo_post
#SBATCH -p gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH -o md/1OKE/apo/logs/slurm_post_%j.out
#SBATCH -e md/1OKE/apo/logs/slurm_post_%j.err
set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR:-$PWD}"; OUT="$ROOT/md/1OKE/apo/out"; cd "$OUT"
CENTER=Protein; FITGRP=Backbone; OUTGRP=Protein; TPR=md.tpr; TRAJ=md.xtc
[[ -f "$TPR" && -f "$TRAJ" ]] || { echo "ERROR: run md.sh first"; exit 1; }
module load Gromacs
printf "System\n"               | gmx trjconv -s "$TPR" -f "$TRAJ"      -o md_whole.xtc  -pbc whole -nobackup
printf "System\n"               | gmx trjconv -s "$TPR" -f md_whole.xtc  -o md_nojump.xtc -pbc nojump -nobackup
printf "%s\nSystem\n" "$CENTER" | gmx trjconv -s "$TPR" -f md_nojump.xtc -o md_center.xtc -pbc mol -center -ur compact -nobackup
printf "%s\n%s\n" "$FITGRP" "$OUTGRP" | gmx trjconv -s "$TPR" -f md_center.xtc -o md_fit.xtc -fit rot+trans -nobackup
printf "%s\n" "$OUTGRP"         | gmx trjconv -s "$TPR" -f md_fit.xtc    -o md_ref.pdb -dump 0 -nobackup
rm -f md_whole.xtc md_nojump.xtc
echo "=== apo Dengue-2 envelope ectodomain (sE; 1OKE chain A) postprocess complete $(date) ==="
