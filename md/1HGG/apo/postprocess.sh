#!/bin/bash
# postprocess.sh — PBC + rot/trans fit for apo Influenza H3 hemagglutinin protomer (X-31; 1HGG chains A=HA1,B=HA2). Submit FROM REPO ROOT: sbatch md/1HGG/apo/postprocess.sh
#SBATCH -J HA_apo_post
#SBATCH -p gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH -o md/1HGG/apo/logs/slurm_post_%j.out
#SBATCH -e md/1HGG/apo/logs/slurm_post_%j.err
set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR:-$PWD}"; OUT="$ROOT/md/1HGG/apo/out"; cd "$OUT"
CENTER=Protein; FITGRP=Backbone; OUTGRP=Protein; TPR=md.tpr; TRAJ=md.xtc
[[ -f "$TPR" && -f "$TRAJ" ]] || { echo "ERROR: run md.sh first"; exit 1; }
module load Gromacs
printf "System\n"               | gmx trjconv -s "$TPR" -f "$TRAJ"      -o md_whole.xtc  -pbc whole -nobackup
printf "System\n"               | gmx trjconv -s "$TPR" -f md_whole.xtc  -o md_nojump.xtc -pbc nojump -nobackup
printf "%s\nSystem\n" "$CENTER" | gmx trjconv -s "$TPR" -f md_nojump.xtc -o md_center.xtc -pbc mol -center -ur compact -nobackup
printf "%s\n%s\n" "$FITGRP" "$OUTGRP" | gmx trjconv -s "$TPR" -f md_center.xtc -o md_fit.xtc -fit rot+trans -nobackup
printf "%s\n" "$OUTGRP"         | gmx trjconv -s "$TPR" -f md_fit.xtc    -o md_ref.pdb -dump 0 -nobackup
rm -f md_whole.xtc md_nojump.xtc
echo "=== apo Influenza H3 hemagglutinin protomer (X-31; 1HGG chains A=HA1,B=HA2) postprocess complete $(date) ==="
