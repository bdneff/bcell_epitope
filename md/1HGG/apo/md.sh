#!/bin/bash
# md.sh — apo production MD (100 ns) for apo Influenza H3 hemagglutinin protomer (X-31; 1HGG chains A=HA1,B=HA2). Submit FROM THE REPO ROOT: sbatch md/1HGG/apo/md.sh
#SBATCH -J HA_apo_prod
#SBATCH -p gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH -o md/1HGG/apo/logs/slurm_prod_%j.out
#SBATCH -e md/1HGG/apo/logs/slurm_prod_%j.err
set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR:-$PWD}"; RUNDIR="$ROOT/md/1HGG/apo"; CFG="$RUNDIR/configs"; OUT="$RUNDIR/out"
[[ -f "$OUT/npt.gro" ]] || { echo "ERROR: $OUT/npt.gro not found — run prep.sh first"; exit 1; }
cd "$OUT"
echo "=== apo Influenza H3 hemagglutinin protomer (X-31; 1HGG chains A=HA1,B=HA2) production | $(date) | node $(hostname) | job ${SLURM_JOB_ID} ==="
module load Gromacs
gmx grompp -f "$CFG/md.mdp" -c npt.gro -t npt.cpt -p topol.top -o md.tpr -nobackup
gmx mdrun -v -deffnm md -ntmpi 1 -ntomp 8 -gpu_id 0 -pin on -nobackup
echo "=== apo Influenza H3 hemagglutinin protomer (X-31; 1HGG chains A=HA1,B=HA2) production finished $(date) — record ns/day, pull md.xtc ==="
