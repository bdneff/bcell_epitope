#!/bin/bash
# =============================================================================
# md_fresean.sh — apo-2JEL FRESEAN production: 5 replicas x 20 ns.
# Submit FROM THE REPO ROOT:  sbatch md/2JEL/apo/md_fresean.sh
# Reuses the equilibrated box (out/npt.gro) with FRESH Maxwell velocities per
# replica via a per-replica gen_seed. Direct production; first 1 ns discarded in
# analysis (user decision 2026-07-01). Protocol: Fast Sampling SI (see
# docs/ENVIRONMENT_fresean.md). Velocities every 20 fs -> .trr for FRESEAN.
# Outputs -> md/2JEL/apo/out_fresean/rep${SLURM_ARRAY_TASK_ID}/ (gitignored).
# =============================================================================
#SBATCH -J 2JEL_fresean
#SBATCH -p gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=1-2
#SBATCH -o md/2JEL/apo/logs/slurm_fresean_%A_%a.out
#SBATCH -e md/2JEL/apo/logs/slurm_fresean_%A_%a.err

set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
RUNDIR="$ROOT/md/2JEL/apo"
CFG="$RUNDIR/configs_fresean"
SRC="$RUNDIR/out"                      # equilibrated box from the 100 ns prep
REP="$RUNDIR/out_fresean/rep${SLURM_ARRAY_TASK_ID}"

[[ -f "$CFG/md_fresean.mdp" ]] || { echo "ERROR: $CFG/md_fresean.mdp missing (submit from repo root)"; exit 1; }
[[ -f "$SRC/npt.gro" ]]        || { echo "ERROR: $SRC/npt.gro missing — run prep first"; exit 1; }
mkdir -p "$REP" "$RUNDIR/logs"
cd "$REP"

echo "=== apo-2JEL FRESEAN rep ${SLURM_ARRAY_TASK_ID} | $(date) | node $(hostname) | job ${SLURM_JOB_ID} ==="
module load Gromacs
gmx --version | grep -i "version" || true

# Reproducible per-replica velocity seed: fixed base + task id, written into a
# per-replica mdp copy (gen_seed is the canonical GROMACS mechanism; grompp has no -seed flag).
SEED=$(( 20260701 + SLURM_ARRAY_TASK_ID ))
sed "s/^gen_seed.*/gen_seed                = $SEED/" "$CFG/md_fresean.mdp" > md_fresean_rep.mdp
echo "gen_seed = $SEED"
grep '^gen_seed\|^gen_vel' md_fresean_rep.mdp

# Fresh velocities from equilibrated coordinates (no -t checkpoint; continuation=no, gen_vel=yes)
gmx grompp -f md_fresean_rep.mdp -c "$SRC/npt.gro" -p "$SRC/topol.top" -o md_fresean.tpr -nobackup

gmx mdrun -v -deffnm md_fresean -ntmpi 1 -ntomp 8 -gpu_id 0 -pin on -nobackup

echo "=== rep ${SLURM_ARRAY_TASK_ID} finished $(date) ==="
echo "Record ns/day (md_fresean.log). .trr (coords+vel @20fs) stays on scratch."
