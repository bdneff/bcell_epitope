#!/bin/bash
# =============================================================================
# md_fresean.sh — apo-2JEL FRESEAN production: 2 replicas x 20 ns.
# Submit FROM THE REPO ROOT:  sbatch md/2JEL/apo/md_fresean.sh
# Reuses the equilibrated box (out/npt.gro) with FRESH Maxwell velocities per
# replica via a per-replica gen_seed. Direct production; first 1 ns discarded in
# analysis (user decision 2026-07-01). Protocol matches the authors' own
# sample-NPT.mdp (Fast Sampling SI): 2 fs, h-bonds, Nose-Hoover tau_t=1.0,
# Parrinello-Rahman tau_p=2.0, coords+velocities to .trr every 20 fs.
#
# DATA HANDLING: the full-system .trr is ~900 GB/replica (1e6 frames x all atoms).
# The FRESEAN mode analysis only uses PROTEIN velocities (2 beads/residue CG), so
# after mdrun we trjconv out a protein-only .trr (velocities preserved, ~40x
# smaller) and DELETE the full .trr — mirrors the authors' 03-CG + clean.sh.
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

# Reproducible per-replica velocity seed
SEED=$(( 20260701 + SLURM_ARRAY_TASK_ID ))
sed "s/^gen_seed.*/gen_seed                = $SEED/" "$CFG/md_fresean.mdp" > md_fresean_rep.mdp
echo "gen_seed = $SEED"; grep '^gen_seed\|^gen_vel' md_fresean_rep.mdp

# grompp: -maxwarn 1 for the benign Parrinello-Rahman+gen_vel note (box is already
# NPT-equilibrated; only velocities are fresh; 1 ns discard covers the startup transient)
gmx grompp -f md_fresean_rep.mdp -c "$SRC/npt.gro" -p "$SRC/topol.top" -o md_fresean.tpr -maxwarn 1 -nobackup

gmx mdrun -v -deffnm md_fresean -ntmpi 1 -ntomp 8 -gpu_id 0 -pin on -nobackup

# --- crash check: did we reach 20 ns? ---
if [[ ! -f md_fresean.gro ]]; then
  echo "ERROR: md_fresean.gro absent — run did not complete; KEEPING full .trr for diagnosis"; exit 1
fi
echo "=== production complete $(date); ns/day: $(grep -A1 'Performance' md_fresean.log | tail -1) ==="

# --- reduce: protein-only .trr (velocities preserved), then delete the full .trr ---
echo "Protein" | gmx trjconv -f md_fresean.trr -s md_fresean.tpr -o md_fresean_prot.trr -pbc mol -nobackup
if [[ -s md_fresean_prot.trr ]]; then
  FULL=$(du -h md_fresean.trr | cut -f1); PROT=$(du -h md_fresean_prot.trr | cut -f1)
  echo "reduced full .trr ($FULL) -> protein-only ($PROT); deleting full .trr"
  rm -f md_fresean.trr
else
  echo "ERROR: protein-only .trr empty — KEEPING full .trr"; exit 1
fi
echo "=== rep ${SLURM_ARRAY_TASK_ID} done $(date). Kept: md_fresean_prot.trr, .tpr, .gro, .log ==="
