#!/bin/bash
# =============================================================================
# md_fresean.sh — apo-1AKI FRESEAN production: 2 replicas x 20 ns.
# Submit FROM THE REPO ROOT:  sbatch md/1AKI/apo/md_fresean.sh
# Reuses the equilibrated box (out/npt.gro) with FRESH Maxwell velocities per
# replica (per-replica gen_seed). Direct production; first 1 ns discarded in
# analysis. Protocol matches the authors' sample-NPT.mdp (Fast Sampling SI):
# 2 fs, h-bonds, Nose-Hoover tau_t=1.0, Parrinello-Rahman tau_p=2.0,
# coords+velocities to .trr every 20 fs.
#
# DATA HANDLING (v2 — velocity-preserving): the full-system .trr is ~900 GB/rep.
# FRESEAN needs PROTEIN VELOCITIES, so we extract a protein-only .trr with a
# PLAIN group selection (NO -pbc/-fit — those DROP velocities; that bug wasted
# the first run). We then VERIFY 'gmx check' reports Velocities in the reduced
# file BEFORE deleting the full .trr. PBC-whole is a coordinate cosmetic that
# FRESEAN's own 03-CG step applies later; it does NOT belong in this reduction.
# =============================================================================
#SBATCH -J 1AKI_fresean
#SBATCH -p gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=1-2
#SBATCH -o md/1AKI/apo/logs/slurm_fresean_%A_%a.out
#SBATCH -e md/1AKI/apo/logs/slurm_fresean_%A_%a.err

set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
RUNDIR="$ROOT/md/1AKI/apo"
CFG="$RUNDIR/configs_fresean"
SRC="$RUNDIR/out"
REP="$RUNDIR/out_fresean/rep${SLURM_ARRAY_TASK_ID}"

[[ -f "$CFG/md_fresean.mdp" ]] || { echo "ERROR: $CFG/md_fresean.mdp missing"; exit 1; }
[[ -f "$SRC/npt.gro" ]]        || { echo "ERROR: $SRC/npt.gro missing"; exit 1; }
mkdir -p "$REP" "$RUNDIR/logs"
cd "$REP"
# clean any velocity-less artifacts from the first (buggy) run
rm -f md_fresean_prot.trr md_fresean.trr

echo "=== apo-1AKI FRESEAN rep ${SLURM_ARRAY_TASK_ID} | $(date) | node $(hostname) | job ${SLURM_JOB_ID} ==="
module load Gromacs
gmx --version | grep -i "version" || true

SEED=$(( 20260701 + SLURM_ARRAY_TASK_ID ))
sed "s/^gen_seed.*/gen_seed                = $SEED/" "$CFG/md_fresean.mdp" > md_fresean_rep.mdp
echo "gen_seed = $SEED"; grep '^gen_seed\|^gen_vel' md_fresean_rep.mdp

gmx grompp -f md_fresean_rep.mdp -c "$SRC/npt.gro" -p "$SRC/topol.top" -o md_fresean.tpr -maxwarn 1 -nobackup
gmx mdrun -v -deffnm md_fresean -ntmpi 1 -ntomp 8 -gpu_id 0 -pin on -nobackup

# --- crash check ---
if [[ ! -f md_fresean.gro ]]; then
  echo "ERROR: md_fresean.gro absent — run did not complete; KEEPING full .trr"; exit 1
fi
echo "=== production complete $(date); perf: $(grep -A1 'Performance' md_fresean.log | tail -1) ==="

# --- reduce: protein-only .trr, PLAIN selection (velocities preserved) ---
echo "Protein" | gmx trjconv -f md_fresean.trr -s md_fresean.tpr -o md_fresean_prot.trr -nobackup

# --- VERIFY velocities present BEFORE deleting the full .trr ---
gmx check -f md_fresean_prot.trr 2>&1 | tee gmxcheck_prot.out | grep -iE 'Velocities|Coords|Last frame' || true
if [[ -s md_fresean_prot.trr ]] && grep -qi 'Velocities' gmxcheck_prot.out; then
  FULL=$(du -h md_fresean.trr | cut -f1); PROT=$(du -h md_fresean_prot.trr | cut -f1)
  echo "VERIFIED velocities in protein-only .trr ($PROT); full was $FULL — deleting full .trr"
  rm -f md_fresean.trr
else
  echo "ERROR: protein-only .trr empty OR missing velocities — KEEPING full .trr for recovery"; exit 1
fi
echo "=== rep ${SLURM_ARRAY_TASK_ID} done $(date). Kept: md_fresean_prot.trr (w/ velocities), .tpr, .gro, .log ==="
