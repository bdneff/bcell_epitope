#!/bin/bash
# =============================================================================
# md_fresean_chunked.sh — apo-1OKE (Dengue E) FRESEAN production, 1 replica,
# CHUNKED into 1 ns pieces to cap the transient full-system .trr footprint.
#
# WHY CHUNKED: 1OKE is 347,139 atoms. A monolithic 20 ns run writing coords+
# velocities every 20 fs = ~8 TB transient full-system .trr. Running 1 ns at a
# time and reducing-then-deleting each chunk caps peak scratch at ~one chunk
# (~416 GB full) + the accumulating protein-only .trr (~7.4 GB/chunk, 6129
# protein atoms; ~147 GB for 20 ns). Fits /scratch (105 TB free) with margin.
#
# RESUMABLE: a state file (fresean_state.txt) records the last COMPLETED chunk.
# Re-submitting the job resumes from there. Each chunk continues from the
# previous chunk's checkpoint (grompp -t prev.cpt, continuation=yes).
#
# VELOCITY-PRESERVING REDUCTION (the bug that wasted an earlier run): extract
# protein-only with a PLAIN group selection (NO -pbc/-fit — those DROP
# velocities). Verify 'gmx check' reports Velocities BEFORE deleting the full
# chunk .trr. PBC-whole is applied later by FRESEAN's own 03-CG step.
#
# FORCE FIELD: AMBER99SB-ILDN / TIP3P, inherited from out/topol.top (matches
# 1AKI and every other benchmark system — required for mode comparability).
#
# Submit FROM THE REPO ROOT:  sbatch md/1OKE/apo/md_fresean_chunked.sh
# Built with Claude Science (session #4, 2026-07-07).
# =============================================================================
#SBATCH -J 1OKE_fresean
#SBATCH -p gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=2-00:00:00
#SBATCH -o md/1OKE/apo/logs/slurm_fresean_%A.out
#SBATCH -e md/1OKE/apo/logs/slurm_fresean_%A.err

set -euo pipefail
export GMX_MAXBACKUP=-1          # module Gromacs mdrun has no -noback flag
ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
RUNDIR="$ROOT/md/1OKE/apo"
CFG="$RUNDIR/configs_fresean"
SRC="$RUNDIR/out"
REP="$RUNDIR/out_fresean/rep1"
NCHUNKS=20                      # 20 x 1 ns = 20 ns
SEED=20260707

mkdir -p "$REP" "$RUNDIR/logs"
cd "$REP"

[[ -f "$CFG/md_fresean_chunk1.mdp" ]] || { echo "ERROR: chunk1 mdp missing"; exit 1; }
[[ -f "$SRC/npt.gro" ]]              || { echo "ERROR: $SRC/npt.gro missing"; exit 1; }

echo "=== apo-1OKE FRESEAN chunked | $(date) | node $(hostname) | job ${SLURM_JOB_ID:-none} ==="
module load Gromacs
gmx --version | grep -i version || true

STATE="$REP/fresean_state.txt"
LAST=0; [[ -f "$STATE" ]] && LAST=$(cat "$STATE")
echo "resuming after chunk $LAST (of $NCHUNKS)"

for i in $(seq $((LAST+1)) $NCHUNKS); do
  echo "----- CHUNK $i / $NCHUNKS  $(date) -----"
  if [[ $i -eq 1 ]]; then
    sed "s/^gen_seed.*/gen_seed = $SEED/" "$CFG/md_fresean_chunk1.mdp" > c${i}.mdp
    gmx grompp -f c${i}.mdp -c "$SRC/npt.gro" -p "$SRC/topol.top" \
        -o chunk${i}.tpr -maxwarn 1 -nobackup
  else
    prev=$((i-1))
    gmx grompp -f "$CFG/md_fresean_cont.mdp" -c chunk${prev}.gro -t chunk${prev}.cpt \
        -p "$SRC/topol.top" -o chunk${i}.tpr -maxwarn 1 -nobackup
  fi

  # mdrun for this 1 ns chunk (checkpoint every 15 min for within-chunk restart)
  gmx mdrun -v -deffnm chunk${i} -ntmpi 1 -ntomp 8 -gpu_id 0 -pin on -cpt 15

  # --- reduce to protein-only (velocities preserved), verify, delete full ---
  echo "Protein" | gmx trjconv -f chunk${i}.trr -s chunk${i}.tpr \
      -o prot_chunk${i}.trr
  gmx check -f prot_chunk${i}.trr 2>&1 | tee check${i}.out | grep -iE 'Velocities|Coords' || true
  if [[ -s prot_chunk${i}.trr ]] && grep -qi 'Velocities' check${i}.out; then
      full=$(du -h chunk${i}.trr | cut -f1); prot=$(du -h prot_chunk${i}.trr | cut -f1)
      echo "chunk $i: verified velocities in protein-only ($prot); full was $full — deleting full .trr"
      rm -f chunk${i}.trr
  else
      echo "ERROR chunk $i: protein-only empty OR no velocities — KEEPING full .trr"; exit 1
  fi

  # keep only the checkpoint+gro needed to seed the next chunk; drop older .gro/.tpr/.cpt-1
  [[ $i -gt 1 ]] && rm -f chunk$((i-1)).tpr chunk$((i-1)).gro chunk$((i-1)).cpt chunk$((i-1))_prev.cpt 2>/dev/null || true
  echo "$i" > "$STATE"
done

# --- concatenate protein-only chunks into the FRESEAN input trajectory ---
echo "----- all chunks done; concatenating -----"
gmx trjcat -f prot_chunk*.trr -o md_fresean_prot.trr -nobackup -settime <<EOF || \
gmx trjcat -f $(ls -v prot_chunk*.trr) -o md_fresean_prot.trr -nobackup
EOF
gmx check -f md_fresean_prot.trr 2>&1 | grep -iE 'Velocities|Coords|Step' | head
# reference gro (first frame, protein-only) for the extraction step
echo "Protein" | gmx trjconv -f md_fresean_prot.trr -s chunk1.tpr -o md_fresean.gro -dump 0 -nobackup 2>/dev/null || \
  cp chunk1.gro md_fresean.gro
cp chunk1.tpr md_fresean.tpr 2>/dev/null || true
echo "=== DONE $(date). Kept: md_fresean_prot.trr (velocities), md_fresean.gro, .tpr, per-chunk prot_chunk*.trr ==="
