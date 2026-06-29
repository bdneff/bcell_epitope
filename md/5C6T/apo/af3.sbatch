#!/usr/bin/env bash
# =============================================================================
# af3.sbatch — AlphaFold3 prediction of the HCMV gB ectodomain (5C6T construct), to
# supply coordinates for the unresolved loops (gB is absent from AlphaFold DB).
# Mirrors the Altin-lab AF3 recipe (episcaf episcaf_pipeline/hpc/sbatch/af3_array.sbatch):
# singularity container + installed weights/DBs on Gemini. Submit FROM THE REPO ROOT:
#   sbatch md/5C6T/apo/af3.sbatch
# Input  : md/5C6T/apo/af3/gB_5C6T.json  (single-sequence; analysis/make_af3_json.py)
# Output : md/5C6T/apo/af3/out/gb_5c6t/  -> *_model.cif (AF3 lowercases the name)
# Next   : pull back, then `bash analysis/repair_gB.sh` (see md/5C6T/apo/prep.sh header).
# =============================================================================
#SBATCH -J gB_af3
#SBATCH -p gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -t 00:30:00
#SBATCH -o md/5C6T/apo/logs/slurm_af3_%j.out
#SBATCH -e md/5C6T/apo/logs/slurm_af3_%j.err

set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
RUNDIR="$ROOT/md/5C6T/apo"
JSON="$RUNDIR/af3/gB_5C6T.json"
OUT="$RUNDIR/af3/out"
mkdir -p "$OUT" "$RUNDIR/logs"
[[ -f "$JSON" ]] || { echo "ERROR: $JSON not found (run analysis/make_af3_json.py)"; exit 1; }

module load singularity

# Lab-installed AF3 3.0.1 (paths from episcaf af3_array.sbatch).
SIF=/tgen_labs/altin/alphafold3/containers/alphafold_3.0.1.sif
MODELS=/ref_genomes/alphafold/alphafold3/models
DBS=/ref_genomes/alphafold/alphafold3/

# Proven recipe: single-sequence (JSON carries a self-MSA) + --norun_data_pipeline.
# To instead build a real MSA (better for this natural protein, IF the genetic DBs are
# installed): regenerate the JSON with `--mode pipeline` and DROP --norun_data_pipeline
# below (and raise -t to a few hours).
echo "=== gB AF3 | $(date) | node $(hostname) | job ${SLURM_JOB_ID} ==="
singularity exec --nv -B /home,/scratch,/tgen_labs,/ref_genomes "$SIF" \
    python /app/alphafold/run_alphafold.py \
        --json_path "$JSON" \
        --model_dir "$MODELS" \
        --db_dir "$DBS" \
        --output_dir "$OUT" \
        --norun_data_pipeline \
        --num_diffusion_samples 1

# Surface the produced model path for convenience.
echo "=== AF3 done $(date). Model(s): ==="
find "$OUT" -name '*_model.cif' -printf '  %p\n' || true
echo "Next: pull md/5C6T/apo/af3/out back, then bash analysis/repair_gB.sh"
