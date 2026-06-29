#!/bin/bash
# submit_apo.sh — submit prep -> production chains for one or more apo systems on Gemini.
# Production is held until its prep finishes cleanly (--dependency=afterok).
#
# Run FROM THE REPO ROOT (scripts anchor to $SLURM_SUBMIT_DIR; submitting elsewhere breaks paths):
#   bash md/submit_apo.sh 1JRH 1OKE 1HGG 1HGU 1AHW 2JEL 1BJ1
#   bash md/submit_apo.sh 1AKI            # re-run HEL alone
#
# Each <PDB> must have md/<PDB>/apo/{prep.sh,md.sh}. Skips (does not abort) any that don't,
# or whose required input is missing. gB (5C6T) is intentionally NOT ready — it needs its
# AF model first (see md/5C6T/apo/prep.sh header).
[[ $# -ge 1 ]] || { echo "usage: bash md/submit_apo.sh <PDB> [<PDB> ...]"; exit 1; }
ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
cd "$ROOT" || { echo "ERROR: cannot cd to repo root $ROOT"; exit 1; }

command -v sbatch >/dev/null || { echo "ERROR: sbatch not found (run on a Gemini submit host)"; exit 1; }

ok=0; skip=0
for s in "$@"; do
  prep="md/$s/apo/prep.sh"; prod="md/$s/apo/md.sh"
  if [[ ! -f "$prep" || ! -f "$prod" ]]; then
    echo "SKIP $s: missing $prep or $prod"; skip=$((skip+1)); continue
  fi
  pj=$(sbatch --parsable "$prep") || { echo "SKIP $s: prep submit failed"; skip=$((skip+1)); continue; }
  qj=$(sbatch --parsable --dependency=afterok:"$pj" "$prod") \
      || { echo "WARN $s: prep=$pj submitted but prod submit failed"; skip=$((skip+1)); continue; }
  echo "$s: prep=$pj  prod=$qj (afterok:$pj)"
  ok=$((ok+1))
done
echo "---"
echo "submitted $ok chain(s), skipped $skip. Watch: squeue -u \$USER   (prod shows PD/(Dependency))"
echo "After: check disulfide counts in md/<PDB>/apo/logs/slurm_prep_*.out before trusting production."
