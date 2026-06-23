#!/bin/bash
# render_ddg_figures.sh — render per-antibody ddG-colored structure panels for the manuscript.
# Reads manuscript/figures/ddg/manifest.tsv (from analysis/make_ddg_coloring.py) and writes
# manuscript/figures/<key>_<pdb>_ddg.png. Runs LOCALLY (VMD + external Tachyon + macOS sips).
#   VMD=/path/to/vmd analysis/render_ddg_figures.sh    (if vmd not on PATH)
# Then rebuild: cd manuscript && tectonic main.tex
set -euo pipefail
cd "$(dirname "$0")/.."
VMD="${VMD:-vmd}"
TACHYON="${TACHYON:-/Applications/VMD.app/Contents/vmd2/lib/tachyon_MACOSXARM64}"
RES="${RES:-1600}"
TCL=analysis/vmd/render_ddg.tcl
MAN="${MAN:-manuscript/figures/ddg/manifest.tsv}"

command -v "$VMD" >/dev/null || { echo "ERROR: VMD not found ('$VMD'). Set VMD=/path/to/vmd."; exit 1; }
[[ -f "$MAN" ]] || { echo "ERROR: $MAN missing — run python3 analysis/make_ddg_coloring.py first"; exit 1; }

n=0
while IFS=$'\t' read -r key pdb ab struct datafile vmin vmax; do
  [[ "$key" == "key" ]] && continue            # header
  [[ -f "$struct" && -f "$datafile" ]] || { echo "skip $key/$pdb (missing struct/data)"; continue; }
  out="manuscript/figures/${key}_${pdb}_ddg.png"
  echo ">>> $key/$pdb ($ab) scale [$vmin,$vmax] -> $out"
  "$VMD" -dispdev text -e "$TCL" -args "$struct" "$datafile" "$vmin" "$vmax" "/tmp/${key}_${pdb}.dat" </dev/null >/dev/null 2>&1
  "$TACHYON" "/tmp/${key}_${pdb}.dat" -res "$RES" "$RES" -aasamples 12 \
      -format TARGA -o "/tmp/${key}_${pdb}_ddg.tga" >/dev/null 2>&1
  sips -s format png "/tmp/${key}_${pdb}_ddg.tga" --out "$out" >/dev/null
  magick "$out" -trim +repage -bordercolor white -border 4% "$out"   # crop whitespace
  n=$((n+1))
done < "$MAN"
echo "rendered $n ddG panel(s). Now: cd manuscript && tectonic main.tex"
