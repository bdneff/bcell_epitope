#!/bin/bash
# render_structures.sh — generate per-antigen structure figures for the manuscript.
# Runs LOCALLY (needs VMD + macOS `sips`). For each system whose postprocess output
# (md_ref.pdb) is present locally, renders a clean cartoon to
# manuscript/figures/<key>_structure.png — which \figorbox then picks up automatically.
#
# VMD not on PATH? run:  VMD=/path/to/vmd analysis/render_structures.sh
# Then rebuild the PDF:  cd manuscript && tectonic main.tex
set -euo pipefail
cd "$(dirname "$0")/.."          # repo root
VMD="${VMD:-vmd}"
TCL=analysis/vmd/render_structure.tcl
FIG=manuscript/figures

# key (matches \figorbox name) -> structure file from that run's postprocess output
MAP=(
  "lysozyme md/1AKI/apo/out/md_ref.pdb"
  "hpr      md/2JEL/apo/out/md_ref.pdb"
  "vegf     md/1BJ1/apo/out/md_ref.pdb"
  "bonta1   md/2NYY/apo/out/md_ref.pdb"
  "mtsp1    md/3NPS/apo/out/md_ref.pdb"
)

command -v "$VMD" >/dev/null || { echo "ERROR: VMD not found ('$VMD'). Set VMD=/path/to/vmd."; exit 1; }
rendered=0
for entry in "${MAP[@]}"; do
  set -- $entry; key=$1; pdb=$2
  if [[ ! -f "$pdb" ]]; then echo "skip $key  ($pdb not present — download/run it first)"; continue; fi
  echo ">>> rendering $key from $pdb"
  "$VMD" -dispdev text -e "$TCL" -args "$pdb" "/tmp/${key}_structure.tga" >/dev/null 2>&1
  sips -s format png "/tmp/${key}_structure.tga" --out "$FIG/${key}_structure.png" >/dev/null
  echo "    -> $FIG/${key}_structure.png"
  rendered=$((rendered+1))
done
echo "rendered $rendered figure(s). Now: cd manuscript && tectonic main.tex (then commit figures + main.pdf)."
