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
# external Tachyon ray-tracer (ships with VMD); override with TACHYON=/path if needed
TACHYON="${TACHYON:-/Applications/VMD.app/Contents/vmd2/lib/tachyon_MACOSXARM64}"
RES="${RES:-1600}"              # output is RES x RES px
TCL=analysis/vmd/render_structure.tcl
FIG=manuscript/figures

# key (matches \figorbox name) -> preferred (relaxed apo frame) | fallback (repo input structure)
MAP=(
  "hel_apo  - structures/1AKI.pdb"
  "lysozyme md/1AKI/apo/out/md_ref.pdb structures/1AKI.pdb"
  "hpr      md/2JEL/apo/out/md_ref.pdb structures/HPr_2JEL.pdb"
  "vegf     md/1BJ1/apo/out/md_ref.pdb structures/VEGF_1BJ1.pdb"
  "bonta1   md/2NYY/apo/out/md_ref.pdb structures/BontA1-Hc_2NYY.pdb"
  "mtsp1    md/3NPS/apo/out/md_ref.pdb -"
)

command -v "$VMD" >/dev/null || { echo "ERROR: VMD not found ('$VMD'). Set VMD=/path/to/vmd."; exit 1; }
rendered=0
for entry in "${MAP[@]}"; do
  set -- $entry; key=$1; ref=$2; fallback=$3
  if   [[ -f "$ref" ]];      then pdb=$ref;      src="relaxed apo frame"
  elif [[ "$fallback" != "-" && -f "$fallback" ]]; then pdb=$fallback; src="input structure"
  else echo "skip $key  (no md_ref.pdb and no input structure yet)"; continue; fi
  echo ">>> rendering $key from $pdb  ($src)"
  "$VMD" -dispdev text -e "$TCL" -args "$pdb" "/tmp/${key}.dat" >/dev/null 2>&1
  "$TACHYON" "/tmp/${key}.dat" -res "$RES" "$RES" -aasamples 12 \
      -format TARGA -o "/tmp/${key}_structure.tga" >/dev/null 2>&1
  sips -s format png "/tmp/${key}_structure.tga" --out "$FIG/${key}_structure.png" >/dev/null
  magick "$FIG/${key}_structure.png" -trim +repage -bordercolor white -border 4% "$FIG/${key}_structure.png"
  echo "    -> $FIG/${key}_structure.png"
  rendered=$((rendered+1))
done
echo "rendered $rendered figure(s). Now: cd manuscript && tectonic main.tex (then commit figures + main.pdf)."
