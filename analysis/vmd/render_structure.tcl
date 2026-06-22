# render_structure.tcl — clean cartoon render of a structure to a TGA (headless VMD).
# Usage:
#   vmd -dispdev text -e render_structure.tcl -args <structure.pdb|.gro> <out.tga> [resid ...]
# Optional trailing resids are drawn as red licorice over the cartoon (e.g. epitope residues).
# Pair with analysis/render_structures.sh, which converts the TGA to PNG (macOS `sips`).

if {[llength $argv] < 2} {
    puts "usage: vmd -dispdev text -e render_structure.tcl -args <structure> <out.tga> \[resid...\]"
    quit
}
set struct [lindex $argv 0]
set outimg [lindex $argv 1]
set hilite [lrange $argv 2 end]

mol new $struct waitfor all
mol delrep 0 top

# base cartoon, colored by secondary structure, matte (ambient-occlusion) finish
mol representation NewCartoon 0.30 12.0 4.5
mol color Structure
mol selection {protein}
mol material AOChalky
mol addrep top

# optional: highlight alanine-scan / epitope residues as red licorice
if {[llength $hilite] > 0} {
    mol representation Licorice 0.30 12.0 12.0
    mol color ColorID 1
    mol selection "protein and resid $hilite"
    mol material AOShiny
    mol addrep top
}

# display / lighting for a publication-style still
display projection Orthographic
display depthcue off
display ambientocclusion on
display shadows on
color Display Background white
axes location Off
catch {display resize 1600 1600}
display resetview
rotate x by -75
rotate y by 20

render TachyonInternal $outimg
quit
