# render_ddg.tcl — cartoon colored per-residue by ddG -> Tachyon scene (headless VMD).
# Usage:
#   vmd -dispdev text -e render_ddg.tcl -args <structure> <datafile> <vmin> <vmax> <out.dat>
# <datafile> has "<resid> <ddg>" lines. Labelled residues are colored by ddG on a
# blue->white->red scale (white at the mid of [vmin,vmax]); all other residues are gray.
# Pair with analysis/render_ddg_figures.sh (external tachyon -> PNG). See [[vmd-render-pipeline]].

set struct   [lindex $argv 0]
set datafile [lindex $argv 1]
set vmin     [lindex $argv 2]
set vmax     [lindex $argv 3]
set outscene [lindex $argv 4]

mol new $struct waitfor all

# read resid -> ddg, stamp ddg into the beta field
set resids {}
set fp [open $datafile r]
while {[gets $fp line] >= 0} {
    if {[string trim $line] eq ""} continue
    lassign $line resid ddg
    set sel [atomselect top "resid $resid"]
    $sel set beta $ddg
    $sel delete
    lappend resids $resid
}
close $fp

mol delrep 0 top
color scale method BWR
catch {color scale midpoint 0.5}

# rep 0: labelled residues, colored by ddG (Beta) on the fixed [vmin,vmax] scale
mol representation NewCartoon 0.30 12.0 4.5
mol color Beta
mol selection "protein and resid $resids"
catch {mol material AOChalky}
mol addrep top
catch {mol scaleminmax top 0 $vmin $vmax}

# rep 1: everything else, gray (distinct from the white midpoint)
mol representation NewCartoon 0.30 12.0 4.5
mol color ColorID 2
mol selection "protein and not (resid $resids)"
catch {mol material AOChalky}
mol addrep top

catch {display projection Orthographic}
catch {display depthcue off}
catch {display ambientocclusion on}
catch {display shadows on}
catch {color Display Background white}
catch {axes location Off}
catch {display resetview}
catch {rotate x by -75}
catch {rotate y by 20}

render Tachyon $outscene
quit
