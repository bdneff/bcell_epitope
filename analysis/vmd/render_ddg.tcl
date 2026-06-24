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

set mid [mol new $struct waitfor all]

# Stamp values into beta and build an EXPLICIT chain-aware selection of labelled
# residues (the sentinel/"beta>x" trick does not colour reliably in VMD 2.0a7).
# Data lines: "<resid> <value>" or "<chain> <resid> <value>" (chain-aware, so HA1/HA2
# numbered from 1 don't collide).
set clauses {}
set fp [open $datafile r]
while {[gets $fp line] >= 0} {
    if {[string trim $line] eq ""} continue
    if {[llength $line] >= 3} {
        lassign $line chain resid val
        set sel [atomselect top "chain $chain and resid $resid"]
        lappend clauses "(chain $chain and resid $resid)"
    } else {
        lassign $line resid val
        set sel [atomselect top "resid $resid"]
        lappend clauses "(resid $resid)"
    }
    $sel set beta $val
    $sel delete
}
close $fp
set labelled [join $clauses " or "]

mol delrep 0 top
color scale method BWR
catch {color scale midpoint 0.5}

# rep 0: labelled residues, colored by value (Beta) on the fixed [vmin,vmax] scale
mol representation NewCartoon 0.30 12.0 4.5
mol color Beta
mol selection "protein and ($labelled)"
catch {mol material AOChalky}
mol addrep $mid
set rep0 [expr {[molinfo $mid get numreps] - 1}]
catch {mol scaleminmax $mid $rep0 $vmin $vmax}

# rep 1: everything else, gray (distinct from the white midpoint)
mol representation NewCartoon 0.30 12.0 4.5
mol color ColorID 2
mol selection "protein and not ($labelled)"
catch {mol material AOChalky}
mol addrep $mid

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
