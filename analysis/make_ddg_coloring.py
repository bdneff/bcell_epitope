#!/usr/bin/env python3
"""Prepare per-antibody ddG residue-coloring data + colorbars for the structure figures.

For each simulated antigen with an available structure, and each antibody scanned against
it, writes a `<resid> <ddg>` data file (consumed by analysis/vmd/render_ddg.tcl) plus a
manifest, and renders a shared per-antigen colorbar PNG. ddG is colored blue->white->red
(white at the mid of the per-antigen range); residues with no label are left gray by the
renderer. Numbering is QC'd against the actual structure before writing.

Outputs:
  manuscript/figures/ddg/<key>_<pdb>.dat     resid ddg  (one panel = one antibody)
  manuscript/figures/ddg/manifest.tsv        key pdb ab struct datafile vmin vmax
  manuscript/figures/<key>_ddgbar.png        shared colorbar for that antigen
Run: python3 analysis/make_ddg_coloring.py
"""
import csv, collections, pathlib
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm, colors

CSV = pathlib.Path("benchmark/antigen_alanine_scan_extracted_SIMPLE_v1.csv")
DDGDIR = pathlib.Path("manuscript/figures/ddg")
FIGDIR = pathlib.Path("manuscript/figures")

# antigen -> (key, input structure with benchmark numbering). mtsp1 has no structure yet.
ANTIGEN = {
    "Lysozyme": ("lysozyme", "structures/1AKI.pdb"),
    "HPr":      ("hpr",      "structures/HPr_2JEL.pdb"),
    "VEGF":     ("vegf",     "structures/VEGF_1BJ1.pdb"),
    "Bont/A1":  ("bonta1",   "structures/BontA1-Hc_2NYY.pdb"),
}
ABBREV = {
    "Neutralizing mAb (fab-12)": "fab-12", "affinity matured mAb (Y0317)": "Y0317",
    "mAb CR1": "CR1", "mAb AR2": "AR2", "Fab Inhibitor E2": "E2", "Fab Inhibitor S4": "S4",
    "Jel42 antibody": "Jel42", "HYHEL-63": "HyHEL-63",
}

def struct_resids(path):
    rs = set()
    for line in open(path):
        if line[:4] == "ATOM":
            try: rs.add(int(line[22:26]))
            except ValueError: pass
    return rs

rows = [r for r in csv.DictReader(open(CSV))
        if r["Use for MD correlation?"].startswith("YES") and r["Antigen"].strip() in ANTIGEN]

# antigen -> pdb -> {resid: ddg}; and antigen -> ordered (pdb, ab)
data = collections.defaultdict(lambda: collections.defaultdict(dict))
complexes = collections.defaultdict(list); seen = collections.defaultdict(set)
for r in rows:
    ag = r["Antigen"].strip(); pdb = r["PDB/model"].strip()
    if pdb not in seen[ag]:
        seen[ag].add(pdb); complexes[ag].append((pdb, r["Antibody"].strip()))
    try:
        data[ag][pdb][int(r["Residue #"])] = float(r["Best ΔΔG kcal/mol"])
    except ValueError:
        pass

DDGDIR.mkdir(parents=True, exist_ok=True)
# ONE global symmetric scale for ALL antigens: +/- the largest |ddG| in the entire
# benchmark, so white=0 and colours are comparable across every figure.
GM = 0.0
for r in csv.DictReader(open(CSV)):
    if not r["Use for MD correlation?"].startswith("YES"): continue
    try: GM = max(GM, abs(float(r["Best ΔΔG kcal/mol"])))
    except ValueError: pass
vmin, vmax = -GM, GM
print(f"global symmetric scale +/-{GM:.2f} (white at 0), shared across all colorbars")

manifest = ["\t".join(["key","pdb","ab","struct","datafile","vmin","vmax"])]
for ag, (key, struct) in ANTIGEN.items():
    if not pathlib.Path(struct).exists():
        print(f"skip {ag}: {struct} missing"); continue
    present = struct_resids(struct)
    for pdb, ab in complexes[ag]:
        d = data[ag][pdb]
        miss = sorted(set(d) - present)
        if miss: print(f"  WARN {key}/{pdb}: resids not in structure (skipped): {miss}")
        df = DDGDIR / f"{key}_{pdb}.dat"
        with open(df, "w") as f:
            for resid in sorted(d):
                if resid in present: f.write(f"{resid} {d[resid]:.3f}\n")
        manifest.append("\t".join([key, pdb, ABBREV.get(ab, ab), struct, str(df), f"{vmin:.3f}", f"{vmax:.3f}"]))
        print(f"  {key}/{pdb} ({ABBREV.get(ab,ab)}): {sum(1 for x in d if x in present)} residues")
(DDGDIR / "manifest.tsv").write_text("\n".join(manifest) + "\n")
print(f"wrote manifest with {len(manifest)-1} panels")

# ONE shared colorbar (identical dimensions/scale/labels/ticks for every figure)
fig, ax = plt.subplots(figsize=(0.7, 3.0))
cb = matplotlib.colorbar.ColorbarBase(ax, cmap=matplotlib.colormaps["bwr"],
        norm=colors.Normalize(vmin=vmin, vmax=vmax), orientation="vertical")
cb.set_label(r"$\Delta\Delta G_\mathrm{bind}$ (kcal/mol)", fontsize=18)
cb.ax.tick_params(labelsize=16)
fig.savefig(FIGDIR / "ddgbar.png", dpi=200, bbox_inches="tight"); plt.close(fig)
print("wrote shared colorbar figures/ddgbar.png")
