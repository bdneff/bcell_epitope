#!/usr/bin/env python
"""Per-residue crystallographic B-factor for each antigen -- the STATIC, experimental proxy
for residue mobility (the Westhof 1984 flexibility hypothesis), and the static analogue of the
RMSF we will get from MD. A B-factor signal here is the teaser for the dynamics features.

Uses the RAW crystal structures (NOT the repaired _fixed ones: AF-grafted loops carry pLDDT,
not crystallographic B). Raw B-factors are not comparable across crystals (resolution,
refinement, overall scaling), so we z-normalise PER STRUCTURE; that z-score is the feature.

CONFOUND: 6 of the 8 antigens here are deposited as antibody-antigen COMPLEXES (only 1AKI and
1HGU are apo). With the antibody bound, the epitope residues' B-factors are damped by the
interface contacts -- which can by itself produce a low-B-at-hotspots signal. Apo-MD RMSF is the
clean unbound measurement; this static B-factor is a confounded teaser (see manuscript sec:baseline).

Output: benchmark/features/bfactor_crystal.csv
  columns: antigen,key,pdb,chain,resid,resname,bfactor,bfactor_z
B per residue = CA B-factor (fallback: mean over residue heavy atoms).

  micromamba run -n bcell-repair python analysis/compute_bfactor.py
"""
import csv
import pathlib
import numpy as np
from Bio.PDB import PDBParser

# antigen -> (key, RAW crystal structure with benchmark numbering). Crystal B-factors only.
ANTIGEN = {
    "Lysozyme": ("lysozyme", "structures/1AKI.pdb"),
    "HPr":      ("hpr",      "structures/HPr_2JEL.pdb"),
    "VEGF":     ("vegf",     "structures/VEGF_1BJ1.pdb"),
    "Bont/A1":  ("bonta1",   "structures/BontA1-Hc_2NYY.pdb"),
    "IFN-gamma receptor":  ("ifngr",  "structures/IFNgR_1JRH.pdb"),
    "Tissue factor":       ("tf",     "structures/TissueFactor_1AHW.pdb"),
    "HCMV glycoprotein B": ("hcmvgb", "structures/HCMVgB_5C6T.pdb"),
    "human growth hormone":("hgh",    "structures/hGH_1HGU.pdb"),
}
AA = {"ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","ILE","LEU","LYS",
      "MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}

OUT = pathlib.Path("benchmark/features/bfactor_crystal.csv")
OUT.parent.mkdir(parents=True, exist_ok=True)
parser = PDBParser(QUIET=True)

rows = []
for antigen, (key, struct) in ANTIGEN.items():
    p = pathlib.Path(struct)
    if not p.exists():
        print(f"skip {antigen}: {struct} missing")
        continue
    pdb = p.stem.split("_")[-1] if "_" in p.stem else p.stem
    model = parser.get_structure(key, struct)[0]
    recs = []  # (chain, resid, resname, B)
    for chain in model:
        for res in chain:
            if res.id[0] != " " or res.resname not in AA:
                continue
            if res.has_id("CA"):
                b = float(res["CA"].get_bfactor())
            else:
                b = float(np.mean([a.get_bfactor() for a in res]))
            recs.append((chain.id, res.id[1], res.resname, b))
    bvals = np.array([r[3] for r in recs])
    mu, sd = bvals.mean(), bvals.std()
    sd = sd if sd > 1e-6 else 1.0
    for ch, resid, resn, b in recs:
        rows.append([antigen, key, pdb, ch, resid, resn, f"{b:.2f}", f"{(b - mu) / sd:.4f}"])
    print(f"  {key:9s} {struct:34s} {len(recs)} res, B {bvals.min():.1f}-{bvals.max():.1f} (mu {mu:.1f})")

with open(OUT, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["antigen", "key", "pdb", "chain", "resid", "resname", "bfactor", "bfactor_z"])
    w.writerows(rows)
print(f"wrote {OUT} ({len(rows)} residues)")
