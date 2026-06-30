#!/usr/bin/env python
"""Per-residue solvent-accessible surface area (SASA) for each apo antigen structure --
the STATIC-exposure baseline feature (the control our dynamics features must beat).

Single-frame SASA (Biopython Shrake-Rupley) on the apo input structures. This is the
pre-trajectory baseline; once the Gemini runs land, swap in trajectory-AVERAGED SASA
(gmx sasa / MDAnalysis over the production xtc) at the same residue keys.

Output: benchmark/features/sasa_apo_singleframe.csv
  columns: antigen,key,pdb,chain,resid,resname,sasa,rel_sasa
rel_sasa = sasa / max-ASA(resname) (Tien et al. 2013 theoretical), so it is comparable
across residue types (0 = buried, ~1 = fully exposed).

Run in the bcell-repair env:
  micromamba run -n bcell-repair python analysis/compute_sasa.py
"""
import csv
import pathlib
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

# antigen -> (key, apo structure with benchmark numbering). Mirrors analysis/make_ddg_coloring.py.
ANTIGEN = {
    "Lysozyme": ("lysozyme", "structures/1AKI.pdb"),
    "HPr":      ("hpr",      "structures/HPr_2JEL.pdb"),
    "VEGF":     ("vegf",     "structures/VEGF_1BJ1.pdb"),
    "Bont/A1":  ("bonta1",   "structures/BontA1-Hc_2NYY.pdb"),
    "IFN-gamma receptor":  ("ifngr",  "structures/IFNgR_1JRH.pdb"),
    "Tissue factor":       ("tf",     "structures/TissueFactor_1AHW_fixed.pdb"),
    "HCMV glycoprotein B": ("hcmvgb", "structures/HCMVgB_5C6T.pdb"),
    "human growth hormone":("hgh",    "structures/hGH_1HGU_fixed.pdb"),
}

# Tien et al. 2013 (PLoS ONE) theoretical max ASA (A^2), for relative SASA normalisation.
MAXASA = {
    "ALA":129.0,"ARG":274.0,"ASN":195.0,"ASP":193.0,"CYS":167.0,"GLU":223.0,"GLN":225.0,
    "GLY":104.0,"HIS":224.0,"ILE":197.0,"LEU":201.0,"LYS":236.0,"MET":224.0,"PHE":240.0,
    "PRO":159.0,"SER":155.0,"THR":172.0,"TRP":285.0,"TYR":263.0,"VAL":174.0,
}

OUT = pathlib.Path("benchmark/features/sasa_apo_singleframe.csv")
OUT.parent.mkdir(parents=True, exist_ok=True)
parser = PDBParser(QUIET=True)
sr = ShrakeRupley()

rows = []
for antigen, (key, struct) in ANTIGEN.items():
    p = pathlib.Path(struct)
    if not p.exists():
        print(f"skip {antigen}: {struct} missing")
        continue
    pdb = p.stem.replace("_fixed", "").split("_")[-1]  # e.g. hGH_1HGU_fixed -> 1HGU; 1AKI -> 1AKI
    model = parser.get_structure(key, struct)[0]
    sr.compute(model, level="R")
    n = 0
    for chain in model:
        for res in chain:
            if res.id[0] != " " or res.resname not in MAXASA:
                continue
            sasa = float(res.sasa)
            rel = sasa / MAXASA[res.resname]
            rows.append([antigen, key, pdb, chain.id, res.id[1], res.resname,
                         f"{sasa:.2f}", f"{rel:.4f}"])
            n += 1
    print(f"  {key:9s} {struct:40s} {n} residues")

with open(OUT, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["antigen", "key", "pdb", "chain", "resid", "resname", "sasa", "rel_sasa"])
    w.writerows(rows)
print(f"wrote {OUT} ({len(rows)} residues)")
