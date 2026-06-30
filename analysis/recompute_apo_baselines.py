#!/usr/bin/env python
"""Recompute the SASA and B-factor baselines on genuine APO (antibody-free) crystal structures,
to remove the apo-vs-complex confound (manuscript sec:baseline): in the original baselines 6/8
structures were antibody-antigen complexes, so the antibody damped the epitope residues' B-factors
(and the SASA used the bound conformation).

Apo replacements (verified to cover the alanine-scan residues; see session notes):
  lysozyme 1AKI*, hGH 1HGU*  (*already apo, local files)
  HPr 1OPD, VEGF 4KZN (homodimer), tissue factor 2HFT, BoNT/A1 Hc 3FUO (+1 numbering),
  HCMV gB 5CXF (apo ectodomain, one protomer)
IFN-gamma receptor omitted (no ligand-free apo structure exists).

Per-residue feature on the antigen chain(s) only, in BENCHMARK numbering (resid + delta):
  - B-factor: CA B, z-normalised per structure
  - rel. SASA: Shrake-Rupley on the kept chains, normalised by Tien max-ASA
Writes benchmark/features/{bfactor,sasa}_apocrystal.csv and saves the prepped chains to structures/apo/.

  micromamba run -n bcell-repair python analysis/recompute_apo_baselines.py
"""
import csv, io, pathlib, urllib.request, warnings
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select
warnings.filterwarnings("ignore")

MAXASA = {"ALA":129.0,"ARG":274.0,"ASN":195.0,"ASP":193.0,"CYS":167.0,"GLU":223.0,"GLN":225.0,
          "GLY":104.0,"HIS":224.0,"ILE":197.0,"LEU":201.0,"LYS":236.0,"MET":224.0,"PHE":240.0,
          "PRO":159.0,"SER":155.0,"THR":172.0,"TRP":285.0,"TYR":263.0,"VAL":174.0}
from Bio.PDB.SASA import ShrakeRupley

# key: (antigen-display, source, chains, delta)  delta added to apo resid -> benchmark numbering
APO = {
    "lysozyme": ("Lysozyme",             "local:structures/1AKI.pdb",      ["A"], 0),
    "hgh":      ("human growth hormone", "local:structures/hGH_1HGU.pdb",  ["A"], 0),
    "hpr":      ("HPr",                  "dl:1OPD", ["A"], 0),
    "vegf":     ("VEGF",                 "dl:4KZN", ["A", "B"], 0),
    "tf":       ("Tissue factor",        "dl:2HFT", ["A"], 0),
    "bonta1":   ("Bont/A1",              "dl:3FUO", ["A"], 1),
    "hcmvgb":   ("HCMV glycoprotein B",  "dl:5CXF", ["A"], 0),
}
parser = PDBParser(QUIET=True)
sr = ShrakeRupley()
apodir = pathlib.Path("structures/apo"); apodir.mkdir(parents=True, exist_ok=True)


def load(src):
    if src.startswith("local:"):
        return parser.get_structure("x", src[6:])[0]
    pdb = src.split(":")[1]
    txt = urllib.request.urlopen(f"https://files.rcsb.org/download/{pdb}.pdb", timeout=40).read().decode()
    return parser.get_structure(pdb, io.StringIO(txt))[0], pdb


brows, srows = [], []
for key, (antigen, src, chains, delta) in APO.items():
    res = load(src)
    if isinstance(res, tuple):
        model, pdb = res
    else:
        model, pdb = res, src.split("/")[-1].split(".")[0].split("_")[-1]
    # keep only the antigen chains; strip hetero/water
    for ch_id in [c.id for c in model]:
        if ch_id not in chains:
            model.detach_child(ch_id)
    for ch in model:
        for r in [r for r in ch if r.id[0] != " " or r.resname not in MAXASA]:
            ch.detach_child(r.id)
    sr.compute(model, level="R")
    recs = []
    for ch in model:
        for r in ch:
            if not r.has_id("CA"):
                continue
            b = float(r["CA"].get_bfactor())
            sasa = float(r.sasa)
            recs.append((ch.id, r.id[1] + delta, r.resname, b, sasa))
    bvals = np.array([x[3] for x in recs]); mu, sd = bvals.mean(), bvals.std() or 1.0
    for ch, resid, resn, b, sasa in recs:
        brows.append([antigen, key, pdb, ch, resid, resn, f"{b:.2f}", f"{(b-mu)/sd:.4f}"])
        srows.append([antigen, key, pdb, ch, resid, resn, f"{sasa:.2f}", f"{sasa/MAXASA[resn]:.4f}"])
    io_ = PDBIO(); io_.set_structure(model); io_.save(str(apodir / f"{key}_{pdb}_apo.pdb"))
    print(f"  {key:9s} {pdb}  chains {chains} delta {delta:+d}  {len(recs)} res, B {bvals.min():.1f}-{bvals.max():.1f}")

with open("benchmark/features/bfactor_apocrystal.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["antigen","key","pdb","chain","resid","resname","bfactor","bfactor_z"]); w.writerows(brows)
with open("benchmark/features/sasa_apocrystal.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["antigen","key","pdb","chain","resid","resname","sasa","rel_sasa"]); w.writerows(srows)
print(f"wrote benchmark/features/bfactor_apocrystal.csv + sasa_apocrystal.csv ({len(brows)} residues, 7 apo antigens)")
