#!/usr/bin/env python
"""
Dengue E (1OKE) single-system combination model — independent reproduction.

Brandon asked: the single features are weak; can a linear combination predict the
binary shotgun epitope label (96 critical / ~298 non-critical, one well-sampled system)?

This is deliberately INDEPENDENT of the Claude Science pipeline: it recomputes SASA and
crystal B-factor locally from structures/DengueE_1OKE.pdb (neither is in any committed CSV),
verifies numbering by WT identity against the label AND the MD feature files, then scores
single features + incremental logistic combinations under two CV schemes:
  - stratified repeated k-fold (Science's scheme; residues i.i.d. -- OPTIMISTIC, spatial leak)
  - spatially-blocked k-fold (contiguous sequence blocks -- honest for one autocorrelated chain)

Run:  micromamba run -n bcell-repair python analysis/dengue_combination.py

NOTE: single system, so this is a within-protein detection test, NOT the leave-one-antigen-out
protocol used for Tier-1. No p-value is claimed for any single feature; the question is only
whether a combination beats the SASA exposure control on held-out residues of THIS protein.
"""
import numpy as np, pandas as pd, pathlib
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import RepeatedStratifiedKFold, cross_val_score
from sklearn.metrics import roc_auc_score

ROOT = pathlib.Path(__file__).resolve().parent.parent
FEAT = ROOT / "benchmark/features"
THREE2ONE = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
             'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
             'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

# ---------------------------------------------------------------- labels + MD numbering
lab = pd.read_csv(ROOT / "benchmark/tier2_labels_dengueE.csv")     # 96 critical residues only
crit_resids = set(lab.resid); wt = dict(zip(lab.resid, lab.WT))
dccm = pd.read_csv(FEAT / "dccm_apomd.csv"); dccm = dccm[dccm.pdb == "1OKE"].copy()
UNIVERSE = sorted(dccm.resid)                                       # 394 residues = MD numbering
res3 = dict(zip(dccm.resid, dccm.resname))                         # MD resname per resid (1-letter check below)
m = sum(THREE2ONE.get(res3[r]) == wt[r] for r in crit_resids)
print(f"[align] MD (dccm) vs label WT identity: {m}/{len(crit_resids)}  (universe {len(UNIVERSE)} residues)")
assert m == len(crit_resids), "MD numbering does not match label WT -- STOP"

# ---------------------------------------------------------------- recompute SASA + B-factor
struct = PDBParser(QUIET=True).get_structure("d", ROOT / "structures/DengueE_1OKE.pdb")
model = struct[0]
# pick the chain whose residue count matches the MD universe (the modelled monomer)
chains = [(c.id, [r for r in c if r.id[0] == ' ']) for c in model]
chain_id, residues = max(chains, key=lambda kv: len(kv[1]))
print(f"[struct] chains={[c for c,_ in chains]}; using chain {chain_id} ({len(residues)} residues) for SASA")
# offset search: crystal author numbering may differ from MD 1..394; find offset by WT identity
struct_res = {r.id[1]: r for r in residues}
struct_one = {rid: THREE2ONE.get(r.resname) for rid, r in struct_res.items()}
best = None
for off in range(-30, 31):
    hit = sum(struct_one.get(r + off) == wt[r] for r in crit_resids)
    if best is None or hit > best[1]:
        best = (off, hit)
offset, hit = best
print(f"[struct] best crystal->MD offset = {offset:+d}  (WT identity {hit}/{len(crit_resids)})")
assert hit == len(crit_resids), "no crystal offset recovers the label WT -- SASA/Bfactor numbering unsafe"

ShrakeRupley().compute(model, level="R")
sasa = {}; bfac = {}; prov = []
for r in residues:
    md_rid = r.id[1] - offset                                      # map crystal -> MD numbering
    sasa[md_rid] = r.sasa
    ca = r["CA"] if "CA" in r else None
    bfac[md_rid] = ca.get_bfactor() if ca is not None else np.nan
    prov.append(("1OKE", chain_id, md_rid, r.resname, r.sasa, bfac[md_rid]))
# persist the two crystal features so the provenance is IN-REPO (they were in no committed CSV)
prov_df = pd.DataFrame(prov, columns=["pdb","chain","resid","resname","sasa","bfactor"]).sort_values("resid")
prov_df.to_csv(FEAT / "dengue_1OKE_crystal_sasa_bfactor.csv", index=False)
print(f"[persist] wrote {FEAT/'dengue_1OKE_crystal_sasa_bfactor.csv'} "
      f"(single-frame Shrake-Rupley SASA + CA B-factor, HOLO 1OKE crystal, MD numbering)")

# ---------------------------------------------------------------- assemble feature matrix
def oke(f, cols):
    d = pd.read_csv(FEAT / f); d = d[d.pdb == "1OKE"]
    return d.set_index("resid")[cols]

df = pd.DataFrame(index=UNIVERSE)
df["SASA"]     = pd.Series(sasa)
df["Bfactor"]  = pd.Series(bfac)
hb = oke("hbonds_apomd.csv", ["hb_intra", "hb_water"])
df["HBintra"], df["HBwater"] = hb.hb_intra, hb.hb_water
g = oke("grinn_intenergy_apomd.csv", ["intEn_vdw", "intEn_elec", "intEn_total"])
df["vdW"], df["elec"], df["Etot"] = g.intEn_vdw, g.intEn_elec, g.intEn_total
df["DCCM"] = dccm.set_index("resid")["coupling_abs"]
dih = oke("dihedral_flex_apomd.csv", ["bb_dih_std", "chi1_std"])
df["bbflex"] = dih.bb_dih_std   # chi1 dropped: undefined for Gly/Ala (165/394 NaN), not a universal feature
df["y"] = [1 if r in crit_resids else 0 for r in df.index]

FEATURES = ["SASA","Bfactor","HBintra","HBwater","vdW","elec","Etot","DCCM","bbflex"]
dense = df[FEATURES + ["y"]].dropna()   # ~392 residues; bbflex has 2 termini NaN
print(f"[matrix] {len(dense)} residues (dense set); critical={int(dense.y.sum())}, "
      f"non-critical={int((dense.y==0).sum())}")

# ---------------------------------------------------------------- single-feature signed AUROC
print("\n[single-feature signed detection AUROC  (P(critical>non-crit); <0.5 => critical are LOWER)]")
for f in FEATURES:
    s = df[[f, "y"]].dropna()
    a = roc_auc_score(s.y.values, s[f].values)
    print(f"   {f:9s} {a:.3f}  n={len(s):3d}   {'(critical HIGHER)' if a>=0.5 else '(critical LOWER)'}")

# ---------------------------------------------------------------- CV helpers
def blocked_folds(n, k=5, block=8):
    """contiguous sequence blocks assigned round-robin to folds -> spatial CV."""
    fold = np.empty(n, int)
    for i, start in enumerate(range(0, n, block)):
        fold[start:start+block] = i % k
    return fold

def cv_auc(cols, scheme):
    clf = make_pipeline(StandardScaler(), LogisticRegression(C=1.0, max_iter=2000))
    sub = df[cols + ["y"]].dropna()          # per-model residue set (keeps N maximal)
    Xc, yy = sub[cols].values, sub.y.values
    if scheme == "strat":
        cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=0)
        s = cross_val_score(clf, Xc, yy, cv=cv, scoring="roc_auc")
        return s.mean(), s.std(), len(sub)
    fold = blocked_folds(len(yy))            # spatially-blocked (residues already in sequence order)
    aucs = []
    for f in range(5):
        tr, te = fold != f, fold == f
        if yy[te].sum() == 0 or yy[te].sum() == te.sum():
            continue
        clf.fit(Xc[tr], yy[tr])
        aucs.append(roc_auc_score(yy[te], clf.predict_proba(Xc[te])[:, 1]))
    return float(np.mean(aucs)), float(np.std(aucs)), len(sub)

# ---------------------------------------------------------------- incremental combinations
INCR = [("SASA (control)", ["SASA"]),
        ("+ Bfactor",      ["SASA","Bfactor"]),
        ("+ DCCM",         ["SASA","Bfactor","DCCM"]),
        ("+ bbflex",       ["SASA","Bfactor","DCCM","bbflex"]),
        ("all 10 features", FEATURES)]
print("\n[incremental logistic AUROC]   (Science ref: SASA 0.628 -> +Bfac 0.684 -> +DCCM 0.696)")
print(f"   {'model':18s} {'strat-5x5':>14s}   {'blocked-5':>14s}   n")
for name, cols in INCR:
    ms, ss, n = cv_auc(cols, "strat"); mb, sb, _ = cv_auc(cols, "blocked")
    print(f"   {name:18s} {ms:.3f}+/-{ss:.3f}   {mb:.3f}+/-{sb:.3f}   {n}")

# full-model coefficient signs (standardized)
clf = make_pipeline(StandardScaler(), LogisticRegression(C=1.0, max_iter=2000)).fit(
    dense[FEATURES].values, dense.y.values)
coef = clf.named_steps["logisticregression"].coef_[0]
print("\n[full-model standardized coefficients  (sign = direction of critical association)]")
for f, c in sorted(zip(FEATURES, coef), key=lambda kv: -abs(kv[1])):
    print(f"   {f:9s} {c:+.3f}")
