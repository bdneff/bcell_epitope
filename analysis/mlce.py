#!/usr/bin/env python
"""mlce.py -- Colombo MLCE (matrix of low coupling energies) on gRINN pairwise
interaction-energy matrices, scored against epitope labels. Tests the energetic-coupling
thesis (Serapian/Colombo 2020, doi:10.1021/acs.jpclett.0c02341) on our systems.

MLCE: eigendecompose the symmetric residue-residue interaction-energy matrix; reconstruct
the coupling map from the most stabilizing (most negative eigenvalue) modes; per-residue
local coupling energy = row-sum of |reconstructed coupling|. LOW coupling surface residues
= predicted antigenic.

Inputs : md/<sys>/apo/grinn/out/average_interaction_energies.csv  (on Gemini scratch;
         pull to mlce_in/<sys>_pairs.csv first)
         benchmark/labels/ddg_by_md_residue.csv
         benchmark/features/dynamics_sasa_modes.csv (relSASA for surface restriction)
Output : figures + mlce_coupling z-score per (pdb,resid)
"""
import pandas as pd, numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score

def build_energy_matrix(csv):
    df = pd.read_csv(csv)
    resids = sorted(set(df["Res1_ResNum"]).union(df["Res2_ResNum"]))
    idx = {r: i for i, r in enumerate(resids)}
    n = len(resids); M = np.zeros((n, n))
    for r in df.itertuples():
        i, j = idx[r.Res1_ResNum], idx[r.Res2_ResNum]
        M[i, j] = r.Avg_Total_Energy; M[j, i] = r.Avg_Total_Energy
    return np.array(resids), M

def mlce(resids, M, n_modes=None):
    w, V = np.linalg.eigh(M)                 # symmetric interaction-energy matrix
    order = np.argsort(w)                    # most stabilizing (negative) first
    if n_modes is None:
        n_modes = max(2, int(0.1 * len(resids)))
    sel = order[:n_modes]
    C = (V[:, sel] * w[sel]) @ V[:, sel].T   # reconstruct coupling from top stabilizing modes
    return np.abs(C).sum(1)                  # per-residue local coupling energy (LOW = antigenic)

if __name__ == "__main__":
    systems = ["1AKI","2JEL","1AHW","1BJ1","1HGU","1JRH","1OKE"]
    rows = []
    for sys in systems:
        resids, M = build_energy_matrix(f"mlce_in/{sys}_pairs.csv")
        lce = mlce(resids, M)
        z = (lce - lce.mean()) / lce.std()
        for r, l in zip(resids, z):
            rows.append(dict(pdb=sys, resid=int(r), mlce_coupling_z=l))
    pd.DataFrame(rows).to_csv("benchmark/features/mlce_coupling.csv", index=False)
    print("wrote benchmark/features/mlce_coupling.csv")
