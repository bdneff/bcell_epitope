#!/usr/bin/env python
"""external_dengue.py -- the strict test: train on 6 Tier-1 antigens (alanine ddG),
predict an unseen antigen (Dengue envelope, 1OKE) scored against the orthogonal shotgun
n/N label. External hold-out + cross-label-type. The capstone of the physics-feature program.

Transferable feature set (exists for BOTH Tier-1 and 1OKE -- NOT rmsf/slow-modes, which 1OKE
lacks until its FRESEAN run lands): intEn_{total,elec,vdw}, hb_{intra,water}, bb_dih_std,
chi1_std, coupling_{abs,pos} (DCCM), mlce_coupling.

Result (2026-07-07): every physics-feature set transfers at chance (0.43-0.52); only exposure
(crystal SASA) transfers (AUROC 0.58, graded rho +0.15).
"""
import pandas as pd, numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
from scipy.stats import spearmanr

FEAT = ["intEn_total_z","intEn_elec_z","intEn_vdw_z","hb_intra_z","hb_water_z",
        "bb_dih_std_z","chi1_std_z","coupling_abs_z","coupling_pos_z","mlce_coupling_z"]

def load(fd, lbl):
    en   = pd.read_csv(f"{fd}/grinn_intenergy_apomd.csv"); en["pdb"]=en["pdb"].str.upper()
    hb   = pd.read_csv(f"{fd}/hbonds_apomd.csv");          hb["pdb"]=hb["pdb"].str.upper()
    dih  = pd.read_csv(f"{fd}/dihedral_flex_apomd.csv");   dih["pdb"]=dih["pdb"].str.upper()
    dccm = pd.read_csv(f"{fd}/dccm_apomd.csv");            dccm["pdb"]=dccm["pdb"].str.upper()
    mlce = pd.read_csv(f"{fd}/mlce_coupling.csv");         mlce["pdb"]=mlce["pdb"].str.upper()
    m = en[["pdb","resid","intEn_total_z","intEn_elec_z","intEn_vdw_z"]]
    for df,cols in [(hb,["hb_intra_z","hb_water_z"]),(dih,["bb_dih_std_z","chi1_std_z"]),
                    (dccm,["coupling_abs_z","coupling_pos_z"]),(mlce,["mlce_coupling_z"])]:
        m = m.merge(df[["pdb","resid"]+cols], on=["pdb","resid"], how="left")
    return m

# ... (surface restriction via relSASA/crystal SASA; fit Tier-1, predict 1OKE; see HANDOFF.md)
if __name__ == "__main__":
    print("See dynamics_modes.tex Fig. external_dengue_test.png for the result.")
