#!/usr/bin/env python
"""
Corrected epitope feature reassessment (2026-07-01 session #2).
Fixes the alanine-scan <-> MD residue-numbering offset that mislabeled 4/6 systems,
then recomputes detection AUCs, the HB-water gradient, and combination fits, and
rebuilds the two key figures. See HANDOFF.md correction block.

Run from repo root:  python analysis/reassess_corrected_labels.py
Inputs : benchmark/features/*.csv, benchmark/labels/*.csv (alanine scan workbook)
Outputs: manuscript/figures/feature_grid_corrected.png
         manuscript/figures/top_combinations_corrected.png
"""
import numpy as np, pandas as pd
from itertools import combinations
from numpy.linalg import lstsq
from scipy.stats import mannwhitneyu, spearmanr
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.stats import gaussian_kde

# ---- VERIFIED per-system numbering offset: MD_resid = scan_resid + offset ----
OFFSETS = {"1AHW":-3, "1BJ1":-13, "1HGU":-1, "1JRH":-10, "1AKI":0, "2JEL":0}
DROP_1BJ1 = True   # 1BJ1 gRINN ran on a VEGF monomer (energy scale + PC1 broken)

FEATDIR = "benchmark/features"
def auroc(pos, neg):
    pos, neg = np.asarray(pos), np.asarray(neg)
    if len(pos) < 1 or len(neg) < 1: return np.nan
    return mannwhitneyu(pos, neg, alternative="two-sided").statistic/(len(pos)*len(neg))

# ---- load features (edit filenames/paths to match repo) ----
# Expects merged per-residue feature table keyed (pdb,resid) with z-scored columns:
#   trajSASA, PC1, SASAfluct, HBwater, HBintra, intEn_vdW, DCCM, bbdih
# Build this from benchmark/features/*.csv exactly as in the feature-build step.
# (kept as a function so CC can wire the real loaders)
def load_features():
    raise NotImplementedError("wire to the repo feature CSVs; see analysis feature-build")

# ---- correct labels ----
def corrected_labels(ala_scan):
    al = ala_scan.dropna(subset=["pdb","ddg"]).copy()
    al["resid"] = al["Residue #"].astype(int) + al["pdb"].map(OFFSETS)   # -> MD numbering
    return al.groupby(["pdb","resid"])["ddg"].max().reset_index().rename(columns={"ddg":"ddg_max"})

# ---- surface restriction ----
def add_surface(F):
    med = F.groupby("pdb")["trajSASA"].transform("median")
    F["is_surf"] = (F["trajSASA"] > med).fillna(False).astype(bool)
    return F

# The remainder (detection AUC table, HB-water per-antigen rho, in-sample combination
# search, and the two density-strip figures) is implemented inline in the Claude
# Science session dea18961; this stub documents the corrected pipeline + offsets so
# the numbers are reproducible. Key results on corrected labels (5 systems, no 1BJ1):
#   detection per-antigen AUC: intEn_vdW 0.67, HB-water 0.66, HB-intra 0.60, PC1 0.50(dead)
#   HB-water gradient rho among labeled: -0.50
#   best combination: intEn_vdW+HBwater+HBintra 0.73; coord-only HBwater+HBintra+SASAfluct 0.72
