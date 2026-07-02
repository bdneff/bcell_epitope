#!/usr/bin/env python
"""
Corrected epitope-feature reassessment — RUNNABLE (replaces the 2026-07-01 stub).

Regenerates the detection-AUROC leaderboard, the HB-water gradient, the in-sample
combination search, and the two figures directly from the committed feature CSVs, so
the numbers are reproducible and locked with asserts.

Fixes baked in (see HANDOFF.md CRITICAL CORRECTION block):
  1. label numbering offset   (numbering_offsets.json; MD_resid = scan_resid + offset)
  2. gRINN 0-index off-by-one (grinn resid += 1, applied in exactly ONE place, at load)
  3. 1BJ1 dropped             (gRINN ran on a VEGF monomer -> energy scale + PC1 broken)
Validation: WT identity of every labeled residue is asserted against the MD structure
(fails loud if inputs change). Surface restriction uses nan-safe per-pdb median.

Run from repo root:  micromamba run -n bcell-repair python analysis/reassess_corrected_labels.py
Inputs : benchmark/features/*.csv, benchmark/labels/numbering_offsets.json,
         benchmark/antigen_alanine_scan_extracted_SIMPLE_v3.csv
Outputs: manuscript/figures/feature_grid_corrected.png, top_combinations_corrected.png
"""
import json, itertools
import numpy as np, pandas as pd
from numpy.linalg import lstsq
from scipy.stats import spearmanr
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

FEATDIR = "benchmark/features"
OFFS = json.load(open("benchmark/labels/numbering_offsets.json"))["offsets"]
ANT2PDB = {"Lysozyme":"1AKI","HPr":"2JEL","VEGF":"1BJ1","IFN-gamma receptor":"1JRH",
           "Tissue factor":"1AHW","human growth hormone":"1HGU"}
DROP = {"1BJ1"}                       # VEGF monomer gRINN bug
AA3TO1 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G',
          'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
          'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
# feature col -> display name (leaderboard order filled after compute)
FEATURES = [("hb_water_z","H-bonds to water"), ("hb_intra_z","H-bonds within protein"),
            ("relsasa_md_z","Trajectory SASA"), ("intEn_vdw_z","Interaction energy (vdW)"),
            ("coupling_abs_z","DCCM coupling"), ("bb_dih_std_z","Backbone flexibility"),
            ("mode1_z","PC1 participation"), ("sasa_sd_z","SASA fluctuation")]


def _read(fname, cols, shift=0):
    df = pd.read_csv(f"{FEATDIR}/{fname}")
    keep = ["pdb", "resid"] + [c for c in cols if c in df.columns]
    df = df[keep].copy()
    df["resid"] = df["resid"].astype(int) + shift        # <-- gRINN +1 lives here, nowhere else
    return df.groupby(["pdb", "resid"], as_index=False).mean(numeric_only=True)


def load_features():
    """Merge the committed per-residue feature CSVs to one (pdb,resid) table."""
    F = _read("grinn_intenergy_apomd.csv", ["intEn_vdw_z", "intEn_total_z"], shift=1)  # 0->1 indexed
    for fn, cols in [("hbonds_apomd.csv", ["hb_water_z", "hb_intra_z"]),
                     ("dynamics_sasa_modes.csv", ["relsasa_md_z", "sasa_sd_z", "mode1_z"]),
                     ("dccm_apomd.csv", ["coupling_abs_z", "coupling_pos_z"]),
                     ("dihedral_flex_apomd.csv", ["bb_dih_std_z", "chi1_std_z"])]:
        add = _read(fn, cols)
        F = F.merge(add, on=["pdb", "resid"], how="outer")
    return F


def resname_map():
    """(pdb,resid) -> 3-letter resname, from a feature CSV that carries it (MD 1-indexed)."""
    d = pd.read_csv(f"{FEATDIR}/dccm_apomd.csv")[["pdb", "resid", "resname"]].drop_duplicates(["pdb", "resid"])
    return {(r.pdb, int(r.resid)): r.resname for r in d.itertuples()}


def corrected_labels():
    """Alanine-scan ddG -> MD numbering via the verified offsets; returns per-(pdb,resid) max ddG
    plus the WT-identity check rows (scan WT vs MD resname)."""
    s = pd.read_csv("benchmark/antigen_alanine_scan_extracted_SIMPLE_v3.csv")
    s = s[s["Use for MD correlation?"].astype(str).str.startswith("YES")].copy()
    s = s[s["Antigen"].str.strip().isin(ANT2PDB)]
    s["pdb"] = s["Antigen"].str.strip().map(ANT2PDB)
    s["resid"] = s["Residue #"].astype(int) + s["pdb"].map(OFFS)
    s["ddg"] = pd.to_numeric(s["Best ΔΔG kcal/mol"], errors="coerce")
    s = s.dropna(subset=["ddg"])
    s["wt"] = s["WT"].astype(str).str.strip()
    ddg = s.groupby(["pdb", "resid"])["ddg"].max().reset_index().rename(columns={"ddg": "ddg_max"})
    return s[["pdb", "resid", "wt"]], ddg


def assert_wt_identity(scan_rows, res):
    """Every labeled residue's scan WT must equal the MD resname. Fail loud."""
    ok = tot = 0
    for _, r in scan_rows.iterrows():
        rn = res.get((r["pdb"], int(r["resid"])))
        if rn is None:
            continue
        tot += 1
        ok += (AA3TO1.get(rn, "?") == r["wt"])
    assert tot > 0 and ok == tot, f"WT-identity check FAILED: {ok}/{tot} match (inputs changed?)"
    print(f"  WT-identity check: {ok}/{tot} labeled residues match MD structure  [OK]")


def auroc(pos, neg):
    pos = np.asarray([x for x in pos if np.isfinite(x)])
    neg = np.asarray([x for x in neg if np.isfinite(x)])
    if not len(pos) or not len(neg):
        return np.nan
    a = np.concatenate([pos, neg]); o = np.argsort(a, kind="mergesort"); sv = a[o]
    rk = np.empty(len(a)); i = 0
    while i < len(sv):
        j = i
        while j + 1 < len(sv) and sv[j + 1] == sv[i]:
            j += 1
        rk[o[i:j + 1]] = (i + j) / 2 + 1; i = j + 1
    return (rk[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def add_surface(F):
    """Surface = trajSASA above the nan-safe per-pdb median."""
    med = F.groupby("pdb")["relsasa_md_z"].transform("median")   # pandas median skips NaN
    F["is_surf"] = F["relsasa_md_z"] > med
    return F


def detection_auc(F, col):
    """Per-antigen surface-restricted directional AUROC = P(feat_hot > feat_rest), averaged over
    the systems where it is defined (5 systems, 1BJ1 dropped)."""
    per = {}
    for pdb, g in F.groupby("pdb"):
        if pdb in DROP:
            continue
        g = g[g["is_surf"]]
        a = auroc(g.loc[g["hot"], col], g.loc[~g["hot"], col])
        if np.isfinite(a):
            per[pdb] = a
    return (np.mean(list(per.values())) if per else np.nan), per


def main():
    F = load_features()
    scan_rows, ddg = corrected_labels()
    assert_wt_identity(scan_rows, resname_map())
    F = F.merge(ddg, on=["pdb", "resid"], how="left")
    F["hot"] = (F["ddg_max"] >= 1.0).fillna(False)
    F = add_surface(F)

    # ---- detection leaderboard ----
    board = []
    for col, name in FEATURES:
        if col not in F.columns:
            continue
        m, per = detection_auc(F, col)
        board.append((name, col, m, per))
    board.sort(key=lambda r: -r[2])
    print("\nDetection per-antigen AUROC (surface, 5 systems, 1BJ1 dropped):")
    for name, col, m, per in board:
        print(f"  {name:26s} {m:.2f}   per-antigen={ {k:round(v,2) for k,v in per.items()} }")

    # ---- HB-water gradient among ALL labeled epitope residues (no gRINN, NO surface
    #      restriction -- surface-restricting it wrongly drops 2JEL, the +0.27 dissenter,
    #      and inflated the earlier -0.50; the honest per-antigen-mean value is -0.33). ----
    lab = F[F["ddg_max"].notna() & (~F["pdb"].isin(DROP))].dropna(subset=["hb_water_z", "ddg_max"])
    per_rho = {p: spearmanr(gg["hb_water_z"], gg["ddg_max"]).correlation
               for p, gg in lab.groupby("pdb") if len(gg) > 3}
    rho = float(np.nanmean(list(per_rho.values())))
    print(f"\nHB-water gradient (per-antigen mean, all labeled): rho = {rho:.2f}"
          f"   per-antigen={ {k: round(v, 2) for k, v in per_rho.items()} }")

    # ---- in-sample combination search (lstsq fit of hot ~ features on surface residues) ----
    feat_cols = [c for c, _ in FEATURES if c in F.columns]
    S = F[F["is_surf"] & (~F["pdb"].isin(DROP))].dropna(subset=feat_cols + ["hot"])
    y = S["hot"].astype(float).values
    combos = []
    for k in (2, 3):
        for sub in itertools.combinations(feat_cols, k):
            X = np.column_stack([S[c].values for c in sub] + [np.ones(len(S))])
            w, *_ = lstsq(X, y, rcond=None)
            score = X @ w
            a = auroc(score[y == 1], score[y == 0])
            combos.append((a, sub))
    combos.sort(key=lambda t: -t[0])
    print("\nTop in-sample feature combinations (surface pool):")
    for a, sub in combos[:6]:
        print(f"  {a:.2f}  {' + '.join(dict(FEATURES)[c] if c in dict(FEATURES) else c for c in sub)}")

    make_figures(F, board, combos)

    # ---- LOCK the reconciled numbers ----
    d = {name: m for name, _, m, _ in board}
    assert abs(d.get("Interaction energy (vdW)", np.nan) - 0.57) <= 0.02, \
        f"intEn_vdW {d.get('Interaction energy (vdW)'):.2f} != 0.57 (is gRINN +1 applied?)"
    assert abs(rho - (-0.33)) <= 0.04, f"HB-water gradient {rho:.2f} != -0.33 (all labeled, per-antigen)"
    # HB-water DETECTION: Science kernel reports 0.656 (per-antigen surface AUROC, 33 on-surface
    # labeled, 5 systems, no 1BJ1). This script's surface set matches (33 on-surface labeled) and
    # intEn reproduces at 0.57, but HB-water detection here is ~0.52 (1JRH the divergent system) --
    # left un-hard-asserted pending a direct diff of the AUROC/negative-set vs the Science kernel.
    print(f"\nLocked: intEn_vdW 0.57, HB-water gradient -0.33. "
          f"HB-water detection={d.get('H-bonds to water', float('nan')):.2f} [reconcile vs Science 0.656].")


def make_figures(F, board, combos):
    # feature grid: per-feature hotspot vs rest on surface residues
    n = len(board); ncol = 4; nrow = -(-n // ncol)
    fig, ax = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.2 * nrow), squeeze=False)
    for i, (name, col, m, per) in enumerate(board):
        a = ax[i // ncol][i % ncol]
        S = F[F["is_surf"] & (~F["pdb"].isin(DROP))]
        hot = S.loc[S["hot"], col].dropna(); rest = S.loc[~S["hot"], col].dropna()
        a.hist(rest, bins=20, density=True, alpha=0.5, color="#bbb", label="rest")
        a.hist(hot, bins=12, density=True, alpha=0.7, color="#c1272d", label="hotspot")
        a.set_title(f"{name}\nAUROC {m:.2f}", fontsize=11, fontweight="bold")
        a.tick_params(labelsize=8)
        if i == 0:
            a.legend(fontsize=8, frameon=False)
    for j in range(n, nrow * ncol):
        ax[j // ncol][j % ncol].axis("off")
    fig.suptitle("Corrected feature reassessment (surface residues, 5 systems, 1BJ1 dropped)\n"
                 "red = alanine-scan hotspot ($\\Delta\\Delta G\\geq1$), grey = rest", fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig("manuscript/figures/feature_grid_corrected.png", dpi=180); plt.close(fig)

    fig, a = plt.subplots(figsize=(8, 4.5))
    top = combos[:8][::-1]
    a.barh([" + ".join(dict(FEATURES).get(c, c) for c in sub) for _, sub in top],
           [x for x, _ in top], color="#4878b8")
    a.set_xlabel("in-sample detection AUROC"); a.set_xlim(0.5, max(0.75, max(x for x, _ in top) + 0.02))
    a.set_title("Top feature combinations (in-sample, surface pool)", fontweight="bold")
    fig.tight_layout(); fig.savefig("manuscript/figures/top_combinations_corrected.png", dpi=180); plt.close(fig)
    print("  wrote feature_grid_corrected.png + top_combinations_corrected.png")


if __name__ == "__main__":
    main()
