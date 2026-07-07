#!/usr/bin/env python
"""
fresean_involvement.py — FRESEAN collective-mode involvement vs epitope labels.

Built with Claude Science (session #4, 2026-07-07).

Computes per-residue involvement in FRESEAN internal modes (7..k) and tests
whether it detects epitope residues, with RMSF and SASA controls, plus a
cumulative-band sweep. Currently wired for 1AKI (the only system with extracted
FRESEAN modes as of this writing); 1OKE/Dengue-E will slot in once its modes
land (see md/1OKE/apo/md_fresean.sh).

INPUTS (per system SYS):
  fresean_extract/<rep>/evec_freq1_mode1-30_cg.xyz   30 modes, 2-bead CG (BACK+SIDE)
  fresean_extract/<rep>/ref-cg.gro                   CG reference coords (nm)
  benchmark/labels/ddg_by_md_residue.csv             per-MD-residue graded ddG (corrected numbering)
  benchmark/features/rmsf_apomd.csv                  apo-MD RMSF (rmsf_nm)
  benchmark/features/dynamics_sasa_modes.csv         trajectory-avg rel SASA (relsasa_md_z)

METHOD:
  involvement_i(7..k) = sum_{m=7}^{k} sum_{beads b in residue i} |e_{mb}|^2
  (eigenvectors are unit-normalized, so total involvement over residues = #modes.)
  Detection = Mann-Whitney AUROC, labeled epitope residues vs the rest.
  Controls: regress RMSF (or SASA) out of involvement, re-test the residual.

KEY 1AKI RESULT (23 alanine-scan labels, single system, N=1 protein):
  involvement(7-30) detection AUROC 0.745 (p=2e-4); SASA control 0.671; RMSF 0.664.
  Survives SASA control (residual 0.66) and RMSF control (residual 0.73).
  Cumulative band peaks at 7-21 (0.777) on a broad 7-14..7-24 plateau (~0.75-0.78);
  do NOT hard-tune the endpoint — it's a within-sample argmax on N=23.
  Epitopes are MORE involved in collective motion (opposite of the old Cartesian-PCA
  'isolation' reading). Graded gradient vs ddG is flat (rho -0.10).
CAVEAT: 1AKI modes vs 1AKI labels only. NOT the Dengue-E epitopeness score.
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr, pearsonr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FEAT = os.path.join(REPO, "benchmark", "features")
LABELS = os.path.join(REPO, "benchmark", "labels")

# ---- mode / structure parsers -------------------------------------------------
def read_modes(path):
    """Return (nmodes, nbeads, 3) eigenvector array from a FRESEAN .xyz mode file."""
    with open(path) as f:
        L = f.readlines()
    n = int(L[0].split()[0]); block = n + 2
    out = []
    for m in range(len(L) // block):
        coords = [list(map(float, L[m*block + 2 + i].split()[1:4])) for i in range(n)]
        out.append(coords)
    return np.array(out)

def parse_cg_gro(path):
    L = open(path).readlines(); n = int(L[1])
    beads = []
    for i in range(n):
        ln = L[2+i]
        beads.append(dict(resid=int(ln[0:5]), resname=ln[5:10].strip(),
                          atom=ln[10:15].strip(),
                          xyz=[float(ln[20:28]), float(ln[28:36]), float(ln[36:44])]))
    return beads

# ---- involvement --------------------------------------------------------------
def involvement_band(modes, beads, lo, hi):
    """Per-residue sum of |e|^2 over modes lo..hi (inclusive, 1-indexed)."""
    resids = np.array([b["resid"] for b in beads]); uniq = sorted(set(resids))
    inv = {r: 0.0 for r in uniq}
    for mi in range(lo-1, hi):
        m2 = (modes[mi]**2).sum(1)
        for k, r in enumerate(resids):
            inv[r] += m2[k]
    return np.array(uniq), np.array([inv[r] for r in uniq])

def det_auroc(vec, labeled_mask):
    l = vec[labeled_mask]; r = vec[~labeled_mask]
    U, p = mannwhitneyu(l, r, alternative="two-sided")
    return U / (len(l) * len(r)), p

def residual_auroc(vec, ctrl, labeled_mask):
    """Detection AUROC of vec after linearly regressing out ctrl."""
    ok = ~np.isnan(ctrl)
    A = np.vstack([ctrl[ok], np.ones(ok.sum())]).T
    coef, *_ = np.linalg.lstsq(A, vec[ok], rcond=None)
    resid = vec[ok] - A @ coef
    lok = labeled_mask[ok]
    U, p = mannwhitneyu(resid[lok], resid[~lok], alternative="two-sided")
    return U / (lok.sum() * (~lok).sum()), p, coef

# ---- rigid-body verification --------------------------------------------------
def rigid_fractions(modes, ref_xyz):
    """Fraction of each mode lying in the 6-dim rigid (trans+rot) subspace."""
    N = ref_xyz.shape[0]; r = ref_xyz - ref_xyz.mean(0)
    T = np.zeros((3, N, 3))
    for a in range(3): T[a, :, a] = 1.0
    R = np.zeros((3, N, 3)); ax = np.eye(3)
    for a in range(3): R[a] = np.cross(ax[a][None, :], r)
    Q, _ = np.linalg.qr(np.concatenate([T, R]).reshape(6, -1).T)
    out = []
    for m in range(len(modes)):
        v = modes[m].reshape(-1); v = v / np.linalg.norm(v)
        out.append(np.linalg.norm(Q @ (Q.T @ v))**2)
    return np.array(out)

def trans_rot_split(modes, ref_xyz, k=6):
    """Per-mode translation-only and rotation-only fractions (separate bases)."""
    N = ref_xyz.shape[0]; r = ref_xyz - ref_xyz.mean(0)
    T = np.zeros((3, N, 3))
    for a in range(3): T[a, :, a] = 1.0
    R = np.zeros((3, N, 3)); ax = np.eye(3)
    for a in range(3): R[a] = np.cross(ax[a][None, :], r)
    QT, _ = np.linalg.qr(T.reshape(3, -1).T)
    QR, _ = np.linalg.qr(R.reshape(3, -1).T)
    rows = []
    for m in range(k):
        v = modes[m].reshape(-1); v = v / np.linalg.norm(v)
        tf = np.linalg.norm(QT @ (QT.T @ v))**2
        rf = np.linalg.norm(QR @ (QR.T @ v))**2
        rows.append((m+1, tf, rf))
    return rows

# ---- label / feature loaders --------------------------------------------------
def load_labels(pdb):
    lab = pd.read_csv(os.path.join(LABELS, "ddg_by_md_residue.csv"))
    a = lab[lab["pdb"].astype(str).str.upper() == pdb.upper()].copy()
    return set(a["resid"]), dict(zip(a["resid"], a["ddg_max"]))

def feature_vec(csv, col, pdb, resid_arr):
    df = pd.read_csv(os.path.join(FEAT, csv))
    d = df[df["pdb"].astype(str).str.upper() == pdb.upper()]
    m = dict(zip(d["resid"], d[col]))
    return np.array([m.get(r, np.nan) for r in resid_arr])

# ---- main analysis for one system --------------------------------------------
def analyze(pdb, mode_xyz, ref_gro, outdir="."):
    modes = read_modes(mode_xyz)
    beads = parse_cg_gro(ref_gro)
    resid_arr, invol = involvement_band(modes, beads, 7, 30)

    lab_resids, ddg_by = load_labels(pdb)
    labeled = np.array([r in lab_resids for r in resid_arr])
    inv_lab, inv_rest = invol[labeled], invol[~labeled]
    auroc, p_mwu = det_auroc(invol, labeled)

    rmsf = feature_vec("rmsf_apomd.csv", "rmsf_nm", pdb, resid_arr)
    sasa = feature_vec("dynamics_sasa_modes.csv", "relsasa_md_z", pdb, resid_arr)
    a_rmsf, p_rmsf = det_auroc(rmsf[~np.isnan(rmsf)],
                               labeled[~np.isnan(rmsf)]) if not np.all(np.isnan(rmsf)) else (np.nan, np.nan)
    a_sasa, p_sasa = det_auroc(sasa[~np.isnan(sasa)],
                               labeled[~np.isnan(sasa)]) if not np.all(np.isnan(sasa)) else (np.nan, np.nan)
    res_r_auroc, res_r_p, _ = residual_auroc(invol, rmsf, labeled)
    res_s_auroc, res_s_p, sasa_coef = residual_auroc(invol, sasa, labeled)

    # graded gradient among labeled
    xg = np.array([invol[list(resid_arr).index(r)] for r in lab_resids if r in resid_arr])
    yg = np.array([ddg_by[r] for r in lab_resids if r in resid_arr])
    rho_g, p_g = spearmanr(xg, yg)

    # cumulative band sweep + single-mode + rigid check
    cum = pd.DataFrame(
        [(f"7-{hi}", hi-6, *det_auroc(involvement_band(modes, beads, 7, hi)[1], labeled))
         for hi in range(8, 31)],
        columns=["band", "nmodes", "auroc", "p"])
    single = pd.DataFrame(
        [(mi, *det_auroc(involvement_band(modes, beads, mi, mi)[1], labeled)) for mi in range(1, 31)],
        columns=["mode", "auroc", "p"])
    ref_cg = np.array([b["xyz"] for b in beads])  # nm; rigid check is scale-free
    rf = rigid_fractions(modes, ref_cg)
    tr = trans_rot_split(modes, ref_cg)

    stats = dict(pdb=pdb, n_labeled=int(labeled.sum()), n_total=len(resid_arr),
                 auroc=auroc, p=p_mwu, sasa_auroc=a_sasa, rmsf_auroc=a_rmsf,
                 invol_sasa_r=float(pearsonr(invol[~np.isnan(sasa)], sasa[~np.isnan(sasa)])[0]) if not np.all(np.isnan(sasa)) else np.nan,
                 invol_rmsf_r=float(pearsonr(invol[~np.isnan(rmsf)], rmsf[~np.isnan(rmsf)])[0]) if not np.all(np.isnan(rmsf)) else np.nan,
                 residual_sasa_auroc=res_s_auroc, residual_rmsf_auroc=res_r_auroc,
                 graded_rho=rho_g, best_band=cum.loc[cum.auroc.idxmax(), "band"],
                 best_band_auroc=float(cum.auroc.max()),
                 rigid_first6=float(rf[:6].mean()), rigid_rest=float(rf[6:].mean()))
    return dict(modes=modes, beads=beads, resid_arr=resid_arr, invol=invol,
                labeled=labeled, inv_lab=inv_lab, inv_rest=inv_rest,
                rmsf=rmsf, sasa=sasa, sasa_coef=sasa_coef, cum=cum, single=single,
                rf=rf, tr=tr, stats=stats,
                a_sasa=a_sasa, p_sasa=p_sasa, res_s_auroc=res_s_auroc)


if __name__ == "__main__":
    import json, glob
    # locate 1AKI modes (rep1) — adjust if run elsewhere
    cand = glob.glob(os.path.join(REPO, "md", "1AKI", "apo", "fresean_extract", "rep1",
                                  "evec_freq1_mode1-30_cg.xyz"))
    if not cand:
        # fall back to a local staging dir
        cand = glob.glob("fresean_1aki/rep1_evec_freq1_mode1-30_cg.xyz")
    mode_xyz = cand[0]
    ref_gro = mode_xyz.replace("evec_freq1_mode1-30_cg.xyz", "ref-cg.gro").replace(
        "rep1_evec_freq1_mode1-30_cg.xyz", "rep1_ref-cg.gro")
    res = analyze("1AKI", mode_xyz, ref_gro)
    print(json.dumps(res["stats"], indent=2, default=float))
