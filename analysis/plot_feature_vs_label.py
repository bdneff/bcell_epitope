#!/usr/bin/env python
"""Correlate a per-residue feature with the epitope labels -- the plot family Jacob used
(feature on X, label on Y, prediction threshold, confusion counts). General over features:
today SASA (static-exposure baseline), later RMSF / non-bonded energy / hydration.

Two questions, two metrics:
  (1) Among SCANNED residues, does the feature track graded importance?  -> Spearman(feature, ddG).
  (2) Over the WHOLE antigen, does the feature rank the energetic hotspots above the rest?
      -> ROC-AUC with positives = scanned residues with ddG >= cutoff, negatives = every other
      residue (caveat: unscanned residues are assumed non-epitope, which alanine scans largely
      justify since they target the known epitope).
ddG per residue = max over all antibodies scanned against that antigen (is it important for ANY mAb).

  micromamba run -n bcell-repair python analysis/plot_feature_vs_label.py \
      --feature-csv benchmark/features/sasa_apo_singleframe.csv --feature-col rel_sasa \
      --label "rel. SASA (apo)" --cutoff 1.0
"""
import argparse, csv, collections, pathlib
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

LABEL_CSV = "benchmark/antigen_alanine_scan_extracted_SIMPLE_v3.csv"
# antigen display key (matches feature 'key' column)
KEYS = {"Lysozyme":"lysozyme","HPr":"hpr","VEGF":"vegf","Bont/A1":"bonta1",
        "IFN-gamma receptor":"ifngr","Tissue factor":"tf",
        "HCMV glycoprotein B":"hcmvgb","human growth hormone":"hgh"}


def load_labels():
    """antigen-key -> {resid: max ddG over antibodies}."""
    lab = collections.defaultdict(dict)
    for r in csv.DictReader(open(LABEL_CSV)):
        if not r["Use for MD correlation?"].startswith("YES"):
            continue
        ag = r["Antigen"].strip()
        if ag not in KEYS:
            continue
        try:
            resid = int(r["Residue #"]); ddg = float(r["Best ΔΔG kcal/mol"])
        except ValueError:
            continue
        k = KEYS[ag]
        lab[k][resid] = max(lab[k].get(resid, -1e9), ddg)
    return lab


def load_feature(path, col):
    """antigen-key -> {resid: feature} (averaged over chains)."""
    acc = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in csv.DictReader(open(path)):
        try:
            acc[r["key"]][int(r["resid"])].append(float(r[col]))
        except (ValueError, KeyError):
            continue
    return {k: {res: float(np.mean(v)) for res, v in d.items()} for k, d in acc.items()}


def auc(pos, neg):
    pos, neg = np.asarray(pos), np.asarray(neg)
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    # ROC-AUC = P(pos > neg) via rank-sum
    allv = np.concatenate([pos, neg]); ranks = allv.argsort().argsort().astype(float)
    # average ranks for ties
    order = np.argsort(allv)
    sv = allv[order]; r = np.empty(len(allv)); i = 0
    while i < len(sv):
        j = i
        while j + 1 < len(sv) and sv[j + 1] == sv[i]:
            j += 1
        r[order[i:j + 1]] = (i + j) / 2.0 + 1
        i = j + 1
    rp = r[:len(pos)].sum()
    return (rp - len(pos) * (len(pos) + 1) / 2.0) / (len(pos) * len(neg))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--feature-csv", required=True)
    ap.add_argument("--feature-col", required=True)
    ap.add_argument("--label", default="feature")
    ap.add_argument("--cutoff", type=float, default=1.0, help="ddG (kcal/mol) for hotspot-positive")
    ap.add_argument("--outdir", default="analysis/figures")
    args = ap.parse_args()

    labels = load_labels()
    feat = load_feature(args.feature_csv, args.feature_col)
    outdir = pathlib.Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    keys = [k for k in KEYS.values() if k in labels and k in feat]
    ncol = 4; nrow = -(-len(keys) // ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 3.2 * nrow), squeeze=False)
    pooled_x, pooled_y = [], []
    print(f"\n{'antigen':10s} {'n_scan':>6s} {'n_hot':>5s} {'Spearman':>9s} {'AUC(whole)':>10s}")
    print("-" * 48)
    summary = []
    for i, k in enumerate(keys):
        ax = axes[i // ncol][i % ncol]
        lab = labels[k]; fe = feat[k]
        # scanned residues with both feature and label
        xs = np.array([fe[res] for res in lab if res in fe])
        ys = np.array([lab[res] for res in lab if res in fe])
        rho = spearmanr(xs, ys).correlation if len(xs) > 2 else float("nan")
        # whole-protein AUC: positives = scanned hotspots, negatives = all other residues
        hot = {res for res, d in lab.items() if d >= args.cutoff}
        pos = [fe[res] for res in hot if res in fe]
        neg = [fe[res] for res in fe if res not in hot]
        A = auc(pos, neg)
        summary.append((k, len(xs), len(pos), rho, A))
        print(f"{k:10s} {len(xs):6d} {len(pos):5d} {rho:9.3f} {A:10.3f}")
        pooled_x += list(xs); pooled_y += list(ys)
        # plot: feature (X) vs ddG (Y), colored by hotspot
        col = ["#1f3fff" if y >= args.cutoff else "#222" for y in ys]
        ax.scatter(xs, ys, c=col, s=18, alpha=0.8, edgecolors="none")
        ax.axhline(args.cutoff, color="orange", ls="--", lw=1)
        ax.set_title(f"{k}: $\\rho$={rho:.2f}, AUC={A:.2f}\n(n={len(xs)}, hot={len(pos)})", fontsize=9)
        ax.set_xlabel(args.label, fontsize=8); ax.set_ylabel(r"$\Delta\Delta G$ (max over mAb)", fontsize=8)
        ax.tick_params(labelsize=7)
    for j in range(len(keys), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    fig.suptitle(f"Static-exposure control: {args.label} vs alanine-scan $\\Delta\\Delta G$ "
                 f"(blue = hotspot $\\geq${args.cutoff} kcal/mol)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    p1 = outdir / "sasa_vs_ddg_per_antigen.png"
    fig.savefig(p1, dpi=150); plt.close(fig)

    # pooled
    px, py = np.array(pooled_x), np.array(pooled_y)
    rho_all = spearmanr(px, py).correlation
    fig2, ax = plt.subplots(figsize=(5, 4))
    col = ["#1f3fff" if y >= args.cutoff else "#222" for y in py]
    ax.scatter(px, py, c=col, s=14, alpha=0.6, edgecolors="none")
    ax.axhline(args.cutoff, color="orange", ls="--", lw=1)
    ax.set_xlabel(args.label); ax.set_ylabel(r"$\Delta\Delta G$ (max over mAb, kcal/mol)")
    ax.set_title(f"Pooled scanned residues (n={len(px)}): Spearman $\\rho$={rho_all:.3f}")
    fig2.tight_layout(); p2 = outdir / "sasa_vs_ddg_pooled.png"
    fig2.savefig(p2, dpi=150); plt.close(fig2)

    print("-" * 48)
    print(f"POOLED scanned-residue Spearman rho = {rho_all:.3f} (n={len(px)})")
    print(f"wrote {p1}\n      {p2}")


if __name__ == "__main__":
    main()
