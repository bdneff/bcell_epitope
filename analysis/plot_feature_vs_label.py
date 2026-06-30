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

# Manuscript-grade defaults: large fonts, clean axes.
plt.rcParams.update({
    "font.size": 15, "axes.titlesize": 15, "axes.labelsize": 17,
    "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 14,
    "axes.linewidth": 1.1, "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 150, "savefig.bbox": "tight", "font.family": "DejaVu Sans",
})
GREY = "#c2c2c2"      # not important
RED = "#c1272d"       # experimentally important
SHADE = "#ffe2b0"     # "SASA would call epitope" region

LABEL_CSV = "benchmark/antigen_alanine_scan_extracted_SIMPLE_v3.csv"
# antigen display key (matches feature 'key' column)
KEYS = {"Lysozyme":"lysozyme","HPr":"hpr","VEGF":"vegf","Bont/A1":"bonta1",
        "IFN-gamma receptor":"ifngr","Tissue factor":"tf",
        "HCMV glycoprotein B":"hcmvgb","human growth hormone":"hgh"}
# antigen key -> PDB code (shown on the plots instead of the antigen name)
KEY2PDB = {"lysozyme":"1AKI","hpr":"2JEL","vegf":"1BJ1","bonta1":"2NYY",
           "ifngr":"1JRH","tf":"1AHW","hcmvgb":"5C6T","hgh":"1HGU"}


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
    ap.add_argument("--feature-threshold", type=float, default=0.25,
                    help="feature value of the dashed decision line (default 0.25 = standard "
                         "relative-SASA 'exposed' cutoff: what exposure alone would call epitope)")
    ap.add_argument("--outdir", default="analysis/figures")
    args = ap.parse_args()

    labels = load_labels()
    feat = load_feature(args.feature_csv, args.feature_col)
    outdir = pathlib.Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    thr = args.feature_threshold
    keys = [k for k in KEYS.values() if k in labels and k in feat]
    ncol = 4; nrow = -(-len(keys) // ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 3.8 * nrow), squeeze=False)
    # pooled confusion + scanned-only correlation accumulators
    P_pred = P_tp = P_imp = P_res = 0
    pooled_xs, pooled_ys = [], []
    print(f"\n{'PDB':6s} {'#res':>5s} {'pred+':>6s} {'TP':>4s} {'#imp':>5s} "
          f"{'prec':>6s} {'recall':>7s} {'Spearman':>9s} {'AUC':>6s}")
    print("-" * 70)
    for i, k in enumerate(keys):
        ax = axes[i // ncol][i % ncol]
        lab = labels[k]; fe = feat[k]
        resids = sorted(fe)                                   # ALL residues of the antigen
        x = np.array([fe[r] for r in resids])
        # y = experimental importance where known, else 0 (unscanned -> assumed not important)
        y = np.array([lab.get(r, 0.0) for r in resids])
        important = np.array([(r in lab and lab[r] >= args.cutoff) for r in resids])
        pred_pos = x >= thr                                  # what SASA alone would label epitope
        n_res = len(resids); n_imp = int(important.sum())
        n_pred = int(pred_pos.sum()); n_tp = int((pred_pos & important).sum())
        prec = n_tp / n_pred if n_pred else float("nan")
        recall = n_tp / n_imp if n_imp else float("nan")
        # scanned-only correlation (does the feature rank graded importance?)
        sx = np.array([fe[r] for r in lab if r in fe]); sy = np.array([lab[r] for r in lab if r in fe])
        rho = spearmanr(sx, sy).correlation if len(sx) > 2 else float("nan")
        A = auc([fe[r] for r in lab if lab[r] >= args.cutoff and r in fe],
                [fe[r] for r in fe if not (r in lab and lab[r] >= args.cutoff)])
        pooled_xs += list(sx); pooled_ys += list(sy)
        P_pred += n_pred; P_tp += n_tp; P_imp += n_imp; P_res += n_res
        print(f"{KEY2PDB[k]:6s} {n_res:5d} {n_pred:6d} {n_tp:4d} {n_imp:5d} "
              f"{prec:6.2f} {recall:7.2f} {rho:9.3f} {A:6.2f}")
        # plot ALL residues: grey = not important, red = experimentally important
        ax.axvspan(thr, max(x.max(), thr) + 0.02, color=SHADE, alpha=0.5, lw=0)
        ax.scatter(x[~important], y[~important], c=GREY, s=20, alpha=0.7, edgecolors="none",
                   label="not important")
        ax.scatter(x[important], y[important], c=RED, s=46, alpha=0.95, edgecolors="k",
                   linewidths=0.4, label="important (binding)", zorder=3)
        ax.axvline(thr, color="#e8820c", ls="--", lw=1.8)
        ax.set_title(f"{KEY2PDB[k]}   precision {prec:.0%} · recall {recall:.0%}",
                     fontsize=14, fontweight="bold")
        ax.set_xlabel(args.label, fontsize=13)
        ax.set_ylabel(r"$\Delta\Delta G$ (kcal/mol)", fontsize=13)
        ax.tick_params(labelsize=12)
    for j in range(len(keys), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    # single figure-level legend BELOW the grid (outside all panels -> hides no data)
    handles, labs = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labs, loc="lower center", ncol=2, fontsize=14, frameon=True,
               markerscale=1.4, bbox_to_anchor=(0.5, 0.005))
    fig.suptitle(f"Static surface exposure as an epitope label (apo, single-frame): every residue by "
                 f"{args.label};\nred = experimentally important, shaded = SASA$\\geq${thr} "
                 f"(what exposure alone would call epitope)", fontsize=16, fontweight="bold")
    fig.tight_layout(rect=[0, 0.05, 1, 0.95])
    p1 = outdir / "sasa_allres_per_antigen.png"
    fig.savefig(p1, dpi=200); plt.close(fig)

    # pooled: all residues, same decision line -- the manuscript figure
    fig2, ax = plt.subplots(figsize=(8.5, 6.0))
    allx, ally, allimp = [], [], []
    for k in keys:
        for r in sorted(feat[k]):
            allx.append(feat[k][r]); ally.append(labels[k].get(r, 0.0))
            allimp.append(r in labels[k] and labels[k][r] >= args.cutoff)
    allx, ally, allimp = np.array(allx), np.array(ally), np.array(allimp)
    ax.axvspan(thr, allx.max() + 0.03, color=SHADE, alpha=0.55, lw=0)
    ax.scatter(allx[~allimp], ally[~allimp], c=GREY, s=16, alpha=0.55, edgecolors="none",
               label="not important for binding")
    ax.scatter(allx[allimp], ally[allimp], c=RED, s=42, alpha=0.95, edgecolors="k",
               linewidths=0.4, label="important for binding", zorder=3)
    ax.axvline(thr, color="#e8820c", ls="--", lw=2.2)
    prec_all = P_tp / P_pred if P_pred else float("nan")
    recall_all = P_tp / P_imp if P_imp else float("nan")
    rho_all = spearmanr(np.array(pooled_xs), np.array(pooled_ys)).correlation
    ax.set_xlabel(f"{args.label}  →  more solvent-exposed", fontsize=18)
    ax.set_ylabel(r"$\Delta\Delta G_\mathrm{bind}$ (kcal/mol)", fontsize=18)
    ax.set_title("If surface exposure were the epitope label", fontsize=19, fontweight="bold", pad=12)
    # annotation box with the confusion summary (top-left, clear of data)
    txt = (f"{P_pred:,} / {P_res:,} residues flagged epitope\n"
           f"by SASA $\\geq$ {thr}; only {P_tp} truly important\n"
           f"$\\Rightarrow$ precision {prec_all:.0%}    recall {recall_all:.0%}")
    ax.text(0.035, 0.97, txt, transform=ax.transAxes, ha="left", va="top", fontsize=14.5,
            bbox=dict(boxstyle="round,pad=0.5", fc="white", ec="#999", alpha=0.97))
    ax.legend(fontsize=13.5, loc="upper left", bbox_to_anchor=(1.01, 1.0), framealpha=1.0,
              markerscale=1.3, title="shaded = SASA $\\geq$ 0.25\n(exposure's epitope call)",
              title_fontsize=12, borderaxespad=0)
    ax.margins(x=0.02)
    fig2.tight_layout(); p2 = outdir / "sasa_allres_pooled.png"
    fig2.savefig(p2, dpi=220); plt.close(fig2)

    print("-" * 70)
    print(f"POOLED: {P_res} residues; SASA>={thr} labels {P_pred} epitope, {P_tp} truly important "
          f"-> precision {prec_all:.1%}, recall {recall_all:.1%}; scanned Spearman rho={rho_all:.3f}")
    print(f"wrote {p1}\n      {p2}")


if __name__ == "__main__":
    main()
