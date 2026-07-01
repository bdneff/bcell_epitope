#!/usr/bin/env python
"""First dynamics feature: per-residue apo-MD RMSF vs alanine-scan ddG_bind.

The trajectory analogue of the static B-factor baseline (recompute_apo_baselines.py),
but from real 100 ns apo-antigen MD instead of crystallographic B-factors. Asks the
same question the B-factor teaser did -- does mobility track binding-hotspot ddG? --
now with genuine dynamics.

Design decisions honored:
  - Per antibody-EPITOPE labels (Design #1 override): ddG kept as a 2D [residue x antibody]
    structure, one column per (antigen, antibody). Never averaged over antibodies -- HEL's
    D1.3 sits on the opposite face from HyHEL-63/10 and a mean would smear two real epitopes.
  - Leakage control (Design #3): correlation reported per-epitope, pooled, and leave-one-
    ANTIGEN-out (LOAO). Split is by antigen, not antibody, because HEL's 3 antibodies share
    one apo trajectory -- an antibody split would leak features across train/test.
  - Apo features (Design #4): features come from the ligand-free antigen only.

Numbering alignment: the MD structure (apo PDB) is numbered differently from the
antibody-complex PDBs the labels come from. We do NOT assume a fixed offset -- we verify
alignment by matching each scanned residue's WT identity against the residue actually present
in the MD structure at that position, and only accept a system if identities agree.
The six Tier-1 systems below all align at delta=0. Two antigens are deliberately absent (and are
NOT the same protein): MT-SP1 (labels 3BN9/3NPS) has Tier-1 ddG labels but no apo MD; Dengue E
(1OKE) has an apo MD trajectory but Tier-2 binary labels -- scored separately on the binary track,
not here.

Feature extraction (run on Gemini, gpu-v100-dev, GROMACS 2023.2-Container):
  gmx trjconv -s md.tpr -f md.xtc -o whole.xtc -pbc mol -center -b 10000
  gmx rmsf   -s md.tpr -f whole.xtc -res -b 10000 -o rmsf_<PDB>.xvg     # CA, per residue
First 10 ns discarded as equilibration; RMSF is least-squares-fit (removes global tumbling);
z-normalised per structure to match bfactor_z.

Run:  python analysis/rmsf_dynamics_baseline.py \
          --xvg-dir hpc/<jobid> \
          --md-resnames handoff/mdres.txt
Writes benchmark/features/rmsf_apomd.csv and manuscript/figures/rmsf_baseline.png.
Recorded in manuscript/figures/FIGURES.md.
"""
import argparse, csv, collections, glob, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

AA3TO1 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G',
          'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
          'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
# antigen display name -> apo-MD PDB code (the simulated structure)
ANTIGEN2PDB = {"Lysozyme":"1AKI", "HPr":"2JEL", "VEGF":"1BJ1", "IFN-gamma receptor":"1JRH",
               "Tissue factor":"1AHW", "human growth hormone":"1HGU"}
DISP = {"1AKI":"Lysozyme","2JEL":"HPr","1BJ1":"VEGF","1JRH":"IFN-\u03b3R",
        "1AHW":"Tissue factor","1HGU":"hGH"}
LABEL_CSV = "benchmark/antigen_alanine_scan_extracted_SIMPLE_v3.csv"
DDG_IMPORTANT = 1.0     # kcal/mol; matches baseline "important" threshold


def read_xvg(fn):
    d = {}
    for line in open(fn):
        if line[:1] in "@#":
            continue
        p = line.split()
        if len(p) >= 2:
            d[int(p[0])] = float(p[1])
    return d


def load_md_resnames(path):
    """handoff file of '===PDB===' blocks then 'chain:resid:RESNAME' tokens -> {pdb:{resid:1letter}}."""
    out = collections.defaultdict(dict)
    cur = None
    for tok in open(path).read().split():
        if tok.startswith("==="):
            cur = re.match(r"===(\w+)", tok).group(1)
            continue
        parts = tok.split(":")
        if len(parts) == 3 and cur:
            try:
                resid = int(parts[1])
            except ValueError:
                continue
            if parts[2] in AA3TO1 and resid not in out[cur]:
                out[cur][resid] = AA3TO1[parts[2]]
    return out


def load_epitopes():
    """(pdb, antibody) -> {resid: ddG}, and (pdb) -> [(resid, WT1)] for alignment check."""
    epi = collections.defaultdict(dict)
    wt = collections.defaultdict(list)
    for r in csv.DictReader(open(LABEL_CSV)):
        if not r["Use for MD correlation?"].startswith("YES"):
            continue
        ag = r["Antigen"].strip()
        if ag not in ANTIGEN2PDB:
            continue
        pdb = ANTIGEN2PDB[ag]
        try:
            resid = int(r["Residue #"]); ddg = float(r["Best \u0394\u0394G kcal/mol"])
        except ValueError:
            continue
        epi[(pdb, r["Antibody"].strip())][resid] = ddg
        w = r["WT"].strip().upper()
        if len(w) == 1:
            wt[pdb].append((resid, w))
    return epi, wt


def aligned_systems(wt, mdres, min_frac=0.95):
    """Return set of PDBs whose scanned WT identities match the MD structure (delta=0)."""
    ok = set()
    for pdb, pairs in wt.items():
        seq = mdres.get(pdb, {})
        if not seq or not pairs:
            continue
        m = sum(1 for resid, w in pairs if seq.get(resid) == w)
        if m >= min_frac * len(pairs):
            ok.add(pdb)
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xvg-dir", required=True, help="dir with rmsf_<PDB>.xvg from Gemini")
    ap.add_argument("--md-resnames", required=True, help="handoff mdres.txt (resid->resname per system)")
    ap.add_argument("--out-csv", default="benchmark/features/rmsf_apomd.csv")
    ap.add_argument("--out-fig", default="manuscript/figures/rmsf_baseline.png")
    args = ap.parse_args()

    rmsf = {os.path.basename(f).split("_")[1].split(".")[0]: read_xvg(f)
            for f in glob.glob(os.path.join(args.xvg_dir, "rmsf_*.xvg"))}
    mdres = load_md_resnames(args.md_resnames)
    epi, wt = load_epitopes()
    ok = aligned_systems(wt, mdres)
    dropped = sorted(set(ANTIGEN2PDB.values()) & set(rmsf) - ok)

    # z-score RMSF per structure; write feature CSV (schema matches bfactor_apocrystal.csv)
    zz = {}
    rows = []
    for pdb in sorted(set(rmsf) & ok):
        ids = sorted(rmsf[pdb]); v = np.array([rmsf[pdb][i] for i in ids])
        z = (v - v.mean()) / v.std(ddof=0)
        zz[pdb] = {i: float(x) for i, x in zip(ids, z)}
        for i, zval, vval in zip(ids, z, v):
            rows.append({"antigen": DISP[pdb], "key": pdb.lower(), "pdb": pdb, "chain": "A",
                         "resid": i, "resname": mdres[pdb].get(i, "UNK"),
                         "rmsf_nm": round(vval, 4), "rmsf_z": round(float(zval), 4)})
    os.makedirs(os.path.dirname(args.out_csv), exist_ok=True)
    with open(args.out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["antigen","key","pdb","chain","resid","resname","rmsf_nm","rmsf_z"])
        w.writeheader(); w.writerows(rows)

    # correlations: per-epitope, pooled, LOAO (by antigen)
    records = []   # (pdb, ab, resid, rmsf_z, ddG)
    per = []
    for (pdb, ab), dd in epi.items():
        if pdb not in ok:
            continue
        xs, ys = [], []
        for resid, ddg in dd.items():
            if resid in zz[pdb]:
                xs.append(zz[pdb][resid]); ys.append(ddg)
                records.append((pdb, ab, resid, zz[pdb][resid], ddg))
        if len(xs) >= 4:
            per.append({"pdb": pdb, "ag": DISP[pdb], "ab": ab, "n": len(xs), "rho": spearmanr(xs, ys)[0]})
    per.sort(key=lambda d: d["rho"])
    X = np.array([r[3] for r in records]); Y = np.array([r[4] for r in records])
    pdbs = [r[0] for r in records]
    pooled_rho, pooled_p = spearmanr(X, Y)
    ag_rho = {}
    for pdb in ok:
        xs = [X[i] for i in range(len(X)) if pdbs[i] == pdb]
        ys = [Y[i] for i in range(len(X)) if pdbs[i] == pdb]
        if len(xs) >= 5:
            ag_rho[pdb] = spearmanr(xs, ys)[0]
    loao = float(np.nanmean(list(ag_rho.values())))

    _plot(per, X, Y, pooled_rho, pooled_p, loao, dropped, args.out_fig)
    print(f"systems used: {sorted(ok)}  dropped(numbering): {dropped}")
    print(f"pooled rho={pooled_rho:.3f} p={pooled_p:.3f}  LOAO mean rho={loao:.3f}")
    print(f"wrote {args.out_csv} ({len(rows)} residues) and {args.out_fig}")


def _plot(per, X, Y, pooled_rho, pooled_p, loao, dropped, out_fig):
    plt.rcParams.update({"font.size": 12, "axes.spines.top": False, "axes.spines.right": False,
                         "savefig.bbox": "tight", "figure.dpi": 150})
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12, 5.4), gridspec_kw={"width_ratios": [1.15, 1]})
    ags = sorted(set(d["ag"] for d in per))
    cmap = dict(zip(ags, plt.cm.tab10(np.linspace(0, 1, 10))))
    for i, d in enumerate(per):
        c = cmap[d["ag"]]
        axL.plot([0, d["rho"]], [i, i], color=c, lw=2, zorder=1)
        axL.scatter(d["rho"], i, s=90, color=c, zorder=2, edgecolor="white", lw=0.8)
        axL.text(d["rho"] + (0.03 if d["rho"] >= 0 else -0.03), i,
                 f'{d["ag"]} \u00b7 {d["ab"].split()[0][:8]} (n={d["n"]})',
                 va="center", ha="left" if d["rho"] >= 0 else "right", fontsize=8)
    axL.axvline(0, color="0.5", lw=1)
    axL.axvline(pooled_rho, color="#c1272d", ls="--", lw=1.6, zorder=0)
    axL.set_xlim(-1.25, 1.25); axL.set_ylim(-0.7, len(per) - 0.3); axL.set_yticks([])
    axL.set_xlabel("Spearman \u03c1  (apo RMSF vs \u0394\u0394G$_{bind}$)")
    axL.set_title("Per antibody-epitope: flexibility does not track hotspots", fontsize=11, loc="left")
    axL.text(pooled_rho - 0.03, len(per) - 0.5, f"pooled \u03c1={pooled_rho:.2f}",
             color="#c1272d", ha="right", fontsize=8.5, style="italic")

    imp = Y >= DDG_IMPORTANT
    axR.axvspan(0, X.max() + 0.3, color="#ffe2b0", alpha=0.6, zorder=0, label="RMSF \u2265 0 (predicted epitope)")
    axR.scatter(X[~imp], Y[~imp], s=42, color="#4878b8", edgecolor="white", lw=0.5,
                label="scanned, weak (\u0394\u0394G<1)", zorder=2)
    axR.scatter(X[imp], Y[imp], s=52, color="#c1272d", edgecolor="white", lw=0.5,
                label="important (\u0394\u0394G\u22651)", zorder=3)
    axR.axvline(0, color="orange", ls="--", lw=1.6, zorder=1)
    axR.set_xlim(X.min() - 0.3, X.max() + 0.3)
    axR.set_xlabel("apo MD RMSF (z)  \u2192  more mobile"); axR.set_ylabel("\u0394\u0394G$_{bind}$ (kcal/mol)")
    n_flag = int((X >= 0).sum()); tp = int((imp & (X >= 0)).sum())
    axR.set_title(f"If MD mobility were the epitope label\nRMSF\u22650 flags {n_flag}/{len(X)} scanned; "
                  f"{tp}/{int(imp.sum())} important caught", fontsize=10.5, loc="left")
    axR.legend(frameon=False, fontsize=7.5, loc="upper right")

    drop_note = f"   [{'/'.join(dropped)} excluded: numbering unaligned]" if dropped else ""
    fig.suptitle("Real apo-MD flexibility (RMSF) mirrors the B-factor baseline: no signal beyond exposure",
                 fontsize=12.5, fontweight="bold", y=1.02)
    fig.text(0.5, 0.965,
             f"pooled \u03c1={pooled_rho:.2f} (p={pooled_p:.2f}) \u00b7 LOAO mean \u03c1={loao:.2f}{drop_note}",
             ha="center", fontsize=8.5, color="0.35")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out_fig, dpi=150, bbox_inches="tight")


if __name__ == "__main__":
    main()
