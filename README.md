# bcell_epitope

**Predicting B-cell epitope residues from antigen-only information** — specifically, whether
**apo-antigen dynamics, solvation, and energetics** predict which residues matter for antibody
binding, *beyond* what static surface exposure already explains.

Residue labels are experimental **alanine-scan ΔΔG_bind** (how much each residue contributes to
binding). We simulate the unbound (apo) antigen and ask whether features from the trajectory —
flexibility (RMSF), intra-protein energetics, hydration thermodynamics, collective-mode
participation — rank the binding-important residues better than static baselines.

> The living record is `manuscript/main.pdf` — read it top to bottom; each section is one
> iteration (what we tried, what we found, what it changed). This README is the orientation map.

## The idea

Most B-cell epitope predictors are trained on **proximity/contact labels**: every residue near a
bound antibody is "positive." But proximity ≠ importance — a few hotspots carry most of the binding
free energy while many contacting residues contribute little. (In IEDB, >99% of structure-linked
conformational epitopes are crystal/EM *contact*-derived, not binding-measured.) We instead use
**energetic labels** (ΔΔG) and ask what antigen-intrinsic signal predicts them. The recurring
theme — following Serapian & Colombo, *"The Answer Lies in the Energy"* — is that the predictive
variable is a residue's **energetic/dynamic coupling** to the rest of the protein, not raw
exposure or raw mobility.

## Data — the label benchmark

| Tier | Labels | Source | Coverage |
|------|--------|--------|----------|
| **Tier 1** (graded) | per-residue ΔΔG_bind | AB-Bind + SKEMPI 2.0 (primary papers cited) | **10 antigens** |
| **Tier 2** (binary) | critical / non-critical | IEDB shotgun-mutagenesis (Doranz/Integral Molecular) | flagships: Dengue E, Influenza HA |

Tier-1 antigens (PDB): lysozyme `1AKI`, HPr `2JEL`, VEGF `1BJ1`, MT-SP1 `3BN9/3NPS`, Bont/A1
`2NYY`, Bont/A2 (homology), IFN-γ receptor `1JRH`, tissue factor `1AHW`, HCMV glycoprotein B
`5C6T`, human growth hormone `1HGU`. Tier-2 flagships: Dengue envelope `1OKE`, influenza HA `1HGG`.
Curated data and provenance live in `benchmark/`.

## Repository layout

| Path | Contents |
|------|----------|
| `manuscript/` | Living LaTeX log: `main.tex` + `sections/` + `figures/` (+ `FIGURES.md` provenance, `refs.bib`) |
| `benchmark/` | Curated alanine-scan labels (ΔΔG CSVs + xlsx), IEDB Tier-2 labels, extracted per-residue `features/` |
| `structures/` | Input + cleaned/repaired antigen structures; `structures/repair/` = AlphaFold models + intermediates |
| `md/<PDB>/apo/` | Per-system simulation dirs: frozen `configs/` (mdp) + `prep.sh`/`md.sh`; outputs to `out/` (gitignored) |
| `analysis/` | Feature extraction + correlation/plotting (see below) |
| `env/` | `repair-environment.yml` — host-side conda/micromamba env for structure repair + analysis |
| `docs/` | `METHODS.md` (frozen stack), `ENVIRONMENT_gemini.md` (cluster), `examples/` |
| `CLAUDE.md` | Operating contract (reproducibility discipline, how to work in this repo) |

## Workflow

**Stage and analyze locally (git-tracked); run production MD on Gemini (Slurm); trajectories live
on Gemini and are gitignored.**

1. **Simulate (Gemini).** `git pull`, then submit prep→production chains from the repo root:
   ```bash
   bash md/submit_apo.sh 1JRH 1OKE 1HGG 1HGU 1AHW 2JEL 1BJ1
   ```
   Stack: AMBER99SB-ILDN / TIP3P / cubic box (1.2 nm) / 0.15 M NaCl / 100 ns single trajectory
   (GROMACS, partition `gpu-a100`). Each prep checks the expected disulfide count.

2. **Repair incomplete structures (local).** Crystals with missing loops are completed by
   AlphaFold-loop grafting in the `bcell-repair` env (`env/repair-environment.yml`):
   `analysis/repair_structure.py` → `finalize_structure.py` → `close_loops.py`. Viral gB
   (not in AlphaFold DB) is modeled with **AlphaFold3 on Gemini** (`md/5C6T/apo/af3.sbatch`,
   the lab container) then grafted (`analysis/repair_gB.sh`).

3. **Correlate features with labels (local).** Per-residue feature → epitope-label plots:
   ```bash
   micromamba run -n bcell-repair python analysis/compute_sasa.py
   micromamba run -n bcell-repair python analysis/plot_feature_vs_label.py \
       --feature-csv benchmark/features/sasa_apo_singleframe.csv --feature-col rel_sasa \
       --label "rel. SASA (apo)" --feature-threshold 0.25
   ```
   `plot_feature_vs_label.py` is general over features (SASA today; RMSF / energetics / hydration
   when trajectories land) and reports per-antigen + pooled precision/recall, Spearman, and AUC.

4. **Build the manuscript.**
   ```bash
   cd manuscript && tectonic main.tex     # self-contained; no TeX Live install needed
   ```

## Static baselines (first results)

Before any dynamics, two single-crystal proxies are tested as controls the dynamics features must
beat (manuscript §"Static baselines"):

- **Surface exposure (SASA).** Flags ~half of all residues as epitope but ~3% are truly important:
  **precision ≈ 3%, recall ≈ 63%** pooled over the 8 antigens with apo structures. Exposure locates
  the surface but has no precision and cannot rank energetic importance.
- **Crystallographic mobility (B-factor).** The energetic hotspots are *rigid*, not mobile (pooled
  Spearman **−0.29**) — opposite the naive flexibility-antigenicity expectation, consistent with
  hotspots being pre-organized, weakly-coupled anchors.

Neither static proxy recovers the hotspots from one structure — the motivation for apo-MD.

## Reproducibility

Every figure is regenerated by a named script recorded in `manuscript/figures/FIGURES.md`; every
number traces to a file or command. Labels are read from the actual benchmark files, never invented.
ML evaluation is **leakage-controlled**: leave-one-antigen-out splits (effective N ≈ number of
antigens, not residues). See `CLAUDE.md` for the full operating discipline and `docs/METHODS.md`
for the frozen methods stack.

---
PI: John Altin. Built with reproducibility-first MD discipline; the PDF is the source of truth.
