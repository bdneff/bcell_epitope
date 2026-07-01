# HANDOFF — bcell_epitope dynamics features (read this first)

**Last updated:** 2026-07-01 LATE PM (Claude Science session, Comp Chem MD specialist) — energetic features + directional result; see top section.
**Purpose:** catch a fresh session (Claude Code or otherwise) up on the dynamics-feature
line of work without re-reading the whole history. Pairs with `CLAUDE.md` (operating
contract) and `manuscript/main.pdf` (the living notebook / record of record).

---

## SESSION UPDATE (2026-07-01 LATE PM) — energetic features landed; the key result is directional

**TL;DR.** All the geometry/exposure/flexibility features are weak (AUROC 0.5–0.68 on the honest
test). The one genuinely new axis — gRINN intra-protein **interaction energy** — is now computed on
all 6 Tier-1 systems + Dengue E. The headline is a *direction*, not a magnitude: **epitope hotspots
sit at the LOW-participation nodes of the slowest collective mode (decoupled from the soft motion) and
are more internally bonded (higher DCCM + interaction energy).** Nothing is a strong predictor at N=6.

### The honest evaluation protocol (use this, not raw pooled ρ)
- **Surface-restricted directional AUROC.** Restrict to surface residues (trajSASA z > per-antigen
  median), then AUROC = P(feature_hotspot > feature_rest). This removes the surface-vs-core confound.
  Report the **direction** (>0.5 hotspots higher, <0.5 lower) — the sign IS the result. Never report
  the 2-sided max(a, 1−a) alone (it hides direction).
- **Why it matters:** trajSASA/SASA-fluct **collapse to 0.51 (chance)** under this test — their apparent
  signal was just "detecting the outside." This is the clean confirmation that raw exposure carries no
  hotspot signal (the baseline the project is built to beat).

### Feature leaderboard (Tier-1, surface-restricted directional AUROC, 6 antigens)
| feature | AUROC | hotspot direction |
|---|---|---|
| PC1 participation (positional PCA) | **0.32** | LOWER — decoupled node |
| backbone dihedral flexibility | **0.33** | LOWER — stiffer backbone |
| χ1 flexibility | 0.53 | higher (weak) |
| SASA fluctuation | 0.51 | at chance |
| trajectory SASA | 0.51 | at chance |
| RMSF | 0.58 | higher |
| intEn electrostatic | 0.61 | higher |
| intEn total | 0.61 | higher |
| DCCM coupling | 0.62 | higher — more wired-in |
| **intEn van der Waals** | **0.65** | higher — top energy term |
Coherent physical picture: **anchored, pre-organized surface patches** — quiet in the soft mode, stiff
backbone, but internally well-bonded. All effects weak; treat as hypotheses for FRESEAN + StaB-ΔΔG.

### gRINN interaction energy — the substantive new feature (§12.5, Eq. 38)
- **INTRAMOLECULAR** residue↔residue energy on the **apo antigen alone** (no antibody). Per pair (i,j):
  real-space Coulomb (`intEn_elec`) + Lennard-Jones (`intEn_vdw`), summed over partners → per-residue.
  PME caveat: long-range reciprocal Ewald is not per-pair decomposable; elec is short-range Coulomb only.
- **Not a surface proxy:** intEn_vdw vs trajSASA Pearson = −0.09. On surface-only it STAYS 0.65 while
  trajSASA collapses to 0.51 — carries near-independent information.
- **On Tier-2 Dengue E (verified 96/96 numbering):** energy leads clearly — total 0.63, elec 0.62,
  vdw 0.61 vs coordinate features 0.51–0.53. (RMSF/SASA/slow-mode NOT yet on Tier-2 → comparison incomplete.)
- Output: `benchmark/features/grinn_intenergy_apomd.csv` (pdb,resid,intEn_total/elec/vdw + _z; 7 systems).

### Learning step (single supervised combination — user asked "can weighting combos do better?")
- **No, at N=6.** Unsupervised PC1 (38% var, exposure-dominated loadings) → 0.56; equal-weight z-sum → 0.57.
  Supervised **L2-logistic under leave-one-ANTIGEN-out = 0.33–0.45 (below chance)** — classic
  overfit/anti-generalization: weights fit on 5 antigens point the wrong way on the held-out 6th.
- **Conclusion:** no transferable linear combination beats the best single feature at this N. This is the
  rigorous "nothing works" — motivates the energetic features and the StaB-ΔΔG label scale-up.
  **Do NOT use random forests/boosting at N=6** (CV unreliable; each fold's test set = 1 antigen).

### CORRECTION made this session (important for honesty)
Agent briefly claimed "interaction energy is the ONLY feature above chance on the surface test" — this is
**FALSE and was retracted.** On Tier-1 surface-only, by SEPARATION STRENGTH (folded |AUROC−0.5|, i.e. the
2-sided value, used here only to rank magnitude), slow-mode PCA (0.68 = complement of the 0.32 directional
value in the table) and backbone dihedral (0.67 = complement of 0.33) are ABOVE intEn_vdw (0.65). Note the
sign flip: PCA/dihedral separate hotspots by being LOWER (directional 0.32/0.33), energy/DCCM by being
HIGHER (0.65/0.62) — the table is directional, this ranking is magnitude-only. The true statement: several features cluster in a weak 0.5–0.68 band, NO
clear winner on Tier-1; energy only stands out on the dense Dengue E label. NOTE slow-mode/backbone-flex/
RMSF are NOT independent (all ~re-measure flexibility; PC1-participation vs RMSF Pearson +0.52).

### HA (1HGG) held out — numbering mismatch
Dengue E (1OKE) shotgun labels align 96/96 at delta=0 (trusted). **HA aligns only 62/84 (74%), no global
offset fixes it** — mismatches cluster at resids ~18–60 and 153+ (HA1 head + HA2 stem; cleaved trimer,
per-chain register issue). **HA is HELD OUT of all Tier-2 scoring/gRINN until per-chain numbering resolved.**

### New files this session
- `benchmark/features/grinn_intenergy_apomd.csv` — per-residue intra-protein interaction energy (7 systems).
- `benchmark/features/dccm_apomd.csv`, `dihedral_flex_apomd.csv` — DCCM coupling + φ/ψ/χ1 flexibility (8 systems).
  (`ss_stability_apomd.csv` is EMPTY — MDAnalysis DSSP failed on every system, "unequal N,CA,C,O atoms";
  needs mdtraj.compute_dssp or complete-backbone selection as a follow-up.)
- `benchmark/features/hbonds_<SYS>.csv` — per-residue H-bond counts, intra-protein + protein-water,
  **explicit AMBER selections** (the default MDAnalysis guesser undercounts ~3–4× — it misses the amber
  backbone amide `H` and polar side-chain H names; fixed 2JEL 0.14→0.55 HB/res). Merge to hbonds_apomd.csv.
- Figures/artifacts (Claude Science store, not repo): `clear_leaderboard.png`, `epitope_feature_direction.png`
  (diverging left/right by hotspot direction — the presentation figure), `feature_cheatsheet.pdf`
  (3-page per-feature math + description + direction reference; built with matplotlib mathtext, NOT LaTeX).

### gRINN operational recipe (for reproducing / extending)
- Env `grinn` ships its OWN GROMACS 2026.0 (CPU AVX2_256) — `conda activate grinn`, do NOT `module load
  Gromacs` (that's a GPU container, fatal nvml on CPU nodes). Trajectory mode REQUIRES both `--top` AND
  `--traj`; the topology must be atom-consistent with the structure → pass a DRY topology (strip SOL/NA/CL
  from `[ molecules ]`) + protein-only PDB. Verify residue numbering with the WT-identity check afterward.

### Still running on Gemini at session end
- **FRESEAN** velocity-mode reruns (`34201968` 1AKI, `34201969` 2JEL) — the PROPER version of the PC1 bar.
  Key re-check when they land: do hotspots sit at the nodes of the TRUE low-freq FRESEAN modes (not the
  Cartesian-PCA placeholder)? If it holds with the frequency-selective operator, it's a real mechanistic result.
- **H-bond** parallel jobs (one per system) — 2JEL/1AKI/1JRH done and sane; 1AHW/1BJ1/1HGU/1OKE finishing.

---

## SESSION UPDATE (2026-07-01 PM) — four parallel tracks running; physics-feature expansion

### What's running on Gemini RIGHT NOW (check before assuming anything is done)
None of these use velocities-on-login; all are Slurm jobs. FRESEAN is GPU (gpu-a100);
the feature jobs are CPU (compute partition), so they don't touch the a100/v100 budget.

| track | job id | partition | what | check |
|---|---|---|---|---|
| FRESEAN 1AKI | Slurm `34201968_[1-2]` | gpu-a100 | 20 ns ×2 rep, velocity-preserving rerun | `sacct -j 34201968` |
| FRESEAN 2JEL | Slurm `34201969_[1-2]` | gpu-a100 | 20 ns ×2 rep | `sacct -j 34201969` |
| coord features | CS job `93eda7b1-8399-460f-89b5-077f27b135f0` | compute | DCCM + dihedral flex + SS stability, 6 systems | `.claude-science/jobs/<id>/slurm-*.out` |
| gRINN env build | CS job `9592b8c4-81f0-4642-a53f-7b5e1aac37ec` | compute | clone gRINN + build conda env `grinn` | same |

**The FRESEAN reruns REPLACE the first run** (34188829/34188830), which produced
velocity-LESS trajectories: the reduction step used `gmx trjconv -pbc mol`, and PBC/fit
transforms DISCARD velocities. FRESEAN's whole observable is the velocity
cross-correlation, so that run was unusable and the velocity-bearing full .trr had
already been auto-deleted → forced a 20 ns re-run. Fixed in commit `3569a56`: reduce with
a PLAIN Protein selection (no -pbc/-fit), and gate the full-.trr delete on `gmx check`
reporting `Velocities`, not on file-non-empty. **GROMACS gotcha to remember: never -pbc/-fit
a .trr you need velocities from.** The analysis-chain validation job `7b0f97e0` (Slurm
34201793) caught this at the velocity check before any FRESEAN math ran — that's expected,
not a new failure.

### The physics-feature roadmap (from `b_cell_eptiope_prediction_update_3rd.pdf` §12)
The update doc's §12 is the antigen-only feature taxonomy. Mapping to what the EXISTING apo
trajectories (all 6 intact WITH water + ions, md.tpr + md.xtc, ~10k frames at 10 ps) can yield
— NO new MD, NO velocities needed except where noted:
- **§12.3 flexibility** — RMSF ✅ done (rmsf_apomd.csv). NOW COMPUTING: backbone φ/ψ
  dihedral fluctuation, side-chain χ1 variability, DSSP secondary-structure persistence
  (coord-feature job above → `dihedral_flex_apomd.csv`, `ss_stability_apomd.csv`).
- **§12.4 correlation/modes** — PCA participation ✅ (placeholder, re-measures RMSF at
  +0.59, invalid at low freq — that's why FRESEAN). FRESEAN 🔄 running. NOW COMPUTING:
  residue-residue displacement correlation matrix / DCCM (Eq. 35) → per-residue coupling
  score (`dccm_apomd.csv`, cols coupling_abs/coupling_pos + z).
- **§12.5 intra-protein energetics (Eq. 38)** — the strongest UNTAPPED feature, the true
  energetic-coupling axis. Computing via **gRINN** (`osercinoglu/grinn`, actively developed
  2026, supports GROMACS 2020.7–2025.2, defaults to amber99sb-ildn = our FF). Env building
  now. Outputs per-residue-pair intEn Total/Elec/VdW → sum to per-residue totals (§12.5) +
  the IEM (§12.4). PME caveat: energygrp decomposition gives short-range Coulomb+LJ within
  cutoff (the physically-local interaction) — reciprocal-space PME can't be per-group
  decomposed; standard and appropriate for "internal energetic environment."
- **§12.6 protein-solvent energetics (Eq. 39)** — residue-water/ion interaction energies,
  same `gmx -rerun` machinery, water is present. NOT started.
- **§12.2 hydration / 3D-2PT** — Matthias's 3D-2PT (like FRESEAN, separate toolchain, needs
  velocities → rides the FRESEAN reruns). Water residence/density needs water saved more
  often than 10 ps. DEFERRED.

**MM-PBSA is dropped.** The update doc's actual label plan is **StaB-ddG** pseudo-labels
(benchmark it vs experimental alanine scans first; full alchemical is "its own PhD, don't
reinvent"). MM-PBSA needs the complex and fights the apo thesis — parked for good.

### New feature CSVs this session will add (schema matches rmsf_apomd.csv)
- `benchmark/features/dccm_apomd.csv` — displacement-correlation coupling (§12.4)
- `benchmark/features/dihedral_flex_apomd.csv` — φ/ψ + χ1 fluctuation (§12.3)
- `benchmark/features/ss_stability_apomd.csv` — DSSP helix/sheet persistence (§12.3)
- `benchmark/features/` gRINN outputs (§12.5) — pending env build + validation on 2JEL
All z-scored per structure; join to labels via the identity-verified numbering in
`analysis/rmsf_dynamics_baseline.py` (6 Tier-1 systems align at delta=0; the driver scripts
still carry the analysis harness). Score each on the leaderboard (pooled ρ + LOAO-by-antigen)
alongside RMSF/SASA/mode-participation. NONE of these are committed yet.

### Driver scripts staged (not yet in repo)
- `handoff/coord_features.py` (in the CS job workdir) — the DCCM/flex/SS analysis. Once it
  lands cleanly, promote to `analysis/coord_features.py` and commit.
- gRINN run: `python grinn_workflow.py <protein.pdb> <out> --top <top> --traj <xtc>
  --nofixpdb --skip 100 --nt 8 --initpairfiltercutoff 10 --no_pen` — validate on 2JEL first
  (smallest, 85 res), then the other 5. `--skip 100` = ~100 frames (avg interaction energy
  converges fine; full 10k frames × pairwise reruns is intractable).

---

## ⚠️ CORRECTION (2026-07-01) — 1OKE is Dengue Envelope (a Tier-2 antigen), NOT MT-SP1
A labeling error propagated through the RMSF iteration: the RMSF driver treats `1OKE` as
"MT-SP1" and excludes it for "protease numbering." Both halves are wrong. Ground truth
(verified against the structure files + all label CSVs on disk):
- **1OKE = Dengue Envelope protein** (`structures/DengueE_1OKE.pdb`). Has an apo MD trajectory.
  Its label is **Tier-2 binary**, in `benchmark/tier2_labels_dengueE.csv` (96 residues; schema
  `chain,resid,WT,n_antibodies,critical` — critical-residue calls from shotgun/flow-cytometry
  scans, 781 mapped antibodies). It is NOT in the Tier-1 graded-ΔΔG CSV, which is why my first
  pass wrongly called it label-less. Dengue E is a legitimate Tier-2 system (see
  `manuscript/sections/binary_labels.tex`), paired with HA (1HGG), which is also Tier-2.
- **MT-SP1 = 3BN9 / 3NPS** (27 + 27 rows in the Tier-1 CSV, protease insertion-code numbering
  e.g. 221a). Has Tier-1 labels but **NO apo MD** (no `md/3BN9`, no `md/3NPS`).
- The driver's `KEY2PDB={"MT-SP1":"1OKE"}` / `PDB2KEY={"1OKE":"MT-SP1"}` is a **false
  equivalence between two unrelated proteins**. The "5/47 WT-identity match → excluded for
  protease numbering" reasoning is WRONG: it compared Dengue Envelope's sequence to MT-SP1's
  scan positions. Wrong-protein, not wrong-numbering.

**Effect on the Tier-1 leaderboard: NONE.** The 6-antigen / 9-epitope / 97-residue Tier-1 set
(1AKI, 2JEL, 1JRH, 1BJ1, 1AHW, 1HGU) never included either — Dengue E is Tier-2, MT-SP1 has no
MD. Every Tier-1 ρ / LOAO number stands. The error is confined to comments, the two dicts, and
prose. Two real facts it surfaces:
- **1OKE (Dengue E) is a Tier-2 antigen with both MD and labels** — it should be scored on the
  **binary** track (critical vs not, AUROC / precision-recall), NOT dropped, and NOT joined to
  the Tier-1 graded set. Same for HA (1HGG). The dynamics features (RMSF, FRESEAN, etc.) apply
  to it directly; only the label type and metric differ.
- **MT-SP1 (3BN9/3NPS) is a candidate for a FUTURE Tier-1 MD run** — good graded labels (27+27),
  no trajectory yet. Correct framing, not "excluded due to numbering."

**Files still carrying the error (fix these):**
- `analysis/rmsf_dynamics_baseline.py` — lines ~22, 49, 51: the `MT-SP1↔1OKE` dict entries +
  the exclusion comment. Drop the false mapping. 1OKE = Dengue E (Tier-2, score on binary
  track separately); MT-SP1 = 3BN9/3NPS (Tier-1 label, no MD). Neither belongs in the Tier-1
  join, so no Tier-1 numbers change.
- `manuscript/sections/dynamics_rmsf.tex` — lines ~27–28, 43, 72: "MT-SP1 (1OKE)" phrasing +
  the "excluded (numbering)" caveat. Reword: MT-SP1 is Tier-1-label-only (no MD); Dengue E
  (1OKE) is Tier-2 and scored on the binary track.
- `manuscript/figures/FIGURES.md` — reconcile the `mtsp1` vs `dengueE` key rows.
- Tier-2 scoring of the dynamics features (RMSF now; FRESEAN later) on Dengue E + HA is a
  genuine TODO — it was never done, and 1OKE's apo trajectory has been sitting unused because
  of the mislabel.

---

## The question this project is testing
Do antigen-only (apo) dynamics/solvation/energetics features predict **B-cell epitope
hotspots** — where "hotspot" = alanine-scan **ΔΔG_bind** (an *energy* label, not a
distance/exposure label). The insight: static surface exposure correlates with
distance-based epitope labels but NOT with the energy label; a better (dynamics/energetics)
feature should. Label is **per-antibody-epitope** (a [residue × antibody] array), never
averaged over antibodies. Evaluation is **leave-one-antigen-out (LOAO)** — effective N =
#antigens (~6), split by antigen (HEL's 3 antibodies share one trajectory → antibody split
would leak).

## Feature leaderboard so far (pooled Spearman ρ vs ΔΔG; LOAO ρ), 6 antigens / 9 epitopes / 97 residues
| feature | pooled ρ | LOAO ρ | file |
|---|---|---|---|
| slow-mode participation (gmx covar PCA, mode1) | −0.10 | −0.26 | benchmark/features/dynamics_sasa_modes.csv |
| 5-slowest-mode participation (PCA) | −0.17 | −0.22 | (same) |
| RMSF | −0.12 | −0.15 | benchmark/features/rmsf_apomd.csv |
| SASA fluctuation | −0.07 | −0.10 | (dynamics_sasa_modes.csv) |
| trajectory SASA (exposure) | −0.16 | −0.08 | (same) |
| apo B-factor (static baseline) | −0.23 | — | benchmark/features/bfactor_apocrystal.csv |

**All Tier-1 features are weak on the energy label** — confirms the insight (exposure fails
on the energy label). Collective-mode participation has the strongest honest LOAO lean, sign
negative → hotspots are LOW-participation residues = **nodes of the soft modes** = decoupled
from the global collective motion (geometric decoupling; NOT the same as energetic coupling
à la Serapian/Colombo, which is a separate feature not yet computed).

**IMPORTANT caveat on the PCA mode feature:** `gmx covar` gives *quasi-harmonic positional
PCA* modes, which are invalid at the lowest frequencies (this is the whole motivation for
FRESEAN). Per-residue PCA participation correlates +0.59 with RMSF — it is largely
re-measuring flexibility, not an independent observable. It is a **placeholder**. The real
quantity is FRESEAN low-frequency mode participation (below).

---

## WHAT IS RUNNING RIGHT NOW (Gemini, Slurm) — as of 2026-07-01
FRESEAN low-frequency collective-mode analysis, the physically-correct replacement for the
PCA placeholder. See `docs/ENVIRONMENT_fresean.md` for the full protocol + build recipe.

- **Jobs:** `34188829_[1-2]` (1AKI / lysozyme), `34188830_[1-2]` (2JEL / HPr) — 2 replicas
  × 20 ns each on gpu-a100. Submitted via raw `sbatch` from repo root (NOT the Claude
  Science compute layer → no auto-notifications; poll the run dirs).
- **Systems:** the two smallest antigens, chosen to validate the pipeline before the large ones.
- **Protocol (frozen):** Fast Sampling SI settings — 2 fs, h-bonds/LINCS, Nosé–Hoover
  τ_t=1.0, Parrinello–Rahman τ_p=2.0, AMBER99SB-ILDN/TIP3P, coords+velocities to `.trr`
  every 20 fs, fresh Maxwell velocities per replica (reproducible gen_seed = 20260702/03/…).
  Direct production from the equilibrated NPT box; **first 1 ns discarded** in analysis.
- **Replica count = 2** (not the paper's 5): we are NOT re-proving FRESEAN reproducibility —
  2 gives matrix-averaging + one rep1-vs-rep2 mode-overlap sanity check.
- **Disk handling:** the full-system `.trr` is ~900 GB/replica (real; matches the authors'
  own mdp). Each job strips to a **protein-only `.trr`** (velocities preserved, ~40× smaller,
  ~47 GB) with `gmx trjconv`, then deletes the full `.trr` — mirrors the repo's 03-CG +
  02-MD/clean.sh. Scratch had 127 TB free.

### To check progress
```
ssh gemini
squeue -u bneff | grep fresean
sacct -j 34188829,34188830 --format=JobID,State,ExitCode,Elapsed | grep -v '\.'
# per-replica MD progress:
tail /scratch/bneff/bcell_epitope/md/1AKI/apo/out_fresean/rep1/md_fresean.log
```
Kept outputs per replica: `md_fresean_prot.trr` (protein, w/ velocities), `.tpr`, `.gro`, `.log`.

---

## NEXT STEPS (when trajectories land) — the FRESEAN analysis chain
FRESEAN toolbox is BUILT on Gemini: `/home/bneff/FRESEAN-metadynamics/bin/` (10 binaries),
conda env `fresean` (gcc14/fftw/gsl). Reference example = the repo's HEWL protocol
(`scripts/protocol_FRESEAN+convergedWTMetad+DCCM/`), which is literally lysozyme (our 1AKI).

Per replica (see repo steps 03-CG → 04-FRESEAN):
1. **Sanity: velocity-ACF settling.** Confirm the velocity autocorrelation has settled by
   1 ns (Nosé–Hoover ringing from the cold velocity start). If not, extend the discard.
2. **Coarse-grain** protein `.trr` → 2 beads/residue (1 for Gly): `fresean mtop` (topol_prot.top)
   then `traj_coarse` / `fresean coarse`. Needs a protein-only topology `topol_prot.top`.
3. **Velocity cross-correlation matrix**, max correlation time 2 ps: `fresean matrix` / `gen-modes_omp`.
4. **Eigendecompose the zero-frequency matrix**: `fresean eigen`. Modes 1–3 = translation,
   4–6 = rotation, **7+ = the anharmonic low-frequency vibrations** (authors extract modes 7,8).
5. **Extract** modes 7+ (`fresean extract`) → per-residue participation.
6. **Average** participation over the 2 replicas; **rep1-vs-rep2 mode overlap** = the stability
   check (high ⇒ feature real, stop; low ⇒ 20 ns insufficient for that system, extend).
7. Map bead participation → per-residue → z-score per structure → write
   `benchmark/features/fresean_lowfreq.csv` (schema: antigen,key,pdb,chain,resid,resname,value,value_z
   — match rmsf_apomd.csv).
8. **Score against ΔΔG** on the same leaderboard (pooled + LOAO). This REPLACES the PCA
   placeholder row. Numbering: reuse the identity-verified alignment in
   `analysis/rmsf_dynamics_baseline.py` (6 systems align at delta=0; 1OKE/MT-SP1 are not in the graded set --- see the CORRECTION block above).
9. Build the figure (figure-style), write `manuscript/sections/dynamics_modes.tex` (notebook
   voice — see dynamics_rmsf.tex), fill the `[verify]` placeholders in
   `docs/ENVIRONMENT_fresean.md` (ns/day, GB/replica, mode indices), commit + push + pull.

## Deferred threads
- **High-frequency (mid-IR) participation** — FRESEAN-2023 spectroscopic settings (0.5 fs,
  unconstrained, 4 fs output). Tier-2b; deliberately deferred to explore slow modes first.
  Could carry hotspot info from local chemical environment.
- **Energetic coupling** (Serapian/Colombo, gRINN-style intra-protein interaction-energy
  matrix) — the *energetic* (not geometric) coupling feature. Not yet computed. This is the
  distinct axis that the mode-participation result gestures at but does not measure.
- **1OKE / MT-SP1** — see the CORRECTION block at the top. 1OKE (Dengue E) has MD but no
  label; MT-SP1 (3BN9/3NPS) has labels but no MD. Neither is in the leaderboard. MT-SP1 is a
  candidate for a future apo MD run (it has 27+27 labels); 1OKE is an orphan trajectory.
- **1HGG (HA)** and other large antigens — FRESEAN on these needs the CG + likely more memory;
  large systems are why the 900 GB `.trr` handling matters. Do after the two small systems validate.
- `benchmark/features/dynamics_sasa_modes.csv` is written but not yet committed / not yet
  written into the manuscript (the SASA + PCA-modes leaderboard section is unwritten).

## Repo conventions (from CLAUDE.md / docs/)
- Stage + analyze locally; run MD on Gemini (Slurm). Trajectories stay on scratch (gitignored).
- Submit from repo root (`/scratch/bneff/bcell_epitope`); scripts anchor to `$SLURM_SUBMIT_DIR`.
- Module is `Gromacs` (capital G) = 2023.2-dev; `gmx` only on GPU nodes.
- Git identity: Brandon Neff <bdneff@asu.edu>; commits carry a "Built with Claude Science" trailer.
- Manuscript builds LOCALLY: `cd manuscript && tectonic main.tex` (self-contained, no TeX Live;
  matches CLAUDE.md/README — NOT latexmk, which needs a full TeX install). Not on Gemini, not in sandbox.
- Feature CSV schema: `antigen,key,pdb,chain,resid,resname,<value>,<value>_z` (z-scored per structure).
