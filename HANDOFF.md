# HANDOFF — bcell_epitope dynamics features (read this first)

**Last updated:** 2026-07-01 (Claude Science session, Comp Chem MD specialist).
**Purpose:** catch a fresh session (Claude Code or otherwise) up on the dynamics-feature
line of work without re-reading the whole history. Pairs with `CLAUDE.md` (operating
contract) and `manuscript/main.pdf` (the living notebook / record of record).

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
