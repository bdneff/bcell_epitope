# HANDOFF — bcell_epitope dynamics features (read this first)

**Last updated:** 2026-07-07 (Claude Science session #5) — **Physics-feature program CLOSED at N=6: FRESEAN collective modes, energetic MLCE coupling, and every combination reduce to exposure under strict evaluation. The strict train-on-6/test-on-Dengue test is the capstone. New manuscript section `dynamics_modes.tex` written. Scale-up (SKEMPI labels + 1OKE dense run) is now the path. 1OKE chunk 1 validated; chunks 2–20 HELD pending user go.**

---

## ⚠️ SESSION #5 (2026-07-07) — READ FIRST

### Headline: the physics-feature program is closed at N=6 (this is a RESULT, not a failure)
Every MD-derived physics feature — FRESEAN collective-mode involvement, dynamical isolation, gRINN
energetics, Colombo's MLCE energetic coupling, and all combinations (linear + MLP) — reduces to
**solvent exposure** once evaluated honestly (surface-restricted, cross-antigen, or held-out antigen).
This was established three independent ways this session. The ceiling is the **label set (N=6 antigens)**,
not the model. The two moves that attack it (SKEMPI expansion, 1OKE dense labels) are staged.

**1. FRESEAN involvement 0.75 was a buried/exposed confound.** On 1AKI's 23 labels, mode involvement
detects epitopes at AUROC 0.745 (p=2e-4) whole-protein — but surface-restricted it falls to **0.638 (n.s.)**,
and the isolation reading falls to **0.552 (chance) and flips direction** (integration-vs-SASA pearson −0.68).
The mode signal mostly encoded which residues are on the outside. (`analysis/fresean_involvement.py`)

**2. "Co-moving decoupled patch" is GENERIC, not epitope-specific.** All 12 dynamical communities of 1AKI
show within-patch coupling ≫ between (mean +0.56 vs −0.04); epitope patches are no more coherent than
non-epitope (0.55 vs 0.56). The property is real and visualisable but has zero discriminative value.

**3. Per-antibody lysozyme decomposition = the one clean POSITIVE (zero new MD).** SKEMPI single-mutant
records for 1DQJ (HyHEL-63, 11 resids) and 1VFB (D1.3, 12 resids) EXACTLY reproduce our pooled 23-residue
1AKI set. Their footprints each co-move internally (+0.20, +0.24) and are decoupled from each other (−0.02)
→ two antibodies target two dynamically distinct patches. (1YQV = single R45K on a non-epitope isolated
module; dropped.) Saved `lysozyme_per_antibody_labels.csv` (24 labels / 3 antibodies, offset 0 to 1AKI).

**4. Colombo MLCE tested on its own terms — negative all three ways.** MLCE (matrix of low coupling
energies) implemented on gRINN's pairwise interaction-energy matrix (`average_interaction_energies.csv`,
on scratch for all 7 systems). (i) Paper's qualitative standard: bottom-15% MLCE surface patch overlaps
epitopes on 2/6 systems (2JEL, 1JRH), fails on 4. (ii) Strict surface per-residue AUROC = **0.49 (chance)**.
(iii) As a feature: MLCE correlates **−0.84** with per-residue interaction energy we already have; adding it
to the surface LOAO model moves 0.40→0.42 (noise). Note: our MLCE uses gRINN in-vacuo Coulomb+LJ; Colombo's
adds GB solvation — a [verify] gap, but the −0.84 redundancy suggests solvation won't rescue it.

**5. The decisive test — train Tier-1, predict unseen Dengue E.** Fit on all 6 Tier-1 antigens (alanine ΔΔG),
predict Dengue envelope scored against the orthogonal shotgun n/N label. **Every physics-feature set lands at
chance (0.43–0.52); coupling features below it. Only exposure transfers: crystal SASA AUROC 0.58, graded
ρ=+0.15.** Exhaustive 2047-subset search (MLCE included) maxes at 0.61 inside the shuffled null; MLP memorises
(train→1.0, test 0.58±0.03). This is the strongest test in the project and it's unambiguous.
- Transferable feature set (exists for BOTH Tier-1 and 1OKE): intEn_total/elec/vdw, hb_intra/water,
  bb_dih_std, chi1_std, coupling_abs/pos (DCCM), mlce_coupling. **RMSF and slow-modes are NOT in the external
  test** because 1OKE has no trajectory-dynamics features yet (that's the pending FRESEAN run) — rerun the
  external test with them once 1OKE modes land.

### Manuscript: NEW section `manuscript/sections/dynamics_modes.tex` (written this session)
- `\section{Collective modes and the strict-evaluation ceiling}\label{sec:modes}`, wired into `main.tex`
  after `dynamics_energetics`. Notebook voice, 7 paragraphs, 5 figures. This is the rewrite the runlog's
  "Note for a later revision" flagged as pending (energetics narrative rested on the retracted energy result).
- Figures staged in `manuscript/figures/`: `surface_isolation_control_1aki.png`, `per_antibody_patches_1aki.png`,
  `patch_generality_1aki.png`, `mlce_as_feature_1aki.png`, `external_dengue_test.png`.
- **New bib entries** `fresean2023` (JCTC 10.1021/acs.jctc.2c01309) and `fastsampling2024` (Sci Adv
  10.1126/sciadv.aea4617, arXiv 2411.08154) — **author lists carry `note = {... [verify]}`; Brandon must
  confirm against the published versions** (Brandon is a co-author on fastsampling2024, NOT on fresean2023).
- `runlog.tex` got a session-#5 paragraph summarising the closure.
- **Manuscript NOT built** — tectonic isn't in the Claude Science sandbox. Structure validated (all cite keys
  defined, all figorbox figures present, begin/end balanced). **Build locally: `cd manuscript && tectonic main.tex`.**

### SKEMPI label expansion (the antigen scale-up path)
- Downloaded SKEMPI 2.0 (`https://life.bsc.es/pid/skempi2/database/download/skempi_v2.csv`; life.bsc.es
  is network-allowlisted). 7085 records / 345 PDBs; ab-antigen subset 897 records / 52 complexes.
- **CRITICAL curation fix:** first pass (594 labels/48 complexes) was WRONG — it counted antibody-CDR
  mutations as epitope labels. After checking every chain against RCSB entities, mutations TRULY on the
  ANTIGEN chain = **211 labels / 19 complexes**. Saved `skempi_antibody_antigen_labels.csv` (⚠️ still contains
  the over-broad 594 set; the true-211 subset is computed in-session, not yet re-saved — REGENERATE before use).
- Usable NEW soluble antigens for MD: **MT-SP1 (3BN9+3NPS, same 227-res construct, 54 alanine labels, 2
  antibodies), enterotoxin B (3W2D, 4), gp120-OD (4JPK, 8)**. Most SKEMPI "antigens" (integrin/EGFR/HER2/IL-13/
  gp120-3NGB/HA) had Fab-side mutations → NOT usable.

### 1OKE (Dengue E) FRESEAN — chunk 1 VALIDATED, chunks 2–20 HELD
- Chunk-1 job `8b8092e2-23e8-43d4-bbec-5cbbb50f9f17` SUCCESS. Velocity gate PASSED: `prot_chunk1.trr` has
  50,001 frames of coords+velocities (the −pbc-mol velocity-drop bug did NOT recur). Full 404 GB `chunk1.trr`
  deleted; 6.9 GB protein-only retained. `fresean_state.txt`=1, resumable. Perf 1.34 h/ns → chunks 2–20 ≈ 25 h.
- **User HELD the launch** to settle MLCE/MT-SP1 first. To resume: `sbatch md/1OKE/apo/md_fresean_chunked.sh`
  (auto-resumes from chunk 2).

### MT-SP1 build — BLOCKED on numbering, do NOT attach labels yet
- 3BN9 + 3NPS fetched. Both clean to build (HETATM only EDO/HOH/SO4/CL/NA — no active-site inhibitor;
  7 SS-bonds/copy). Antigen is the IDENTICAL 227-res construct (resid 16–244, 100% identical) → **ONE apo MD
  serves both antibodies.**
- **BLOCKER:** SKEMPI renumbered the PDB's chymotrypsin insertion-code numbering (raw 60,60a,60b,60c,60e,60f)
  sequentially; NO integer offset maps them (best +5 = 10/27 WT match). **Labels require proper sequence
  alignment before any MD build — same class as the earlier project label bugs.** Not started.

---

## ⚠️ SESSION #4 (2026-07-07) — READ FIRST

### 1. FRESEAN 1AKI mode extraction COMPLETE + rigid-body verified
- Both 1AKI replicas extracted 30 zero-freq modes (`md/1AKI/apo/fresean_extract/rep{1,2}/evec_freq1_mode1-30_cg.xyz`, `eval_matrix_cg.mmat.dat`, `ref-cg.gro`). CG = 246 beads / 129 residues (2-bead BACK+SIDE; Gly has BACK only).
- **Rigid-body check PASSED:** modes 1-6 rigid (combined trans+rot subspace overlap 0.91; reps 0.907/0.911), modes 7-30 internal (~0.015). Verified split: modes 1-3 translation-dominated (0.87-0.92), 4-6 rotation-dominated (0.86-0.96). NOTE FRESEAN sorts eigenvalues large→small so this ordering was computed, not assumed.
- **The "failed" compute-job status is a FALSE ALARM** (known): the job's final `ls eigenval*cg*` matches nothing (file is `eval_matrix_cg.mmat.dat`) → nonzero exit under `set -eo pipefail`. Outputs are valid. Harvest globs also failed; pull with `c.download(remote, local=None)`.

### 2. NEW RESULT — collective-mode involvement DETECTS 1AKI epitopes, beats exposure
Feature = per-residue involvement in internal modes: `involvement_i = Σ_{m=7..k} Σ_{beads∈i} |e_{mb}|²` (eigenvectors unit-normalized → total over residues = #modes). Script: `analysis/fresean_involvement.py`.
- **Detection AUROC 0.745 (p=2e-4)** on 23 alanine-scan labels — epitopes are MORE involved in collective motion (OPPOSITE of the old Cartesian-PCA "isolation/decoupling" reading).
- **Beats the exposure control:** SASA (relsasa_md_z) 0.671, RMSF 0.664. Involvement survives SASA control (residual 0.66) AND RMSF control (residual 0.73) — not just exposure, not just flexibility.
- **Band sweep:** cumulative 7→k peaks at **7-21 (0.777)** on a broad 7-14..7-24 plateau (~0.75-0.78). DO NOT hard-tune the endpoint — it's a within-sample argmax on N=23. Slowest modes alone are weak (7-8: 0.62, 7-9: 0.655 — below SASA); the signal needs the mid-frequency internal band. No single mode dominates (strongest: 21, 27, 23).
- **Graded gradient vs ΔΔG is flat** (ρ=-0.10) — detects epitope-vs-not, doesn't rank hotspot strength.
- **CAVEAT (critical): N=1 protein.** This is 1AKI's modes vs 1AKI's labels — NOT the Dengue-E epitopeness score. Promising single-system signal, not a result. The whole point of the 1OKE run below is to test it with dense labels + real power.
- Figures (artifacts, not committed): `fresean_involvement_1aki.png`, `involvement_vs_sasa_1aki.png`, `mode_band_sweep_1aki.png`, `fresean_mode_verification.png`. Animations: `cartoon_mode{01,04,07,08,09}.gif` (PyMOL), `fresean_modes_grid.gif` (3x3 labeled), `mode_XX.gif` (matplotlib backbone). NMWiz file for local VMD: `~/Desktop/1AKI_rep1_FRESEAN_modes.nmd` (`vmd -e <file>`).

### 3. 1OKE (Dengue E) FRESEAN run — STAGED, NOT LAUNCHED
- **Why it needs new MD:** 1OKE's existing trajectory is xtc-only (no velocities); FRESEAN needs coords+velocities in a full-precision .trr. Reuses the equilibrated box (`out/npt.gro`) + `out/topol.top` (AMBER99SB-ILDN/TIP3P — same FF as 1AKI, REQUIRED for mode comparability) with fresh Maxwell velocities.
- **Chunked to cap storage:** 1OKE is 347,139 atoms (6,129 protein). A monolithic 20 ns velocity trr = ~8 TB transient. Chunked 1 ns × 20 with reduce-then-delete caps peak scratch at ~one full chunk (~416 GB) + accumulating protein-only (~7.4 GB/chunk → ~147 GB). Fits /scratch (105 TB free).
- **Files (on scratch + repo `md/1OKE/apo/`):** `configs_fresean/md_fresean_chunk1.mdp` (gen_vel), `md_fresean_cont.mdp` (continuation), `md_fresean_chunked.sh` (resumable via `fresean_state.txt`, velocity-preserving reduction with gmx-check gate, trjcat at end). 1 replica.
- **Submit:** `sbatch md/1OKE/apo/md_fresean_chunked.sh` (gpu-a100, runs on compute node — NEVER login). ~20 chunks over many hours. NOT yet submitted — awaiting Brandon's go.
- **Tamarind ruled out:** their GROMACS only emits xtc+pdb (no velocities) → can't feed FRESEAN.

### 4. Reconciliation asks from Claude Code (local session) — OWED
- **Blocked-CV disagreement:** CC's Dengue combination under spatially-blocked CV (contiguous sequence-block holdout, since epitopes are patches) finds **B-factor is the ONLY robust increment over SASA (+0.05); DCCM's +0.007 is spatial leakage** (vanishes under blocking), vdW coef collapses to ~0.01. Honest headline: exposure + crystal rigidity ~0.66, MD features don't add under an honest split. → re-run our combination under blocked/grouped CV to confirm.
- **1OKE apo-MD features OWED (quick, from existing 1OKE trajectory, no new MD):** (a) per-residue apo-MD RMSF, (b) trajectory-averaged apo SASA — both schema-matched to `rmsf_apomd.csv`, MD 1..394 numbering, 96/96 WT-identity check. Commit to `benchmark/features/`.
- **Provenance gap:** commit the Dengue feature-build script (1OKE SASA + crystal B-factor live only in the Science workspace, in NO committed CSV).

---

## ⚠️ SESSION #3 (2026-07-06) — READ FIRST

### 1. gRINN energy feature was INVALID (position ramp) — now REBUILT
- **Symptom Brandon caught:** per-residue `intEn_vdw` trended upward with residue *number* in the Dengue profile.
- **Diagnosis:** the committed `grinn_intenergy_apomd.csv` per-residue values ramp with resid (Pearson **+0.75..+0.86 in ALL 7 systems**) and were all-POSITIVE (unphysical). They have **~0 correlation** (pearson −0.05..−0.11) with the true pairwise energies. The per-residue **aggregation step** was broken (one-sided / mis-summed), NOT gRINN itself.
- **Raw pairwise file is fine:** `md/<SYS>/apo/grinn/out/average_interaction_energies.csv` (cols `Res1_ResNum,Res2_ResNum,Avg_VdW/Elec/Total_Energy`; each unordered pair stored ONCE).
- **FIX (applied):** symmetric re-aggregation — sum each pairwise E_ij into BOTH residues i and j. Verified: symmetric per-residue sum = 2× pairwise total. No MD/gRINN re-run needed.
- **Rebuilt `benchmark/features/grinn_intenergy_apomd.csv` for all 7 systems** (UNCOMMITTED as of this writing — `git status` shows it modified). After fix: position corr −0.05..+0.17 (ramp gone), vdw non-positive. `resid` is gRINN author numbering (1AHW 4-211, 1BJ1 14-107, 1OKE 1-394, etc.) — **must verify alignment to dccm/hbonds (1..N MD numbering) before Tier-1 use.**
- **Corrected 1OKE energy detection:** intEn_vdw AUROC **0.567** (real now: epitope-critical = less-favorable vdW packing / protrusion), intEn_elec **0.495** (dead — was riding the ramp), intEn_total **0.515**. Only vdW carries signal. NOTE the broken vdw also read ~0.57 by positional coincidence (epitopes cluster at high resid).

### 2. Dengue E (1OKE) single-system analysis — the statistically honest track
Pivot away from N=5 Tier-1 (3-17 labeled residues/system, fragile means) to the ONE well-sampled shotgun system: **96 critical / 298 non-critical residues, one system, real p-values.**
- **Per-feature detection AUROC (raw, direction-agnostic):** SASA **0.63** (exposure control), B-factor 0.40 (rigid→critical), HB-water 0.57, vdW 0.57 (corrected), HB-intra 0.41, DCCM 0.53, bb-flex 0.51.
- **Honest 5-fold CV combination (L2 logistic, 5 repeats, 392 residues):** SASA-only **0.628** → +B-factor **0.684** → +DCCM **0.696** → +bb-flex **0.697** (nominal best; the +bb-flex step is +0.001, within CV noise, so SASA+B-factor+DCCM is the parsimonious model). **Dynamics/rigidity add +0.07 on top of the exposure control.**
- **B-factor is the key additive feature** (single-feature runner-up to SASA at 0.59). vdW and HB-water do NOT add on top of SASA+B-factor+DCCM (redundant with exposure).
- **B-factor provenance caveat** (`benchmark/features/BFACTOR_1OKE_NOTE.md`): 1OKE crystal is HOLO (β-OG detergent bound, not antibody). Verified apo counterpart **1OAN** gives B-factor AUROC 0.386 ≈ holo 0.40 — signal is NOT a ligand artifact. B-factor+SASA are crystal/structural (single-frame); the rest are apo-MD (source mismatch to keep in mind).
- Figures: `dengue_profile_raw.png` (all 7 metrics vs epitope score, raw units), `dengue_combination.png` (single-feature bars + incremental CV model).
- **OWED:** apo-MD RMSF for 1OKE (offered, not yet run — the MD-native flexibility metric to pair with/replace crystal B-factor).

### 3. FRESEAN extraction — MD done, extraction pipeline underway (1AKI first)
- **All 4 replicas complete with velocities:** 1AKI rep1/rep2, 2JEL rep1/rep2 — each 1,000,001 frames, Coords+Velocities+Box, 0.02 ps, protein-only trr (1AKI 1960 atoms).
- **Pipeline (Heyden `protocol_FRESEAN+convergedWTMetad+DCCM`):** 03-CG (mtop → coarse-grain → matrix) → 04-FRESEAN (matrix → eigen → extract top 30 zero-freq modes).
- **Gotchas resolved this session:** (a) FRESEAN bins need `LD_LIBRARY_PATH=/home/bneff/.conda/envs/fresean/lib` for libfftw3; (b) protein-only `topol_prot.top` built by cutting full topol.top through the `#endif` of the POSRES block (cutting before it leaves `#ifdef POSRES` open → "No molecules defined"); (c) the `-pbc mol` velocity-drop bug does NOT bite when protein-only trr is matched to a protein-only tpr (1960=1960) — verified both `-pbc whole` and `-pbc mol` preserve velocities in gmx2026 with matched atom counts.
- **Working dir:** `md/1AKI/apo/fresean_extract/` (has `topol_prot.top`, `topol_prot-aa.mtop`, `prot.tpr`, `posre.itp`).
- **Status:** diagnostic jobs passed; full CG→FRESEAN run being submitted for 1AKI both reps. 2JEL after 1AKI verifies.

### FRESEAN is a GATE, not a result — decided 2026-07-06 (Brandon)
1AKI/2JEL are Tier-1 (3–17 labeled residues) → cannot be scored for a hotspot p-value. Treat these
two runs as a **methods/implementation gate + convergence test**, NOT a hypothesis test. Do NOT quote a
FRESEAN AUROC off the fragile Tier-1 labels (that is the exact trap the Dengue pivot escaped). Two
go/no-go checks when the extraction lands:
1. **Convergence:** rep1-vs-rep2 mode overlap on the top zero-freq modes (7,8,…). High ⇒ 20 ns is
   enough; low ⇒ the run length is insufficient and any Dengue run must be longer.
2. **Independence from RMSF:** per-residue FRESEAN participation vs RMSF Pearson. The PCA placeholder
   FAILED here (+0.52 = just re-measured flexibility). Near-zero ⇒ genuine new axis; ~+0.5 ⇒ a fancier
   flexibility proxy that adds nothing over B-factor/RMSF — park it.
   (1AKI bonus: the repo reference FRESEAN protocol is HEWL = same protein → a ground-truth mode-7/8
   comparison to validate the eigenvectors independent of any epitope question.)

**Decision the gate drives:** if BOTH pass (clean implementation + converges + RMSF-independent) → the
pipeline works, so run a **velocity-bearing 20 ns FRESEAN on Dengue E (1OKE)** and score participation as
a 4th increment in the honest CV (SASA → +B-factor → +DCCM → **+FRESEAN?**). This needs a DEDICATED run —
the standard apo Dengue trajectory saves coords every 10 ps with no velocities; FRESEAN needs 20 fs
velocity output. If either check fails → park FRESEAN, note why in the notebook, don't spend a Dengue GPU
run on it. Do NOT write `dynamics_modes.tex` until the gate passes (leave it unwritten).

---

## ⚠️ CRITICAL CORRECTION (2026-07-01 late, session #2) — READ BEFORE TRUSTING ANY ΔΔG RESULT

**Bug found:** the alanine-scan ΔΔG table (`Residue #`) uses the antibody-**complex** PDB author
numbering. The apo-MD structures are numbered differently. All feature↔ΔΔG joins this project made
on raw residue number were therefore attaching ΔΔG to the **wrong residue** for 4 of 6 systems.

**Verification method:** compare scan WT one-letter against the residue actually present at that number
in the MD structure. Before fix: **39/101** labeled residues matched. After per-system offset: **101/101**.

**Offset map (MD_resid = scan_resid + offset):**

| system | offset | notes |
|---|---|---|
| 1AHW | −3 | |
| 1BJ1 | −13 | (also has separate energy-scale bug, see below) |
| 1HGU | −1 | |
| 1JRH | −10 | |
| 1AKI | 0 | already correct |
| 2JEL | 0 | already correct |

Persisted: `benchmark/labels/numbering_offsets.json` + corrected `benchmark/labels/ddg_by_md_residue.csv`.
**Any feature-build code that joins ΔΔG on residue number MUST apply this offset first.**

**THIRD bug (found by Claude Code repro, 2026-07-01 late):** `grinn_intenergy_apomd.csv` is
**0-indexed** (resid 0…N−1); every other feature file + corrected labels are **1-indexed**. Joining
gRINN on raw resid attaches energies **one residue off**. FIX: `grinn["resid"] += 1` (verified clean
+1 shift, identical per-system counts). Any build using gRINN MUST apply +1 first.

**Impact:** correcting it drops interaction-energy detection AUC **0.67 → 0.57** — intEn_vdW was NEVER
a real feature; its apparent signal was the misalignment. **The "internal energy / protrusion /
under-satisfied side chain" mechanism is DEAD — do not repeat it.**

**FULLY-CORRECTED leaderboard (both label offsets + gRINN +1; 5 systems, 1BJ1 excluded, surface, per-antigen AUROC):**
| feature | detection AUC |
|---|---|
| **H-bonds to water** | **0.66** ← single best, survived every correction |
| H-bonds within protein | 0.60 |
| Trajectory SASA | 0.57 |
| Interaction energy (vdW) | 0.57 (was 0.67, off-by-one artifact) |
| DCCM coupling | 0.56 |
| Backbone flexibility | 0.52 |
| PC1 participation | 0.50 (dead) |
| SASA fluctuation | 0.40 |

- Gradient (HB-water ρ among labeled): **−0.50** (unchanged; uses no gRINN).
- HB-water is the one feature that survived every correction (detects 0.66, grades −0.50, no gRINN, no 1BJ1 dependence). It is the honest headline.
- **Claude Code is writing the runnable eval script** (`analysis/reassess_corrected_labels.py`, currently a STUB) to regenerate both figures from committed CSVs and lock these numbers with asserts. The two committed figures used PRE-gRINN-fix energies → their intEn_vdW panels are STALE, superseded by CC's regen.
- Combination fits may move: intEn_vdW was in the top combos (0.73); at 0.57 the H-bond-only combo (HBwater+HBintra+SASAfluct, 0.72, no gRINN) is likely the honest winner — CC to confirm.

**Second, separate 1BJ1 bug:** gRINN ran on a single VEGF **monomer** (94 residues, chain 0–93) but
VEGF is an obligate homodimer. Missing partner-chain contacts → energy scale broken (vdW median 126 vs
43 elsewhere; elec +36 vs −2) AND PCA mode-1 = 90% variance rigid-body monomer rocking. 1BJ1 inflated
BOTH the vdW-energy and PC1 detection results. Exclude 1BJ1 or re-run gRINN on the dimer before trusting energy features.

**Status of prior figures:** every feature-vs-ΔΔG figure and combination fit from session #1 and earlier
in session #2 used the mislabeled join → **stale**. Rebuild on corrected labels next.

---

**(prior header, retained):** energetic features + directional result; see section below.
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
