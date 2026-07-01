# Compute environment — FRESEAN low-frequency mode analysis (Gemini)

Tier-2 dynamics feature: per-residue participation in **anharmonic low-frequency collective modes**,
computed with the Heyden-lab FRESEAN toolbox. This is the physically-correct replacement for the
`gmx covar` quasi-harmonic PCA placeholder (see `manuscript/sections/dynamics_modes.tex`): positional
PCA is invalid at the lowest frequencies, which is the entire motivation for FRESEAN
(Sauer & Heyden, *J. Chem. Theory Comput.* 2023, DOI 10.1021/acs.jctc.2c01309; PMID 37515568).
Protocol follows the co-authored *Fast Sampling of Protein Conformational Dynamics*
(Sauer, Mondal, **Neff**, Maiti, Heyden; *Sci. Adv.* 2025, DOI 10.1126/sciadv.aea4617; arXiv 2411.08154).

## FRESEAN code (confirmed on Gemini 2026-07-01)
- Repo: **`HeydenLabASU-collab/FRESEAN-metadynamics`** (the collaborator fork the author uses),
  cloned to **`/home/bneff/FRESEAN-metadynamics`**.
- Build toolchain env (conda-forge, matched to the repo's tested versions to avoid ABI mismatch with
  Gemini's system gcc 11.5 / GSL 2.4–2.5 modules): **`fresean`** conda env at
  `/home/bneff/.conda/envs/fresean` —
  `gcc_linux-64=14 gxx_linux-64=14 gsl fftw make python=3.12 numpy=2.2.6 matplotlib=3.10`.
  (Repo `environment.yml` only supplies the Python analysis deps; the C build deps FFTW+GSL are added here.)
- Build: from the repo, `make && make install && make clean` — compiles the C tools
  (`traj_coarse`, `gen-modes_omp`, `eigen`, `extract`, …) into `bin/`, linked against `-lfftw3 -lgsl -lgslcblas`
  with `-fopenmp`. The `fresean` wrapper dispatches subroutines (`fresean <tool> -h` for usage).
- [verify] confirm `make install` PATH export (repo appends to `~/.bashrc`); record `fresean` version/commit.

## MD protocol — **Fast Sampling settings** (frozen, 2026-07-01)
Decision (user): explore the **slow collective-motion** picture first; the high-frequency
(0.5 fs unconstrained, mid-IR) participation is a deliberately deferred Tier-2b thread.
Per-system frozen mdp: `md/<PDB>/apo/configs_fresean/md_fresean.mdp`. Departures from the house
`configs/md.mdp` are marked `[FRESEAN]` in the file. The four physics-relevant differences:

| knob | house md.mdp | **FRESEAN (Fast Sampling SI)** | why |
|------|--------------|--------------------------------|-----|
| coords+velocities | `nstxout=0 nstvout=0` (xtc only) | **`nstxout=10 nstvout=10`** → `.trr` every 20 fs | velocity TCF is the FRESEAN observable |
| thermostat | V-rescale, τ_t=0.1 | **Nosé–Hoover, τ_t=1.0** | deterministic; gentler on the velocity autocorrelation |
| production length | 100 ns × 1 | **20 ns × 2 replicas** | low-freq modes converge from short reps (no rare-event dependence) |
| velocities | `gen_vel=no` (continue) | **`gen_vel=yes`, unique seed/replica** | independent Maxwell resamples |

Unchanged and kept comparable to the rest of the apo work: **AMBER99SB-ILDN / TIP3P**, dt = 2 fs,
h-bonds/LINCS, PME (rc=rvdw=1.0), Parrinello–Rahman (τ_p=2.0), 300 K / 1 bar. First systems:
**HPr (2JEL)** and **lysozyme (1AKI)** — the two smallest, to validate the pipeline before the large antigens.

## FRESEAN analysis chain (per replica, then averaged)
Following the Fast Sampling SI:
1. Store coords+velocities every **20 fs** (done by production mdp).
2. Rotate trajectory into the crystal reference frame (remove global rotation).
3. Coarse-grain to **2 beads / residue** (1 bead for Gly) — `fresean coarse` / `traj_coarse`.
4. Build the mass-weighted **velocity cross-correlation matrix**, max correlation time **2 ps** — `fresean matrix`.
5. Eigen-decompose the **zero-frequency** matrix — `fresean eigen`; **modes 1–3 = translation, 4–6 = rotation,
   modes 7+ = anharmonic low-frequency vibrations**. Extract with `fresean extract`.
6. Per-residue participation in modes 7+ → the feature, scored against ΔΔG on the same leaderboard as
   RMSF / SASA / (deprecated) PCA participation. Average the participation over the 2 replicas.
   NOTE: replica count is **2**, not the paper's 5 — we are not re-proving FRESEAN's reproducibility
   (Matthias's Fast Sampling paper established that); 2 replicas give matrix-averaging plus a single
   rep1-vs-rep2 mode-overlap sanity check. High overlap ⇒ feature is stable, stop; low overlap ⇒ the
   20 ns low-frequency mode isn't well-defined for that system and we extend.

## Job geometry
- MD reruns: same as the house apo runs (`--gres=gpu:1`, GPU offload `-nb gpu -pme gpu -bonded gpu`);
  20 ns × 2 replicas per system. [verify] record ns/day and walltime.
- FRESEAN analysis: CPU + OpenMP (`gen-modes_omp`, `eigen` are `-fopenmp`); runs where the `.trr` live on scratch.
- [verify] `.trr` at 20 fs cadence for 20 ns is large (~[verify] GB/replica) — keep on scratch, harvest only
  the small eigenvector/participation outputs back to the repo.

## Sync workflow
Same bridge as the rest of the project: stage frozen `configs_fresean/` + submit scripts here (git-tracked),
commit + push, `git pull` on Gemini, submit from repo root. Trajectories stay on scratch (gitignored);
only per-residue participation CSVs come back. Record job ids + ns/day + `fresean` commit in
`manuscript/sections/runlog.tex`.
