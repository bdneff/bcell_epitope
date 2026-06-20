# Compute environment — Gemini cluster

Stage and analyze in this repo; run production MD on **Gemini** (Slurm). Trajectories live on
Gemini and are gitignored. Settings below are inherited from the working TGen scripts in
`../../md/scripts/` and **must be confirmed on Gemini** before the first run.

## Confirmed on Gemini (2026-06-18)
- [x] Module: **`Gromacs`** (capital G; lowercase `gromacs` errors on a missing `shared` dep).
      GROMACS **2023.2-dev** (`2023.2-dev-20230712-a3d6917-dirty-unknown`, avx2_256 build).
      `gmx` only runs on a GPU node, not the login node.
- [x] GPU partition: **`gpu-a100`** (4-day walltime limit, A100 = fastest; also `gpu-v100`
      8-day, and `gpu-v100-dev` 2-day-limit interactive). Using `gpu-a100` for prep + production.
- [x] Repo cloned to **`/scratch/bneff/bcell_epitope`**; trajectories stay there (gitignored).
- [ ] Scratch quota/retention policy — still to verify.

## Assumed job geometry (from md/scripts, adjust after confirming)
- `--gres=gpu:1`, `--cpus-per-task=8`, `--mem=32G`, `-pin on`, `-ntmpi 1 -ntomp 8`.
- Prep walltime ~4 h; 100 ns production walltime ~10--20 h on one A100 (record actual ns/day).
- GPU offload depends on the integrator: **EM (`steep`) = `-nb gpu` only** (`-bonded gpu`
  errors — needs a dynamical integrator); NVT/NPT (`md`) = `-nb gpu -pme gpu -bonded gpu`;
  production `-gpu_id 0 -pin on` (auto-assign).

## Run-dir layout
Per-PDB, per-state run dirs: `md/<PDB>/<state>/` (e.g. `md/1AKI/apo/`). Each owns its frozen
`configs/` (mdps, git-tracked) + `prep.sh`/`md.sh` (git-tracked). Scripts write **all** GROMACS
output into that run's `out/` (gitignored) and slurm logs into `logs/` — nothing lands in the
repo root. Structures stay shared in `structures/<PDB>.pdb`.

## Sync workflow
1. Stage inputs + frozen configs + sbatch here (git-tracked); commit + push.
2. `git pull` on Gemini. Do **not** `rsync --delete` toward Gemini data.
3. **Submit from the repo root** (`/scratch/bneff/bcell_epitope`) — scripts anchor paths to
   `$SLURM_SUBMIT_DIR`; Slurm spools the script, so `$0`-relative paths fail.
4. `sbatch md/1AKI/apo/prep.sh` → equilibrated system in `md/1AKI/apo/out/` (stale outputs from a
   failed run are overwritten; `prep.sh` uses `-nobackup`). Check `pdb2gmx.log` reports 4 SS bonds.
5. `sbatch md/1AKI/apo/md.sh` → production; pull back `out/md.xtc` (compressed) for analysis.
6. **Move logs/outputs off scratch promptly** — scratch is purged (logs were lost once). Keep the
   small text logs (`logs/slurm_*.out`, `out/pdb2gmx.log`, `out/md.log`).
7. Record job ids + ns/day + commands in `manuscript/sections/runlog.tex`.
