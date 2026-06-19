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
- GPU offload: EM/equil `-nb gpu -pme gpu -bonded gpu`; production `-gpu_id 0 -pin on`.

## Sync workflow
1. Stage inputs + frozen configs + sbatch here (git-tracked).
2. Sync repo → Gemini folder (rsync/git). Do **not** `rsync --delete` toward Gemini data.
3. `sbatch md/apo_hel/prep.sh structures/1AKI.pdb` → equilibrated system.
4. `sbatch md/apo_hel/md.sh` → production; pull back `.xtc` (compressed) for analysis.
5. Record job ids + commands in `manuscript/sections/runlog.tex`.
