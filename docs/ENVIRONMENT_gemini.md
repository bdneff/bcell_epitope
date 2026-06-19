# Compute environment — Gemini cluster

Stage and analyze in this repo; run production MD on **Gemini** (Slurm). Trajectories live on
Gemini and are gitignored. Settings below are inherited from the working TGen scripts in
`../../md/scripts/` and **must be confirmed on Gemini** before the first run.

## To confirm on Gemini (then pin in the run log)
- [ ] GROMACS module name + **version** (`module avail gromacs`; `gmx --version`).
      Working scripts use `module load gromacs` (also seen as `Gromacs` — confirm exact case).
- [ ] GPU partition name + QOS (TGen scripts use `-p gpu-a100`, `--gres=gpu:1`).
- [ ] Scratch/quota policy for trajectory storage; where run dirs should live.
- [ ] Path to sync this repo into (the "Gemini folder" this project syncs with).

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
