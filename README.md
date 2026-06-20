# bcell_epitope

Predict B-cell epitope residues from **antigen-only** information. Residue labels are
alanine-scan ΔΔG_bind (binding importance); the question is whether antigen dynamics /
solvation / energetics predict it **beyond static surface exposure**.

- **Living record:** `manuscript/main.pdf` (build: `cd manuscript && latexmk -pdf main.tex`). Read top to bottom.
- **How to work here:** `CLAUDE.md` (operating contract). Methods: `docs/METHODS.md`. Compute: `docs/ENVIRONMENT_gemini.md`.
- **Labels:** `benchmark/antigen_alanine_scanning_benchmark_v2.xlsx` (audited; 163-row clean core).
- **First system:** apo hen egg-white lysozyme (`structures/1AKI.pdb`); run staged in `md/1AKI/apo/`.

Stage + analyze here (git); run on Gemini (Slurm); trajectories live on Gemini (gitignored).
