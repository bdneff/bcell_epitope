#!/bin/bash
# =============================================================================
# prep.sh — apo-HEL system preparation (bcell_epitope). Submit with: sbatch prep.sh ...
#   pdb2gmx -> editconf -> solvate -> genion -> EM -> NVT -> NPT
# Frozen stack: AMBER99SB-ILDN / TIP3P / dodecahedron 1.2 nm / 0.15 M NaCl
#               (see ../../docs/METHODS.md, signed off 2026-06-18)
#
# Usage (from this directory, after syncing repo to Gemini):
#   sbatch prep.sh ../../structures/1AKI.pdb
# Output: npt.gro, npt.cpt, topol.top  (inputs for md.sh)
# =============================================================================
#SBATCH -J hel_prep
#SBATCH -p gpu-a100              # confirmed on Gemini 2026-06-18 (4-day limit, A100 fastest)
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH -o slurm_prep_%j.out
#SBATCH -e slurm_prep_%j.err

set -euo pipefail

INPUT_PDB="${1:?Usage: sbatch prep.sh <input.pdb>}"
[[ -f "$INPUT_PDB" ]] || { echo "ERROR: $INPUT_PDB not found"; exit 1; }
# Slurm copies the batch script to a private spool dir, so $0/dirname won't find configs/.
# Anchor to the submit dir (repo root when you run `sbatch md/apo_hel/prep.sh ...`).
CFG="${SLURM_SUBMIT_DIR:-$PWD}/md/apo_hel/configs"
[[ -d "$CFG" ]] || { echo "ERROR: configs dir not found at $CFG (submit from repo root)"; exit 1; }

echo "=== apo-HEL prep | $(date) | node $(hostname) | job ${SLURM_JOB_ID} ==="
module load Gromacs              # Gemini module = "Gromacs" (capital G); GROMACS 2023.2-dev
gmx --version | grep -i version || true

check() { [[ -f "$1" ]] || { echo "ERROR: missing expected output $1 ($2)"; exit 1; }; echo "  -> OK: $1"; }

# --- 0. Clean: drop crystal waters (78 HOH in 1AKI); keep protein only -------
grep -v '^HETATM' "$INPUT_PDB" | grep -vE '^ATOM.{13}HOH' > clean.pdb || true
# (1AKI has no non-water HETATM; clean.pdb = protein chain A)

# --- 1. pdb2gmx --------------------------------------------------------------
# Disulfides: HEL has 4 (6-127, 30-115, 64-80, 76-94). GROMACS forms SS bonds
# automatically by SG-SG distance. VERIFY the log below reports 4 SS bonds.
gmx pdb2gmx -f clean.pdb -o protein.gro -p topol.top -i posre.itp \
    -water tip3p -ff amber99sb-ildn -ignh -nobackup 2>&1 | tee pdb2gmx.log
echo "=== Disulfide check (expect 4) ==="; grep -iE "disulfide|S-S|Linking .*CYS" pdb2gmx.log || \
    echo "  WARNING: no SS-bond lines matched — inspect pdb2gmx.log manually"
check protein.gro pdb2gmx; check topol.top pdb2gmx

# --- 2. Box: dodecahedron, 1.2 nm clearance ---------------------------------
gmx editconf -f protein.gro -o box.gro -c -d 1.2 -bt dodecahedron -nobackup
check box.gro editconf

# --- 3. Solvate (TIP3P) ------------------------------------------------------
gmx solvate -cp box.gro -cs spc216.gro -o solvated.gro -p topol.top -nobackup
check solvated.gro solvate

# --- 4. Ions: neutralize + 0.15 M NaCl --------------------------------------
gmx grompp -f "$CFG/ions.mdp" -c solvated.gro -p topol.top -o ions.tpr -maxwarn 1 -nobackup
printf "SOL\n" | gmx genion -s ions.tpr -o ions.gro -p topol.top \
    -pname NA -nname CL -neutral -conc 0.15 -nobackup
check ions.gro genion

# --- 5. Energy minimization --------------------------------------------------
gmx grompp -f "$CFG/em.mdp" -c ions.gro -p topol.top -o em.tpr -nobackup
# EM uses the 'steep' integrator; GPU bonded offload requires a dynamical integrator (md/sd),
# so NO -bonded gpu here (it errors). -nb gpu is fine. NVT/NPT/prod use md, so they keep it.
gmx mdrun -v -deffnm em -ntmpi 1 -ntomp 8 -nb gpu -pin on -nobackup
check em.gro "mdrun EM"

# --- 6. NVT (100 ps, restrained) --------------------------------------------
gmx grompp -f "$CFG/nvt.mdp" -c em.gro -r em.gro -p topol.top -o nvt.tpr -nobackup
gmx mdrun -v -deffnm nvt -ntmpi 1 -ntomp 8 -nb gpu -pme gpu -bonded gpu -pin on -nobackup
check nvt.gro "mdrun NVT"; check nvt.cpt "mdrun NVT"

# --- 7. NPT (100 ps, restrained, C-rescale) ---------------------------------
gmx grompp -f "$CFG/npt.mdp" -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -nobackup
gmx mdrun -v -deffnm npt -ntmpi 1 -ntomp 8 -nb gpu -pme gpu -bonded gpu -pin on -nobackup
check npt.gro "mdrun NPT"; check npt.cpt "mdrun NPT"

echo "=== prep complete $(date) — ready for: sbatch md.sh ==="
