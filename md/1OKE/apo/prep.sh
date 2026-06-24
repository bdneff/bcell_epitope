#!/bin/bash
# =============================================================================
# prep.sh — apo Dengue-2 envelope ectodomain (sE; 1OKE chain A) system prep. Submit FROM THE REPO ROOT:  sbatch md/1OKE/apo/prep.sh
# Tier-2 (binary critical-residue labels). DENV2 sE monomer; native E is a homodimer.
# Frozen stack: AMBER99SB-ILDN / TIP3P / cubic 1.2 nm / 0.15 M NaCl. Expected disulfides: 6.
# Outputs -> md/1OKE/apo/out/ (gitignored); slurm logs -> md/1OKE/apo/logs/.
# =============================================================================
#SBATCH -J DengueE_apo_prep
#SBATCH -p gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH -o md/1OKE/apo/logs/slurm_prep_%j.out
#SBATCH -e md/1OKE/apo/logs/slurm_prep_%j.err

set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
RUNDIR="$ROOT/md/1OKE/apo"; CFG="$RUNDIR/configs"; OUT="$RUNDIR/out"
INPUT_PDB="$ROOT/structures/DengueE_1OKE.pdb"
[[ -d "$CFG" ]]       || { echo "ERROR: configs not found at $CFG (submit from repo root)"; exit 1; }
[[ -f "$INPUT_PDB" ]] || { echo "ERROR: $INPUT_PDB not found"; exit 1; }
mkdir -p "$OUT"; cd "$OUT"

echo "=== apo Dengue-2 envelope ectodomain (sE; 1OKE chain A) prep | $(date) | node $(hostname) | job ${SLURM_JOB_ID} ==="
module load Gromacs
gmx --version | grep -i version || true
check() { [[ -f "$1" ]] || { echo "ERROR: missing $1 ($2)"; exit 1; }; echo "  -> OK: $1"; }

cp "$INPUT_PDB" start.pdb
gmx pdb2gmx -f start.pdb -o protein.gro -p topol.top -i posre.itp \
    -water tip3p -ff amber99sb-ildn -ignh -nobackup 2>&1 | tee pdb2gmx.log
echo "=== Disulfide check (expect 6) ==="; grep -iE "disulfide|S-S|Linking .*CYS" pdb2gmx.log || echo "  (none found)"
check protein.gro pdb2gmx; check topol.top pdb2gmx

gmx editconf -f protein.gro -o box.gro -c -d 1.2 -bt cubic -nobackup; check box.gro editconf
gmx solvate -cp box.gro -cs spc216.gro -o solvated.gro -p topol.top -nobackup; check solvated.gro solvate
gmx grompp -f "$CFG/ions.mdp" -c solvated.gro -p topol.top -o ions.tpr -maxwarn 1 -nobackup
printf "SOL\n" | gmx genion -s ions.tpr -o ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15 -nobackup; check ions.gro genion
gmx grompp -f "$CFG/em.mdp"  -c ions.gro -p topol.top -o em.tpr -nobackup
gmx mdrun -v -deffnm em -ntmpi 1 -ntomp 8 -nb gpu -pin on -nobackup; check em.gro "mdrun EM"
gmx grompp -f "$CFG/nvt.mdp" -c em.gro  -r em.gro  -p topol.top -o nvt.tpr -nobackup
gmx mdrun -v -deffnm nvt -ntmpi 1 -ntomp 8 -nb gpu -pme gpu -bonded gpu -pin on -nobackup; check nvt.gro "mdrun NVT"; check nvt.cpt "mdrun NVT"
gmx grompp -f "$CFG/npt.mdp" -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -nobackup
gmx mdrun -v -deffnm npt -ntmpi 1 -ntomp 8 -nb gpu -pme gpu -bonded gpu -pin on -nobackup; check npt.gro "mdrun NPT"; check npt.cpt "mdrun NPT"
echo "=== apo Dengue-2 envelope ectodomain (sE; 1OKE chain A) prep complete $(date) — next: sbatch md/1OKE/apo/md.sh ==="
