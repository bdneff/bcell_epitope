# Frozen methods — apo-HEL run #1 (SIGNED OFF 2026-06-18)

Project-local, frozen record of the methods stack for the first simulation. Purpose of run #1:
**pipeline validation on Gemini + first flexibility/RMSF feature** (the Westhof re-test). The
house default stack (`../../md/docs/METHODS_DEFAULTS.md`) was signed off as the starting point;
deviations and their rationale are recorded here. Specialized FRESEAN / 3D-2PT runs are a
*separate* protocol (deferred — see "Deferred").

## System
- Antigen: hen egg-white lysozyme (HEL), **apo**, from **RCSB 1AKI** (1.5 Å X-ray, chain A,
  residues 1–129). Downloaded 2026-06-18 from `https://files.rcsb.org/download/1AKI.pdb`.
- Verified: 1AKI sequence matches the benchmark WT identity at all 23 labeled positions →
  labels map 1:1 onto 1–129 numbering. (Hotspots K96/K97/Y20 confirmed.)
- Apo crystal chosen (not antigen-stripped-from-complex) to match the prediction scenario;
  HEL is rigid (apo≈bound) so the choice is low-risk.

## Signed-off stack (run #1)
| Knob | Choice | vs house default | Rationale |
|------|--------|------------------|-----------|
| Protein FF | AMBER99SB-ILDN | same | signed-off default; well-validated for HEL |
| Water | TIP3P | same | OK for pipeline validation + relative RMSF; **flagged** — upgrade to TIP4P/2005 for dynamics-grade features |
| Box | dodecahedron, 1.2 nm clearance | reconciled (house scripts disagreed cubic/dodec) | efficient for globular HEL |
| Salt | neutralize + 0.15 M NaCl | same | physiological-ish |
| Integrator/dt | leapfrog md, 2 fs | same | constraints on h-bonds permit 2 fs |
| Constraints | h-bonds, LINCS | same | standard; **flagged** — remove for FRESEAN/2PT runs |
| Electrostatics | PME, rcoulomb 1.0, fourier 0.16 | same | |
| vdW | rvdw 1.0, potential-shift-Verlet, DispCorr EnerPres | same | |
| Thermostat | V-rescale, 300 K, tau_t 0.1, groups Protein/Non-Protein | same | single bath OK (equilibrium run, not non-eq) |
| Barostat | C-rescale for NPT equil; Parrinello-Rahman for production | **deviation** | PR rings if started far from equilibrium; C-rescale is the modern equilibration choice (fixes a house-script issue) |
| Equilibration | EM → 100 ps NVT (posre) → 100 ps NPT (posre) | same | |
| Production length | 100 ns, single trajectory | documented | one long run for RMSF statistics (Brandon, 2026-06-18); replicas deferred |
| Seeds | gen_seed = -1, ld_seed = -1 (random) | same | Brandon: bitwise repro is not the goal. Reproducibility = git-tracked inputs + documented protocol + saved derived outputs; random seeds also give genuinely independent replicas. If a trr is purged from scratch, re-pull and re-run. |
| Output | xtc + energy + log every 10 ps; no trr | same | sufficient for RMSF |

## HEL-specific prep checkpoints (do not gloss over)
- **Disulfides:** HEL has 4 (Cys6–Cys127, Cys30–Cys115, Cys64–Cys80, Cys76–Cys94). pdb2gmx must
  form all 4 (CYS2). **Verify the pdb2gmx log reports 4 SS bonds** before proceeding — the fold
  depends on them.
- **Protonation:** default pdb2gmx at pH 7. HEL His15 → check the assigned state (HID/HIE);
  default is acceptable but record what was assigned.
- **Termini:** default charged (NH3+ / COO-).
- **Crystal waters:** the 78 HOH in 1AKI are removed; fresh TIP3P added by solvate. (Record.)
- **`-maxwarn`:** read every grompp warning; do not let maxwarn hide a real one.

## Deferred (separate protocol, not run #1)
- **FRESEAN modes / 3D-2PT:** require high-frequency velocity output (nstvout every few fs over
  short windows), constraint removal / sub-fs dt, and a water model validated for dynamics
  (TIP4P/2005). Pull exact settings from the FRESEAN + energy-transfer manuscripts when we set
  this up. This is the run that powers the §12.4 mode-participation and §12.2/12.6 hydration
  features — the project's distinctive contribution.
- Replicas (≥3) + RMSF convergence test for the production flexibility feature.
