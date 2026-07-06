# Dengue E (1OKE) crystal B-factor — provenance caveat

The crystal B-factor track in the 1OKE epitope-feature profile is from the
**deposited 1OKE crystal structure**, which is the DENV-2 soluble E ectodomain
crystallized WITH the detergent **n-octyl-beta-D-glucoside (beta-OG)** bound in
the fusion-loop hydrophobic pocket. It is NOT an antibody complex and NOT a
ligand-free apo structure.

Implications for interpreting B-factor as "flexibility":
- Crystal-packing contacts damp surface B-factors (all crystallographic B-factors).
- The bound beta-OG ligand rigidifies the fusion-loop / kl-hairpin region
  (residues ~98-110), which overlaps the epitope-critical patch — so low
  B-factor there may be ligand-induced, not intrinsic.
- Every OTHER feature in the profile is from the APO MD trajectory; the B-factor
  row is the only crystallographic (holo/ligand-bound) quantity — an apo-vs-holo
  mismatch analogous to the Tier-1 apo-crystal issue.

Cleaner apples-to-apples flexibility metric = per-residue RMSF from the apo MD
trajectory (same simulation as the other features).
