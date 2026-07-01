# Apo-antigen features for antibody-hotspot prediction — a note for the lab

The question behind all of this: can antibody-free (apo) molecular dynamics of an antigen tell us
which of its surface residues an antibody will care about — before we know the antibody, and beyond
what surface exposure alone already says? We have apo trajectories for six antigens (lysozyme, HPr,
tissue factor, VEGF, human growth hormone, IFN-gamma receptor) and per-residue alanine-scan binding
energies (ddG) as the label. This note is the state of that search, and one result I want to flag.

**How we score, and why it matters.** The honest number here is not "does the feature correlate over
all residues" — that is dominated by the trivial fact that buried residues and surface residues differ
in everything. We restrict to *surface* residues first, which removes the surface-versus-core confound,
and then ask a directional question: in that surface population, do hotspots score higher or lower on
the feature? And we compute it per antigen — each antigen its own baseline — because residues within a
fold share a trajectory, so our effective N is the number of antigens (six), not residues. This is the
same leave-one-antigen-out logic as the rest of the project.

**Exposure is confirmed useless on this label — which is the point.** Once we are inside the surface
population, trajectory SASA separates hotspots from other surface residues at chance (AUROC ~0.51). Raw
exposure only ever "found the outside." That is exactly the premise we set out to test: on an
energy-based label, the answer does not reduce to solvent accessibility, the way it does for the
distance-based labels everyone else uses.

**The physical picture that emerges is coherent, even though every piece is weak.** Across features that
are more-or-less independent, the hotspots come out as *anchored, pre-organized, desolvated* patches.
They are less hydrated (fewer hydrogen bonds to water), they sit at the quiet nodes of the slowest
collective motion (low PC1 participation — decoupled from the soft mode), their backbones are stiffer,
and yet they are more tightly bonded *internally* (more intra-protein hydrogen bonds, higher dynamic
cross-correlation, more favorable intra-protein interaction energy). These are not four independent
votes — being internally well-bonded is what makes a residue rigid, which is what keeps it off the soft
mode. And it lines up with the entropy argument already in the manuscript: a pre-organized, rigid anchor
pays little conformational entropy when the antibody arrives, so it can carry more of the binding energy.

**Reading the label as a gradient instead of a threshold changes the ranking — and this is worth
internalizing.** If we score against the *continuous* ddG (Spearman) rather than a binary hotspot cutoff,
the hydration signal is the strongest single feature (per-antigen rho = -0.48), and the intra-protein
interaction energy — which had *led* the binary hotspot test at AUROC 0.65 — goes essentially flat on the
gradient (rho ~ +0.08). The interpretation is clean: energy flags the handful of extreme hotspots at the
top, but it does not track the whole binding-energy gradient; hydration does. Binary labeling would have
hidden this from us entirely.

**And it is hydration, not hydrophobicity.** Static residue-type hydrophobicity (Kyte-Doolittle) is flat
(rho ~ 0.00). The signal is not "hotspots are hydrophobic amino acids" — it is that hotspot residues
actually make fewer hydrogen bonds to water *in their structural context in the trajectory*, which the
sequence cannot see. It correlates with SASA at only +0.19, so it is genuine context-dependent
desolvation, not exposure wearing a different hat.

**The result I want to flag: combining hydration with collective motion does better than either alone.**
The two strongest low-direction features — hydrogen bonds to water (desolvation) and PC1 participation
(collective rigidity) — carry partly independent information. Their sign-aligned, equal-weight sum, with
no fitting at all, tracks the graded ddG at per-antigen rho = -0.60, above hydration alone (-0.48) and
PC1 alone (-0.32), and it is negative in all four antigens we can score. Because the weights are fixed —
not fit to the data — this is not the overfitting that sent a fitted logistic *below* chance at N=6; it
is just averaging two real, complementary signals. Adding a third feature makes it worse, which tells me
these two axes (desolvation and rigidity) are specifically the ones carrying the information.

**Where this actually sits.** Everything above is weak. N is six antigens, and the per-antigen numbers
rest on four of them; the ordering among the middling features is noise. The honest bottleneck is not the
model — no fitted combination beat the best hand-chosen fixed one — it is the size of the label set, which
is why scaling up the labels (the StaB-ddG pseudo-label plan) matters more than any modeling cleverness
right now. Two caveats specifically: PC1 participation is a Cartesian-PCA placeholder that partly
re-measures flexibility (it correlates with RMSF at +0.52), so the "hotspots at the nodes" story only
becomes real if it survives the FRESEAN velocity-mode analysis (running now) — the proper operator, and
our lab's method. And the interaction energies are intramolecular: how tightly a residue is held by its
own fold, with no antibody in the calculation at all.

---

**The three things attached:**
- *Leaderboard* (epitope_feature_direction.png, epitope_feature_distributions.png) — every feature scored
  on surface-restricted directional separation of hotspots; the distributions version shows the actual
  hotspot-vs-rest densities so you can see the (weak) separation rather than trust a number.
- *Cheat sheet* (feature_cheatsheet.pdf) — the reference for what each feature is: the exact math, a
  plain-language description, how it is computed from the trajectory, the hotspot direction, and the
  caveats. Includes how MD geometrically defines a hydrogen bond, and the binary-vs-graded scoring point.
- *Combination* (combined_cv_ddg.png) — the hydration + collective-motion collective variable.

*Built with Claude Science; features computed on Gemini, scoring leave-one-antigen-out.*
