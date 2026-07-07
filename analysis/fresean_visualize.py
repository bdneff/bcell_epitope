#!/usr/bin/env python
"""
fresean_visualize.py — regenerate FRESEAN mode figures and animations.

Built with Claude Science (session #4, 2026-07-07).

Produces, for one system's extracted modes:
  1. fresean_mode_verification.png  — rigid-body fraction + eigenvalue spectrum
  2. cartoon_mode{NN}.gif           — PyMOL cartoon animations (needs pymol env)
  3. fresean_modes_grid.gif         — 3x3 labeled grid (green=trans, red=rot, blue=internal)
  4. <system>_modes.nmd             — NMWiz file for local VMD (all 30 modes)

Two-env workflow (kernels don't share state):
  - PyMOL rendering (steps 2-3 frames): env 'pymol' (pymol-open-source); Pillow NOT there.
  - GIF assembly + matplotlib figures: env 'python' (has PIL via matplotlib).

INPUTS: fresean_extract/<rep>/evec_freq1_mode1-30_cg.xyz, ref-cg.gro,
        eval_matrix_cg.mmat.dat, and the all-atom ref.pdb (1960 atoms for 1AKI).

MODE DISPLACEMENT (cartoon): displace each residue's all-atom coords RIGIDLY by
that residue's backbone-bead eigenvector (modes live in 2-bead CG space, so side
chains ride with their residue). Amplitude scaled for visibility, not physical.

See analysis/fresean_involvement.py for the parsers (read_modes, parse_cg_gro,
rigid_fractions, trans_rot_split) — imported here to avoid duplication.
"""
import os, sys, glob
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fresean_involvement import read_modes, parse_cg_gro, rigid_fractions, trans_rot_split

def gro_xyz(path):
    L = open(path).readlines(); n = int(L[1])
    return np.array([[float(L[2+i][20:28]), float(L[2+i][28:36]), float(L[2+i][36:44])]
                     for i in range(n)])

def back_indices(beads):
    return np.array([i for i, b in enumerate(beads) if b["atom"] == "BACK"])

def verification_figure(mode_xyz, ref_gro, eig_dat, out="fresean_mode_verification.png", title=""):
    """Rigid-body fraction (per mode) + eigenvalue spectrum. Needs figure-style loaded."""
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    modes = read_modes(mode_xyz); beads = parse_cg_gro(ref_gro); ref = gro_xyz(ref_gro)
    rf = rigid_fractions(modes, ref)
    ev = np.array([float(x) for x in open(eig_dat).read().split()])
    ndof = len(beads) * 3; ev_block = ev[:ndof]
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11, 4.2))
    x = np.arange(1, len(modes)+1)
    axL.bar(x[:6], rf[:6], color="#c0392b", label="rigid-body (1-6)", zorder=3)
    axL.bar(x[6:], rf[6:], color="#4a78b5", label="internal (7+)", zorder=3)
    axL.axvline(6.5, color="#888", ls="--", lw=1)
    axL.set_xlabel("FRESEAN mode index"); axL.set_ylabel("rigid-body fraction")
    axL.set_ylim(0, 1.02); axL.legend(fontsize=8.5, frameon=False)
    axL.set_title("Modes 1-6 are translation + rotation", fontsize=11, loc="left")
    axR.semilogy(np.arange(1, 31), ev_block[:30], "o-", color="#16408a", ms=5)
    axR.axvline(6.5, color="#888", ls="--", lw=1)
    axR.set_xlabel("eigenvalue index (large->small)"); axR.set_ylabel("eigenvalue")
    axR.set_title("Eigenvalue spectrum", fontsize=11, loc="left")
    fig.suptitle(title or "FRESEAN mode verification", fontsize=12, y=1.0)
    fig.tight_layout(rect=[0, 0, 1, 0.95]); fig.savefig(out, dpi=200, bbox_inches="tight")
    return dict(rigid_first6=float(rf[:6].mean()), rigid_rest=float(rf[6:].mean()),
                trans_rot=trans_rot_split(modes, ref))

def write_nmd(mode_xyz, ref_gro, out, title="FRESEAN_modes"):
    """NMWiz .nmd (one pseudo-CA/residue = BACK bead) for local VMD `vmd -e out.nmd`."""
    modes = read_modes(mode_xyz); beads = parse_cg_gro(ref_gro)
    bi = back_indices(beads); ref = gro_xyz(ref_gro)
    coords = ref[bi] * 10.0  # nm -> Angstrom
    resids = [beads[i]["resid"] for i in bi]; resn = [beads[i]["resname"] for i in bi]
    with open(out, "w") as f:
        f.write("nmwiz_load %s\n" % os.path.basename(out))
        f.write("title %s\n" % title)
        f.write("names " + " ".join(["CA"]*len(resids)) + "\n")
        f.write("resnames " + " ".join(resn) + "\n")
        f.write("resids " + " ".join(map(str, resids)) + "\n")
        f.write("chids " + " ".join(["A"]*len(resids)) + "\n")
        f.write("coordinates " + " ".join("%.3f" % c for c in coords.reshape(-1)) + "\n")
        for mi in range(len(modes)):
            v = modes[mi][bi]
            f.write("mode %d %s\n" % (mi+1, " ".join("%.4f" % c for c in v.reshape(-1))))
    return out

def write_multistate_pdb(mode_idx, ref_pdb, mode_xyz, ref_gro, out, amp=9.0, nframes=20):
    """All-atom multi-state PDB sweeping one mode (for PyMOL cartoon animation)."""
    modes = read_modes(mode_xyz); beads = parse_cg_gro(ref_gro); bi = back_indices(beads)
    bres = [beads[i]["resid"] for i in bi]
    res2idx = {r: k for k, r in enumerate(bres)}
    lines = [l for l in open(ref_pdb) if l.startswith("ATOM")]
    ares = np.array([int(l[22:26]) for l in lines])
    axyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in lines])
    vb = modes[mode_idx][bi]
    disp = np.array([vb[res2idx[r]] for r in ares])
    with open(out, "w") as f:
        for fr in range(nframes):
            s = amp * np.sin(2*np.pi*fr/nframes)
            f.write("MODEL     %d\n" % (fr+1))
            for k, l in enumerate(lines):
                nx, ny, nz = axyz[k] + s*disp[k]
                f.write("%s%8.3f%8.3f%8.3f%s" % (l[:30], nx, ny, nz, l[54:]))
            f.write("ENDMDL\n")
    return out

# PyMOL render + GIF assembly are driven by render_modes.py (pymol env) and
# assemble_gifs (python env) — see the __main__ block for the 1AKI recipe.

if __name__ == "__main__":
    print("Import this module; see docstring. 1AKI recipe:")
    print("  from fresean_visualize import verification_figure, write_nmd, write_multistate_pdb")
