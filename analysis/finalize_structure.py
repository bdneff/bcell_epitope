#!/usr/bin/env python
"""Complete partial side chains on a grafted antigen chain with pdbfixer and write
the frozen MD input (heavy atoms; pdb2gmx -ignh rebuilds H). Disordered termini are
NOT extended -- we explicitly clear missingResidues so pdbfixer only fills missing
*atoms* on residues that already exist (manuscript sec:repair: termini are capped,
not modelled).

  micromamba run -n bcell-repair python analysis/finalize_structure.py \
      structures/repair/TissueFactor_1AHW_grafted.pdb structures/TissueFactor_1AHW_fixed.pdb
"""
import sys
from pdbfixer import PDBFixer
from openmm.app import PDBFile

inp, out = sys.argv[1], sys.argv[2]
fixer = PDBFixer(filename=inp)
fixer.findMissingResidues()
if fixer.missingResidues:
    print(f"  pdbfixer still sees missing residues {fixer.missingResidues} -- clearing "
          "(internal gaps were grafted from AF; termini are capped, not modelled)")
fixer.missingResidues = {}                      # do not extend termini / re-add grafted loops
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()
fixer.findMissingAtoms()
n_atoms = sum(len(v) for v in fixer.missingAtoms.values())
print(f"  completing {n_atoms} missing heavy atoms on {len(fixer.missingAtoms)} residues")
fixer.addMissingAtoms()
with open(out, 'w') as fh:
    PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
print(f"  wrote {out}")
