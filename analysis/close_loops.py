#!/usr/bin/env python
"""Close grafted-loop junctions by restrained energy minimization (OpenMM, implicit
solvent). All atoms outside the freed loop window are held by a stiff position
restraint, so every crystallographic coordinate and disulfide is preserved; only the
loop and a couple of flanking residues relax to close the spliced peptide bonds.
A strong harmonic C(i)-N(i+1) restraint is added across any junction left open by the
graft, so the loop is pulled to a real peptide-bond length. This is local geometry
repair, not sampling.

NOTE: OpenMM's minimizeEnergy ignores particle mass (mass only fixes atoms in
dynamics), so we restrain positions explicitly and use constraints=None.

  micromamba run -n bcell-repair python analysis/close_loops.py \
      structures/TissueFactor_1AHW_fixed.pdb structures/repair/TF_closed.pdb 81-92
"""
import sys
import numpy as np
from openmm.app import PDBFile, ForceField, Modeller, CutoffNonPeriodic, Simulation
from openmm import LangevinIntegrator, CustomExternalForce, CustomBondForce
from openmm.unit import kelvin, picosecond, femtosecond, nanometer, kilojoule_per_mole

inp, out, free_spec = sys.argv[1], sys.argv[2], sys.argv[3]
free = set()
for part in free_spec.split(','):
    a, b = part.split('-')
    free.update(range(int(a), int(b) + 1))

pdb = PDBFile(inp)
ff = ForceField('amber99sbildn.xml', 'amber99_obc.xml')
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(ff)
system = ff.createSystem(modeller.topology, nonbondedMethod=CutoffNonPeriodic,
                         nonbondedCutoff=1.5 * nanometer, constraints=None)
pos_nm = np.array(modeller.positions.value_in_unit(nanometer))

# (1) hold every atom outside the freed loop window with a stiff position restraint
restraint = CustomExternalForce('0.5*kr*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
restraint.addGlobalParameter('kr', 500000.0)     # kJ/mol/nm^2 -- effectively immobile
for p in ('x0', 'y0', 'z0'):
    restraint.addPerParticleParameter(p)
n_free = 0
for atom in modeller.topology.atoms():
    if int(atom.residue.id) in free:
        n_free += 1
    else:
        x, y, z = pos_nm[atom.index]
        restraint.addParticle(atom.index, [x, y, z])
system.addForce(restraint)
print(f"  freed {n_free} atoms in residues {sorted(free)[0]}..{sorted(free)[-1]}; rest position-restrained")

# (2) close any junction the graft left open (>1.6 A)
closer = CustomBondForce('0.5*kc*(r-r0)^2')
closer.addGlobalParameter('kc', 300000.0)
closer.addGlobalParameter('r0', 0.133)
residues = list(modeller.topology.residues())
n_closed = 0
for r1, r2 in zip(residues, residues[1:]):
    if r1.chain.id == r2.chain.id and int(r2.id) - int(r1.id) == 1:
        c = next((a for a in r1.atoms() if a.name == 'C'), None)
        n = next((a for a in r2.atoms() if a.name == 'N'), None)
        if c and n:
            d = np.linalg.norm(pos_nm[c.index] - pos_nm[n.index])
            if d > 0.16:
                closer.addBond(c.index, n.index, [])
                n_closed += 1
                print(f"  closing open junction {r1.id}-{r2.id} ({d*10:.2f} A -> 1.33 A)")
if n_closed:
    system.addForce(closer)

integ = LangevinIntegrator(300 * kelvin, 1 / picosecond, 1 * femtosecond)
sim = Simulation(modeller.topology, system, integ)
sim.context.setPositions(modeller.positions)
e0 = sim.context.getState(getEnergy=True).getPotentialEnergy()
sim.minimizeEnergy(maxIterations=10000, tolerance=1.0 * kilojoule_per_mole / nanometer)
e1 = sim.context.getState(getEnergy=True).getPotentialEnergy()
print(f"  energy {e0} -> {e1}")
pos = sim.context.getState(getPositions=True).getPositions()
with open(out, 'w') as fh:
    PDBFile.writeFile(modeller.topology, pos, fh, keepIds=True)
print(f"  wrote {out}")
