#!/usr/bin/env python
"""Repair an incomplete crystal antigen chain for apo-MD by grafting missing
internal loops from an AlphaFold model, then completing side chains with pdbfixer.

Strategy (manuscript sec:repair): the crystal supplies all experimentally resolved
coordinates; AlphaFold supplies ONLY the unresolved internal residues. We first
check the AF model agrees with the crystal where both exist (CA-RMSD over matched,
identical residues) -- the "similar enough in matched regions" test -- then for each
internal gap we locally superpose AF onto the crystal on the flanking resolved
residues and splice the AF loop in at crystal numbering. Disordered N-/C-termini
are NOT modelled (they bear no labels); they are left for capping.

Run inside the bcell-repair env:
  micromamba run -n bcell-repair python analysis/repair_structure.py \
      --crystal structures/hGH_1HGU.pdb --chain A \
      --af structures/repair/AF-P01241.pdb --out structures/hGH_1HGU_fixed.pdb
"""
import argparse
import sys
import numpy as np
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Superimposer
from Bio.PDB.Polypeptide import is_aa
from Bio.Align import PairwiseAligner


def load_first_model(path, sid):
    """Parse PDB or mmCIF (AlphaFold Server returns CIF) and return model 0."""
    parser = MMCIFParser(QUIET=True) if path.lower().endswith('.cif') else PDBParser(QUIET=True)
    return parser.get_structure(sid, path)[0]

THREE2ONE = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
    'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
    'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
}


def chain_residues(chain):
    """Ordered list of standard amino-acid residues (skip hetero/water)."""
    out = []
    for res in chain:
        if res.id[0] == ' ' and res.resname in THREE2ONE:
            out.append(res)
    return out


def seq_of(residues):
    return ''.join(THREE2ONE[r.resname] for r in residues)


def map_offset(crys_res, af_res):
    """Align crystal vs AF sequence; return integer offset k s.t. af_resnum ~= crys_resnum + k,
    plus the list of (crys_res, af_res) matched (identical) residue pairs."""
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    cseq, aseq = seq_of(crys_res), seq_of(af_res)
    aln = aligner.align(cseq, aseq)[0]
    pairs = []
    offsets = []
    # aln.aligned gives matched block index ranges in each sequence
    for (c0, c1), (a0, a1) in zip(aln.aligned[0], aln.aligned[1]):
        for ci, ai in zip(range(c0, c1), range(a0, a1)):
            cr, ar = crys_res[ci], af_res[ai]
            if cr.resname == ar.resname:
                pairs.append((cr, ar))
                offsets.append(ar.id[1] - cr.id[1])
    if not offsets:
        sys.exit("ERROR: no identical aligned residues between crystal and AF model")
    # robust constant offset = mode
    vals, counts = np.unique(offsets, return_counts=True)
    k = int(vals[np.argmax(counts)])
    frac = counts.max() / len(offsets)
    print(f"  numbering offset af = crys + {k}  (consistent for {frac*100:.0f}% of matched residues)")
    return k, pairs


def ca(res):
    return res['CA'] if res.has_id('CA') else None


def superpose_rms(crys_pairs):
    fixed = [ca(c) for c, a in crys_pairs if ca(c) and ca(a)]
    moving = [ca(a) for c, a in crys_pairs if ca(c) and ca(a)]
    sup = Superimposer()
    sup.set_atoms(fixed, moving)
    return sup, sup.rms, len(fixed)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--crystal', required=True)
    ap.add_argument('--chain', required=True)
    ap.add_argument('--af', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--flank', type=int, default=6, help='resolved residues each side for local fit')
    ap.add_argument('--max-gap', type=int, default=20, help='refuse to graft gaps longer than this')
    args = ap.parse_args()

    crys_model = load_first_model(args.crystal, 'x')
    af_chain = chain_residues(list(load_first_model(args.af, 'a'))[0])
    crys_chain_obj = crys_model[args.chain]
    crys_res = chain_residues(crys_chain_obj)

    print(f"crystal {args.crystal} chain {args.chain}: {len(crys_res)} residues "
          f"({crys_res[0].id[1]}..{crys_res[-1].id[1]})")
    print(f"AF model {args.af}: {len(af_chain)} residues")

    k, pairs = map_offset(crys_res, af_chain)
    af_by_num = {r.id[1]: r for r in af_chain}

    # --- global agreement over matched residues (the "similar enough" test) ---
    _, global_rms, n = superpose_rms(pairs)
    print(f"  GLOBAL CA-RMSD (AF vs crystal, {n} matched residues): {global_rms:.2f} A")
    if global_rms > 3.0:
        print("  WARNING: AF and crystal disagree globally (>3 A) -- AF loops are NOT trustworthy here.")

    # --- find internal gaps (missing residues between resolved ones) ---
    present = sorted(r.id[1] for r in crys_res)
    gaps = []
    for a, b in zip(present, present[1:]):
        if b - a > 1:
            gaps.append((a + 1, b - 1))
    if not gaps:
        print("  no internal gaps; only side-chain completion needed")
    grafted = []
    for (g0, g1) in gaps:
        n_missing = g1 - g0 + 1
        if n_missing > args.max_gap:
            print(f"  gap {g0}-{g1} ({n_missing} res): TOO LONG (> {args.max_gap}) -- not grafting "
                  f"(needs a full predicted model, not a spliced loop)")
            continue
        # flank residues present in both crystal and AF
        flank_nums = [g0 - i for i in range(1, args.flank + 1)] + [g1 + i for i in range(1, args.flank + 1)]
        fl_pairs = []
        for cn in flank_nums:
            af = af_by_num.get(cn + k)
            cr = next((r for r in crys_res if r.id[1] == cn), None)
            if cr and af and cr.resname == af.resname and ca(cr) and ca(af):
                fl_pairs.append((cr, af))
        if len(fl_pairs) < 4:
            print(f"  gap {g0}-{g1}: too few flank anchors ({len(fl_pairs)}) -- skipped")
            continue
        sup, lrms, nfl = superpose_rms(fl_pairs)
        # transform AF gap residues into crystal frame and splice in
        added = 0
        for gn in range(g0, g1 + 1):
            af_res = af_by_num.get(gn + k)
            if af_res is None:
                print(f"    AF missing residue for crystal {gn} -- gap incompletely modelled")
                continue
            sup.apply(list(af_res.get_atoms()))
            af_res.detach_parent()               # detach before renumber (avoids id-collision warning)
            af_res.id = (' ', gn, ' ')           # renumber to crystal numbering
            crys_chain_obj.add(af_res)
            added += 1
        print(f"  gap {g0}-{g1} ({n_missing} res): local CA-RMSD {lrms:.2f} A on {nfl} flanks; grafted {added}")
        grafted.append((g0, g1, added, lrms))

    # reorder residues by number so the chain is sequential
    crys_chain_obj.child_list.sort(key=lambda r: r.id[1])
    crys_chain_obj.child_dict = {r.id: r for r in crys_chain_obj.child_list}

    io = PDBIO()
    io.set_structure(crys_model)
    io.save(args.out)
    print(f"  wrote {args.out}  ({len(chain_residues(crys_chain_obj))} residues)")
    print("  NOTE: run pdbfixer --add-atoms=heavy on this file to complete partial side chains.")


if __name__ == '__main__':
    main()
