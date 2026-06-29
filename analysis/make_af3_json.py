#!/usr/bin/env python
"""Build an AlphaFold3 input JSON from a single-chain FASTA, matching the format the
Altin-lab AF3 pipeline uses on Gemini (see episcaf episcaf_pipeline/stages/stage04).

Two modes:
  --mode pipeline   (default) emit sequence only -> AF3 builds a real MSA from the
                    installed databases (run WITHOUT --norun_data_pipeline). Best for a
                    NATURAL protein like gB.
  --mode single     emit a self-sequence MSA + empty templates -> single-sequence
                    prediction (run WITH --norun_data_pipeline). Mirrors the lab's
                    de-novo-design recipe exactly; faster, lower accuracy for naturals.

  python analysis/make_af3_json.py structures/repair/gB_5C6T_ectodomain.fasta \
      --name gB_5C6T --out md/5C6T/apo/af3/gB_5C6T.json --mode pipeline
"""
import argparse
import json
from pathlib import Path


def read_fasta_seq(path):
    seq = []
    for line in Path(path).read_text().splitlines():
        if line.startswith('>') or not line.strip():
            continue
        seq.append(line.strip())
    return ''.join(seq)


def payload(name, seq, seed, mode):
    if mode == 'single':
        msa = f">query\n{seq}\n"
        protein = {"id": "A", "sequence": seq,
                   "unpairedMsa": msa, "pairedMsa": msa, "templates": []}
    else:  # pipeline: let AF3 search the databases for a real MSA
        protein = {"id": "A", "sequence": seq}
    return {"name": name, "modelSeeds": [seed],
            "sequences": [{"protein": protein}],
            "dialect": "alphafold3", "version": 1}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('fasta')
    ap.add_argument('--name', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--mode', choices=['pipeline', 'single'], default='pipeline')
    args = ap.parse_args()
    seq = read_fasta_seq(args.fasta)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload(args.name, seq, args.seed, args.mode), indent=2))
    print(f"wrote {out}  (name={args.name}, {len(seq)} residues, mode={args.mode}, seed={args.seed})")


if __name__ == '__main__':
    main()
