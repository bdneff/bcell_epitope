#!/usr/bin/env python3
"""Mine IEDB for antigen-side, mutation-derived (NON-proximity) epitope labels.

Why: standard B-cell epitope labels are proximity/contact-defined (a residue is
"epitope" if it is near the antibody in a structure). This project wants
*energetic* importance instead. IEDB's structure-linked conformational epitopes
are ~entirely contact-derived (of 4,164 discontinuous epitopes with a pdb_id,
>99% are x-ray/cryo-EM contact), so we deliberately do NOT use the pdb_id-linked
records. Instead we pull *functional* conformational epitopes (defined by binding
assays) and, within them, the FLOW-CYTOMETRY subset = Integral Molecular / Doranz
"shotgun mutagenesis" alanine scans, whose residue sets ARE the critical residues
(binary critical/non-critical labels, energetically grounded).

Source: IEDB next-gen query API (PostgREST, no auth), https://query-api.iedb.org
  table bcell_search; structure_type='Discontinuous peptide'. Residues are in
  `structure_description` (UniProt numbering -> needs SIFTS mapping to a PDB for MD).

Output: benchmark/iedb_shotgun_candidates_v1.csv  (one row per antigen)
Run:    python3 analysis/fetch_iedb_shotgun.py
"""
import urllib.request, urllib.parse, json, collections, re, csv, pathlib

API = "https://query-api.iedb.org/"
OUT = pathlib.Path("benchmark/iedb_shotgun_candidates_v1.csv")
DISC = "structure_type=eq." + urllib.parse.quote("Discontinuous peptide")
COLS = ("bcell_id,assay_names,parent_source_antigen_name,"
        "parent_source_antigen_source_org_name,structure_description,pubmed_id")
# assays that mean the epitope was defined by STRUCTURAL CONTACT (proximity) -> reject
STRUCTURAL = ("x-ray", "electron microscopy", "3d structure", "crystallog",
              "angstrom", "nmr ", "modeling", "docking")
# antigens already in our benchmark -> flag as known
HAVE = ("lysozyme", "hpr", "vegf", "matriptase", "mt-sp1", "botulinum",
        "growth hormone", "interferon gamma", "tissue factor", "glycoprotein b")


def get(path):
    req = urllib.request.Request(API + path, headers={"Accept": "application/json"})
    return json.load(urllib.request.urlopen(req, timeout=120))


def pull_all_discontinuous():
    rows, off = [], 0
    while True:
        page = get(f"bcell_search?select={COLS}&{DISC}&order=bcell_id&limit=1000&offset={off}")
        rows += page
        if len(page) < 1000:
            return rows
        off += 1000


def sim_category(name):
    n = name.lower()
    if any(k in n for k in ("spike", "hemagglutinin")): return "viral glycoprotein (large/trimeric)"
    if "polyprotein" in n: return "flavi/alphavirus E-protein (soluble ectodomain avail)"
    if "envelope" in n: return "viral envelope glycoprotein"
    if any(k in n for k in ("chemokine", "gro-", "cxcl", "interleukin", "growth-regulated")):
        return "soluble chemokine/cytokine (EASY)"
    if any(k in n for k in ("hla", "histocompatibility", "mhc")): return "MHC/HLA (simulatable)"
    if any(k in n for k in ("claudin", "thyrotropin", "gpcr", "tetraspanin", "cd20", "cd81")):
        return "membrane protein (hard; check)"
    if "receptor" in n: return "receptor (check soluble ectodomain)"
    return "other / check per-antigen"


def main():
    rows = pull_all_discontinuous()
    functional = [r for r in rows
                  if r["assay_names"] and not any(k in r["assay_names"].lower() for k in STRUCTURAL)]
    shot = [r for r in functional if "flow cytometry" in r["assay_names"].lower()]
    print(f"discontinuous={len(rows)}  functional={len(functional)}  shotgun(flow cytometry)={len(shot)}")

    ag = collections.defaultdict(lambda: {"epi": 0, "refs": set(), "res": set(), "org": "", "uni": ""})
    for r in shot:
        name = (r["parent_source_antigen_name"] or "").strip()
        if not name:
            continue
        m = re.search(r"UniProt:(\w+)", name)
        base = re.sub(r"\s*\(UniProt:\w+\)", "", name).strip()
        d = ag[base]
        d["epi"] += 1
        d["org"] = (r["parent_source_antigen_source_org_name"] or "").strip()
        d["uni"] = m.group(1) if m else ""
        if r["pubmed_id"]:
            d["refs"].add(r["pubmed_id"])
        for x in (r["structure_description"] or "").split(","):
            if x.strip():
                d["res"].add(x.strip())

    out = []
    for name, d in sorted(ag.items(), key=lambda kv: -kv[1]["epi"]):
        known = any(h in name.lower() for h in HAVE)
        out.append(dict(antigen=name, uniprot=d["uni"], organism=d["org"], n_epitopes=d["epi"],
                        n_refs=len(d["refs"]), n_critical_residues=len(d["res"]),
                        already_in_benchmark="yes" if known else "", sim_category=sim_category(name)))
    OUT.parent.mkdir(exist_ok=True)
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)
    print(f"wrote {len(out)} candidate antigens -> {OUT}")


if __name__ == "__main__":
    main()
