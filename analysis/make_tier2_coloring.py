#!/usr/bin/env python3
"""Build Tier-2 (binary + consensus-frequency) critical-residue coloring for IEDB shotgun antigens.

For each configured antigen, pulls its flow-cytometry (shotgun-mutagenesis) conformational
epitopes from IEDB, parses the curated critical residues (WT+position), and WT-VERIFIES each
against the chosen simulation structure under the allowed (chain, numbering-offset) candidates
-- discarding any residue that does not match (other serotype/protein/numbering). Surviving
residues are accumulated into a consensus: how many antibodies call each (chain,resid) critical.

Writes, per antigen <key>:
  manuscript/figures/ddg/<key>_freq.dat     "chain resid n_antibodies"  (consensus frequency)
  manuscript/figures/ddg/<key>_binary.dat   "chain resid 1"             (binary critical)
  manuscript/figures/<key>_freqbar.png       white->red colorbar (0..max antibodies)
  benchmark/tier2_labels_<key>.csv           provenance (chain,resid,WT,n_antibodies)
  manuscript/figures/ddg/tier2_manifest.tsv  key panel struct datafile vmin vmax

Render with analysis/render_ddg.tcl (now chain-aware) via render_ddg_figures.sh.
Run: python3 analysis/make_tier2_coloring.py
"""
import urllib.request, urllib.parse, json, collections, re, csv, pathlib
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors

API = "https://query-api.iedb.org/"
DISC = "structure_type=eq." + urllib.parse.quote("Discontinuous peptide")
DDGDIR = pathlib.Path("manuscript/figures/ddg")
FIGDIR = pathlib.Path("manuscript/figures")
BENCH = pathlib.Path("benchmark")
T2O = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
       'TYR':'Y','VAL':'V'}

# key -> (display, uniprot, structure, [(chain, offset)] verification candidates: resid = pos+offset)
ANTIGENS = {
    "dengueE": ("Dengue E (DENV2)", "P17763", "structures/DengueE_1OKE.pdb", [("A", 0), ("A", -280)]),
    "ha":      ("Influenza HA (H3)", "P03452", "structures/HA_1HGG.pdb",      [("A", -16), ("B", -345)]),
}


def get(path):
    req = urllib.request.Request(API + path, headers={"Accept": "application/json"})
    return json.load(urllib.request.urlopen(req, timeout=120))


def struct_seq(path):
    ch = collections.defaultdict(dict)
    for line in open(path):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            ch[line[21]][int(line[22:26])] = T2O.get(line[17:20].strip(), 'X')
    return ch


def pull(uni):
    rows, off = [], 0
    filt = "&parent_source_antigen_name=ilike." + urllib.parse.quote(f"*{uni}*")
    while True:
        p = get(f"bcell_search?select=structure_description,assay_names{filt}&{DISC}&order=bcell_id&limit=1000&offset={off}")
        rows += p
        if len(p) < 1000:
            break
        off += 1000
    return [r for r in rows if "flow cytometry" in (r["assay_names"] or "").lower()]


def main():
    DDGDIR.mkdir(parents=True, exist_ok=True)
    manifest = ["\t".join(["key", "panel", "struct", "datafile", "vmin", "vmax"])]
    for key, (disp, uni, struct, cands) in ANTIGENS.items():
        seq = struct_seq(struct)
        recs = pull(uni)
        freq = collections.Counter()      # (chain, resid) -> n_antibodies
        wtmap = {}
        for r in recs:
            seen = set()
            for tok in (r["structure_description"] or "").split(","):
                m = re.match(r"^([A-Z])(\d+)$", tok.strip())
                if not m:
                    continue
                aa, pos = m.group(1), int(m.group(2))
                for ch, offv in cands:
                    if seq.get(ch, {}).get(pos + offv) == aa:
                        k = (ch, pos + offv)
                        if k not in seen:
                            freq[k] += 1
                            seen.add(k)
                            wtmap[k] = aa
                        break
        nmax = max(freq.values())
        # data files (chain-aware)
        with open(DDGDIR / f"{key}_freq.dat", "w") as f:
            for (ch, r), n in sorted(freq.items()):
                f.write(f"{ch} {r} {n}\n")
        with open(DDGDIR / f"{key}_binary.dat", "w") as f:
            for (ch, r) in sorted(freq):
                f.write(f"{ch} {r} 1\n")
        # provenance CSV
        with open(BENCH / f"tier2_labels_{key}.csv", "w", newline="") as f:
            w = csv.writer(f); w.writerow(["chain", "resid", "WT", "n_antibodies", "critical"])
            for (ch, r), n in sorted(freq.items()):
                w.writerow([ch, r, wtmap[(ch, r)], n, 1])
        # frequency uses symmetric scale so the white->red half of BWR shows 0..nmax
        manifest.append("\t".join([key, "freq", struct, str(DDGDIR / f"{key}_freq.dat"), f"{-nmax}", f"{nmax}"]))
        manifest.append("\t".join([key, "binary", struct, str(DDGDIR / f"{key}_binary.dat"), "-1", "1"]))
        # per-antigen frequency colorbar (white->red, 0..nmax)
        fig, ax = plt.subplots(figsize=(0.7, 3.0))
        cb = matplotlib.colorbar.ColorbarBase(ax, cmap=matplotlib.colormaps["Reds"],
                norm=colors.Normalize(vmin=0, vmax=nmax), orientation="vertical")
        cb.set_label("antibodies (critical)", fontsize=16)
        cb.ax.tick_params(labelsize=14)
        fig.savefig(FIGDIR / f"{key}_freqbar.png", dpi=200, bbox_inches="tight"); plt.close(fig)
        print(f"{key}: {len(recs)} antibody records -> {len(freq)} WT-verified critical residues "
              f"(max {nmax} antibodies); wrote freq+binary data, colorbar, provenance CSV")
    (DDGDIR / "tier2_manifest.tsv").write_text("\n".join(manifest) + "\n")
    print(f"wrote {DDGDIR/'tier2_manifest.tsv'}")


if __name__ == "__main__":
    main()
