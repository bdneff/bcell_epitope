#!/usr/bin/env python3
"""Sync the flat SIMPLE benchmark CSV from the curated v2 workbook.

The table/figure scripts read benchmark/antigen_alanine_scan_extracted_SIMPLE_v1.csv,
but that file was a partial extraction that stopped after the first complex of some
antigens (e.g. Lysozyme had only HyHEL-63/1DQJ, missing D1.3/1VFB and HyHEL-10/3HFM,
both present in the v2 workbook's Simple_Primary sheet).

This script is idempotent: it reads Simple_Primary, and APPENDS to the CSV any primary
single-residue alanine row (keyed by Antigen+PDB+Mutation) that is not already present.
Existing rows are left byte-for-byte untouched. The auxiliary columns the CSV carries but
the workbook does not (fold-change ratio, temperature, formatted source) are filled to
match the existing curation: ratio = exp(ddG/RT) at T=298.15 K (RT=0.5925 kcal/mol).

Run:    python3 analysis/sync_benchmark_csv.py
Then:   python3 analysis/make_label_tables.py && python3 analysis/make_ddg_coloring.py
Recorded in manuscript/figures/FIGURES.md.
"""
import csv, math, pathlib
import openpyxl

CSV  = pathlib.Path("benchmark/antigen_alanine_scan_extracted_SIMPLE_v1.csv")
XLSX = pathlib.Path("benchmark/antigen_alanine_scanning_benchmark_v2.xlsx")
SHEET = "Simple_Primary"
RT_298 = 0.0019872 * 298.15  # kcal/mol, matches existing fold-change column

# Only sync antigens we actually simulate (avoids pulling the whole benchmark in).
SIMULATED = {"Lysozyme", "HPr", "VEGF", "MT-SP1", "Bont/A1"}

def keyof(antigen, pdb, mut):
    return (antigen.strip(), pdb.strip(), mut.strip())

# --- existing CSV (preserve order + exact columns) ---
with open(CSV, newline="") as f:
    reader = csv.DictReader(f)
    fields = reader.fieldnames
    existing_rows = list(reader)
have = {keyof(r["Antigen"], r["PDB/model"], r["Mutation"]) for r in existing_rows}

# --- workbook ---
wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
ws = wb[SHEET]
xrows = list(ws.iter_rows(values_only=True))
H = {h: i for i, h in enumerate(xrows[0])}
def gx(r, col):
    v = r[H[col]]
    return "" if v is None else str(v).strip()

def src_fmt(source_row):
    # workbook "AB-Bind line 623" -> existing CSV style "AB-Bind row line 623; AB-Bind line 623"
    n = source_row.replace("AB-Bind line", "").strip()
    return f"AB-Bind row line {n}; AB-Bind line {n}" if n else source_row

appended = []
for r in xrows[1:]:
    ag = gx(r, "Antigen")
    if ag not in SIMULATED:
        continue
    if not gx(r, "Use for MD correlation?").startswith("YES"):
        continue
    k = keyof(ag, gx(r, "PDB"), gx(r, "Mutation"))
    if k in have:
        continue
    have.add(k)
    try:
        ddg = float(gx(r, "Best ΔΔG kcal/mol"))
        ratio = f"{math.exp(ddg / RT_298):.3f}"
    except ValueError:
        ratio = ""
    appended.append({
        "Use for MD correlation?": gx(r, "Use for MD correlation?"),
        "Dataset": gx(r, "Dataset"),
        "Antigen": ag,
        "Antibody": gx(r, "Antibody"),
        "PDB/model": gx(r, "PDB"),
        "Mutation": gx(r, "Mutation"),
        "Residue #": gx(r, "Residue #"),
        "WT": gx(r, "WT"),
        "Mutant": gx(r, "Mutant"),
        "Best ΔΔG kcal/mol": gx(r, "Best ΔΔG kcal/mol"),
        "Label type": gx(r, "Label type"),
        "Binding fold-change/ratio": ratio,
        "Assay": gx(r, "Assay"),
        "Temperature K": "298.15",
        "Source": src_fmt(gx(r, "Source row")),
        "DOI": gx(r, "DOI"),
        "Plain-English interpretation": gx(r, "Plain-English interpretation"),
        "Caution / notes": gx(r, "Notes"),
    })

if not appended:
    print("CSV already in sync with workbook; nothing appended.")
else:
    with open(CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(existing_rows)
        w.writerows(appended)
    print(f"appended {len(appended)} row(s):")
    from collections import Counter
    for (ag, pdb, ab), n in Counter(
            (a["Antigen"], a["PDB/model"], a["Antibody"]) for a in appended).items():
        print(f"  +{n:2d}  {ag} / {pdb} / {ab}")
