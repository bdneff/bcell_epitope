Antigen-side alanine scanning benchmark v2

What changed in v2:
- Added missing AB-Bind antibody-antigen PDB groups that were not in v1:
  * 1VFB: Lysozyme / D1.3 Fv, 12 antigen-side single-alanine ddG rows.
  * 3HFM: Lysozyme / HyHEL-10 Fv, 4 antigen-side single-alanine ddG rows.
- Added a secondary hGH Table 2.4 aggregate profile:
  * 75 rows manually transcribed from the last column of Lei Jin thesis Table 2.4 screenshots.
  * This hGH table is NOT antibody-specific. It is the total ΔΔG across 21 anti-hGH MAbs and is included as an antigen-level score only.

Main files:
- antigen_alanine_scan_extracted_SIMPLE_v2.csv: readable table for quick use.
- antigen_alanine_scan_extracted_FULL_v2.csv: full schema table.
- antigen_alanine_scanning_benchmark_v2.xlsx: workbook with summary, readable rows, full rows, hGH aggregate profile, and inventory.
- hgh_table2_4_total_profile_MANUAL_v2.csv: secondary hGH aggregate profile.
- antigen_alanine_scan_paper_inventory_v2.csv: source inventory.

How to use:
- For MD correlation using single-residue antibody-specific alanine labels, use:
  antigen_alanine_scan_extracted_SIMPLE_v2.csv
  and filter:
  Use for MD correlation? = YES - primary single-residue alanine label
- Main label column:
  Best ΔΔG kcal/mol
- Positive ΔΔG = mutation weakens binding.
- Negative ΔΔG = mutation improves binding.

Current primary dataset:
- Total extracted rows: 203
- Primary single-residue labels: 199
- Dataset groups: 13
- PDB groups: 12

Important caveat:
The hGH Table 2.4 aggregate profile is manually transcribed from scanned PDF screenshots and should be used as secondary/provisional data only until QC'd against the original table image.
