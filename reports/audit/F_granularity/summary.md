# hdWGCNA wrapper arm (W3+W4): results summary

Granularity benchmark, CD8 T cells, 99 donors (age 21-75). Three hdWGCNA 0.4.12
arms (metacell / pseudobulk / singlecell granularity) run through the same
pipeline (`hdwgna_pipeline.R`), then joined to the 9 existing EB-family/cNMF
arms in the 12-arm concordance comparison (`compare_F_hdwgna.py`).

## Pipeline corrections applied before final runs

1. **Power-rule fix (CONFIRMED DEFECT in our own D2 code, fixed)**: module
   power was selected by `which.max(SFT.R.sq)`; the hdWGCNA tutorial convention
   (basic_tutorial.Rmd ~line 320, and `ConstructNetwork` auto-selection) is the
   FIRST power with scale-free fit R^2 >= 0.80. Fixed uniformly across arms;
   the singlecell arm went from DEGENERATE (power 30 at cap, 0 modules) to
   2 modules at power 5. Old exports archived under `powmax_archive/`.
2. **Per-cell module eigengenes (design, not bug)**: `ModuleEigengenes` runs on
   the parent object and returns per-CELL MEs even when metacells are built
   (tutorial line 413 confirms this is the designed quantity). Donor scores are
   aggregated adaptively: n_cells-weighted metacell mean for direct-network
   arms, per-cell mean for the metacell-built arm.
3. **min_cells deviation (registered)**: singlecell arm min_cells 25 -> 80.
   `ConstructMetacells` forms metacells from k=25 neighborhoods under
   max_shared=10, feasible only when group n > k^2/max_shared ~= 63 cells;
   hdWGCNA 0.4.12 crashes on per-group failure (empirical: n=27 donor ->
   "Metacell failed" -> NULL dereference). 93/99 donors retained (donor sizes
   27-77 for the six smallest; median 265).
4. **ModuleConnectivity guarded**: wrapped in tryCatch with `conn_ok` flag;
   degenerate arms write module_genes + eigengenes only.

## Arm results (final, convention power rule)

| arm | genes in network | modules (non-grey) | soft power | min_cells | metacells built |
|---|---|---|---|---|---|
| hdwgna_metacell | 709 | 2 | 3 | 25 | yes |
| hdwgna_pseudobulk | 1486 | 2 | 4 | 25 | yes |
| hdwgna_singlecell | 447 | 2 | 5 | 80 (deviation) | yes (93/99 donors) |

All three arms non-degenerate. G5 (descriptive, no pass bar): SVD reconstruction
of the metacell-routed network mean |cor| 0.96 (min 0.95), mean top-50 gene
Jaccard 0.66 (`g5_hdwgna_svd.csv`, `g5_hdwgna_summary.json`).

## 12-arm concordance (mean matched |cor| on gene loadings, Hungarian match)

`comparison_table.csv` (66 pairs), `concordance_heatmap.png` (12x12; genenmf
arms absent = documented K=0 empty-file limitation).

Key numbers:
- hdwgna_metacell vs hdwgna_singlecell: **0.96** - the wrapper's two
  metacell-built granularities agree strongly.
- cebmf_metacell vs cebmf_singlecell: **0.82** - cEBMF route robustness.
- cebmf_singlecell vs hdwgna_singlecell: **0.76**; cebmf_metacell vs
  hdwgna_metacell: **0.72** - cross-method concordance is substantial at
  matched granularity.
- cnmf arms: ~0.00 against everything, including other cnmf arms - consistent
  with the S5 synthetic anchor finding (cNMF k-selection underselects and its
  programs do not align with the EB-family structure on this data).
- ebmf_pseudobulk vs ebmf_singlecell: 0.19 - plain ebmf is much less
  route-robust than cebmf (0.82).

## Interpretation

The hdWGCNA wrapper produces stable, non-degenerate modules at all three
granularities and its programs correlate substantially with the cEBMF programs
(~0.7-0.76 at matched granularity), supporting the plan's premise that the
co-expression module view and the EB factor view capture overlapping biology.
The wrapper's own two metacell-built arms agree at 0.96, i.e. the residual
method differences are not driven by the metacell builder. hdWGCNA's K=2
modules are coarser than the K=10 cEBMF programs, so per-factor concordance is
bounded; the comparison is descriptive (G5, no pass bar) as registered.

## Files

- `hdwgna_pipeline.R`, `hdwgna_wrapper.R` - pipeline + wrapper
- `F_hdwgna_{arm}_loadings.csv`, `F_hdwgna_{arm}_donor_scores.csv`,
  `F_hdwgna_{arm}_scores.csv`, `eigengenes_{arm}.csv`,
  `metacells_per_donor_{arm}.csv` per arm
- `run_stats_{arm}.json` - per-arm stats incl. `degenerate` flag, powers, timing
- `comparison_table.csv`, `concordance_heatmap.png` - 12-arm comparison
- `g5_hdwgna_svd.csv`, `g5_hdwgna_summary.json` - G5 descriptive numbers
- `powmax_archive/` - pre-fix exports (which.max power rule)
