# REPORT: SEACells vs SuperCell 2.0 metacell benchmark

Verdict (PREREG rule): **REPLACE (SuperCell as default R-native path)**

- Mean |ARI gap| (dubious detection): 0.0190 (PREREG tolerance 0.05)
- Mean relative fidelity gap: 5.4% (PREREG tolerance 10%)

## Tables
- metrics_all.csv / metrics_summary.csv (mcRigor ARI/FF/purity, ASW, compactness, runtime)
- recovery_all.csv / recovery_summary.csv (program detection + fidelity)
- permc_*.csv (per-metacell purity/size), timing_*.csv (membership wall-clock)

## Figures
- fig_ari.png, fig_ff.png, fig_asw.png, fig_detection.png, fig_fidelity.png (line + ribbon)
- fig_purity_ecdf.png (purity distributions, mirrors mcRigor Fig 1b)
- fig_size_purity.png (size-vs-purity confound control)
- fig_granularity_control.png (achieved vs target N)
- fig_paired_ari.png, fig_paired_fidelity.png (per-seed paired comparison)
- fig_detection_heatmap.png (per-cell-type program detection)
- fig_runtime.png (membership wall-clock, log scale)
- fig_downstream_cost_time.png, fig_downstream_cost_mem.png (program-fit time and peak memory per target N)

## Limitations
- Synthetic fixtures only (planted types/donors/programs, 5% injected doublets).
- Geometry in PCA space, not diffusion maps (paper's choice) - pure-R stand-in.
- SEACells run pooled (single group) to match SCimplify's pooled partition.
- SuperCell pinned to supercell-2.0 @ 89b34c078dba0b91289a5a23cf75a6db3a46d232.
- Detection = max Spearman(loading, binary signature); the binary target bounds
  attainable correlation below 1, so absolute values are conservative.
- Program-fit memory is a gc() high-water proxy (Max Used), not RSS.

## Deviations from PREREG
- SEACells 15k reduced to 2 seeds (runtime); SuperCell 15k kept 3 seeds.
- SEACells wrapper fix: cells_per_metacell = floor(n_cells/N) must be passed
  alongside target_metacells=N (default 75 capped metacell count at N>66).
- mcRigor calling convention pinned against HEAD c70cd6f: membership labels
  are MC-prefixed strings named by cell id (as.character() drops vector names
  in R 4.4); DETECT receives a raw Seurat object; fields=character(0).
- max_factors in recovery is adaptive: min(6, min donors per type - 1);
  coarse N yields thin dominant-donor designs (part of the granularity signal).

