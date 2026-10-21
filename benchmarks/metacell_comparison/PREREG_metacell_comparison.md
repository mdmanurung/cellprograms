# PREREG: SEACells vs SuperCell 2.0 metacell benchmark

Written before any benchmark run. Pins and analysis choices below are fixed before
results are seen; deviations discovered mid-run get documented in REPORT.md, not
silently amended here.

## Question

Can SuperCell 2.0 (pure R) replace SEACells (Python/basilisk) as the metacell backend
for cellprograms, judged by (1) mcRigor metacell QC, (2) metacell geometry, and
(3) downstream planted-program recovery through `fit_celltype_programs`?

## Pins

- SuperCell: `GfellerLab/SuperCell@supercell-2.0`, branch HEAD `89b34c078dba0b91289a5a23cf75a6db3a46d232`
  at PREREG time (install-time commit recorded in results).
- SEACells 0.3.3 via the package's basilisk env (existing `MetacellsByGroupsSEACells` path).
- mcRigor: JSB-UCLA/mcRigor (already a package dependency; version at run time recorded).
- R 4.4.2, Seurat 5.6.0 (worker-1).

## Fixture (synthetic; no real data in this benchmark)

Per seed: 5 planted cell types with donor structure (10 donors), ~3,000 genes,
planted donor-level programs (one per cell type, shared donor activity driving a
gene set), injected heterotypic doublets at a 5% rate with known cell-of-origin
pairs (ground-truth dubious labels). Scales: 5,000 and 15,000 cells. Seeds: 1, 2, 3.
Fixture code mirrors the mcRigor_reproduce synthetic design (planted clusters +
injected doublets + continuous donor structure) using the package's simulation
machinery where applicable.

## Grid

Granularity-matched at target metacell counts N in {50, 100, 200, 400}:

- SuperCell: `SCimplify` gamma searched per (scale, seed, N) to land within +/-10% of N.
- SEACells: `MetacellsByGroupsSEACells` with n_metacells = N directly.

24 runs per method (2 scales x 4 N x 3 seeds). Runs that fail do so loudly and are
recorded as failures in the results table (a method that cannot cover the grid loses
on coverage).

## Metrics (fixed in advance)

1. mcRigor QC: `mcRigor_DETECT` on (single-cell object, cell_membership), then
   `mcRigor_threshold(pur_metric = 'metacell_purity', test_cutoff = 0.05)`.
   Dubious-cell detection vs ground truth: ARI (aricode), recall, precision,
   FF = 2/(1/recall + 1/precision). Per-metacell purity distributions.
2. Geometry (PCA space, pure R; deviation from the paper's diffusion maps is a
   documented limitation): ASW = `cluster::silhouette(membership, dist(embeddings))`
   per-metacell average width; compactness = mean over PCs of within-metacell
   variance of cell embeddings (paper definitions,
   SuperCellMultiomicsAnalyses R/functions/metacellCompactnessSeparation.R).
3. Program recovery (decision-critical): metacell counts -> `fit_celltype_programs`
   (same settings across methods) -> (a) detection: max correlation of each recovered
   program's loadings with each planted program's signature genes; (b) fidelity:
   correlation of metacell-level program scores with the planted donor-level scores.

## Primary comparison

At each (scale, N): mean +/- sd across seeds of every metric, SuperCell vs SEACells.
Verdict rule (fixed now): SuperCell is recommended as the default R-native path if it
matches SEACells within 0.05 ARI on dubious detection and within 10% relative on
program-recovery fidelity, averaged over the grid, with no coverage failures. Ties
favor SuperCell (dependency reduction). Otherwise: demote SEACells to optional or
keep, per which tier drives the gap.

## Scope

Report + recommendation only (benchmarks/metacell_comparison/REPORT.md). No package
code changes in this arc; a SuperCell backend implementation is a separate follow-up.
