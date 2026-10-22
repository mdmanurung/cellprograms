<div class="page-header">

# Optimize SEACells granularity with mcRigor

<div class="d-none name">

`OptimizeGammaSEACells.Rd`

</div>

</div>

<div class="ref-description section level2">

Runs SEACells at several granularities (cells per metacell) per group and calls \`mcRigor::mcRigor_OPTIMIZE\` to select the granularity maximizing the mcRigor Score (1 - w\*DubRate - (1-w)\*ZeroRate). Expensive: one SEACells fit per group per candidate granularity.

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
OptimizeGammaSEACells(
  seurat_obj,
  group.by = "seurat_clusters",
  gammas = c(25, 50, 100),
  reduction = "pca",
  dims = NULL,
  assay = "RNA",
  slot = "counts",
  layer = "counts",
  n_top_genes = 1500,
  n_pc = 50,
  n_waypoint_eigs = 10,
  waypt_proportion = 0.9,
  n_neighbors = 15,
  convergence_epsilon = 1e-05,
  max_iter = 20,
  use_sparse = TRUE,
  min_cells = 100,
  mcrigor_args = list(),
  seed = 20261016,
  verbose = TRUE
)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

seurat_obj<a href="#arg-seurat-obj" class="anchor" aria-label="anchor"></a><br>
A Seurat object.

group.by<a href="#arg-group-by" class="anchor" aria-label="anchor"></a><br>
Metadata column(s) defining groups (as in \[MetacellsByGroupsSEACells()\]).

gammas<a href="#arg-gammas" class="anchor" aria-label="anchor"></a><br>
Candidate granularities as cells-per-metacell values. Default \`c(25, 50, 100)\` (mcRigor semi-synthetic optimum was 50 for SEACells).

reduction<a href="#arg-reduction" class="anchor" aria-label="anchor"></a><br>
Passed to \[MetacellsByGroupsSEACells()\] internals.

dims<a href="#arg-dims" class="anchor" aria-label="anchor"></a><br>
Passed to \[MetacellsByGroupsSEACells()\] internals.

assay<a href="#arg-assay" class="anchor" aria-label="anchor"></a><br>
Passed to \[MetacellsByGroupsSEACells()\] internals.

slot<a href="#arg-slot" class="anchor" aria-label="anchor"></a><br>
Passed to \[MetacellsByGroupsSEACells()\] internals.

layer<a href="#arg-layer" class="anchor" aria-label="anchor"></a><br>
Passed to \[MetacellsByGroupsSEACells()\] internals.

n_top_genes<a href="#arg-n-top-genes" class="anchor" aria-label="anchor"></a><br>
Passed to \[MetacellsByGroupsSEACells()\] internals.

n_pc<a href="#arg-n-pc" class="anchor" aria-label="anchor"></a><br>
Passed to \[MetacellsByGroupsSEACells()\] internals.

n_waypoint_eigs<a href="#arg-n-waypoint-eigs" class="anchor" aria-label="anchor"></a><br>
Passed to \[MetacellsByGroupsSEACells()\] internals.

waypt_proportion<a href="#arg-waypt-proportion" class="anchor" aria-label="anchor"></a><br>
Passed to \[MetacellsByGroupsSEACells()\] internals.

n_neighbors<a href="#arg-n-neighbors" class="anchor" aria-label="anchor"></a><br>
Passed to \[MetacellsByGroupsSEACells()\] internals.

convergence_epsilon<a href="#arg-convergence-epsilon" class="anchor" aria-label="anchor"></a><br>
Passed to \[MetacellsByGroupsSEACells()\] internals.

max_iter<a href="#arg-max-iter" class="anchor" aria-label="anchor"></a><br>
Passed to \[MetacellsByGroupsSEACells()\] internals.

use_sparse<a href="#arg-use-sparse" class="anchor" aria-label="anchor"></a><br>
Passed to \[MetacellsByGroupsSEACells()\] internals.

min_cells<a href="#arg-min-cells" class="anchor" aria-label="anchor"></a><br>
Passed to \[MetacellsByGroupsSEACells()\] internals.

mcrigor_args<a href="#arg-mcrigor-args" class="anchor" aria-label="anchor"></a><br>
Optional named list overriding mcRigor_OPTIMIZE arguments (\`D_bw\`, \`optim_method\`, \`weight\`, \`Nrep\`, ...).

seed<a href="#arg-seed" class="anchor" aria-label="anchor"></a><br>
Passed to \[MetacellsByGroupsSEACells()\] internals.

verbose<a href="#arg-verbose" class="anchor" aria-label="anchor"></a><br>
Passed to \[MetacellsByGroupsSEACells()\] internals.

</div>

<div class="section level2">

## Value<a href="#value" class="anchor" aria-label="anchor"></a>

A list with \`optimal_gamma\` (named list, per group), \`scores\` (per group x gamma Score table), and \`membership\` (named list of per-group membership vectors for the chosen granularity, one entry per cell; NULL when no group produced a result). Gammas whose membership collapses to a single metacell cannot be scored by mcRigor and are dropped with a warning (\`NA\` in \`scores\`).

</div>
