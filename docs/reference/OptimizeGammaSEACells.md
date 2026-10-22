# Optimize SEACells granularity with mcRigor

Runs SEACells at several granularities (cells per metacell) per group
and calls \`mcRigor::mcRigor_OPTIMIZE\` to select the granularity
maximizing the mcRigor Score (1 - w\*DubRate - (1-w)\*ZeroRate).
Expensive: one SEACells fit per group per candidate granularity.

## Usage

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

## Arguments

- seurat_obj:

  A Seurat object.

- group.by:

  Metadata column(s) defining groups (as in
  \[MetacellsByGroupsSEACells()\]).

- gammas:

  Candidate granularities as cells-per-metacell values. Default \`c(25,
  50, 100)\` (mcRigor semi-synthetic optimum was 50 for SEACells).

- reduction:

  Passed to \[MetacellsByGroupsSEACells()\] internals.

- dims:

  Passed to \[MetacellsByGroupsSEACells()\] internals.

- assay:

  Passed to \[MetacellsByGroupsSEACells()\] internals.

- slot:

  Passed to \[MetacellsByGroupsSEACells()\] internals.

- layer:

  Passed to \[MetacellsByGroupsSEACells()\] internals.

- n_top_genes:

  Passed to \[MetacellsByGroupsSEACells()\] internals.

- n_pc:

  Passed to \[MetacellsByGroupsSEACells()\] internals.

- n_waypoint_eigs:

  Passed to \[MetacellsByGroupsSEACells()\] internals.

- waypt_proportion:

  Passed to \[MetacellsByGroupsSEACells()\] internals.

- n_neighbors:

  Passed to \[MetacellsByGroupsSEACells()\] internals.

- convergence_epsilon:

  Passed to \[MetacellsByGroupsSEACells()\] internals.

- max_iter:

  Passed to \[MetacellsByGroupsSEACells()\] internals.

- use_sparse:

  Passed to \[MetacellsByGroupsSEACells()\] internals.

- min_cells:

  Passed to \[MetacellsByGroupsSEACells()\] internals.

- mcrigor_args:

  Optional named list overriding mcRigor_OPTIMIZE arguments (\`D_bw\`,
  \`optim_method\`, \`weight\`, \`Nrep\`, ...).

- seed:

  Passed to \[MetacellsByGroupsSEACells()\] internals.

- verbose:

  Passed to \[MetacellsByGroupsSEACells()\] internals.

## Value

A list with \`optimal_gamma\` (named list, per group), \`scores\` (per
group x gamma Score table), and \`membership\` (named list of per-group
membership vectors for the chosen granularity, one entry per cell; NULL
when no group produced a result). Gammas whose membership collapses to a
single metacell cannot be scored by mcRigor and are dropped with a
warning (\`NA\` in \`scores\`).
