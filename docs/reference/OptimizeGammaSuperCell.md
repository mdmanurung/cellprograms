<div class="page-header">

# Optimize SuperCell granularity with mcRigor

<div class="d-none name">

`OptimizeGammaSuperCell.Rd`

</div>

</div>

<div class="ref-description section level2">

Runs SuperCell's SCimplify at several granularities (cells per metacell - SuperCell's native \`gamma\`) per group and calls \`mcRigor::mcRigor_OPTIMIZE\` to select the granularity maximizing the mcRigor Score (1 - w\*DubRate - (1-w)\*ZeroRate). Cheap relative to the SEACells analogue: one graph partition per group per candidate granularity (seconds).

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
OptimizeGammaSuperCell(
  seurat_obj,
  group.by = "seurat_clusters",
  gammas = c(25, 50, 100),
  reduction = "pca",
  dims = NULL,
  assay = "RNA",
  slot = "counts",
  layer = "counts",
  min_cells = 100,
  k.knn = 30,
  n_top_genes = 2000,
  n_pc = 30,
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
Metadata column(s) defining groups (as in \[MetacellsByGroupsSuperCell()\]).

gammas<a href="#arg-gammas" class="anchor" aria-label="anchor"></a><br>
Candidate granularities as cells-per-metacell values. Default \`c(25, 50, 100)\`.

reduction<a href="#arg-reduction" class="anchor" aria-label="anchor"></a><br>
Passed to \[MetacellsByGroupsSuperCell()\] internals.

dims<a href="#arg-dims" class="anchor" aria-label="anchor"></a><br>
Passed to \[MetacellsByGroupsSuperCell()\] internals.

assay<a href="#arg-assay" class="anchor" aria-label="anchor"></a><br>
Passed to \[MetacellsByGroupsSuperCell()\] internals.

slot<a href="#arg-slot" class="anchor" aria-label="anchor"></a><br>
Passed to \[MetacellsByGroupsSuperCell()\] internals.

layer<a href="#arg-layer" class="anchor" aria-label="anchor"></a><br>
Passed to \[MetacellsByGroupsSuperCell()\] internals.

min_cells<a href="#arg-min-cells" class="anchor" aria-label="anchor"></a><br>
Passed to \[MetacellsByGroupsSuperCell()\] internals.

k.knn<a href="#arg-k-knn" class="anchor" aria-label="anchor"></a><br>
Passed to \[MetacellsByGroupsSuperCell()\] internals.

n_top_genes<a href="#arg-n-top-genes" class="anchor" aria-label="anchor"></a><br>
Passed to \[MetacellsByGroupsSuperCell()\] internals.

n_pc<a href="#arg-n-pc" class="anchor" aria-label="anchor"></a><br>
Passed to \[MetacellsByGroupsSuperCell()\] internals.

mcrigor_args<a href="#arg-mcrigor-args" class="anchor" aria-label="anchor"></a><br>
Optional named list overriding mcRigor_OPTIMIZE arguments (\`Nrep\`, \`weight\`, ...).

seed<a href="#arg-seed" class="anchor" aria-label="anchor"></a><br>
Passed to \[MetacellsByGroupsSuperCell()\] internals.

verbose<a href="#arg-verbose" class="anchor" aria-label="anchor"></a><br>
Passed to \[MetacellsByGroupsSuperCell()\] internals.

</div>

<div class="section level2">

## Value<a href="#value" class="anchor" aria-label="anchor"></a>

A list with \`optimal_gamma\` (named list, per group), \`scores\` (per group x gamma Score table), and \`membership\` (named list of per-group membership vectors for the chosen granularity, one entry per cell; NULL when no group produced a result). Gammas whose membership collapses to a single metacell cannot be scored by mcRigor and are dropped with a warning (\`NA\` in \`scores\`).

</div>
