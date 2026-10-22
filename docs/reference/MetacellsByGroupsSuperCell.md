<div class="page-header">

# Construct metacells by groups with SuperCell and mcRigor QC

<div class="d-none name">

`MetacellsByGroupsSuperCell.Rd`

</div>

</div>

<div class="ref-description section level2">

Constructs transcriptomically homogeneous metacells from a Seurat object using SuperCell 2.0 (Bilous et al. 2022) as the aggregation backend, computed separately within each level of the supplied groupings (typically donor, cell type, or condition), following the hdWGCNA \`MetacellsByGroups\` interface (Morabito et al. 2023). This is the package's \*\*default\*\*, pure-R metacell path: a preregistered benchmark (\`benchmarks/metacell_comparison/\`) found SuperCell + mcRigor comparable to SEACells on dubious-detection ARI and program-recovery fidelity at ~200x lower construction cost. Metacell quality is governed by mcRigor (Liu & Li 2025): per-metacell mcDiv dubious detection with double-permutation nulls, DubRate/ZeroRate scoring, and optional granularity optimization via \[OptimizeGammaSuperCell()\]. The SEACells backend remains available in \[MetacellsByGroupsSEACells()\].

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
MetacellsByGroupsSuperCell(
  seurat_obj,
  group.by = "seurat_clusters",
  ident.group = NULL,
  reduction = "pca",
  dims = NULL,
  assay = "RNA",
  slot = "counts",
  layer = "counts",
  mode = c("sum", "average"),
  cells_per_metacell = 75,
  target_metacells = NULL,
  min_cells = 100,
  k.knn = 30,
  n_top_genes = 2000,
  n_pc = 30,
  run_mcrigor = TRUE,
  mcrigor_args = list(),
  seed = 20261016,
  name = "SC",
  wgcna_name = "SuperCell",
  verbose = TRUE
)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

seurat_obj<a href="#arg-seurat-obj" class="anchor" aria-label="anchor"></a><br>
A Seurat object with raw counts and a dimensionality reduction (default \`pca\`). Set \`reduction = NULL\` to compute a standard Seurat PCA per group on the fly.

group.by<a href="#arg-group-by" class="anchor" aria-label="anchor"></a><br>
Character vector of metadata column names defining the stratification (e.g. \`c("donor", "celltype")\`). Metacells are constructed separately within each combination of levels.

ident.group<a href="#arg-ident-group" class="anchor" aria-label="anchor"></a><br>
Optional metadata column name used as the metacell-name prefix; defaults to the interaction of \`group.by\`.

reduction<a href="#arg-reduction" class="anchor" aria-label="anchor"></a><br>
Name of the dimensionality reduction SCimplify builds the cell-to-cell graph on. Default \`"pca"\`; \`NULL\` computes a PCA per group.

dims<a href="#arg-dims" class="anchor" aria-label="anchor"></a><br>
Dimensions of \`reduction\` to use (indices); default all computed dimensions of the reduction.

assay<a href="#arg-assay" class="anchor" aria-label="anchor"></a><br>
Assay to aggregate. Default \`RNA\`.

slot<a href="#arg-slot" class="anchor" aria-label="anchor"></a><br>
Slot (Seurat v4) to extract counts from. Default \`counts\`.

layer<a href="#arg-layer" class="anchor" aria-label="anchor"></a><br>
Layer (Seurat v5) to extract counts from. Default \`counts\`.

mode<a href="#arg-mode" class="anchor" aria-label="anchor"></a><br>
How to aggregate constituent cells into metacells: \`"sum"\` (raw-count sums, the SEACells default and the benchmark configuration) or \`"average"\` (divide by metacell size). Default \`"sum"\`.

cells_per_metacell<a href="#arg-cells-per-metacell" class="anchor" aria-label="anchor"></a><br>
Target number of single cells per metacell - SuperCell's native \`gamma\`. Default 75 (SEACells guidance; hdWGCNA k guidance: 20-75).

target_metacells<a href="#arg-target-metacells" class="anchor" aria-label="anchor"></a><br>
Optional cap on the number of metacells per group. When set, gamma is bisected (log scale) to land within +/-10 target; otherwise \`gamma = cells_per_metacell\` is used directly.

min_cells<a href="#arg-min-cells" class="anchor" aria-label="anchor"></a><br>
Groups with fewer than this many cells are skipped with a warning. Default 100.

k.knn<a href="#arg-k-knn" class="anchor" aria-label="anchor"></a><br>
Number of neighbors for the cell-to-cell KNN graph. Default 30 (SCimplify default).

n_top_genes<a href="#arg-n-top-genes" class="anchor" aria-label="anchor"></a><br>
Number of highly variable genes for the on-the-fly PCA when \`reduction = NULL\`. Default 2000.

n_pc<a href="#arg-n-pc" class="anchor" aria-label="anchor"></a><br>
Number of PCs for the on-the-fly PCA when \`reduction = NULL\`. Default 30.

run_mcrigor<a href="#arg-run-mcrigor" class="anchor" aria-label="anchor"></a><br>
If TRUE (default), run mcRigor dubious-metacell detection per group and attach mcDiv/DubRate/ZeroRate/Score to the QC table.

mcrigor_args<a href="#arg-mcrigor-args" class="anchor" aria-label="anchor"></a><br>
Optional named list overriding mcRigor_DETECT arguments: \`Nrep\` (50), \`gene_filter\` (0.1), \`feature_use\` (2000), \`test_cutoff\` (0.01), \`cor_method\` ("pearson"), \`aggregate_method\` (defaults to \`"sum"\`/\`"mean"\` following \`mode\`).

seed<a href="#arg-seed" class="anchor" aria-label="anchor"></a><br>
Random seed (set before each group's SCimplify call; the graph construction is deterministic but the parameter is kept for API parity and on-the-fly PCA reproducibility). Default 20261016.

name<a href="#arg-name" class="anchor" aria-label="anchor"></a><br>
Reserved for interface parity with \[MetacellsByGroupsSEACells()\] (metacell names are \`\<group\>\_\<i\>\`). Default \`"SC"\`.

wgcna_name<a href="#arg-wgcna-name" class="anchor" aria-label="anchor"></a><br>
Name under which results are stored in \`seurat_obj@misc\`. Default \`"SuperCell"\`.

verbose<a href="#arg-verbose" class="anchor" aria-label="anchor"></a><br>
Print progress. Default TRUE.

</div>

<div class="section level2">

## Value<a href="#value" class="anchor" aria-label="anchor"></a>

The Seurat object with results stored in \`seurat_obj@misc\[\[wgcna_name\]\]\`: \`metacell_object\` (Seurat object of metacells, counts = aggregated expression, metadata = Group, n_cells, the individual \`group.by\` columns, and mcRigor QC columns), \`soft_assignments\` (always NULL for this backend: SCimplify returns hard memberships only), \`qc\` (per-metacell QC data.frame), \`fit_stats\` (per-group gamma and metacell count), \`params\` (all call parameters), \`skipped_groups\` (data.frame), and \`call\`. Metacells are named \`\<group\>\_\<i\>\` and column-bound across groups.

</div>

<div class="section level2">

## Details<a href="#details" class="anchor" aria-label="anchor"></a>

Groups are processed independently; each group's cells are subsetted (with their reduction), SCimplify partitions the KNN graph, metacell counts are summarized from raw counts, and mcRigor QC is computed against the group's single-cell counts. The result is hdWGCNA-interoperable in the same way as the SEACells backend.

</div>

<div class="section level2">

## References<a href="#references" class="anchor" aria-label="anchor"></a>

Bilous et al. (2022) Bioinformatics, doi:10.1093/bioinformatics/btac104; Liu & Li (2025) Nat Commun, doi:10.1038/s41467-025-63626-5; Morabito et al. (2023) Cell Rep Methods, doi:10.1016/j.crmeth.2023.100498.

</div>
