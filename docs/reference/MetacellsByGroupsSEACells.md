<div class="page-header">

# Construct metacells by groups with SEACells and mcRigor QC

<div class="d-none name">

`MetacellsByGroupsSEACells.Rd`

</div>

</div>

<div class="ref-description section level2">

Constructs transcriptomically homogeneous metacells from a Seurat object using SEACells (Persad et al. 2023) as the aggregation backend, computed separately within each level of the supplied groupings (typically donor, cell type, or condition), following the hdWGCNA \`MetacellsByGroups\` interface (Morabito et al. 2023). Metacell quality is governed by mcRigor (Liu & Li 2025): per-metacell mcDiv dubious detection with double-permutation nulls, DubRate/ZeroRate scoring, and optional granularity optimization via \[OptimizeGammaSEACells()\].

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
MetacellsByGroupsSEACells(
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
  n_top_genes = 1500,
  n_pc = 50,
  n_waypoint_eigs = 10,
  waypt_proportion = 0.9,
  n_neighbors = 15,
  convergence_epsilon = 1e-05,
  max_iter = 20,
  use_sparse = TRUE,
  use_soft_assignments = TRUE,
  run_mcrigor = TRUE,
  mcrigor_args = list(),
  seed = 20261016,
  name = "SEAC",
  wgcna_name = "SEACells",
  verbose = TRUE
)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

seurat_obj<a href="#arg-seurat-obj" class="anchor" aria-label="anchor"></a><br>
A Seurat object with raw counts and a dimensionality reduction (default \`pca\`).

group.by<a href="#arg-group-by" class="anchor" aria-label="anchor"></a><br>
Character vector of metadata column names defining the stratification (e.g. \`c("donor", "celltype")\`). Metacells are constructed separately within each combination of levels.

ident.group<a href="#arg-ident-group" class="anchor" aria-label="anchor"></a><br>
Optional metadata column name used as the metacell-name prefix; defaults to the interaction of \`group.by\`.

reduction<a href="#arg-reduction" class="anchor" aria-label="anchor"></a><br>
Name of the dimensionality reduction to build the SEACells kernel on. The stored embeddings are exported and used directly (\`build_kernel_on = 'X_pca'\`). Set to \`NULL\` to let the Python side compute a scanpy PCA (50 PCs on HVGs) instead.

dims<a href="#arg-dims" class="anchor" aria-label="anchor"></a><br>
Dimensions of \`reduction\` to use (indices or names); default all.

assay<a href="#arg-assay" class="anchor" aria-label="anchor"></a><br>
Assay to aggregate. Default \`RNA\`.

slot<a href="#arg-slot" class="anchor" aria-label="anchor"></a><br>
Slot (Seurat v4) to extract counts from. Default \`counts\`.

layer<a href="#arg-layer" class="anchor" aria-label="anchor"></a><br>
Layer (Seurat v5) to extract counts from. Default \`counts\`.

mode<a href="#arg-mode" class="anchor" aria-label="anchor"></a><br>
How to aggregate constituent cells into metacells: \`"sum"\` (SEACells default, raw-count sums) or \`"average"\` (divide by metacell size). Default \`"sum"\`.

cells_per_metacell<a href="#arg-cells-per-metacell" class="anchor" aria-label="anchor"></a><br>
Target number of single cells per metacell (SEACells guidance: 75; hdWGCNA k guidance: 20-75). Default 75.

target_metacells<a href="#arg-target-metacells" class="anchor" aria-label="anchor"></a><br>
Optional cap on the number of metacells per group.

min_cells<a href="#arg-min-cells" class="anchor" aria-label="anchor"></a><br>
Groups with fewer than this many cells are skipped with a warning. Default 100.

n_top_genes<a href="#arg-n-top-genes" class="anchor" aria-label="anchor"></a><br>
Number of highly variable genes used for the scanpy preprocessing on the Python side. Default 1500.

n_pc<a href="#arg-n-pc" class="anchor" aria-label="anchor"></a><br>
Number of PCs computed by scanpy when \`reduction = NULL\`. Default 50.

n_waypoint_eigs<a href="#arg-n-waypoint-eigs" class="anchor" aria-label="anchor"></a><br>
Number of eigenvectors for waypoint initialization. Default 10.

waypt_proportion<a href="#arg-waypt-proportion" class="anchor" aria-label="anchor"></a><br>
Proportion of archetypes initialized by waypoint analysis. Default 0.9.

n_neighbors<a href="#arg-n-neighbors" class="anchor" aria-label="anchor"></a><br>
Number of neighbors for the SEACells graph. Default 15.

convergence_epsilon<a href="#arg-convergence-epsilon" class="anchor" aria-label="anchor"></a><br>
Franke-Wolfe convergence threshold. Default 1e-5 (hdWGCNA vignette value).

max_iter<a href="#arg-max-iter" class="anchor" aria-label="anchor"></a><br>
Maximum SEACells fit iterations. Default 20; non-convergence is recorded in fit_stats, not fatal.

use_sparse<a href="#arg-use-sparse" class="anchor" aria-label="anchor"></a><br>
Reserved for API compatibility; currently \*\*ignored\*\* — the pinned SEACells 0.3.3 fork has no \`use_sparse\` kwarg (a warning is issued when set to a non-default value). Default TRUE.

use_soft_assignments<a href="#arg-use-soft-assignments" class="anchor" aria-label="anchor"></a><br>
If TRUE (default), export the cells x metacells soft-assignment matrix (\`A\_.T\`) and store it on the output.

run_mcrigor<a href="#arg-run-mcrigor" class="anchor" aria-label="anchor"></a><br>
If TRUE (default), run mcRigor dubious-metacell detection per group and attach DubRate/ZeroRate/Score to the QC table.

mcrigor_args<a href="#arg-mcrigor-args" class="anchor" aria-label="anchor"></a><br>
Optional named list overriding mcRigor_DETECT arguments: \`Nrep\` (50), \`gene_filter\` (0.1), \`feature_use\` (2000), \`test_cutoff\` (0.01), \`cor_method\` ("pearson"), \`aggregate_method\` (defaults to \`"sum"\`/\`"mean"\` following \`mode\`).

seed<a href="#arg-seed" class="anchor" aria-label="anchor"></a><br>
Random seed for SEACells initialization. Default 20261016.

name<a href="#arg-name" class="anchor" aria-label="anchor"></a><br>
String appended to metacell names. Default \`"SEAC"\`.

wgcna_name<a href="#arg-wgcna-name" class="anchor" aria-label="anchor"></a><br>
Name under which results are stored in \`seurat_obj@misc\`. Default \`"SEACells"\`.

verbose<a href="#arg-verbose" class="anchor" aria-label="anchor"></a><br>
Print progress. Default TRUE.

</div>

<div class="section level2">

## Value<a href="#value" class="anchor" aria-label="anchor"></a>

The Seurat object with results stored in \`seurat_obj@misc\[\[wgcna_name\]\]\`: \`metacell_object\` (Seurat object of metacells, counts = aggregated expression, metadata = Group, n_cells, the individual \`group.by\` columns, and mcRigor QC columns), \`soft_assignments\` (dgCMatrix, cells x metacells, if \`use_soft_assignments\`; dropped with a warning for groups whose fit collapses to fewer metacells than archetypes fitted, since archetype columns are not metacell-aligned), \`qc\` (per-metacell QC data.frame), \`params\` (all call parameters), \`skipped_groups\` (data.frame), and \`call\`.

</div>

<div class="section level2">

## Details<a href="#details" class="anchor" aria-label="anchor"></a>

The Python SEACells backend runs inside an isolated basilisk-managed environment (built on first use; requires the suggested package \`basilisk\`); no Python knowledge is required and the R interface is self-contained. SEACells is GPLv2 and is used strictly as an external dependency, never source-copied.

Groups are processed independently; each group's cells are exported (MatrixMarket + CSV, the hdWGCNA vignette round-trip pattern), SEACells is fit with \`n_SEACells = floor(n_cells / cells_per_metacell)\` (capped by \`target_metacells\`), metacells are summarized from raw counts, and mcRigor QC is computed against the group's single-cell counts. Metacells are named \`\<group\>\_\<i\>\` and column-bound across groups. The result is hdWGCNA-interoperable: pass \`hdWGCNA::SetMetacellObject(obj, GetMetacellSEACells(seurat_obj)\$metacell_object)\` after \`hdWGCNA::SetupForWGCNA\` to use it in a WGCNA network.

</div>

<div class="section level2">

## References<a href="#references" class="anchor" aria-label="anchor"></a>

Persad et al. (2023) Nat Biotechnol, doi:10.1038/s41587-023-01716-9; Liu & Li (2025) Nat Commun, doi:10.1038/s41467-025-63626-5; Morabito et al. (2023) Cell Rep Methods, doi:10.1016/j.crmeth.2023.100498.

</div>
