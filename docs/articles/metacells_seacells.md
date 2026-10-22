<div class="page-header">

# Metacells with SEACells and mcRigor gamma optimization

<div class="d-none name">

`metacells_seacells.Rmd`

</div>

</div>

When the starting point is single-cell data, `cellprograms` can aggregate cells into **metacells** — transcriptomically homogeneous small groups — before pseudobulk factorization. Metacells are built with <a href="https://github.com/dpeerlab/SEACells" class="external-link">SEACells</a> (run in an isolated, pinned basilisk Python environment built on first use) and QC’d with mcRigor, stratified so that **metacells never cross the grouping variables** (donor x cell type by default).

This vignette runs the full workflow on a simulated dataset with planted cell-type identity, then feeds the metacell counts into [`fit_celltype_programs()`](../reference/fit_celltype_programs.html).

<div class="section level2">

## Simulate a small single-cell dataset<a href="#simulate-a-small-single-cell-dataset" class="anchor" aria-label="anchor"></a>

Two donors x two cell types, 80 cells each, 200 genes. Each cell type gets a private marker block; everything else is background Poisson noise.

<div id="cb1" class="sourceCode">

\
<span class="kw"><a href="https://rdrr.io/r/base/library.html" class="external-link"><code class="sourceCode r R">library</code></a></span><span class="op">`(`</span><span class="va">`cellprograms`</span><span class="op">`)`</span>\
<span class="fu"><a href="https://rdrr.io/r/base/Random.html" class="external-link"><code class="sourceCode r R">set.seed</code></a></span><span class="op">`(`</span><span class="fl">`11`</span><span class="op">`)`</span>\
\
<span class="va">`n_per`</span>` `<span class="op">`<-`</span>` `<span class="fl">`80`</span>\
<span class="va">`types`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link"><code class="sourceCode r R">c</code></a></span><span class="op">`(`</span><span class="st">`"Tcells"`</span>`, `<span class="st">`"Bcells"`</span><span class="op">`)`</span>\
<span class="va">`donors`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link"><code class="sourceCode r R">c</code></a></span><span class="op">`(`</span><span class="st">`"donor1"`</span>`, `<span class="st">`"donor2"`</span><span class="op">`)`</span>\
<span class="va">`groups`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/rep.html" class="external-link"><code class="sourceCode r R">rep</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/paste.html" class="external-link"><code class="sourceCode r R">paste0</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/rep.html" class="external-link"><code class="sourceCode r R">rep</code></a></span><span class="op">`(`</span><span class="va">`donors`</span>`, each `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/length.html" class="external-link"><code class="sourceCode r R">length</code></a></span><span class="op">`(`</span><span class="va">`types`</span><span class="op">`)`</span><span class="op">`)`</span>`, `<span class="st">`"_"`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/rep.html" class="external-link"><code class="sourceCode r R">rep</code></a></span><span class="op">`(`</span><span class="va">`types`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/length.html" class="external-link"><code class="sourceCode r R">length</code></a></span><span class="op">`(`</span><span class="va">`donors`</span><span class="op">`)`</span><span class="op">`)`</span><span class="op">`)`</span>`,`\
`              each `<span class="op">`=`</span>` `<span class="va">`n_per`</span><span class="op">`)`</span>\
\
<span class="va">`lam`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/matrix.html" class="external-link"><code class="sourceCode r R">matrix</code></a></span><span class="op">`(`</span><span class="fl">`0.5`</span>`, nrow `<span class="op">`=`</span>` `<span class="fl">`200`</span>`, ncol `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/length.html" class="external-link"><code class="sourceCode r R">length</code></a></span><span class="op">`(`</span><span class="va">`groups`</span><span class="op">`)`</span><span class="op">`)`</span>\
<span class="va">`markers`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/list.html" class="external-link"><code class="sourceCode r R">list</code></a></span><span class="op">`(`</span>`Tcells `<span class="op">`=`</span>` `<span class="fl">`1`</span><span class="op">`:`</span><span class="fl">`12`</span>`, Bcells `<span class="op">`=`</span>` `<span class="fl">`13`</span><span class="op">`:`</span><span class="fl">`24`</span><span class="op">`)`</span>\
<span class="kw">`for`</span>` `<span class="op">`(`</span><span class="va">`g`</span>` `<span class="kw">`in`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/unique.html" class="external-link"><code class="sourceCode r R">unique</code></a></span><span class="op">`(`</span><span class="va">`groups`</span><span class="op">`)`</span><span class="op">`)`</span>` `<span class="op">`{`</span>\
`  `<span class="va">`type`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/grep.html" class="external-link"><code class="sourceCode r R">sub</code></a></span><span class="op">`(`</span><span class="st">`".*_"`</span>`, `<span class="st">`""`</span>`, `<span class="va">`g`</span><span class="op">`)`</span>\
`  `<span class="va">`cols`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/which.html" class="external-link"><code class="sourceCode r R">which</code></a></span><span class="op">`(`</span><span class="va">`groups`</span>` `<span class="op">`==`</span>` `<span class="va">`g`</span><span class="op">`)`</span>\
`  `<span class="va">`lam`</span><span class="op">`[`</span><span class="va">`markers`</span><span class="op">`[[`</span><span class="va">`type`</span><span class="op">`]`</span><span class="op">`]`</span>`, `<span class="va">`cols`</span><span class="op">`]`</span>` `<span class="op">`<-`</span>` `<span class="fl">`8`</span>\
<span class="op">`}`</span>\
<span class="va">`counts`</span>` `<span class="op">`<-`</span>` `<span class="fu">`Matrix`</span><span class="fu">`::`</span><span class="fu"><a href="https://rdrr.io/pkg/Matrix/man/Matrix.html" class="external-link"><code class="sourceCode r R">Matrix</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/stats/Poisson.html" class="external-link"><code class="sourceCode r R">rpois</code></a></span><span class="op">`(`</span><span class="fl">`200`</span>` `<span class="op">`*`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/length.html" class="external-link"><code class="sourceCode r R">length</code></a></span><span class="op">`(`</span><span class="va">`groups`</span><span class="op">`)`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/vector.html" class="external-link"><code class="sourceCode r R">as.vector</code></a></span><span class="op">`(`</span><span class="va">`lam`</span><span class="op">`)`</span><span class="op">`)`</span>`,`\
`                         nrow `<span class="op">`=`</span>` `<span class="fl">`200`</span>`, ncol `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/length.html" class="external-link"><code class="sourceCode r R">length</code></a></span><span class="op">`(`</span><span class="va">`groups`</span><span class="op">`)`</span>`, sparse `<span class="op">`=`</span>` `<span class="cn">`TRUE`</span><span class="op">`)`</span>\
<span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">rownames</code></a></span><span class="op">`(`</span><span class="va">`counts`</span><span class="op">`)`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/paste.html" class="external-link"><code class="sourceCode r R">paste0</code></a></span><span class="op">`(`</span><span class="st">`"gene"`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/seq.html" class="external-link"><code class="sourceCode r R">seq_len</code></a></span><span class="op">`(`</span><span class="fl">`200`</span><span class="op">`)`</span><span class="op">`)`</span>\
<span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">colnames</code></a></span><span class="op">`(`</span><span class="va">`counts`</span><span class="op">`)`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/paste.html" class="external-link"><code class="sourceCode r R">paste0</code></a></span><span class="op">`(`</span><span class="st">`"cell"`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/seq.html" class="external-link"><code class="sourceCode r R">seq_len</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/length.html" class="external-link"><code class="sourceCode r R">length</code></a></span><span class="op">`(`</span><span class="va">`groups`</span><span class="op">`)`</span><span class="op">`)`</span><span class="op">`)`</span>\
\
<span class="va">`meta`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/data.frame.html" class="external-link"><code class="sourceCode r R">data.frame</code></a></span><span class="op">`(`</span>\
`  donor `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/grep.html" class="external-link"><code class="sourceCode r R">sub</code></a></span><span class="op">`(`</span><span class="st">`"_.*"`</span>`, `<span class="st">`""`</span>`, `<span class="va">`groups`</span><span class="op">`)`</span>`,`\
`  celltype `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/grep.html" class="external-link"><code class="sourceCode r R">sub</code></a></span><span class="op">`(`</span><span class="st">`".*_"`</span>`, `<span class="st">`""`</span>`, `<span class="va">`groups`</span><span class="op">`)`</span>`,`\
`  row.names `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">colnames</code></a></span><span class="op">`(`</span><span class="va">`counts`</span><span class="op">`)`</span>\
<span class="op">`)`</span>\
\
<span class="va">`obj`</span>` `<span class="op">`<-`</span>` `<span class="fu">`Seurat`</span><span class="fu">`::`</span><span class="fu"><a href="https://satijalab.github.io/seurat-object/reference/CreateSeuratObject.html" class="external-link"><code class="sourceCode r R">CreateSeuratObject</code></a></span><span class="op">`(`</span>`counts `<span class="op">`=`</span>` `<span class="va">`counts`</span>`, meta.data `<span class="op">`=`</span>` `<span class="va">`meta`</span><span class="op">`)`</span>\
<span class="va">`obj`</span>` `<span class="op">`<-`</span>` `<span class="fu">`Seurat`</span><span class="fu">`::`</span><span class="fu"><a href="https://satijalab.org/seurat/reference/NormalizeData.html" class="external-link"><code class="sourceCode r R">NormalizeData</code></a></span><span class="op">`(`</span><span class="va">`obj`</span>`, verbose `<span class="op">`=`</span>` `<span class="cn">`FALSE`</span><span class="op">`)`</span>\
<span class="va">`obj`</span>` `<span class="op">`<-`</span>` `<span class="fu">`Seurat`</span><span class="fu">`::`</span><span class="fu"><a href="https://satijalab.org/seurat/reference/FindVariableFeatures.html" class="external-link"><code class="sourceCode r R">FindVariableFeatures</code></a></span><span class="op">`(`</span><span class="va">`obj`</span>`, nfeatures `<span class="op">`=`</span>` `<span class="fl">`60`</span>`, verbose `<span class="op">`=`</span>` `<span class="cn">`FALSE`</span><span class="op">`)`</span>\
<span class="va">`obj`</span>` `<span class="op">`<-`</span>` `<span class="fu">`Seurat`</span><span class="fu">`::`</span><span class="fu"><a href="https://satijalab.org/seurat/reference/ScaleData.html" class="external-link"><code class="sourceCode r R">ScaleData</code></a></span><span class="op">`(`</span><span class="va">`obj`</span>`, verbose `<span class="op">`=`</span>` `<span class="cn">`FALSE`</span><span class="op">`)`</span>\
<span class="va">`obj`</span>` `<span class="op">`<-`</span>` `<span class="fu">`Seurat`</span><span class="fu">`::`</span><span class="fu"><a href="https://satijalab.org/seurat/reference/RunPCA.html" class="external-link"><code class="sourceCode r R">RunPCA</code></a></span><span class="op">`(`</span><span class="va">`obj`</span>`, npcs `<span class="op">`=`</span>` `<span class="fl">`10`</span>`, verbose `<span class="op">`=`</span>` `<span class="cn">`FALSE`</span><span class="op">`)`</span>

</div>

</div>

<div class="section level2">

## Build metacells per group<a href="#build-metacells-per-group" class="anchor" aria-label="anchor"></a>

[`MetacellsByGroupsSEACells()`](../reference/MetacellsByGroupsSEACells.html) fits SEACells separately within each `group.by` stratum and stores everything under the Seurat object’s `misc`:

<div id="cb2" class="sourceCode">

\
<span class="va">`obj`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`MetacellsByGroupsSEACells`](../reference/MetacellsByGroupsSEACells.html)</span><span class="op">`(`</span><span class="va">`obj`</span>`,`\
`                                 group.by `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link"><code class="sourceCode r R">c</code></a></span><span class="op">`(`</span><span class="st">`"donor"`</span>`, `<span class="st">`"celltype"`</span><span class="op">`)`</span>`,`\
`                                 cells_per_metacell `<span class="op">`=`</span>` `<span class="fl">`20`</span>`,`\
`                                 min_cells `<span class="op">`=`</span>` `<span class="fl">`50`</span>`,`\
`                                 n_top_genes `<span class="op">`=`</span>` `<span class="fl">`60`</span>`,`\
`                                 n_pc `<span class="op">`=`</span>` `<span class="fl">`10`</span>`,`\
`                                 seed `<span class="op">`=`</span>` `<span class="fl">`11`</span>`, verbose `<span class="op">`=`</span>` `<span class="cn">`FALSE`</span><span class="op">`)`</span>\
<span class="va">`res`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`GetMetacellSEACells`](../reference/GetMetacellSEACells.html)</span><span class="op">`(`</span><span class="va">`obj`</span><span class="op">`)`</span>\
<span class="fu"><a href="https://rdrr.io/r/base/names.html" class="external-link"><code class="sourceCode r R">names</code></a></span><span class="op">`(`</span><span class="va">`res`</span><span class="op">`)`</span>

</div>

The result carries the metacell Seurat object, the soft (cells x archetypes) assignment matrix, per-metacell mcRigor QC, fit statistics, the parameters used, and any skipped groups:

<div id="cb3" class="sourceCode">

\
<span class="fu"><a href="https://rdrr.io/r/utils/head.html" class="external-link"><code class="sourceCode r R">head</code></a></span><span class="op">`(`</span><span class="va">`res`</span><span class="op">`$`</span><span class="va">`qc`</span><span class="op">`)`</span>`                 `<span class="co">`# mcDiv, dubious flag, ZeroRate per metacell`</span>\
<span class="va">`res`</span><span class="op">`$`</span><span class="va">`skipped_groups`</span>`           `<span class="co">`# groups dropped (e.g. too few cells) - usually empty`</span>

</div>

Groups with too few cells (below `min_cells`) are skipped with a warning rather than failing the call, and degenerate SEACells fits (all hard assignments collapsing to one archetype) drop that group’s soft assignments gracefully — check `skipped_groups` when working with small strata.

</div>

<div class="section level2">

## Optimize the number of metacells (gamma)<a href="#optimize-the-number-of-metacells-gamma" class="anchor" aria-label="anchor"></a>

SEACells’ `gamma` is the number of archetypes — the granularity. [`OptimizeGammaSEACells()`](../reference/OptimizeGammaSEACells.html) refits SEACells across candidate gammas per group and picks the winner by the mcRigor optimization score:

<div id="cb4" class="sourceCode">

\
<span class="va">`opt`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`OptimizeGammaSEACells`](../reference/OptimizeGammaSEACells.html)</span><span class="op">`(`</span><span class="va">`obj`</span>`,`\
`                             group.by `<span class="op">`=`</span>` `<span class="st">`"donor"`</span>`,`\
`                             gammas `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link"><code class="sourceCode r R">c</code></a></span><span class="op">`(`</span><span class="fl">`10`</span>`, `<span class="fl">`20`</span><span class="op">`)`</span>`,`\
`                             n_top_genes `<span class="op">`=`</span>` `<span class="fl">`60`</span>`,`\
`                             n_pc `<span class="op">`=`</span>` `<span class="fl">`10`</span>`,`\
`                             min_cells `<span class="op">`=`</span>` `<span class="fl">`50`</span>`,`\
`                             mcrigor_args `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/list.html" class="external-link"><code class="sourceCode r R">list</code></a></span><span class="op">`(`</span>`Nrep `<span class="op">`=`</span>` `<span class="fl">`5`</span><span class="op">`)`</span>`,`\
`                             seed `<span class="op">`=`</span>` `<span class="fl">`11`</span>`, verbose `<span class="op">`=`</span>` `<span class="cn">`FALSE`</span><span class="op">`)`</span>\
<span class="va">`opt`</span><span class="op">`$`</span><span class="va">`optimal_gamma`</span>`            `<span class="co">`# winning gamma per donor group`</span>\
<span class="va">`opt`</span><span class="op">`$`</span><span class="va">`scores`</span>`                   `<span class="co">`# mcRigor score per gamma x group`</span>

</div>

Gammas whose membership labels degenerate (fewer than two distinct labels after Seurat v5 aggregation) are dropped with a warning and excluded from the pick rather than crashing the run.

</div>

<div class="section level2">

## Metacell counts feed the program fit directly<a href="#metacell-counts-feed-the-program-fit-directly" class="anchor" aria-label="anchor"></a>

<div id="cb5" class="sourceCode">

\
<span class="va">`mc`</span>` `<span class="op">`<-`</span>` `<span class="va">`res`</span><span class="op">`$`</span><span class="va">`metacell_object`</span>\
<span class="va">`counts_mc`</span>` `<span class="op">`<-`</span>` `<span class="fu">`SeuratObject`</span><span class="fu">`::`</span><span class="fu"><a href="https://satijalab.github.io/seurat-object/reference/Layers.html" class="external-link"><code class="sourceCode r R">LayerData</code></a></span><span class="op">`(`</span><span class="va">`mc`</span><span class="op">`[[`</span><span class="st">`"RNA"`</span><span class="op">`]`</span><span class="op">`]`</span>`, layer `<span class="op">`=`</span>` `<span class="st">`"counts"`</span><span class="op">`)`</span>`  `<span class="co">`# genes x metacells`</span>\
<span class="va">`by_type`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/split.html" class="external-link"><code class="sourceCode r R">split</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">colnames</code></a></span><span class="op">`(`</span><span class="va">`mc`</span><span class="op">`)`</span>`, `<span class="va">`mc`</span><span class="op">`$`</span><span class="va">`celltype`</span><span class="op">`)`</span>\
<span class="va">`pb`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/lapply.html" class="external-link"><code class="sourceCode r R">lapply</code></a></span><span class="op">`(`</span><span class="va">`by_type`</span>`, `<span class="kw">`function`</span><span class="op">`(`</span><span class="va">`cells`</span><span class="op">`)`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/t.html" class="external-link"><code class="sourceCode r R">t</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/matrix.html" class="external-link"><code class="sourceCode r R">as.matrix</code></a></span><span class="op">`(`</span><span class="va">`counts_mc`</span><span class="op">`[`</span>`, `<span class="va">`cells`</span>`, drop `<span class="op">`=`</span>` `<span class="cn">`FALSE`</span><span class="op">`]`</span><span class="op">`)`</span><span class="op">`)`</span><span class="op">`)`</span>\
\
<span class="va">`x`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`as_cell_program_data`](../reference/as_cell_program_data.html)</span><span class="op">`(`</span><span class="va">`pb`</span><span class="op">`)`</span>\
<span class="va">`fit`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`fit_celltype_programs`](../reference/fit_celltype_programs.html)</span><span class="op">`(`</span><span class="va">`x`</span>`, loading_prior `<span class="op">`=`</span>` `<span class="st">`"point_laplace"`</span>`,`\
`                             max_factors `<span class="op">`=`</span>` `<span class="fl">`4`</span>`, features `<span class="op">`=`</span>` `<span class="st">`"variable"`</span>`,`\
`                             n_variable_genes `<span class="op">`=`</span>` `<span class="fl">`100`</span><span class="op">`)`</span>\
<span class="va">`fit`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`canonicalize_programs`](../reference/canonicalize_programs.html)</span><span class="op">`(`</span><span class="va">`fit`</span><span class="op">`)`</span>\
<span class="fu">[`program_summary`](../reference/program_summary.html)</span><span class="op">`(`</span><span class="va">`fit`</span><span class="op">`)`</span><span class="op">`[`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link"><code class="sourceCode r R">c</code></a></span><span class="op">`(`</span><span class="st">`"cell_type"`</span>`, `<span class="st">`"n_programs"`</span>`, `<span class="st">`"R2"`</span>`, `<span class="st">`"elbo"`</span><span class="op">`)`</span><span class="op">`]`</span>

</div>

The planted structure survives the aggregation: each cell type’s programs are dominated by its private marker block.

<div id="cb6" class="sourceCode">

\
<span class="va">`marker_expr`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/matrix.html" class="external-link"><code class="sourceCode r R">as.matrix</code></a></span><span class="op">`(`</span><span class="va">`counts_mc`</span><span class="op">`[`</span><span class="fu"><a href="https://rdrr.io/r/base/unlist.html" class="external-link"><code class="sourceCode r R">unlist</code></a></span><span class="op">`(`</span><span class="va">`markers`</span><span class="op">`)`</span>`, , drop `<span class="op">`=`</span>` `<span class="cn">`FALSE`</span><span class="op">`]`</span><span class="op">`)`</span>\
<span class="va">`marker_expr`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/t.html" class="external-link"><code class="sourceCode r R">t</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/scale.html" class="external-link"><code class="sourceCode r R">scale</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/t.html" class="external-link"><code class="sourceCode r R">t</code></a></span><span class="op">`(`</span><span class="va">`marker_expr`</span><span class="op">`)`</span><span class="op">`)`</span><span class="op">`)`</span>\
<span class="va">`df`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/data.frame.html" class="external-link"><code class="sourceCode r R">data.frame</code></a></span><span class="op">`(`</span>\
`  metacell `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/factor.html" class="external-link"><code class="sourceCode r R">factor</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">colnames</code></a></span><span class="op">`(`</span><span class="va">`marker_expr`</span><span class="op">`)`</span>`, levels `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">colnames</code></a></span><span class="op">`(`</span><span class="va">`marker_expr`</span><span class="op">`)`</span><span class="op">`)`</span>`,`\
`  gene `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/rep.html" class="external-link"><code class="sourceCode r R">rep</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">rownames</code></a></span><span class="op">`(`</span><span class="va">`marker_expr`</span><span class="op">`)`</span>`, times `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/nrow.html" class="external-link"><code class="sourceCode r R">ncol</code></a></span><span class="op">`(`</span><span class="va">`marker_expr`</span><span class="op">`)`</span><span class="op">`)`</span>`,`\
`  z `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/vector.html" class="external-link"><code class="sourceCode r R">as.vector</code></a></span><span class="op">`(`</span><span class="va">`marker_expr`</span><span class="op">`)`</span>\
<span class="op">`)`</span>\
<span class="fu">`ggplot2`</span><span class="fu">`::`</span><span class="fu"><a href="https://ggplot2.tidyverse.org/reference/ggplot.html" class="external-link"><code class="sourceCode r R">ggplot</code></a></span><span class="op">`(`</span><span class="va">`df`</span>`, `<span class="fu">`ggplot2`</span><span class="fu">`::`</span><span class="fu"><a href="https://ggplot2.tidyverse.org/reference/aes.html" class="external-link"><code class="sourceCode r R">aes</code></a></span><span class="op">`(`</span><span class="va">`metacell`</span>`, `<span class="va">`gene`</span>`, fill `<span class="op">`=`</span>` `<span class="va">`z`</span><span class="op">`)`</span><span class="op">`)`</span>` `<span class="op">`+`</span>\
`  `<span class="fu">`ggplot2`</span><span class="fu">`::`</span><span class="fu"><a href="https://ggplot2.tidyverse.org/reference/geom_tile.html" class="external-link"><code class="sourceCode r R">geom_tile</code></a></span><span class="op">`(`</span><span class="op">`)`</span>` `<span class="op">`+`</span>\
`  `<span class="fu">`ggplot2`</span><span class="fu">`::`</span><span class="fu"><a href="https://ggplot2.tidyverse.org/reference/scale_gradient.html" class="external-link"><code class="sourceCode r R">scale_fill_gradient2</code></a></span><span class="op">`(`</span><span class="op">`)`</span>` `<span class="op">`+`</span>\
`  `<span class="fu">`ggplot2`</span><span class="fu">`::`</span><span class="fu"><a href="https://ggplot2.tidyverse.org/reference/labs.html" class="external-link"><code class="sourceCode r R">labs</code></a></span><span class="op">`(`</span>`x `<span class="op">`=`</span>` `<span class="st">`"metacell"`</span>`, y `<span class="op">`=`</span>` `<span class="st">`"marker gene"`</span><span class="op">`)`</span>` `<span class="op">`+`</span>\
`  `<span class="fu">`ggplot2`</span><span class="fu">`::`</span><span class="fu"><a href="https://ggplot2.tidyverse.org/reference/theme.html" class="external-link"><code class="sourceCode r R">theme</code></a></span><span class="op">`(`</span>`axis.text.x `<span class="op">`=`</span>` `<span class="fu">`ggplot2`</span><span class="fu">`::`</span><span class="fu"><a href="https://ggplot2.tidyverse.org/reference/element.html" class="external-link"><code class="sourceCode r R">element_blank</code></a></span><span class="op">`(`</span><span class="op">`)`</span>`,`\
`                 axis.text.y `<span class="op">`=`</span>` `<span class="fu">`ggplot2`</span><span class="fu">`::`</span><span class="fu"><a href="https://ggplot2.tidyverse.org/reference/element.html" class="external-link"><code class="sourceCode r R">element_text</code></a></span><span class="op">`(`</span>`size `<span class="op">`=`</span>` `<span class="fl">`6`</span><span class="op">`)`</span><span class="op">`)`</span>

</div>

</div>

<div class="section level2">

## Session info<a href="#session-info" class="anchor" aria-label="anchor"></a>

    #> R version 4.4.2 (2024-10-31)
    #> Platform: x86_64-conda-linux-gnu
    #> Running under: Ubuntu 24.04.5 LTS
    #>
    #> Matrix products: default
    #> BLAS/LAPACK: /opt/conda/lib/libopenblasp-r0.3.34.so;  LAPACK version 3.12.0
    #>
    #> locale:
    #>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8<br>
    #>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8<br>
    #>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C<br>
    #> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C<br>
    #>
    #> time zone: Etc/UTC
    #> tzcode source: system (glibc)
    #>
    #> attached base packages:
    #> [1] stats     graphics  grDevices utils     datasets  methods   base<br>
    #>
    #> other attached packages:
    #> [1] cellprograms_0.2.0.9001
    #>
    #> loaded via a namespace (and not attached):
    #>   [1] deldir_2.0-4           pbapply_1.7-4          gridExtra_2.3.1<br>
    #>   [4] rlang_1.3.0            magrittr_2.0.5         RcppAnnoy_0.0.23<br>
    #>   [7] otel_0.2.0             spatstat.geom_3.8-2    matrixStats_1.5.0<br>
    #>  [10] ggridges_0.5.7         compiler_4.4.2         png_0.1-9<br>
    #>  [13] systemfonts_1.3.2      vctrs_0.7.3            reshape2_1.4.5<br>
    #>  [16] stringr_1.6.0          pkgconfig_2.0.3        fastmap_1.2.0<br>
    #>  [19] promises_1.5.0         rmarkdown_2.32         ragg_1.5.2<br>
    #>  [22] purrr_1.2.2            xfun_0.61              cachem_1.1.0<br>
    #>  [25] jsonlite_2.0.0         goftest_1.2-3          later_1.4.8<br>
    #>  [28] spatstat.utils_3.2-5   irlba_2.4.1            parallel_4.4.2<br>
    #>  [31] cluster_2.1.8.3        R6_2.6.1               ica_1.0-3<br>
    #>  [34] spatstat.data_3.1-9    bslib_0.12.0           stringi_1.8.7<br>
    #>  [37] RColorBrewer_1.1-3     reticulate_1.47.0      spatstat.univar_3.2-0
    #>  [40] parallelly_1.48.0      lmtest_0.9-40          jquerylib_0.1.4<br>
    #>  [43] scattermore_1.2        Rcpp_1.1.2             knitr_1.52<br>
    #>  [46] tensor_1.5.1           future.apply_1.21.0    zoo_1.9-0<br>
    #>  [49] sctransform_0.4.3      httpuv_1.6.17          Matrix_1.7-6<br>
    #>  [52] splines_4.4.2          igraph_2.1.4           tidyselect_1.2.1<br>
    #>  [55] abind_1.4-8            yaml_2.3.12            spatstat.random_3.5-1
    #>  [58] spatstat.explore_3.8-2 codetools_0.2-20       miniUI_0.1.2<br>
    #>  [61] listenv_1.1.0          lattice_0.23-1         tibble_3.3.1<br>
    #>  [64] plyr_1.8.9             shiny_1.14.0           S7_0.2.2<br>
    #>  [67] ROCR_1.0-12            evaluate_1.0.5         Rtsne_0.17<br>
    #>  [70] future_1.80.0          fastDummies_1.7.6      desc_1.4.3<br>
    #>  [73] survival_3.8-12        polyclip_1.10-7        fitdistrplus_1.2-6<br>
    #>  [76] pillar_1.11.1          Seurat_5.5.1           KernSmooth_2.23-27<br>
    #>  [79] plotly_4.12.1          generics_0.1.4         RcppHNSW_0.7.0<br>
    #>  [82] sp_2.2-3               ggplot2_4.0.3          scales_1.4.0<br>
    #>  [85] globals_0.20.0         xtable_1.8-8           glue_1.8.1<br>
    #>  [88] tools_4.4.2            data.table_1.18.6.1    RSpectra_0.16-2<br>
    #>  [91] RANN_2.6.2             fs_2.1.0               dotCall64_1.2<br>
    #>  [94] cowplot_1.2.0          grid_4.4.2             tidyr_1.3.2<br>
    #>  [97] nlme_3.1-170           patchwork_1.3.2        cli_3.6.6<br>
    #> [100] spatstat.sparse_3.2-0  textshaping_0.4.0      spam_2.11-4<br>
    #> [103] viridisLite_0.4.3      dplyr_1.2.1            uwot_0.2.5<br>
    #> [106] gtable_0.3.6           sass_0.4.10            digest_0.6.39<br>
    #> [109] progressr_1.1.0        ggrepel_0.9.8          htmlwidgets_1.6.4<br>
    #> [112] SeuratObject_5.4.0     farver_2.1.2           htmltools_0.5.9<br>
    #> [115] pkgdown_2.2.1          lifecycle_1.0.5        httr_1.4.9<br>
    #> [118] mime_0.13              MASS_7.3-66

</div>
