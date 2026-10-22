<div class="page-header">

# Quick start: from pseudobulk to shared programs

<div class="d-none name">

`quickstart.Rmd`

</div>

</div>

`cellprograms` discovers transcriptional programs **de novo within each cell type** from pseudobulk single-cell RNA-seq, then investigates cross-cell-type coordination and sharing afterwards — never assuming it. This vignette runs the whole pipeline on a small simulated dataset.

<div class="section level2">

## Simulate pseudobulk data<a href="#simulate-pseudobulk-data" class="anchor" aria-label="anchor"></a>

[`simulate_pseudobulk()`](../reference/simulate_pseudobulk.html) plants private and shared programs with known ground truth. The result is a named list of observations x genes matrices, one per cell type.

<div id="cb1" class="sourceCode">

\
<span class="kw"><a href="https://rdrr.io/r/base/library.html" class="external-link"><code class="sourceCode r R">library</code></a></span><span class="op">`(`</span><span class="va">`cellprograms`</span><span class="op">`)`</span>\
\
<span class="va">`sim`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`simulate_pseudobulk`](../reference/simulate_pseudobulk.html)</span><span class="op">`(`</span>`scenario `<span class="op">`=`</span>` `<span class="st">`"mixed"`</span>`, n_samples `<span class="op">`=`</span>` `<span class="fl">`60`</span>`,`\
`                           n_genes `<span class="op">`=`</span>` `<span class="fl">`800`</span>`, program_size `<span class="op">`=`</span>` `<span class="fl">`40`</span>`, seed `<span class="op">`=`</span>` `<span class="fl">`1`</span><span class="op">`)`</span>\
<span class="fu"><a href="https://rdrr.io/r/utils/str.html" class="external-link"><code class="sourceCode r R">str</code></a></span><span class="op">`(`</span><span class="va">`sim`</span><span class="op">`$`</span><span class="va">`pseudobulk`</span>`, max.level `<span class="op">`=`</span>` `<span class="fl">`1`</span><span class="op">`)`</span>\
<span class="co">`#> List of 5`</span>\
<span class="co">`#>  $ B   : num [1:60, 1:800] 1.146 -2.403 0.573 0.375 -0.425 ...`</span>\
<span class="co">`#>   ..- attr(*, "dimnames")=List of 2`</span>\
<span class="co">`#>  $ CD4 : num [1:60, 1:800] -1.474 0.59 -0.13 0.372 1.452 ...`</span>\
<span class="co">`#>   ..- attr(*, "dimnames")=List of 2`</span>\
<span class="co">`#>  $ CD8 : num [1:60, 1:800] 0.0679 0.9149 0.8625 -1.1996 -2.2594 ...`</span>\
<span class="co">`#>   ..- attr(*, "dimnames")=List of 2`</span>\
<span class="co">`#>  $ NK  : num [1:60, 1:800] -0.0478 1.5893 1.175 0.0746 -0.9374 ...`</span>\
<span class="co">`#>   ..- attr(*, "dimnames")=List of 2`</span>\
<span class="co">`#>  $ Mono: num [1:60, 1:800] -0.7613 -0.4367 -1.4659 0.0574 -1.4384 ...`</span>\
<span class="co">`#>   ..- attr(*, "dimnames")=List of 2`</span>

</div>

</div>

<div class="section level2">

## Build the data object and fit programs<a href="#build-the-data-object-and-fit-programs" class="anchor" aria-label="anchor"></a>

<div id="cb2" class="sourceCode">

\
<span class="va">`x`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`as_cell_program_data`](../reference/as_cell_program_data.html)</span><span class="op">`(`</span><span class="va">`sim`</span><span class="op">`$`</span><span class="va">`pseudobulk`</span><span class="op">`)`</span>\
<span class="fu">[`validate_cell_program_data`](../reference/validate_cell_program_data.html)</span><span class="op">`(`</span><span class="va">`x`</span><span class="op">`)`</span>\
\
<span class="va">`fit`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`fit_celltype_programs`](../reference/fit_celltype_programs.html)</span><span class="op">`(`</span><span class="va">`x`</span>`, loading_prior `<span class="op">`=`</span>` `<span class="st">`"point_laplace"`</span>`,`\
`                             max_factors `<span class="op">`=`</span>` `<span class="fl">`6`</span>`, features `<span class="op">`=`</span>` `<span class="st">`"variable"`</span>`,`\
`                             n_variable_genes `<span class="op">`=`</span>` `<span class="fl">`400`</span><span class="op">`)`</span>\
<span class="va">`fit`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`canonicalize_programs`](../reference/canonicalize_programs.html)</span><span class="op">`(`</span><span class="va">`fit`</span><span class="op">`)`</span>\
<span class="va">`fit`</span>\
<span class="co">`#> `<span style="font-weight: bold;">`cell_program_data`</span></span>\
<span class="co">`#> Cell types: 5 `<span style="color: #0000BB;">`"B"`</span>`, `<span style="color: #0000BB;">`"CD4"`</span>`, `<span style="color: #0000BB;">`"CD8"`</span>`, `<span style="color: #0000BB;">`"NK"`</span>`, and `<span style="color: #0000BB;">`"Mono"`</span></span>\
<span class="co">`#> B: 60 observations x 800 genes`</span>\
<span class="co">`#> CD4: 60 observations x 800 genes`</span>\
<span class="co">`#> CD8: 60 observations x 800 genes`</span>\
<span class="co">`#> NK: 60 observations x 800 genes`</span>\
<span class="co">`#> Mono: 60 observations x 800 genes`</span>\
<span class="co">`#> Shared genes across all cell types: 800`</span>

</div>

[`fit_celltype_programs()`](../reference/fit_celltype_programs.html) fits one empirical-Bayes matrix factorization per cell type (‘flashier’); [`canonicalize_programs()`](../reference/canonicalize_programs.html) resolves sign and scale ambiguity so programs are comparable downstream. Real-data fits should additionally pass `covariates` (e.g. `c("site", "sex")`) and use `max_factors = 30` when the goal is CMP detection (see the CMP benchmark note in `README`).

</div>

<div class="section level2">

## From single cells: metacells (optional)<a href="#from-single-cells-metacells-optional" class="anchor" aria-label="anchor"></a>

When the starting point is a single-cell object rather than pseudobulk profiles, [`MetacellsByGroupsSEACells()`](../reference/MetacellsByGroupsSEACells.html) aggregates cells into transcriptomically homogeneous metacells with SEACells (run in an isolated basilisk-managed Python environment, built on first use) and mcRigor QC, stratified so metacells never cross donors or cell types. The metacell counts feed [`fit_celltype_programs()`](../reference/fit_celltype_programs.html) directly:

<div id="cb3" class="sourceCode">

\
<span class="co">`# Requires: install.packages("Seurat"); remotes::install_github("JSB-UCLA/mcRigor");`</span>\
<span class="co">`# BiocManager::install("basilisk")`</span>\
<span class="va">`seu`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`MetacellsByGroupsSEACells`](../reference/MetacellsByGroupsSEACells.html)</span><span class="op">`(`</span>\
`  `<span class="va">`seurat_obj`</span>`,                         `<span class="co">`# Seurat object with raw counts + PCA`</span>\
`  group.by `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link"><code class="sourceCode r R">c</code></a></span><span class="op">`(`</span><span class="st">`"donor"`</span>`, `<span class="st">`"celltype"`</span><span class="op">`)`</span>`,  `<span class="co">`# stratification: metacells stay within donor x cell type`</span>\
`  cells_per_metacell `<span class="op">`=`</span>` `<span class="fl">`75`</span>`,`\
`  min_cells `<span class="op">`=`</span>` `<span class="fl">`100`</span>\
<span class="op">`)`</span>\
<span class="va">`res`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`GetMetacellSEACells`](../reference/GetMetacellSEACells.html)</span><span class="op">`(`</span><span class="va">`seu`</span><span class="op">`)`</span>\
<span class="va">`mc`</span>` `<span class="op">`<-`</span>` `<span class="va">`res`</span><span class="op">`$`</span><span class="va">`metacell_object`</span>\
<span class="fu"><a href="https://rdrr.io/r/utils/head.html" class="external-link"><code class="sourceCode r R">head</code></a></span><span class="op">`(`</span><span class="va">`res`</span><span class="op">`$`</span><span class="va">`qc`</span><span class="op">`)`</span>`  `<span class="co">`# per-metacell mcRigor QC (mcDiv, divergence threshold, dubious flag)`</span>\
\
<span class="co">`# Metacell counts -> named list of metacells x genes matrices, one per cell type`</span>\
<span class="va">`counts`</span>` `<span class="op">`<-`</span>` `<span class="fu">`SeuratObject`</span><span class="fu">`::`</span><span class="fu"><a href="https://satijalab.github.io/seurat-object/reference/Layers.html" class="external-link"><code class="sourceCode r R">LayerData</code></a></span><span class="op">`(`</span><span class="va">`mc`</span><span class="op">`[[`</span><span class="st">`"RNA"`</span><span class="op">`]`</span><span class="op">`]`</span>`, layer `<span class="op">`=`</span>` `<span class="st">`"counts"`</span><span class="op">`)`</span>`  `<span class="co">`# genes x metacells`</span>\
<span class="va">`by_type`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/split.html" class="external-link"><code class="sourceCode r R">split</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">colnames</code></a></span><span class="op">`(`</span><span class="va">`mc`</span><span class="op">`)`</span>`, `<span class="va">`mc`</span><span class="op">`$`</span><span class="va">`celltype`</span><span class="op">`)`</span>\
<span class="va">`pb`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/lapply.html" class="external-link"><code class="sourceCode r R">lapply</code></a></span><span class="op">`(`</span><span class="va">`by_type`</span>`, `<span class="kw">`function`</span><span class="op">`(`</span><span class="va">`cells`</span><span class="op">`)`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/t.html" class="external-link"><code class="sourceCode r R">t</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/matrix.html" class="external-link"><code class="sourceCode r R">as.matrix</code></a></span><span class="op">`(`</span><span class="va">`counts`</span><span class="op">`[`</span>`, `<span class="va">`cells`</span>`, drop `<span class="op">`=`</span>` `<span class="cn">`FALSE`</span><span class="op">`]`</span><span class="op">`)`</span><span class="op">`)`</span><span class="op">`)`</span>\
\
<span class="va">`x`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`as_cell_program_data`](../reference/as_cell_program_data.html)</span><span class="op">`(`</span><span class="va">`pb`</span><span class="op">`)`</span>\
<span class="va">`fit`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`fit_celltype_programs`](../reference/fit_celltype_programs.html)</span><span class="op">`(`</span><span class="va">`x`</span>`, loading_prior `<span class="op">`=`</span>` `<span class="st">`"point_laplace"`</span>`,`\
`                             max_factors `<span class="op">`=`</span>` `<span class="fl">`6`</span>`, features `<span class="op">`=`</span>` `<span class="st">`"variable"`</span><span class="op">`)`</span>\
<span class="va">`fit`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`canonicalize_programs`](../reference/canonicalize_programs.html)</span><span class="op">`(`</span><span class="va">`fit`</span><span class="op">`)`</span>

</div>

[`OptimizeGammaSEACells()`](../reference/OptimizeGammaSEACells.html) selects the granularity (`cells_per_metacell`) by the mcRigor score when you prefer data-driven granularity over the default.

</div>

<div class="section level2">

## Inspect programs<a href="#inspect-programs" class="anchor" aria-label="anchor"></a>

<div id="cb4" class="sourceCode">

\
<span class="fu">[`program_summary`](../reference/program_summary.html)</span><span class="op">`(`</span><span class="va">`fit`</span><span class="op">`)`</span>`          `<span class="co">`# one row per cell type: n_programs, R2, elbo`</span>\
<span class="co">`#>   cell_type n_obs n_genes n_programs        R2      elbo`</span>\
<span class="co">`#> 1         B    60     400          2 0.2439779 -36065.09`</span>\
<span class="co">`#> 2       CD4    60     400          1 0.1707640 -35709.99`</span>\
<span class="co">`#> 3       CD8    60     400          1 0.2010540 -35859.83`</span>\
<span class="co">`#> 4        NK    60     400          2 0.2948849 -36081.50`</span>\
<span class="co">`#> 5      Mono    60     400          2 0.2997384 -35859.12`</span>\
<span class="fu">[`top_genes`](../reference/top_genes.html)</span><span class="op">`(`</span><span class="va">`fit`</span>`, `<span class="st">`"B_1"`</span>`, n `<span class="op">`=`</span>` `<span class="fl">`8`</span><span class="op">`)`</span>`  `<span class="co">`# top loadings of program B_1`</span>\
<span class="co">`#>    gene   loading          lfsr`</span>\
<span class="co">`#> 1  G474  4.582835 2.630908e-244`</span>\
<span class="co">`#> 3  G733  3.164218 1.075143e-117`</span>\
<span class="co">`#> 2  G422  2.724563  7.611642e-88`</span>\
<span class="co">`#> 4  G197 -2.688759  0.000000e+00`</span>\
<span class="co">`#> 6  G462  2.555132  1.603963e-77`</span>\
<span class="co">`#> 8  G349 -2.213170  0.000000e+00`</span>\
<span class="co">`#> 9  G514  2.067245  1.752200e-51`</span>\
<span class="co">`#> 10 G720  2.031627  8.535883e-50`</span>

</div>

Scores (donor x program activities) and loadings (gene x program) are one call away:

<div id="cb5" class="sourceCode">

\
<span class="va">`Z`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`program_scores`](../reference/program_scores.html)</span><span class="op">`(`</span><span class="va">`fit`</span><span class="op">`)`</span>`    `<span class="co">`# donors x all programs (wide matrix)`</span>\
<span class="va">`W`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`program_loadings`](../reference/program_loadings.html)</span><span class="op">`(`</span><span class="va">`fit`</span><span class="op">`)`</span>`  `<span class="co">`# genes x all programs`</span>\
<span class="fu"><a href="https://rdrr.io/r/base/dim.html" class="external-link"><code class="sourceCode r R">dim</code></a></span><span class="op">`(`</span><span class="va">`Z`</span><span class="op">`)`</span>`; `<span class="fu"><a href="https://rdrr.io/r/base/dim.html" class="external-link"><code class="sourceCode r R">dim</code></a></span><span class="op">`(`</span><span class="va">`W`</span><span class="op">`)`</span>\
<span class="co">`#> [1] 60  8`</span>\
<span class="co">`#> NULL`</span>

</div>

The same views as figures — donor x program score heatmap (with the planted context as an annotation strip), top loadings, and per-program variance explained:

<div id="cb6" class="sourceCode">

\
<span class="va">`sev`</span>` `<span class="op">`<-`</span>` `<span class="va">`sim`</span><span class="op">`$`</span><span class="va">`truth`</span><span class="op">`$`</span><span class="va">`B`</span><span class="op">`[[`</span><span class="fl">`1`</span><span class="op">`]`</span><span class="op">`]`</span><span class="op">`$`</span><span class="va">`z`</span>` `<span class="op">`>`</span>` `<span class="fu"><a href="https://rdrr.io/r/stats/median.html" class="external-link"><code class="sourceCode r R">median</code></a></span><span class="op">`(`</span><span class="va">`sim`</span><span class="op">`$`</span><span class="va">`truth`</span><span class="op">`$`</span><span class="va">`B`</span><span class="op">`[[`</span><span class="fl">`1`</span><span class="op">`]`</span><span class="op">`]`</span><span class="op">`$`</span><span class="va">`z`</span><span class="op">`)`</span>\
<span class="va">`md`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/data.frame.html" class="external-link"><code class="sourceCode r R">data.frame</code></a></span><span class="op">`(`</span>`observation_id `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/paste.html" class="external-link"><code class="sourceCode r R">paste0</code></a></span><span class="op">`(`</span><span class="st">`"obs"`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/seq.html" class="external-link"><code class="sourceCode r R">seq_along</code></a></span><span class="op">`(`</span><span class="va">`sev`</span><span class="op">`)`</span><span class="op">`)`</span>`,`\
`                 planted `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/factor.html" class="external-link"><code class="sourceCode r R">factor</code></a></span><span class="op">`(`</span><span class="va">`sev`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link"><code class="sourceCode r R">c</code></a></span><span class="op">`(`</span><span class="cn">`FALSE`</span>`, `<span class="cn">`TRUE`</span><span class="op">`)`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link"><code class="sourceCode r R">c</code></a></span><span class="op">`(`</span><span class="st">`"low"`</span>`, `<span class="st">`"high"`</span><span class="op">`)`</span><span class="op">`)`</span><span class="op">`)`</span>\
<span class="va">`fit_ann`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`as_cell_program_data`](../reference/as_cell_program_data.html)</span><span class="op">`(`</span><span class="va">`sim`</span><span class="op">`$`</span><span class="va">`pseudobulk`</span>`, sample_metadata `<span class="op">`=`</span>` `<span class="va">`md`</span><span class="op">`)`</span>\
<span class="va">`fit_ann`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`canonicalize_programs`](../reference/canonicalize_programs.html)</span><span class="op">`(`</span><span class="fu">[`fit_celltype_programs`](../reference/fit_celltype_programs.html)</span><span class="op">`(`</span><span class="va">`fit_ann`</span>`, loading_prior `<span class="op">`=`</span>` `<span class="st">`"point_laplace"`</span>`,`\
`                                                       max_factors `<span class="op">`=`</span>` `<span class="fl">`6`</span>`, features `<span class="op">`=`</span>` `<span class="st">`"variable"`</span>`,`\
`                                                       n_variable_genes `<span class="op">`=`</span>` `<span class="fl">`400`</span><span class="op">`)`</span><span class="op">`)`</span>\
<span class="fu">[`plot_program_scores`](../reference/plot_program_scores.html)</span><span class="op">`(`</span><span class="va">`fit_ann`</span>`, annotation `<span class="op">`=`</span>` `<span class="st">`"planted"`</span><span class="op">`)`</span>

</div>

<img src="quickstart_files/figure-html/viz-inspect-1.png" class="r-plt" width="1200" />

<div id="cb7" class="sourceCode">

\
<span class="fu">[`plot_top_genes`](../reference/plot_top_genes.html)</span><span class="op">`(`</span><span class="va">`fit`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link"><code class="sourceCode r R">c</code></a></span><span class="op">`(`</span><span class="st">`"B_1"`</span>`, `<span class="st">`"CD4_1"`</span><span class="op">`)`</span>`, n `<span class="op">`=`</span>` `<span class="fl">`8`</span><span class="op">`)`</span>

</div>

<img src="quickstart_files/figure-html/viz-inspect-2.png" class="r-plt" width="1200" />

<div id="cb8" class="sourceCode">

\
<span class="fu">[`plot_variance_explained`](../reference/plot_variance_explained.html)</span><span class="op">`(`</span><span class="va">`fit`</span><span class="op">`)`</span>

</div>

<img src="quickstart_files/figure-html/viz-inspect-3.png" class="r-plt" width="1200" />

</div>

<div class="section level2">

## Stability<a href="#stability" class="anchor" aria-label="anchor"></a>

[`assess_program_stability()`](../reference/assess_program_stability.html) refits on donor subsamples and reports how often each program is recovered. Refits that fail are counted, not silently dropped: `fit$stability` carries `valid_refits`/`failed_refits` and a per-replicate `replicate_log`; `recovery_freq` keeps its legacy meaning (hits over requested iterations) while `conditional_recovery_freq` is computed over successful refits (all refits failing yields `NA`, not 0).

<div id="cb9" class="sourceCode">

\
<span class="va">`fit`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`assess_program_stability`](../reference/assess_program_stability.html)</span><span class="op">`(`</span><span class="va">`fit`</span>`, n_boot `<span class="op">`=`</span>` `<span class="fl">`4`</span>`, resample_frac `<span class="op">`=`</span>` `<span class="fl">`0.7`</span>`,`\
`                                seed `<span class="op">`=`</span>` `<span class="fl">`1`</span>`, progress `<span class="op">`=`</span>` `<span class="cn">`FALSE`</span><span class="op">`)`</span>\
<span class="va">`fit`</span><span class="op">`$`</span><span class="va">`stability`</span><span class="op">`$`</span><span class="va">`per_program`</span>\
<span class="co">`#>   program_id cell_type n_hit n_eligible recovery_freq conditional_recovery_freq`</span>\
<span class="co">`#> 1        B_1         B     4          4             1                         1`</span>\
<span class="co">`#> 2        B_2         B     4          4             1                         1`</span>\
<span class="co">`#> 3      CD4_1       CD4     4          4             1                         1`</span>\
<span class="co">`#> 4      CD8_1       CD8     4          4             1                         1`</span>\
<span class="co">`#> 5       NK_1        NK     4          4             1                         1`</span>\
<span class="co">`#> 6       NK_2        NK     4          4             1                         1`</span>\
<span class="co">`#> 7     Mono_1      Mono     4          4             1                         1`</span>\
<span class="co">`#> 8     Mono_2      Mono     4          4             1                         1`</span>\
<span class="co">`#>   median_score_cor median_loading_cos sign_consistency top_gene_jaccard`</span>\
<span class="co">`#> 1        0.9999741          0.9977389                1        0.9523810`</span>\
<span class="co">`#> 2        0.9997932          0.9917401                1        0.9047619`</span>\
<span class="co">`#> 3        0.9999745          0.9969061                1        0.9047619`</span>\
<span class="co">`#> 4        0.9999727          0.9973919                1        1.0000000`</span>\
<span class="co">`#> 5        0.9999817          0.9990403                1        1.0000000`</span>\
<span class="co">`#> 6        0.9999712          0.9972391                1        1.0000000`</span>\
<span class="co">`#> 7        0.9999909          0.9989417                1        0.9047619`</span>\
<span class="co">`#> 8        0.9998692          0.9952278                1        0.9047619`</span>\
<span class="fu">[`plot_stability`](../reference/plot_stability.html)</span><span class="op">`(`</span><span class="va">`fit`</span><span class="op">`)`</span>

</div>

<img src="quickstart_files/figure-html/stability-1.png" class="r-plt" width="900" />

</div>

<div class="section level2">

## Coordinated multicellular programs (CMPs)<a href="#coordinated-multicellular-programs-cmps" class="anchor" aria-label="anchor"></a>

[`find_coordinated_programs()`](../reference/find_coordinated_programs.html) tests whether programs from *different* cell types co-vary across donors, with family-wise error control by permutation. Gene overlap plays no role in the call.

<div id="cb10" class="sourceCode">

\
<span class="va">`cmp`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`find_coordinated_programs`](../reference/find_coordinated_programs.html)</span><span class="op">`(`</span><span class="va">`fit`</span>`, n_perm `<span class="op">`=`</span>` `<span class="fl">`199`</span>`, seed `<span class="op">`=`</span>` `<span class="fl">`1`</span><span class="op">`)`</span>\
<span class="va">`cmp`</span>\
<span class="co">`#> <cmp_result> 1 coordinated multicellular program(s) at FWER 0.05 (mode: marginal, ambient: drop)`</span>\
<span class="co">`#>   CMP 1: B_1, CD4_1 (-), NK_2 (-), Mono_1 (-)`</span>

</div>

The coordination network (node size = pooled evidence A_v, edge width = \|T\|, color = CMP membership) and the per-program node statistics:

<div id="cb11" class="sourceCode">

\
<span class="fu">[`plot_cmp_network`](../reference/plot_cmp_network.html)</span><span class="op">`(`</span><span class="va">`cmp`</span><span class="op">`)`</span>

</div>

<img src="quickstart_files/figure-html/viz-cmp-1.png" class="r-plt" width="1050" />

<div id="cb12" class="sourceCode">

\
<span class="fu">[`plot_cmp_nodes`](../reference/plot_cmp_nodes.html)</span><span class="op">`(`</span><span class="va">`cmp`</span><span class="op">`)`</span>

</div>

<img src="quickstart_files/figure-html/viz-cmp-2.png" class="r-plt" width="1050" />

</div>

<div class="section level2">

## Shared vs specific programs<a href="#shared-vs-specific-programs" class="anchor" aria-label="anchor"></a>

[`program_sharing()`](../reference/program_sharing.html) goes the other way: it matches programs across cell types by loading correlation (permutation-calibrated), then asks with ‘mashr’ whether the matched programs share their gene *effects* — correlated loadings alone do not imply shared effects.

<div id="cb13" class="sourceCode">

\
<span class="kw">`if`</span>` `<span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/ns-load.html" class="external-link"><code class="sourceCode r R">requireNamespace</code></a></span><span class="op">`(`</span><span class="st">`"mashr"`</span>`, quietly `<span class="op">`=`</span>` `<span class="cn">`TRUE`</span><span class="op">`)`</span><span class="op">`)`</span>` `<span class="op">`{`</span>\
`  `<span class="va">`shr`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`program_sharing`](../reference/program_sharing.html)</span><span class="op">`(`</span><span class="va">`fit`</span>`, n_perm `<span class="op">`=`</span>` `<span class="fl">`99`</span>`, seed `<span class="op">`=`</span>` `<span class="fl">`1`</span><span class="op">`)`</span>\
`  `<span class="va">`shr`</span>\
`  `<span class="fu">[`plot_sharing`](../reference/plot_sharing.html)</span><span class="op">`(`</span><span class="va">`shr`</span><span class="op">`)`</span>`   `<span class="co">`# pre-screen correlations x mashr classification`</span>\
<span class="op">`}`</span>

</div>

<img src="quickstart_files/figure-html/sharing-1.png" class="r-plt" width="900" />

</div>

<div class="section level2">

## Artifact genes<a href="#artifact-genes" class="anchor" aria-label="anchor"></a>

Immunoglobulin/TCR V(D)J segments, known donor-genotype markers, and sex-chromosome genes can dominate programs for technical or genetic rather than biological reasons. Flag them in existing fits, or exclude them at fit time (exclusion happens before variable-gene selection):

<div id="cb14" class="sourceCode">

\
<span class="fu">[`plot_artifact_flags`](../reference/plot_artifact_flags.html)</span><span class="op">`(`</span><span class="fu">[`flag_artifact_genes`](../reference/flag_artifact_genes.html)</span><span class="op">`(`</span><span class="va">`fit`</span><span class="op">`)`</span><span class="op">`)`</span>

</div>

<img src="quickstart_files/figure-html/viz-artifacts-1.png" class="r-plt" width="1050" />

<div id="cb15" class="sourceCode">

\
<span class="fu">[`flag_artifact_genes`](../reference/flag_artifact_genes.html)</span><span class="op">`(`</span><span class="va">`fit`</span><span class="op">`)`</span>\
<span class="co">`#>   cell_type program_id vdj_frac vdj_flag genotype_marker_frac`</span>\
<span class="co">`#> 1         B        B_1        0    FALSE                    0`</span>\
<span class="co">`#> 2         B        B_2        0    FALSE                    0`</span>\
<span class="co">`#> 3       CD4      CD4_1        0    FALSE                    0`</span>\
<span class="co">`#> 4       CD8      CD8_1        0    FALSE                    0`</span>\
<span class="co">`#> 5        NK       NK_1        0    FALSE                    0`</span>\
<span class="co">`#> 6        NK       NK_2        0    FALSE                    0`</span>\
<span class="co">`#> 7      Mono     Mono_1        0    FALSE                    0`</span>\
<span class="co">`#> 8      Mono     Mono_2        0    FALSE                    0`</span>\
<span class="co">`#>   genotype_marker_flag`</span>\
<span class="co">`#> 1                FALSE`</span>\
<span class="co">`#> 2                FALSE`</span>\
<span class="co">`#> 3                FALSE`</span>\
<span class="co">`#> 4                FALSE`</span>\
<span class="co">`#> 5                FALSE`</span>\
<span class="co">`#> 6                FALSE`</span>\
<span class="co">`#> 7                FALSE`</span>\
<span class="co">`#> 8                FALSE`</span>\
\
<span class="co">`# fit-time masking:`</span>\
<span class="co">`# excl <- artifact_genes(colnames(sim$pseudobulk$B))`</span>\
<span class="co">`# fit_masked <- fit_celltype_programs(x, exclude_genes = excl, ...)`</span>

</div>

</div>

<div class="section level2">

## Next steps<a href="#next-steps" class="anchor" aria-label="anchor"></a>

- The real-data case-study article (pkgdown site) walks through the Stephenson and COMBAT COVID cohorts, including the artifact-masking sensitivity analysis.
- `benchmarks/` contains the pre-registered simulation benchmarks behind the default settings and the CMP caller.

</div>
