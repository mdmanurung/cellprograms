<div class="page-header">

# Joint and tensor factorization: EB-CP and block EBMF

<div class="d-none name">

`coordination_joint.Rmd`

</div>

</div>

The default `cellprograms` pipeline fits each cell type separately and only then asks whether programs are shared. Two experimental (v0.2.0.9001) extensions invert the ordering and factorize across cell types directly:

- **EB-CP** ([`fit_eb_cp()`](../reference/fit_eb_cp.html) / [`eb_cp_from_fit()`](../reference/eb_cp_from_fit.html)): a tensor-style greedy EBMF over the stacked cell-type x donor matrices. Each rank-1 component is a gene loading (`U`), a donor activity vector (`V`) shared across the cell types the component spans, and a cell-type loading vector (`W`). Components are accepted by a held-out z gate on a donor-synchronous split, so the number of components is selected by out-of-sample predictive evidence, not by a variance cutoff.
- **Block EBMF** ([`fit_joint_block()`](../reference/fit_joint_block.html) / [`block_clusters()`](../reference/block_clusters.html)): a shared-first flashier pass over the stacked expression matrix, followed by per-cell-type private passes on the residuals. [`block_clusters()`](../reference/block_clusters.html) collapses shared factors into gene-set clusters by loading overlap and returns donor activity per cluster.

This vignette runs both on a small simulated dataset and states the scope caveats measured in the coordination benchmark (see `benchmarks/coordination/` in the repository).

<div id="cb1" class="sourceCode">

\
<span class="kw"><a href="https://rdrr.io/r/base/library.html" class="external-link"><code class="sourceCode r R">library</code></a></span><span class="op">`(`</span><span class="va">`cellprograms`</span><span class="op">`)`</span>\
\
<span class="va">`sim`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`simulate_pseudobulk`](../reference/simulate_pseudobulk.html)</span><span class="op">`(`</span>`scenario `<span class="op">`=`</span>` `<span class="st">`"mixed"`</span>`, n_samples `<span class="op">`=`</span>` `<span class="fl">`60`</span>`,`\
`                           n_genes `<span class="op">`=`</span>` `<span class="fl">`800`</span>`, program_size `<span class="op">`=`</span>` `<span class="fl">`40`</span>`, seed `<span class="op">`=`</span>` `<span class="fl">`1`</span><span class="op">`)`</span>\
<span class="va">`x`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`as_cell_program_data`](../reference/as_cell_program_data.html)</span><span class="op">`(`</span><span class="va">`sim`</span><span class="op">`$`</span><span class="va">`pseudobulk`</span><span class="op">`)`</span>\
<span class="va">`fit`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`fit_celltype_programs`](../reference/fit_celltype_programs.html)</span><span class="op">`(`</span><span class="va">`x`</span>`, loading_prior `<span class="op">`=`</span>` `<span class="st">`"point_laplace"`</span>`,`\
`                             max_factors `<span class="op">`=`</span>` `<span class="fl">`10L`</span>`, seed `<span class="op">`=`</span>` `<span class="fl">`1`</span><span class="op">`)`</span>\
<span class="va">`fit`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`canonicalize_programs`](../reference/canonicalize_programs.html)</span><span class="op">`(`</span><span class="va">`fit`</span><span class="op">`)`</span>

</div>

<div class="section level2">

## EB-CP: greedy components with held-out selection<a href="#eb-cp-greedy-components-with-held-out-selection" class="anchor" aria-label="anchor"></a>

[`eb_cp_from_fit()`](../reference/eb_cp_from_fit.html) runs EB-CP on the adjusted per-cell-type matrices of a fitted `cell_program_fit`. If `strata` is not supplied, the donor-synchronous train/held-out split is derived from the covariates used at fit time, so the split respects the adjustment; pass `strata` explicitly when you want a different stratification.

<div id="cb2" class="sourceCode">

\
<span class="va">`cp`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`eb_cp_from_fit`](../reference/eb_cp_from_fit.html)</span><span class="op">`(`</span><span class="va">`fit`</span><span class="op">`)`</span>\
<span class="va">`cp`</span><span class="op">`$`</span><span class="va">`n_components`</span>\
<span class="co">`#> [1] 8`</span>\
<span class="va">`cp`</span><span class="op">`$`</span><span class="va">`kept_heldout_z`</span>`          `<span class="co">`# held-out z of each kept component`</span>\
<span class="co">`#> [1] 141.06089 154.97144 103.76617 107.39676  94.15278  63.45945  93.65694`</span>\
<span class="co">`#> [8]  51.93119`</span>\
<span class="va">`cp`</span><span class="op">`$`</span><span class="va">`trace_z`</span>`                 `<span class="co">`# z trace including the first failed component`</span>\
<span class="co">`#> [1] 141.06089 154.97144 103.76617 107.39676  94.15278  63.45945  93.65694`</span>\
<span class="co">`#> [8]  51.93119`</span>

</div>

[`eb_cp_programs()`](../reference/eb_cp_programs.html) turns components into program calls. Components spanning \>= `min_members` cell types at the membership threshold are flagged as coordination programs; the rest are kept as private programs, flagged rather than dropped.

<div id="cb3" class="sourceCode">

\
<span class="va">`pr`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`eb_cp_programs`](../reference/eb_cp_programs.html)</span><span class="op">`(`</span><span class="va">`cp`</span><span class="op">`)`</span>\
<span class="va">`pr`</span><span class="op">`$`</span><span class="va">`programs`</span>\
<span class="co">`#>   program cell_type n_members is_coordination heldout_z`</span>\
<span class="co">`#> 1    CP01      Mono         1           FALSE 141.06089`</span>\
<span class="co">`#> 2    CP02    CD8+NK         2            TRUE 154.97144`</span>\
<span class="co">`#> 3    CP03         B         1           FALSE 103.76617`</span>\
<span class="co">`#> 4    CP04       CD4         1           FALSE 107.39676`</span>\
<span class="co">`#> 5    CP05        NK         1           FALSE  94.15278`</span>\
<span class="co">`#> 6    CP06      Mono         1           FALSE  63.45945`</span>\
<span class="co">`#> 7    CP07    CD8+NK         2            TRUE  93.65694`</span>\
<span class="co">`#> 8    CP08         B         1           FALSE  51.93119`</span>\
<span class="fu"><a href="https://rdrr.io/r/utils/head.html" class="external-link"><code class="sourceCode r R">head</code></a></span><span class="op">`(`</span><span class="va">`pr`</span><span class="op">`$`</span><span class="va">`activity`</span><span class="op">`[`</span>`, `<span class="fl">`1`</span><span class="op">`:`</span><span class="fl">`2`</span>`, drop `<span class="op">`=`</span>` `<span class="cn">`FALSE`</span><span class="op">`]`</span><span class="op">`)`</span>`   `<span class="co">`# donor activity per coordination component`</span>\
<span class="co">`#>        cluster_02  cluster_07`</span>\
<span class="co">`#> obs1   2.53382698 -2.76698019`</span>\
<span class="co">`#> obs10  0.49417548 -0.54670036`</span>\
<span class="co">`#> obs11 -0.36878957  0.34052851`</span>\
<span class="co">`#> obs12 -1.12515663  1.19204711`</span>\
<span class="co">`#> obs13 -0.07259629  0.05518669`</span>\
<span class="co">`#> obs14 -0.71861307  0.58362302`</span>

</div>

</div>

<div class="section level2">

## Block EBMF: shared-first pass plus private residuals<a href="#block-ebmf-shared-first-pass-plus-private-residuals" class="anchor" aria-label="anchor"></a>

[`fit_joint_block()`](../reference/fit_joint_block.html) stacks the per-cell-type matrices and fits shared factors first, then per-cell-type private factors on the residuals.

<div id="cb4" class="sourceCode">

\
<span class="va">`bfit`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`fit_joint_block`](../reference/fit_joint_block.html)</span><span class="op">`(`</span><span class="va">`fit`</span>`, level `<span class="op">`=`</span>` `<span class="st">`"expression"`</span>`, order `<span class="op">`=`</span>` `<span class="st">`"shared_first"`</span>`,`\
`                        max_factors `<span class="op">`=`</span>` `<span class="fl">`10L`</span>`, seed `<span class="op">`=`</span>` `<span class="fl">`1`</span>`, verbose `<span class="op">`=`</span>` `<span class="fl">`0L`</span><span class="op">`)`</span>\
<span class="va">`bfit`</span><span class="op">`$`</span><span class="va">`shared_fit`</span><span class="op">`$`</span><span class="va">`n_factors`</span>\
<span class="co">`#> [1] 4`</span>\
<span class="va">`bc`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`block_clusters`](../reference/block_clusters.html)</span><span class="op">`(`</span><span class="va">`bfit`</span>`, mass_thresh `<span class="op">`=`</span>` `<span class="fl">`0.1`</span><span class="op">`)`</span>\
<span class="fu"><a href="https://rdrr.io/r/utils/head.html" class="external-link"><code class="sourceCode r R">head</code></a></span><span class="op">`(`</span><span class="va">`bc`</span><span class="op">`$`</span><span class="va">`clusters`</span>`, `<span class="fl">`10`</span><span class="op">`)`</span>\
<span class="co">`#>   cluster cell_type program_id      stat fwer_p`</span>\
<span class="co">`#> 1       1         B       SF01 0.3283760     NA`</span>\
<span class="co">`#> 2       1       CD4       SF01 0.3283760     NA`</span>\
<span class="co">`#> 3       1        NK       SF01 0.3283760     NA`</span>\
<span class="co">`#> 4       1      Mono       SF01 0.3283760     NA`</span>\
<span class="co">`#> 5       2       CD8       SF02 0.5278789     NA`</span>\
<span class="co">`#> 6       2        NK       SF02 0.5278789     NA`</span>

</div>

Each row of `bc$clusters` is a gene-set cluster (a group of shared factors with overlapping loadings); `bc$activity` holds donor activity per cluster. The block-mass rule (`mass_thresh`) keeps clusters whose membership spreads across at least two cell types.

</div>

<div class="section level2">

## Scope caveats from the coordination benchmark<a href="#scope-caveats-from-the-coordination-benchmark" class="anchor" aria-label="anchor"></a>

Both extensions trade robustness for directness, and the preregistered coordination benchmark (`benchmarks/coordination/`, tag `prereg-coord-v1`) quantifies where:

- **Block EBMF (`block_E`) is powerful but over-calls without calibration.** On null data the shared-first pass still returns factors whose loadings spread across blocks, and the block-mass rule alone called spurious shared programs in 16/16 null replicates. A permutation-calibrated significance layer (within-strata donor permutation, max-T family-wise control) is available in the benchmark harness and is being productized; until then, treat uncalibrated [`block_clusters()`](../reference/block_clusters.html) output as a candidate list, not as evidence of sharing.
- **EB-CP is conservative.** The held-out z gate keeps specificity high (near-zero false positives on null data) but loses power when shared programs have cell-type-specific gene loadings, because the model’s shared donor-activity assumption is only approximately met.
- **Low-K regime.** Both methods were validated with few recurring programs (1-3) against many private ones. With dozens of true per-cell-type programs, set `max_factors` above the default and expect the shared-first ordering to matter more.

For the default-by-construction path (per-cell-type fits, then [`sharing_spectrum()`](../reference/sharing_spectrum.html) / [`plot_sharing_spectrum()`](../reference/plot_sharing_spectrum.html)), see [`vignette("quickstart")`](../articles/quickstart.html).

</div>
