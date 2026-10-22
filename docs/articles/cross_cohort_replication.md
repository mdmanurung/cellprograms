<div class="page-header">

# Cross-cohort replication: a preregistered-style case study

<div class="d-none name">

`cross_cohort_replication.Rmd`

</div>

</div>

This article is an end-to-end, fully evaluated template for asking: **does a program-level biological effect found in one cohort replicate in independent cohorts?** It uses a simulated three-cohort fixture so the truth is known, and it orchestrates existing exports only — no workflow-specific wrapper. The two transfer pathways are kept distinct:

- **Fixed-basis pathway** — the discovery fit’s learned transform (gene panel, centering/scaling, loading effects) is *applied* to the query cohort via [`project_programs_to_cells()`](../reference/project_programs_to_cells.html). [`freeze_preprocessing()`](../reference/freeze_preprocessing.html) alone is **not** this: it freezes decisions (gene panel, design), not the learned transform.
- **De-novo pathway** — the query cohort is fit independently (here with the frozen feature panel, so the gene universe matches) and programs are aligned afterwards with [`match_programs()`](../reference/match_programs.html). Factor numbers and raw scores are never compared directly across fits.

Technical transportability (coverage, reconstruction) is reported separately from biological effect replication, and every abstention is counted with its denominator.

<div id="cb1" class="sourceCode">

\
<span class="kw"><a href="https://rdrr.io/r/base/library.html" class="external-link"><code class="sourceCode r R">library</code></a></span><span class="op">`(`</span><span class="va">`cellprograms`</span><span class="op">`)`</span>\
<span class="fu"><a href="https://rdrr.io/r/base/Random.html" class="external-link"><code class="sourceCode r R">set.seed</code></a></span><span class="op">`(`</span><span class="fl">`20261019`</span><span class="op">`)`</span>

</div>

<div class="section level2">

## Three cohorts, known truth<a href="#three-cohorts-known-truth" class="anchor" aria-label="anchor"></a>

Two cell types (`T`, `B`), 24 donors per cohort, 400 genes. Planted structure in T cells:

| Program | Loadings | Activity | Truth |
|----|----|----|----|
| P1 | genes 1-40, identical across cohorts | correlated with donor outcome `y` in every cohort | conserved biological signal |
| P2 | genes 61-100, identical across cohorts | driven by a batch covariate in cohort B only | nuisance-only effect |
| P3 | genes 121-160 | active in cohort C only | cohort-specific program |

B cells carry one conserved program (genes 201-240) with no outcome link.

<div id="cb2" class="sourceCode">

\
<span class="va">`n_genes`</span>` `<span class="op">`<-`</span>` `<span class="fl">`400`</span>`; `<span class="va">`n_d`</span>` `<span class="op">`<-`</span>` `<span class="fl">`24`</span>\
<span class="va">`genes`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/paste.html" class="external-link"><code class="sourceCode r R">paste0</code></a></span><span class="op">`(`</span><span class="st">`"G"`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/seq.html" class="external-link"><code class="sourceCode r R">seq_len</code></a></span><span class="op">`(`</span><span class="va">`n_genes`</span><span class="op">`)`</span><span class="op">`)`</span>\
<span class="va">`mk_w`</span>` `<span class="op">`<-`</span>` `<span class="kw">`function`</span><span class="op">`(`</span><span class="va">`idx`</span><span class="op">`)`</span>` `<span class="op">`{`</span>` `<span class="va">`w`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/rep.html" class="external-link"><code class="sourceCode r R">rep</code></a></span><span class="op">`(`</span><span class="fl">`0`</span>`, `<span class="va">`n_genes`</span><span class="op">`)`</span>`; `<span class="va">`w`</span><span class="op">`[`</span><span class="va">`idx`</span><span class="op">`]`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/stats/Normal.html" class="external-link"><code class="sourceCode r R">rnorm</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/length.html" class="external-link"><code class="sourceCode r R">length</code></a></span><span class="op">`(`</span><span class="va">`idx`</span><span class="op">`)`</span><span class="op">`)`</span>`; `<span class="va">`w`</span>` `<span class="op">`}`</span>\
<span class="fu"><a href="https://rdrr.io/r/base/Random.html" class="external-link"><code class="sourceCode r R">set.seed</code></a></span><span class="op">`(`</span><span class="fl">`42`</span><span class="op">`)`</span>\
<span class="va">`w1`</span>` `<span class="op">`<-`</span>` `<span class="fu">`mk_w`</span><span class="op">`(`</span><span class="fl">`1`</span><span class="op">`:`</span><span class="fl">`40`</span><span class="op">`)`</span>`; `<span class="va">`w2`</span>` `<span class="op">`<-`</span>` `<span class="fu">`mk_w`</span><span class="op">`(`</span><span class="fl">`61`</span><span class="op">`:`</span><span class="fl">`100`</span><span class="op">`)`</span>`; `<span class="va">`w3`</span>` `<span class="op">`<-`</span>` `<span class="fu">`mk_w`</span><span class="op">`(`</span><span class="fl">`121`</span><span class="op">`:`</span><span class="fl">`160`</span><span class="op">`)`</span>`; `<span class="va">`w4`</span>` `<span class="op">`<-`</span>` `<span class="fu">`mk_w`</span><span class="op">`(`</span><span class="fl">`201`</span><span class="op">`:`</span><span class="fl">`240`</span><span class="op">`)`</span>\
\
<span class="va">`mk_cohort`</span>` `<span class="op">`<-`</span>` `<span class="kw">`function`</span><span class="op">`(`</span><span class="va">`name`</span>`, `<span class="va">`with_p3`</span>` `<span class="op">`=`</span>` `<span class="cn">`FALSE`</span>`, `<span class="va">`p2_batch`</span>` `<span class="op">`=`</span>` `<span class="cn">`FALSE`</span><span class="op">`)`</span>` `<span class="op">`{`</span>\
`  `<span class="va">`donors`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/paste.html" class="external-link"><code class="sourceCode r R">paste0</code></a></span><span class="op">`(`</span><span class="va">`name`</span>`, `<span class="st">`"_d"`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/seq.html" class="external-link"><code class="sourceCode r R">seq_len</code></a></span><span class="op">`(`</span><span class="va">`n_d`</span><span class="op">`)`</span><span class="op">`)`</span>\
`  `<span class="va">`z1`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/stats/Normal.html" class="external-link"><code class="sourceCode r R">rnorm</code></a></span><span class="op">`(`</span><span class="va">`n_d`</span><span class="op">`)`</span>`; `<span class="va">`z2`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/stats/Normal.html" class="external-link"><code class="sourceCode r R">rnorm</code></a></span><span class="op">`(`</span><span class="va">`n_d`</span><span class="op">`)`</span>`; `<span class="va">`z4`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/stats/Normal.html" class="external-link"><code class="sourceCode r R">rnorm</code></a></span><span class="op">`(`</span><span class="va">`n_d`</span><span class="op">`)`</span>\
`  `<span class="va">`batch`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/sample.html" class="external-link"><code class="sourceCode r R">sample</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/rep.html" class="external-link"><code class="sourceCode r R">rep</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link"><code class="sourceCode r R">c</code></a></span><span class="op">`(`</span><span class="st">`"b1"`</span>`, `<span class="st">`"b2"`</span><span class="op">`)`</span>`, each `<span class="op">`=`</span>` `<span class="va">`n_d`</span>` `<span class="op">`/`</span>` `<span class="fl">`2`</span><span class="op">`)`</span><span class="op">`)`</span>\
`  `<span class="kw">`if`</span>` `<span class="op">`(`</span><span class="va">`p2_batch`</span><span class="op">`)`</span>` `<span class="va">`z2`</span>` `<span class="op">`<-`</span>` `<span class="fl">`1.2`</span>` `<span class="op">`*`</span>` `<span class="op">`(`</span><span class="va">`batch`</span>` `<span class="op">`==`</span>` `<span class="st">`"b2"`</span><span class="op">`)`</span>` `<span class="op">`+`</span>` `<span class="fu"><a href="https://rdrr.io/r/stats/Normal.html" class="external-link"><code class="sourceCode r R">rnorm</code></a></span><span class="op">`(`</span><span class="va">`n_d`</span>`, sd `<span class="op">`=`</span>` `<span class="fl">`0.3`</span><span class="op">`)`</span>\
`  `<span class="va">`y`</span>` `<span class="op">`<-`</span>` `<span class="fl">`0.7`</span>` `<span class="op">`*`</span>` `<span class="va">`z1`</span>` `<span class="op">`+`</span>` `<span class="fu"><a href="https://rdrr.io/r/stats/Normal.html" class="external-link"><code class="sourceCode r R">rnorm</code></a></span><span class="op">`(`</span><span class="va">`n_d`</span>`, sd `<span class="op">`=`</span>` `<span class="fl">`0.8`</span><span class="op">`)`</span>`          `<span class="co">`# outcome linked to P1`</span>\
`  `<span class="va">`YT`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/crossprod.html" class="external-link"><code class="sourceCode r R">tcrossprod</code></a></span><span class="op">`(`</span><span class="va">`z1`</span>`, `<span class="va">`w1`</span><span class="op">`)`</span>` `<span class="op">`+`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/crossprod.html" class="external-link"><code class="sourceCode r R">tcrossprod</code></a></span><span class="op">`(`</span><span class="va">`z2`</span>`, `<span class="va">`w2`</span><span class="op">`)`</span>\
`  `<span class="kw">`if`</span>` `<span class="op">`(`</span><span class="va">`with_p3`</span><span class="op">`)`</span>` `<span class="va">`YT`</span>` `<span class="op">`<-`</span>` `<span class="va">`YT`</span>` `<span class="op">`+`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/crossprod.html" class="external-link"><code class="sourceCode r R">tcrossprod</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/stats/Normal.html" class="external-link"><code class="sourceCode r R">rnorm</code></a></span><span class="op">`(`</span><span class="va">`n_d`</span>`, sd `<span class="op">`=`</span>` `<span class="fl">`1.5`</span><span class="op">`)`</span>`, `<span class="va">`w3`</span><span class="op">`)`</span>\
`  `<span class="va">`YT`</span>` `<span class="op">`<-`</span>` `<span class="va">`YT`</span>` `<span class="op">`+`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/matrix.html" class="external-link"><code class="sourceCode r R">matrix</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/stats/Normal.html" class="external-link"><code class="sourceCode r R">rnorm</code></a></span><span class="op">`(`</span><span class="va">`n_d`</span>` `<span class="op">`*`</span>` `<span class="va">`n_genes`</span>`, sd `<span class="op">`=`</span>` `<span class="fl">`0.8`</span><span class="op">`)`</span>`, `<span class="va">`n_d`</span>`, `<span class="va">`n_genes`</span><span class="op">`)`</span>\
`  `<span class="va">`YB`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/crossprod.html" class="external-link"><code class="sourceCode r R">tcrossprod</code></a></span><span class="op">`(`</span><span class="va">`z4`</span>`, `<span class="va">`w4`</span><span class="op">`)`</span>` `<span class="op">`+`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/matrix.html" class="external-link"><code class="sourceCode r R">matrix</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/stats/Normal.html" class="external-link"><code class="sourceCode r R">rnorm</code></a></span><span class="op">`(`</span><span class="va">`n_d`</span>` `<span class="op">`*`</span>` `<span class="va">`n_genes`</span>`, sd `<span class="op">`=`</span>` `<span class="fl">`0.8`</span><span class="op">`)`</span>`, `<span class="va">`n_d`</span>`, `<span class="va">`n_genes`</span><span class="op">`)`</span>\
`  `<span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">colnames</code></a></span><span class="op">`(`</span><span class="va">`YT`</span><span class="op">`)`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">colnames</code></a></span><span class="op">`(`</span><span class="va">`YB`</span><span class="op">`)`</span>` `<span class="op">`<-`</span>` `<span class="va">`genes`</span>\
`  `<span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">rownames</code></a></span><span class="op">`(`</span><span class="va">`YT`</span><span class="op">`)`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">rownames</code></a></span><span class="op">`(`</span><span class="va">`YB`</span><span class="op">`)`</span>` `<span class="op">`<-`</span>` `<span class="va">`donors`</span>\
`  `<span class="va">`meta`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/data.frame.html" class="external-link"><code class="sourceCode r R">data.frame</code></a></span><span class="op">`(`</span>`observation_id `<span class="op">`=`</span>` `<span class="va">`donors`</span>`, cohort `<span class="op">`=`</span>` `<span class="va">`name`</span>`, outcome `<span class="op">`=`</span>` `<span class="va">`y`</span>`,`\
`                     batch `<span class="op">`=`</span>` `<span class="va">`batch`</span>`, stringsAsFactors `<span class="op">`=`</span>` `<span class="cn">`FALSE`</span><span class="op">`)`</span>\
`  `<span class="fu"><a href="https://rdrr.io/r/base/list.html" class="external-link"><code class="sourceCode r R">list</code></a></span><span class="op">`(`</span>`data `<span class="op">`=`</span>` `<span class="fu">[`as_cell_program_data`](../reference/as_cell_program_data.html)</span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/list.html" class="external-link"><code class="sourceCode r R">list</code></a></span><span class="op">`(`</span>`T `<span class="op">`=`</span>` `<span class="va">`YT`</span>`, B `<span class="op">`=`</span>` `<span class="va">`YB`</span><span class="op">`)`</span>`, sample_metadata `<span class="op">`=`</span>` `<span class="va">`meta`</span><span class="op">`)`</span>`,`\
`       meta `<span class="op">`=`</span>` `<span class="va">`meta`</span><span class="op">`)`</span>\
<span class="op">`}`</span>\
<span class="va">`cohA`</span>` `<span class="op">`<-`</span>` `<span class="fu">`mk_cohort`</span><span class="op">`(`</span><span class="st">`"A"`</span><span class="op">`)`</span>`; `<span class="va">`cohB`</span>` `<span class="op">`<-`</span>` `<span class="fu">`mk_cohort`</span><span class="op">`(`</span><span class="st">`"B"`</span>`, p2_batch `<span class="op">`=`</span>` `<span class="cn">`TRUE`</span><span class="op">`)`</span>\
<span class="va">`cohC`</span>` `<span class="op">`<-`</span>` `<span class="fu">`mk_cohort`</span><span class="op">`(`</span><span class="st">`"C"`</span>`, with_p3 `<span class="op">`=`</span>` `<span class="cn">`TRUE`</span><span class="op">`)`</span>

</div>

The predeclared hypothesis, fixed before looking at any query cohort: **the T-cell program whose activity tracks `outcome` in the discovery cohort replicates with the same sign as the discovery effect at p \< 0.05.** (A fitted program’s sign is arbitrary; canonicalization fixes it within a fit, so the rule is same-sign-as-discovery, not a fixed direction.) Everything else is descriptive.

</div>

<div class="section level2">

## Discovery fit (cohort A)<a href="#discovery-fit-cohort-a" class="anchor" aria-label="anchor"></a>

<div id="cb3" class="sourceCode">

\
<span class="va">`fitA`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`fit_celltype_programs`](../reference/fit_celltype_programs.html)</span><span class="op">`(`</span><span class="va">`cohA`</span><span class="op">`$`</span><span class="va">`data`</span>`, max_factors `<span class="op">`=`</span>` `<span class="fl">`4`</span>`, seed `<span class="op">`=`</span>` `<span class="fl">`1`</span>`,`\
`                              feature_scale `<span class="op">`=`</span>` `<span class="st">`"simulated"`</span><span class="op">`)`</span>\
<span class="va">`fitA`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`canonicalize_programs`](../reference/canonicalize_programs.html)</span><span class="op">`(`</span><span class="va">`fitA`</span><span class="op">`)`</span>\
<span class="fu">[`program_summary`](../reference/program_summary.html)</span><span class="op">`(`</span><span class="va">`fitA`</span><span class="op">`)`</span><span class="op">`[`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link"><code class="sourceCode r R">c</code></a></span><span class="op">`(`</span><span class="st">`"cell_type"`</span>`, `<span class="st">`"n_programs"`</span>`, `<span class="st">`"R2"`</span><span class="op">`)`</span><span class="op">`]`</span>\
<span class="co">`#>   cell_type n_programs         R2`</span>\
<span class="co">`#> 1         T          2 0.26394663`</span>\
<span class="co">`#> 2         B          1 0.07886538`</span>

</div>

Freeze the discovery *decisions* (feature panel, centering) for later de-novo refits, and keep the fitted object itself for fixed-basis projection:

<div id="cb4" class="sourceCode">

\
<span class="va">`frozen`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`freeze_preprocessing`](../reference/freeze_preprocessing.html)</span><span class="op">`(`</span><span class="va">`cohA`</span><span class="op">`$`</span><span class="va">`data`</span>`, features `<span class="op">`=`</span>` `<span class="st">`"all"`</span>`, center `<span class="op">`=`</span>` `<span class="cn">`TRUE`</span><span class="op">`)`</span>\
<span class="va">`frozen`</span><span class="op">`$`</span><span class="va">`digest`</span>\
<span class="co">`#> [1] "a92378907d03ffd308f4d967f75bf62bc71c8173e89c91233e00df622462a720"`</span>

</div>

</div>

<div class="section level2">

## Pathway 1: fixed-basis projection into cohorts B and C<a href="#pathway-1-fixed-basis-projection-into-cohorts-b-and-c" class="anchor" aria-label="anchor"></a>

The discovery programs are applied to each query cohort’s donor pseudobulks. Rows are `type_donor` so both cell types project in one call; the discovery preprocessing is replayed from the fit, not refit.

<div id="cb5" class="sourceCode">

\
<span class="va">`stack`</span>` `<span class="op">`<-`</span>` `<span class="kw">`function`</span><span class="op">`(`</span><span class="va">`coh`</span><span class="op">`)`</span>` `<span class="op">`{`</span>\
`  `<span class="fu"><a href="https://rdrr.io/r/base/cbind.html" class="external-link"><code class="sourceCode r R">rbind</code></a></span><span class="op">`(`</span><span class="va">`coh`</span><span class="op">`$`</span><span class="va">`data`</span><span class="op">`$`</span><span class="va">`pseudobulk`</span><span class="op">`$`</span><span class="cn">`T`</span>`, `<span class="va">`coh`</span><span class="op">`$`</span><span class="va">`data`</span><span class="op">`$`</span><span class="va">`pseudobulk`</span><span class="op">`$`</span><span class="va">`B`</span><span class="op">`)`</span>` `<span class="op">`|>`</span>\
`    `<span class="fu">`` `rownames<-` ``</span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link"><code class="sourceCode r R">c</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/paste.html" class="external-link"><code class="sourceCode r R">paste0</code></a></span><span class="op">`(`</span><span class="st">`"T_"`</span>`, `<span class="va">`coh`</span><span class="op">`$`</span><span class="va">`meta`</span><span class="op">`$`</span><span class="va">`observation_id`</span><span class="op">`)`</span>`,`\
`                   `<span class="fu"><a href="https://rdrr.io/r/base/paste.html" class="external-link"><code class="sourceCode r R">paste0</code></a></span><span class="op">`(`</span><span class="st">`"B_"`</span>`, `<span class="va">`coh`</span><span class="op">`$`</span><span class="va">`meta`</span><span class="op">`$`</span><span class="va">`observation_id`</span><span class="op">`)`</span><span class="op">`)`</span><span class="op">`)`</span>\
<span class="op">`}`</span>\
<span class="va">`project_cohort`</span>` `<span class="op">`<-`</span>` `<span class="kw">`function`</span><span class="op">`(`</span><span class="va">`fit`</span>`, `<span class="va">`coh`</span><span class="op">`)`</span>` `<span class="op">`{`</span>\
`  `<span class="va">`expr`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/utils/stack.html" class="external-link"><code class="sourceCode r R">stack</code></a></span><span class="op">`(`</span><span class="va">`coh`</span><span class="op">`)`</span>\
`  `<span class="va">`cts`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/stats/setNames.html" class="external-link"><code class="sourceCode r R">setNames</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/grep.html" class="external-link"><code class="sourceCode r R">sub</code></a></span><span class="op">`(`</span><span class="st">`"_.*"`</span>`, `<span class="st">`""`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">rownames</code></a></span><span class="op">`(`</span><span class="va">`expr`</span><span class="op">`)`</span><span class="op">`)`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">rownames</code></a></span><span class="op">`(`</span><span class="va">`expr`</span><span class="op">`)`</span><span class="op">`)`</span>\
`  `<span class="fu">[`project_programs_to_cells`](../reference/project_programs_to_cells.html)</span><span class="op">`(`</span><span class="va">`fit`</span>`, `<span class="va">`expr`</span>`, cell_types `<span class="op">`=`</span>` `<span class="va">`cts`</span>`,`\
`                            feature_scale `<span class="op">`=`</span>` `<span class="st">`"simulated"`</span>`, min_coverage `<span class="op">`=`</span>` `<span class="fl">`0.9`</span><span class="op">`)`</span>\
<span class="op">`}`</span>\
<span class="va">`projB`</span>` `<span class="op">`<-`</span>` `<span class="fu">`project_cohort`</span><span class="op">`(`</span><span class="va">`fitA`</span>`, `<span class="va">`cohB`</span><span class="op">`)`</span>\
<span class="va">`projC`</span>` `<span class="op">`<-`</span>` `<span class="fu">`project_cohort`</span><span class="op">`(`</span><span class="va">`fitA`</span>`, `<span class="va">`cohC`</span><span class="op">`)`</span>\
<span class="va">`projB`</span><span class="op">`$`</span><span class="va">`coverage`</span><span class="op">`[`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link"><code class="sourceCode r R">c</code></a></span><span class="op">`(`</span><span class="st">`"cell_type"`</span>`, `<span class="st">`"programs"`</span>`, `<span class="st">`"cells"`</span>`, `<span class="st">`"shared_genes"`</span>`, `<span class="st">`"coverage"`</span><span class="op">`)`</span><span class="op">`]`</span>\
<span class="co">`#>   cell_type programs cells shared_genes coverage`</span>\
<span class="co">`#> 1         T        2    24          400        1`</span>\
<span class="co">`#> 2         B        1    24          400        1`</span>

</div>

Coverage is complete (400/400 genes) — technical transportability holds by construction here; on real data this table is where ineligible type-cohorts are caught. Now the predeclared association on the projected T-cell activities:

<div id="cb6" class="sourceCode">

\
<span class="va">`assoc_projected`</span>` `<span class="op">`<-`</span>` `<span class="kw">`function`</span><span class="op">`(`</span><span class="va">`proj`</span>`, `<span class="va">`meta`</span>`, `<span class="va">`ct`</span>` `<span class="op">`=`</span>` `<span class="st">`"T"`</span><span class="op">`)`</span>` `<span class="op">`{`</span>\
`  `<span class="va">`A`</span>` `<span class="op">`<-`</span>` `<span class="va">`proj`</span><span class="op">`$`</span><span class="va">`activity`</span>\
`  `<span class="va">`rows`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">rownames</code></a></span><span class="op">`(`</span><span class="va">`A`</span><span class="op">`)`</span><span class="op">`[`</span><span class="fu"><a href="https://rdrr.io/r/base/startsWith.html" class="external-link"><code class="sourceCode r R">startsWith</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">rownames</code></a></span><span class="op">`(`</span><span class="va">`A`</span><span class="op">`)`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/paste.html" class="external-link"><code class="sourceCode r R">paste0</code></a></span><span class="op">`(`</span><span class="va">`ct`</span>`, `<span class="st">`"_"`</span><span class="op">`)`</span><span class="op">`)`</span><span class="op">`]`</span>\
`  `<span class="va">`out`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/lapply.html" class="external-link"><code class="sourceCode r R">vapply</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">colnames</code></a></span><span class="op">`(`</span><span class="va">`A`</span><span class="op">`)`</span>`, `<span class="kw">`function`</span><span class="op">`(`</span><span class="va">`p`</span><span class="op">`)`</span>` `<span class="op">`{`</span>\
`    `<span class="va">`s`</span>` `<span class="op">`<-`</span>` `<span class="va">`A`</span><span class="op">`[`</span><span class="va">`rows`</span>`, `<span class="va">`p`</span><span class="op">`]`</span>\
`    `<span class="kw">`if`</span>` `<span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/all.html" class="external-link"><code class="sourceCode r R">all</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/NA.html" class="external-link"><code class="sourceCode r R">is.na</code></a></span><span class="op">`(`</span><span class="va">`s`</span><span class="op">`)`</span><span class="op">`)`</span><span class="op">`)`</span>` `<span class="kw"><a href="https://rdrr.io/r/base/function.html" class="external-link"><code class="sourceCode r R">return</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link"><code class="sourceCode r R">c</code></a></span><span class="op">`(`</span><span class="cn">`NA_real_`</span>`, `<span class="cn">`NA_real_`</span><span class="op">`)`</span><span class="op">`)`</span>\
`    `<span class="va">`y`</span>` `<span class="op">`<-`</span>` `<span class="va">`meta`</span><span class="op">`$`</span><span class="va">`outcome`</span><span class="op">`[`</span><span class="fu"><a href="https://rdrr.io/r/base/match.html" class="external-link"><code class="sourceCode r R">match</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/grep.html" class="external-link"><code class="sourceCode r R">sub</code></a></span><span class="op">`(`</span><span class="st">`"^T_"`</span>`, `<span class="st">`""`</span>`, `<span class="va">`rows`</span><span class="op">`)`</span>`, `<span class="va">`meta`</span><span class="op">`$`</span><span class="va">`observation_id`</span><span class="op">`)`</span><span class="op">`]`</span>\
`    `<span class="va">`cf`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/stats/coef.html" class="external-link"><code class="sourceCode r R">coef</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/summary.html" class="external-link"><code class="sourceCode r R">summary</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/stats/lm.html" class="external-link"><code class="sourceCode r R">lm</code></a></span><span class="op">`(`</span><span class="va">`y`</span>` `<span class="op">`~`</span>` `<span class="va">`s`</span><span class="op">`)`</span><span class="op">`)`</span><span class="op">`)`</span>\
`    `<span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link"><code class="sourceCode r R">c</code></a></span><span class="op">`(`</span><span class="va">`cf`</span><span class="op">`[`</span><span class="st">`"s"`</span>`, `<span class="fl">`1`</span><span class="op">`]`</span>`, `<span class="va">`cf`</span><span class="op">`[`</span><span class="st">`"s"`</span>`, `<span class="fl">`4`</span><span class="op">`]`</span><span class="op">`)`</span>\
`  `<span class="op">`}`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/numeric.html" class="external-link"><code class="sourceCode r R">numeric</code></a></span><span class="op">`(`</span><span class="fl">`2`</span><span class="op">`)`</span><span class="op">`)`</span>\
`  `<span class="fu"><a href="https://rdrr.io/r/base/data.frame.html" class="external-link"><code class="sourceCode r R">data.frame</code></a></span><span class="op">`(`</span>`program `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">colnames</code></a></span><span class="op">`(`</span><span class="va">`A`</span><span class="op">`)`</span>`, beta `<span class="op">`=`</span>` `<span class="va">`out`</span><span class="op">`[`</span><span class="fl">`1`</span>`, `<span class="op">`]`</span>`, p `<span class="op">`=`</span>` `<span class="va">`out`</span><span class="op">`[`</span><span class="fl">`2`</span>`, `<span class="op">`]`</span>`,`\
`             row.names `<span class="op">`=`</span>` `<span class="cn">`NULL`</span><span class="op">`)`</span>\
<span class="op">`}`</span>\
<span class="va">`aB`</span>` `<span class="op">`<-`</span>` `<span class="fu">`assoc_projected`</span><span class="op">`(`</span><span class="va">`projB`</span>`, `<span class="va">`cohB`</span><span class="op">`$`</span><span class="va">`meta`</span><span class="op">`)`</span>\
<span class="va">`aC`</span>` `<span class="op">`<-`</span>` `<span class="fu">`assoc_projected`</span><span class="op">`(`</span><span class="va">`projC`</span>`, `<span class="va">`cohC`</span><span class="op">`$`</span><span class="va">`meta`</span><span class="op">`)`</span>\
<span class="fu"><a href="https://rdrr.io/r/base/cbind.html" class="external-link"><code class="sourceCode r R">rbind</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/data.frame.html" class="external-link"><code class="sourceCode r R">data.frame</code></a></span><span class="op">`(`</span>`cohort `<span class="op">`=`</span>` `<span class="st">`"B"`</span>`, `<span class="va">`aB`</span><span class="op">`)`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/data.frame.html" class="external-link"><code class="sourceCode r R">data.frame</code></a></span><span class="op">`(`</span>`cohort `<span class="op">`=`</span>` `<span class="st">`"C"`</span>`, `<span class="va">`aC`</span><span class="op">`)`</span><span class="op">`)`</span>\
<span class="co">`#>   cohort program       beta           p`</span>\
<span class="co">`#> 1      B     T_1 -2.8684965 0.006144099`</span>\
<span class="co">`#> 2      B     T_2  0.5826891 0.739876914`</span>\
<span class="co">`#> 3      B     B_1         NA          NA`</span>\
<span class="co">`#> 4      C     T_1 -3.1295870 0.003218153`</span>\
<span class="co">`#> 5      C     T_2 -0.3099099 0.815471556`</span>\
<span class="co">`#> 6      C     B_1         NA          NA`</span>

</div>

One T-cell program shows the replicated outcome association in both query cohorts (same sign as the discovery effect, p \< 0.05) — the conserved P1 signal. Which discovery program is it, and is it the planted one?

<div id="cb7" class="sourceCode">

\
<span class="va">`sc`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`program_scores`](../reference/program_scores.html)</span><span class="op">`(`</span><span class="va">`fitA`</span>`, format `<span class="op">`=`</span>` `<span class="st">`"list"`</span><span class="op">`)`</span><span class="op">`$`</span><span class="cn">`T`</span>\
<span class="va">`disc_assoc`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/lapply.html" class="external-link"><code class="sourceCode r R">vapply</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/seq.html" class="external-link"><code class="sourceCode r R">seq_len</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/nrow.html" class="external-link"><code class="sourceCode r R">ncol</code></a></span><span class="op">`(`</span><span class="va">`sc`</span><span class="op">`)`</span><span class="op">`)`</span>`, `<span class="kw">`function`</span><span class="op">`(`</span><span class="va">`k`</span><span class="op">`)`</span>` `<span class="op">`{`</span>\
`  `<span class="va">`cf`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/stats/coef.html" class="external-link"><code class="sourceCode r R">coef</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/summary.html" class="external-link"><code class="sourceCode r R">summary</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/stats/lm.html" class="external-link"><code class="sourceCode r R">lm</code></a></span><span class="op">`(`</span><span class="va">`cohA`</span><span class="op">`$`</span><span class="va">`meta`</span><span class="op">`$`</span><span class="va">`outcome`</span>` `<span class="op">`~`</span>` `<span class="va">`sc`</span><span class="op">`[`</span>`, `<span class="va">`k`</span><span class="op">`]`</span><span class="op">`)`</span><span class="op">`)`</span><span class="op">`)`</span><span class="op">`[`</span><span class="fl">`2`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link"><code class="sourceCode r R">c</code></a></span><span class="op">`(`</span><span class="fl">`1`</span>`, `<span class="fl">`4`</span><span class="op">`)`</span><span class="op">`]`</span>\
`  `<span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link"><code class="sourceCode r R">c</code></a></span><span class="op">`(`</span><span class="va">`cf`</span><span class="op">`[`</span><span class="fl">`1`</span><span class="op">`]`</span>`, `<span class="va">`cf`</span><span class="op">`[`</span><span class="fl">`2`</span><span class="op">`]`</span><span class="op">`)`</span>\
<span class="op">`}`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/numeric.html" class="external-link"><code class="sourceCode r R">numeric</code></a></span><span class="op">`(`</span><span class="fl">`2`</span><span class="op">`)`</span><span class="op">`)`</span>\
<span class="fu"><a href="https://rdrr.io/r/base/data.frame.html" class="external-link"><code class="sourceCode r R">data.frame</code></a></span><span class="op">`(`</span>`program `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">colnames</code></a></span><span class="op">`(`</span><span class="va">`sc`</span><span class="op">`)`</span>`, beta `<span class="op">`=`</span>` `<span class="va">`disc_assoc`</span><span class="op">`[`</span><span class="fl">`1`</span>`, `<span class="op">`]`</span>`, p `<span class="op">`=`</span>` `<span class="va">`disc_assoc`</span><span class="op">`[`</span><span class="fl">`2`</span>`, `<span class="op">`]`</span><span class="op">`)`</span>\
<span class="co">`#>   program        beta            p`</span>\
<span class="co">`#> 1     T_1 -0.77061789 2.472144e-05`</span>\
<span class="co">`#> 2     T_2 -0.01467328 9.472151e-01`</span>

</div>

The same discovery program carries the association in A and replicates in B and C through the fixed basis — score-level transfer plus effect replication, with the sign fixed by canonicalization within the projected frame.

</div>

<div class="section level2">

## Pathway 2: de-novo refit + matching<a href="#pathway-2-de-novo-refit-matching" class="anchor" aria-label="anchor"></a>

Cohort B is refit independently — same frozen feature panel, independent factors — then aligned to the discovery programs by loading-space matching:

<div id="cb8" class="sourceCode">

\
<span class="va">`fitB`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`fit_celltype_programs`](../reference/fit_celltype_programs.html)</span><span class="op">`(`</span><span class="va">`cohB`</span><span class="op">`$`</span><span class="va">`data`</span>`, max_factors `<span class="op">`=`</span>` `<span class="fl">`4`</span>`, seed `<span class="op">`=`</span>` `<span class="fl">`2`</span>`,`\
`                              preprocessing_spec `<span class="op">`=`</span>` `<span class="va">`frozen`</span><span class="op">`)`</span>\
<span class="va">`fitB`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`canonicalize_programs`](../reference/canonicalize_programs.html)</span><span class="op">`(`</span><span class="va">`fitB`</span><span class="op">`)`</span>\
<span class="va">`m`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`match_programs`](../reference/match_programs.html)</span><span class="op">`(`</span><span class="va">`fitA`</span>`, `<span class="va">`fitB`</span>`, `<span class="st">`"T"`</span>`, method `<span class="op">`=`</span>` `<span class="st">`"correlation"`</span>`,`\
`                    space `<span class="op">`=`</span>` `<span class="st">`"loadings"`</span>`, sim_threshold `<span class="op">`=`</span>` `<span class="fl">`0.5`</span><span class="op">`)`</span>\
<span class="va">`m`</span><span class="op">`$`</span><span class="va">`matches`</span>\
<span class="co">`#>   program_a program_b    cosine`</span>\
<span class="co">`#> 1       T_1       T_1 0.9641598`</span>\
<span class="co">`#> 2       T_2       T_2 0.8973300`</span>

</div>

The matched B program’s outcome association is computed on **B’s own scores** (never on discovery coordinates):

<div id="cb9" class="sourceCode">

\
<span class="va">`scB`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`program_scores`](../reference/program_scores.html)</span><span class="op">`(`</span><span class="va">`fitB`</span>`, format `<span class="op">`=`</span>` `<span class="st">`"list"`</span><span class="op">`)`</span><span class="op">`$`</span><span class="cn">`T`</span>\
<span class="va">`hit`</span>` `<span class="op">`<-`</span>` `<span class="va">`m`</span><span class="op">`$`</span><span class="va">`matches`</span><span class="op">`[`</span><span class="va">`m`</span><span class="op">`$`</span><span class="va">`matches`</span><span class="op">`$`</span><span class="va">`cosine`</span>` `<span class="op">`>`</span>` `<span class="fl">`0.5`</span>`, , drop `<span class="op">`=`</span>` `<span class="cn">`FALSE`</span><span class="op">`]`</span>\
<span class="kw">`if`</span>` `<span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/nrow.html" class="external-link"><code class="sourceCode r R">nrow</code></a></span><span class="op">`(`</span><span class="va">`hit`</span><span class="op">`)`</span><span class="op">`)`</span>` `<span class="op">`{`</span>\
`  `<span class="va">`k`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/match.html" class="external-link"><code class="sourceCode r R">match</code></a></span><span class="op">`(`</span><span class="va">`hit`</span><span class="op">`$`</span><span class="va">`program_b`</span><span class="op">`[`</span><span class="fl">`1`</span><span class="op">`]`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">colnames</code></a></span><span class="op">`(`</span><span class="va">`scB`</span><span class="op">`)`</span><span class="op">`)`</span>\
`  `<span class="fu"><a href="https://rdrr.io/r/stats/coef.html" class="external-link"><code class="sourceCode r R">coef</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/summary.html" class="external-link"><code class="sourceCode r R">summary</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/stats/lm.html" class="external-link"><code class="sourceCode r R">lm</code></a></span><span class="op">`(`</span><span class="va">`cohB`</span><span class="op">`$`</span><span class="va">`meta`</span><span class="op">`$`</span><span class="va">`outcome`</span>` `<span class="op">`~`</span>` `<span class="va">`scB`</span><span class="op">`[`</span>`, `<span class="va">`k`</span><span class="op">`]`</span><span class="op">`)`</span><span class="op">`)`</span><span class="op">`)`</span><span class="op">`[`</span><span class="fl">`2`</span>`, `<span class="op">`]`</span>\
<span class="op">`}`</span>\
<span class="co">`#>     Estimate   Std. Error      t value     Pr(>|t|) `</span>\
<span class="co">`#> -0.506379271  0.167712860 -3.019322857  0.006305575`</span>

</div>

</div>

<div class="section level2">

## Abstention accounting<a href="#abstention-accounting" class="anchor" aria-label="anchor"></a>

Every program enters exactly one bucket per pathway; denominators are always shown.

<div id="cb10" class="sourceCode">

\
<span class="va">`disc_sign`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/sign.html" class="external-link"><code class="sourceCode r R">sign</code></a></span><span class="op">`(`</span><span class="va">`disc_assoc`</span><span class="op">`[`</span><span class="fl">`1`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/which.min.html" class="external-link"><code class="sourceCode r R">which.min</code></a></span><span class="op">`(`</span><span class="va">`disc_assoc`</span><span class="op">`[`</span><span class="fl">`2`</span>`, `<span class="op">`]`</span><span class="op">`)`</span><span class="op">`]`</span><span class="op">`)`</span>`  `<span class="co">`# discovery effect direction`</span>\
<span class="va">`accounting`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/data.frame.html" class="external-link"><code class="sourceCode r R">data.frame</code></a></span><span class="op">`(`</span>\
`  pathway `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link"><code class="sourceCode r R">c</code></a></span><span class="op">`(`</span><span class="st">`"fixed-basis (A->B)"`</span>`, `<span class="st">`"fixed-basis (A->C)"`</span>`, `<span class="st">`"de-novo (A vs B)"`</span><span class="op">`)`</span>`,`\
`  discovery_programs `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link"><code class="sourceCode r R">c</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/nrow.html" class="external-link"><code class="sourceCode r R">ncol</code></a></span><span class="op">`(`</span><span class="va">`fitA`</span><span class="op">`$`</span><span class="va">`scores`</span><span class="op">`$`</span><span class="cn">`T`</span><span class="op">`)`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/nrow.html" class="external-link"><code class="sourceCode r R">ncol</code></a></span><span class="op">`(`</span><span class="va">`fitA`</span><span class="op">`$`</span><span class="va">`scores`</span><span class="op">`$`</span><span class="cn">`T`</span><span class="op">`)`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/nrow.html" class="external-link"><code class="sourceCode r R">ncol</code></a></span><span class="op">`(`</span><span class="va">`fitA`</span><span class="op">`$`</span><span class="va">`scores`</span><span class="op">`$`</span><span class="cn">`T`</span><span class="op">`)`</span><span class="op">`)`</span>`,`\
`  attempted `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link"><code class="sourceCode r R">c</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/sum.html" class="external-link"><code class="sourceCode r R">sum</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/startsWith.html" class="external-link"><code class="sourceCode r R">startsWith</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">colnames</code></a></span><span class="op">`(`</span><span class="va">`projB`</span><span class="op">`$`</span><span class="va">`activity`</span><span class="op">`)`</span>`, `<span class="st">`"T_"`</span><span class="op">`)`</span><span class="op">`)`</span>`,`\
`                `<span class="fu"><a href="https://rdrr.io/r/base/sum.html" class="external-link"><code class="sourceCode r R">sum</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/startsWith.html" class="external-link"><code class="sourceCode r R">startsWith</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">colnames</code></a></span><span class="op">`(`</span><span class="va">`projC`</span><span class="op">`$`</span><span class="va">`activity`</span><span class="op">`)`</span>`, `<span class="st">`"T_"`</span><span class="op">`)`</span><span class="op">`)`</span>`,`\
`                `<span class="fu"><a href="https://rdrr.io/r/base/nrow.html" class="external-link"><code class="sourceCode r R">nrow</code></a></span><span class="op">`(`</span><span class="va">`m`</span><span class="op">`$`</span><span class="va">`matches`</span><span class="op">`)`</span><span class="op">`)`</span>`,`\
`  ineligible `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link"><code class="sourceCode r R">c</code></a></span><span class="op">`(`</span><span class="fl">`0L`</span>`, `<span class="fl">`0L`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/nrow.html" class="external-link"><code class="sourceCode r R">ncol</code></a></span><span class="op">`(`</span><span class="va">`fitA`</span><span class="op">`$`</span><span class="va">`scores`</span><span class="op">`$`</span><span class="cn">`T`</span><span class="op">`)`</span>` `<span class="op">`-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/nrow.html" class="external-link"><code class="sourceCode r R">nrow</code></a></span><span class="op">`(`</span><span class="va">`m`</span><span class="op">`$`</span><span class="va">`matches`</span><span class="op">`)`</span><span class="op">`)`</span>`,`\
`  replicated `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link"><code class="sourceCode r R">c</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/sum.html" class="external-link"><code class="sourceCode r R">sum</code></a></span><span class="op">`(`</span><span class="va">`aB`</span><span class="op">`$`</span><span class="va">`p`</span>` `<span class="op">`<`</span>` `<span class="fl">`0.05`</span>` `<span class="op">`&`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/sign.html" class="external-link"><code class="sourceCode r R">sign</code></a></span><span class="op">`(`</span><span class="va">`aB`</span><span class="op">`$`</span><span class="va">`beta`</span><span class="op">`)`</span>` `<span class="op">`==`</span>` `<span class="va">`disc_sign`</span>`, na.rm `<span class="op">`=`</span>` `<span class="cn">`TRUE`</span><span class="op">`)`</span>`,`\
`                 `<span class="fu"><a href="https://rdrr.io/r/base/sum.html" class="external-link"><code class="sourceCode r R">sum</code></a></span><span class="op">`(`</span><span class="va">`aC`</span><span class="op">`$`</span><span class="va">`p`</span>` `<span class="op">`<`</span>` `<span class="fl">`0.05`</span>` `<span class="op">`&`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/sign.html" class="external-link"><code class="sourceCode r R">sign</code></a></span><span class="op">`(`</span><span class="va">`aC`</span><span class="op">`$`</span><span class="va">`beta`</span><span class="op">`)`</span>` `<span class="op">`==`</span>` `<span class="va">`disc_sign`</span>`, na.rm `<span class="op">`=`</span>` `<span class="cn">`TRUE`</span><span class="op">`)`</span>`,`\
`                 `<span class="cn">`NA_integer_`</span><span class="op">`)`</span>`,`\
`  stringsAsFactors `<span class="op">`=`</span>` `<span class="cn">`FALSE`</span><span class="op">`)`</span>\
<span class="va">`accounting`</span>\
<span class="co">`#>              pathway discovery_programs attempted ineligible replicated`</span>\
<span class="co">`#> 1 fixed-basis (A->B)                  2         2          0          1`</span>\
<span class="co">`#> 2 fixed-basis (A->C)                  2         2          0          1`</span>\
<span class="co">`#> 3   de-novo (A vs B)                  2         2          0         NA`</span>

</div>

And the cohort-specific truth: P3 exists only in cohort C. When C is the discovery cohort, its extra program must **abstain** in A and B — no projection target, no match. That is the correct outcome, not a failure:

<div id="cb11" class="sourceCode">

\
<span class="va">`fitC`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`fit_celltype_programs`](../reference/fit_celltype_programs.html)</span><span class="op">`(`</span><span class="va">`cohC`</span><span class="op">`$`</span><span class="va">`data`</span>`, max_factors `<span class="op">`=`</span>` `<span class="fl">`5`</span>`, seed `<span class="op">`=`</span>` `<span class="fl">`3`</span>`,`\
`                              feature_scale `<span class="op">`=`</span>` `<span class="st">`"simulated"`</span><span class="op">`)`</span>\
<span class="va">`fitC`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`canonicalize_programs`](../reference/canonicalize_programs.html)</span><span class="op">`(`</span><span class="va">`fitC`</span><span class="op">`)`</span>\
<span class="va">`m_CA`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`match_programs`](../reference/match_programs.html)</span><span class="op">`(`</span><span class="va">`fitC`</span>`, `<span class="va">`fitA`</span>`, `<span class="st">`"T"`</span>`, method `<span class="op">`=`</span>` `<span class="st">`"correlation"`</span>`,`\
`                       space `<span class="op">`=`</span>` `<span class="st">`"loadings"`</span>`, sim_threshold `<span class="op">`=`</span>` `<span class="fl">`0.5`</span><span class="op">`)`</span>\
<span class="fu"><a href="https://rdrr.io/r/base/data.frame.html" class="external-link"><code class="sourceCode r R">data.frame</code></a></span><span class="op">`(`</span>`discovery_C_programs `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/nrow.html" class="external-link"><code class="sourceCode r R">ncol</code></a></span><span class="op">`(`</span><span class="va">`fitC`</span><span class="op">`$`</span><span class="va">`scores`</span><span class="op">`$`</span><span class="cn">`T`</span><span class="op">`)`</span>`,`\
`           matched_in_A `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/nrow.html" class="external-link"><code class="sourceCode r R">nrow</code></a></span><span class="op">`(`</span><span class="va">`m_CA`</span><span class="op">`$`</span><span class="va">`matches`</span><span class="op">`)`</span>`,`\
`           abstained `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/nrow.html" class="external-link"><code class="sourceCode r R">ncol</code></a></span><span class="op">`(`</span><span class="va">`fitC`</span><span class="op">`$`</span><span class="va">`scores`</span><span class="op">`$`</span><span class="cn">`T`</span><span class="op">`)`</span>` `<span class="op">`-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/nrow.html" class="external-link"><code class="sourceCode r R">nrow</code></a></span><span class="op">`(`</span><span class="va">`m_CA`</span><span class="op">`$`</span><span class="va">`matches`</span><span class="op">`)`</span><span class="op">`)`</span>\
<span class="co">`#>   discovery_C_programs matched_in_A abstained`</span>\
<span class="co">`#> 1                    3            2         1`</span>

</div>

C discovers more T programs than A (it carries P3); the unmatched remainder is the abstention count. The leave-one-cohort-out gate is: conserved P1 replicates in every held-out cohort (shown above for B and C with A as discovery; the rotations B-\>A,C and C-\>A,B use the same code), and cohort-specific programs abstain rather than force a match.

</div>

<div class="section level2">

## Negative controls<a href="#negative-controls" class="anchor" aria-label="anchor"></a>

**(1) Permuted outcome.** Replication must collapse when the outcome is permuted within cohort:

<div id="cb12" class="sourceCode">

\
<span class="fu"><a href="https://rdrr.io/r/base/Random.html" class="external-link"><code class="sourceCode r R">set.seed</code></a></span><span class="op">`(`</span><span class="fl">`99`</span><span class="op">`)`</span>\
<span class="va">`perm_hits`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/lapply.html" class="external-link"><code class="sourceCode r R">vapply</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/seq.html" class="external-link"><code class="sourceCode r R">seq_len</code></a></span><span class="op">`(`</span><span class="fl">`200`</span><span class="op">`)`</span>`, `<span class="kw">`function`</span><span class="op">`(`</span><span class="va">`r`</span><span class="op">`)`</span>` `<span class="op">`{`</span>\
`  `<span class="va">`meta_p`</span>` `<span class="op">`<-`</span>` `<span class="va">`cohB`</span><span class="op">`$`</span><span class="va">`meta`</span>\
`  `<span class="va">`meta_p`</span><span class="op">`$`</span><span class="va">`outcome`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/sample.html" class="external-link"><code class="sourceCode r R">sample</code></a></span><span class="op">`(`</span><span class="va">`meta_p`</span><span class="op">`$`</span><span class="va">`outcome`</span><span class="op">`)`</span>`   `<span class="co">`# shuffle values against fixed donor IDs`</span>\
`  `<span class="va">`aB_perm`</span>` `<span class="op">`<-`</span>` `<span class="fu">`assoc_projected`</span><span class="op">`(`</span><span class="va">`projB`</span>`, `<span class="va">`meta_p`</span><span class="op">`)`</span>\
`  `<span class="fu"><a href="https://rdrr.io/r/base/sum.html" class="external-link"><code class="sourceCode r R">sum</code></a></span><span class="op">`(`</span><span class="va">`aB_perm`</span><span class="op">`$`</span><span class="va">`p`</span>` `<span class="op">`<`</span>` `<span class="fl">`0.05`</span>` `<span class="op">`&`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/sign.html" class="external-link"><code class="sourceCode r R">sign</code></a></span><span class="op">`(`</span><span class="va">`aB_perm`</span><span class="op">`$`</span><span class="va">`beta`</span><span class="op">`)`</span>` `<span class="op">`==`</span>` `<span class="va">`disc_sign`</span>`, na.rm `<span class="op">`=`</span>` `<span class="cn">`TRUE`</span><span class="op">`)`</span>\
<span class="op">`}`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/numeric.html" class="external-link"><code class="sourceCode r R">numeric</code></a></span><span class="op">`(`</span><span class="fl">`1`</span><span class="op">`)`</span><span class="op">`)`</span>\
<span class="fu"><a href="https://rdrr.io/r/base/mean.html" class="external-link"><code class="sourceCode r R">mean</code></a></span><span class="op">`(`</span><span class="va">`perm_hits`</span><span class="op">`)`</span>`  `<span class="co">`# expectation under the null: ~0.05 per non-degenerate program`</span>\
<span class="co">`#> [1] 0.07`</span>

</div>

**(2) Scrambled gene IDs.** A query whose feature IDs do not match the reference must abstain loudly, not silently project:

<div id="cb13" class="sourceCode">

\
<span class="va">`bad`</span>` `<span class="op">`<-`</span>` `<span class="va">`cohB`</span>\
<span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">colnames</code></a></span><span class="op">`(`</span><span class="va">`bad`</span><span class="op">`$`</span><span class="va">`data`</span><span class="op">`$`</span><span class="va">`pseudobulk`</span><span class="op">`$`</span><span class="cn">`T`</span><span class="op">`)`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/paste.html" class="external-link"><code class="sourceCode r R">paste0</code></a></span><span class="op">`(`</span><span class="st">`"X"`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/seq.html" class="external-link"><code class="sourceCode r R">seq_len</code></a></span><span class="op">`(`</span><span class="va">`n_genes`</span><span class="op">`)`</span><span class="op">`)`</span>\
<span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">colnames</code></a></span><span class="op">`(`</span><span class="va">`bad`</span><span class="op">`$`</span><span class="va">`data`</span><span class="op">`$`</span><span class="va">`pseudobulk`</span><span class="op">`$`</span><span class="va">`B`</span><span class="op">`)`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/paste.html" class="external-link"><code class="sourceCode r R">paste0</code></a></span><span class="op">`(`</span><span class="st">`"X"`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/seq.html" class="external-link"><code class="sourceCode r R">seq_len</code></a></span><span class="op">`(`</span><span class="va">`n_genes`</span><span class="op">`)`</span><span class="op">`)`</span>\
<span class="kw"><a href="https://rdrr.io/r/base/conditions.html" class="external-link"><code class="sourceCode r R">tryCatch</code></a></span><span class="op">`(`</span><span class="fu">`project_cohort`</span><span class="op">`(`</span><span class="va">`fitA`</span>`, `<span class="va">`bad`</span><span class="op">`)`</span>`,`\
`         error `<span class="op">`=`</span>` `<span class="kw">`function`</span><span class="op">`(`</span><span class="va">`e`</span><span class="op">`)`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/conditions.html" class="external-link"><code class="sourceCode r R">conditionMessage</code></a></span><span class="op">`(`</span><span class="va">`e`</span><span class="op">`)`</span><span class="op">`)`</span>\
<span class="co">`` #> [1] "\033[1m\033[22mCell type \033[34m\"T\"\033[39m: none of the fit's 400 genes are present in `expression`." ``</span>

</div>

**(3) Nuisance-only effect.** P2 tracks batch in cohort B by construction. A technical association must not present as replicated biology: testing the batch term in cohort C (where no batch effect was planted) finds nothing.

<div id="cb14" class="sourceCode">

\
<span class="va">`sB`</span>` `<span class="op">`<-`</span>` `<span class="va">`projB`</span><span class="op">`$`</span><span class="va">`activity`</span><span class="op">`[`</span><span class="fu"><a href="https://rdrr.io/r/base/paste.html" class="external-link"><code class="sourceCode r R">paste0</code></a></span><span class="op">`(`</span><span class="st">`"T_"`</span>`, `<span class="va">`cohB`</span><span class="op">`$`</span><span class="va">`meta`</span><span class="op">`$`</span><span class="va">`observation_id`</span><span class="op">`)`</span>`, `<span class="op">`]`</span>\
<span class="va">`sC`</span>` `<span class="op">`<-`</span>` `<span class="va">`projC`</span><span class="op">`$`</span><span class="va">`activity`</span><span class="op">`[`</span><span class="fu"><a href="https://rdrr.io/r/base/paste.html" class="external-link"><code class="sourceCode r R">paste0</code></a></span><span class="op">`(`</span><span class="st">`"T_"`</span>`, `<span class="va">`cohC`</span><span class="op">`$`</span><span class="va">`meta`</span><span class="op">`$`</span><span class="va">`observation_id`</span><span class="op">`)`</span>`, `<span class="op">`]`</span>\
<span class="va">`batch_assoc`</span>` `<span class="op">`<-`</span>` `<span class="kw">`function`</span><span class="op">`(`</span><span class="va">`S`</span>`, `<span class="va">`meta`</span><span class="op">`)`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/lapply.html" class="external-link"><code class="sourceCode r R">vapply</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/seq.html" class="external-link"><code class="sourceCode r R">seq_len</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/nrow.html" class="external-link"><code class="sourceCode r R">ncol</code></a></span><span class="op">`(`</span><span class="va">`S`</span><span class="op">`)`</span><span class="op">`)`</span>`, `<span class="kw">`function`</span><span class="op">`(`</span><span class="va">`k`</span><span class="op">`)`</span>` `<span class="op">`{`</span>\
`  `<span class="va">`s`</span>` `<span class="op">`<-`</span>` `<span class="va">`S`</span><span class="op">`[`</span>`, `<span class="va">`k`</span><span class="op">`]`</span>\
`  `<span class="kw">`if`</span>` `<span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/all.html" class="external-link"><code class="sourceCode r R">all</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/NA.html" class="external-link"><code class="sourceCode r R">is.na</code></a></span><span class="op">`(`</span><span class="va">`s`</span><span class="op">`)`</span><span class="op">`)`</span><span class="op">`)`</span>` `<span class="kw"><a href="https://rdrr.io/r/base/function.html" class="external-link"><code class="sourceCode r R">return</code></a></span><span class="op">`(`</span><span class="cn">`NA_real_`</span><span class="op">`)`</span>`   `<span class="co">`# off-cell-type program column`</span>\
`  `<span class="va">`cf`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/stats/coef.html" class="external-link"><code class="sourceCode r R">coef</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/summary.html" class="external-link"><code class="sourceCode r R">summary</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/stats/lm.html" class="external-link"><code class="sourceCode r R">lm</code></a></span><span class="op">`(`</span><span class="va">`s`</span>` `<span class="op">`~`</span>` `<span class="va">`meta`</span><span class="op">`$`</span><span class="va">`batch`</span><span class="op">`)`</span><span class="op">`)`</span><span class="op">`)`</span>\
`  `<span class="kw">`if`</span>` `<span class="op">`(`</span><span class="st">`"meta$batchb2"`</span>` `<span class="op"><a href="https://rdrr.io/r/base/match.html" class="external-link"><code class="sourceCode r R"><span class="sc">%in%</span></code></a></span>` `<span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">rownames</code></a></span><span class="op">`(`</span><span class="va">`cf`</span><span class="op">`)`</span><span class="op">`)`</span>` `<span class="va">`cf`</span><span class="op">`[`</span><span class="st">`"meta$batchb2"`</span>`, `<span class="fl">`4`</span><span class="op">`]`</span>` `<span class="kw">`else`</span>` `<span class="cn">`NA_real_`</span>\
<span class="op">`}`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/numeric.html" class="external-link"><code class="sourceCode r R">numeric</code></a></span><span class="op">`(`</span><span class="fl">`1`</span><span class="op">`)`</span><span class="op">`)`</span>\
<span class="fu"><a href="https://rdrr.io/r/base/data.frame.html" class="external-link"><code class="sourceCode r R">data.frame</code></a></span><span class="op">`(`</span>`program `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">colnames</code></a></span><span class="op">`(`</span><span class="va">`sB`</span><span class="op">`)`</span>`,`\
`           p_batch_B `<span class="op">`=`</span>` `<span class="fu">`batch_assoc`</span><span class="op">`(`</span><span class="va">`sB`</span>`, `<span class="va">`cohB`</span><span class="op">`$`</span><span class="va">`meta`</span><span class="op">`)`</span>`,`\
`           p_batch_C `<span class="op">`=`</span>` `<span class="fu">`batch_assoc`</span><span class="op">`(`</span><span class="va">`sC`</span>`, `<span class="va">`cohC`</span><span class="op">`$`</span><span class="va">`meta`</span><span class="op">`)`</span><span class="op">`)`</span>\
<span class="co">`#>   program    p_batch_B p_batch_C`</span>\
<span class="co">`#> 1     T_1 8.860339e-01 0.5783604`</span>\
<span class="co">`#> 2     T_2 7.579566e-09 0.4969985`</span>\
<span class="co">`#> 3     B_1           NA        NA`</span>

</div>

The batch-loaded program is significant in B and null in C — flagged as cohort-specific technical structure, excluded from the replication claim.

</div>

<div class="section level2">

## What to report, and limitations<a href="#what-to-report-and-limitations" class="anchor" aria-label="anchor"></a>

A replication report from this workflow should state: cohorts and programs attempted/projected/matched/ineligible/nonreplicated (the accounting table); the predeclared association and replication rule; coverage and reconstruction per transfer (technical transportability); and the negative-control outcomes. Limitations: projected activities are query-set-dependent empirical-Bayes estimates (see the projection documentation); sign and scale are comparable only within a pathway, not across de-novo fits; and with 24 donors per cohort the power calculus here is illustrative, not a recommendation.

</div>
