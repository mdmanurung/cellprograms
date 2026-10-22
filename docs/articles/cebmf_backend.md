<div class="page-header">

# The cebmf_torch backend: covariate-EBMF in torch

<div class="d-none name">

`cebmf_backend.Rmd`

</div>

</div>

[`fit_celltype_programs()`](../reference/fit_celltype_programs.html) ships two interchangeable EBMF backends:

|  | `backend = "flashier"` (default) | `backend = "cebmf_torch"` |
|----|----|----|
| Implementation | R (flashier.utils) | Python (torch, covariate-EBMF) |
| Failure mode | nullcheck fallback ladder | fails loudly, no fallback |
| ELBO | reported | `NA` (not produced; the trace is) |
| Convergence record | ELBO trace | `objective_trace` |

Everything around the factorization — preprocessing, canonicalization, program extraction, plots, sharing — is backend-agnostic. This vignette runs the torch backend on a small simulated count dataset and shows where the two contracts differ.

<div class="section level2">

## Python environment<a href="#python-environment" class="anchor" aria-label="anchor"></a>

The torch backend needs a Python environment with `torch`, `numpy`, `scipy`, `pandas`, and `cebmf_torch`. By default the package uses a **basilisk-managed environment** with pinned versions (Python 3.11, `torch==2.7.1` CPU build, `cebmf_torch` at a pinned commit), built automatically on first use and cached across sessions. To point at your own environment instead:

<div id="cb1" class="sourceCode">

\
<span class="fu"><a href="https://rdrr.io/r/base/options.html" class="external-link"><code class="sourceCode r R">options</code></a></span><span class="op">`(`</span>`cellprograms.cebmf_python_env `<span class="op">`=`</span>` `<span class="st">`"/path/to/venv"`</span><span class="op">`)`</span>`  `<span class="co">`# venv dir or python binary`</span>

</div>

</div>

<div class="section level2">

## Input regime: normalize counts first<a href="#input-regime-normalize-counts-first" class="anchor" aria-label="anchor"></a>

The torch backend estimates column-wise noise and **diverges immediately on raw unnormalized Gaussian-scale matrices** (the objective goes non-finite at iteration 0; the backend fails loudly rather than returning a garbage fit — this is by design and is the documented reason the test suite’s registered fixture normalizes counts). Feed it what the single-cell pipeline produces: counts, size-factor normalization, log1p, gene z-scaling.

<div id="cb2" class="sourceCode">

\
<span class="kw"><a href="https://rdrr.io/r/base/library.html" class="external-link"><code class="sourceCode r R">library</code></a></span><span class="op">`(`</span><span class="va">`cellprograms`</span><span class="op">`)`</span>\
<span class="fu"><a href="https://rdrr.io/r/base/Random.html" class="external-link"><code class="sourceCode r R">set.seed</code></a></span><span class="op">`(`</span><span class="fl">`3`</span><span class="op">`)`</span>\
<span class="va">`n`</span>` `<span class="op">`<-`</span>` `<span class="fl">`60`</span>`; `<span class="va">`g`</span>` `<span class="op">`<-`</span>` `<span class="fl">`400`</span>`; `<span class="va">`p`</span>` `<span class="op">`<-`</span>` `<span class="fl">`2`</span>\
\
<span class="co">`# block-structured counts: 2 planted programs over private gene sets`</span>\
<span class="va">`W`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/matrix.html" class="external-link"><code class="sourceCode r R">matrix</code></a></span><span class="op">`(`</span><span class="fl">`0`</span>`, `<span class="va">`g`</span>`, `<span class="va">`p`</span><span class="op">`)`</span>\
<span class="kw">`for`</span>` `<span class="op">`(`</span><span class="va">`j`</span>` `<span class="kw">`in`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/seq.html" class="external-link"><code class="sourceCode r R">seq_len</code></a></span><span class="op">`(`</span><span class="va">`p`</span><span class="op">`)`</span><span class="op">`)`</span>` `<span class="va">`W`</span><span class="op">`[`</span><span class="fu"><a href="https://rdrr.io/r/base/sample.html" class="external-link"><code class="sourceCode r R">sample.int</code></a></span><span class="op">`(`</span><span class="va">`g`</span>`, `<span class="fl">`60`</span><span class="op">`)`</span>`, `<span class="va">`j`</span><span class="op">`]`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/stats/Uniform.html" class="external-link"><code class="sourceCode r R">runif</code></a></span><span class="op">`(`</span><span class="fl">`60`</span>`, `<span class="fl">`0.5`</span>`, `<span class="fl">`1.5`</span><span class="op">`)`</span>\
<span class="va">`S`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/matrix.html" class="external-link"><code class="sourceCode r R">matrix</code></a></span><span class="op">`(`</span><span class="fl">`0`</span>`, `<span class="va">`n`</span>`, `<span class="va">`p`</span><span class="op">`)`</span>\
<span class="kw">`for`</span>` `<span class="op">`(`</span><span class="va">`j`</span>` `<span class="kw">`in`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/seq.html" class="external-link"><code class="sourceCode r R">seq_len</code></a></span><span class="op">`(`</span><span class="va">`p`</span><span class="op">`)`</span><span class="op">`)`</span>` `<span class="va">`S`</span><span class="op">`[`</span><span class="fu"><a href="https://rdrr.io/r/base/sample.html" class="external-link"><code class="sourceCode r R">sample.int</code></a></span><span class="op">`(`</span><span class="va">`n`</span>`, `<span class="fl">`20`</span><span class="op">`)`</span>`, `<span class="va">`j`</span><span class="op">`]`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/stats/Uniform.html" class="external-link"><code class="sourceCode r R">runif</code></a></span><span class="op">`(`</span><span class="fl">`20`</span>`, `<span class="fl">`1`</span>`, `<span class="fl">`2.5`</span><span class="op">`)`</span>\
<span class="va">`lam`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/Log.html" class="external-link"><code class="sourceCode r R">exp</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/Extremes.html" class="external-link"><code class="sourceCode r R">pmin</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/Extremes.html" class="external-link"><code class="sourceCode r R">pmax</code></a></span><span class="op">`(`</span><span class="va">`S`</span>` `<span class="op"><a href="https://rdrr.io/r/base/matmult.html" class="external-link"><code class="sourceCode r R"><span class="sc">%*%</span></code></a></span>` `<span class="fu"><a href="https://rdrr.io/r/base/t.html" class="external-link"><code class="sourceCode r R">t</code></a></span><span class="op">`(`</span><span class="va">`W`</span><span class="op">`)`</span>`, `<span class="fl">`0`</span><span class="op">`)`</span>`, `<span class="fl">`4`</span><span class="op">`)`</span><span class="op">`)`</span>` `<span class="op">`+`</span>` `<span class="fl">`0.05`</span>\
<span class="va">`counts`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/matrix.html" class="external-link"><code class="sourceCode r R">matrix</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/stats/Poisson.html" class="external-link"><code class="sourceCode r R">rpois</code></a></span><span class="op">`(`</span><span class="va">`n`</span>` `<span class="op">`*`</span>` `<span class="va">`g`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/vector.html" class="external-link"><code class="sourceCode r R">as.vector</code></a></span><span class="op">`(`</span><span class="va">`lam`</span><span class="op">`)`</span><span class="op">`)`</span>`, `<span class="va">`n`</span>`, `<span class="va">`g`</span><span class="op">`)`</span>\
\
<span class="co">`# size-factor log1p + gene z-scale: the backend's input regime`</span>\
<span class="va">`sf`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/colSums.html" class="external-link"><code class="sourceCode r R">rowSums</code></a></span><span class="op">`(`</span><span class="va">`counts`</span><span class="op">`)`</span>`; `<span class="va">`sf`</span>` `<span class="op">`<-`</span>` `<span class="va">`sf`</span>` `<span class="op">`/`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/mean.html" class="external-link"><code class="sourceCode r R">mean</code></a></span><span class="op">`(`</span><span class="va">`sf`</span><span class="op">`)`</span>\
<span class="va">`Y`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/scale.html" class="external-link"><code class="sourceCode r R">scale</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/Log.html" class="external-link"><code class="sourceCode r R">log1p</code></a></span><span class="op">`(`</span><span class="va">`counts`</span>` `<span class="op">`/`</span>` `<span class="va">`sf`</span><span class="op">`)`</span><span class="op">`)`</span>\
<span class="fu"><a href="https://rdrr.io/r/base/dimnames.html" class="external-link"><code class="sourceCode r R">dimnames</code></a></span><span class="op">`(`</span><span class="va">`Y`</span><span class="op">`)`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/list.html" class="external-link"><code class="sourceCode r R">list</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/paste.html" class="external-link"><code class="sourceCode r R">paste0</code></a></span><span class="op">`(`</span><span class="st">`"s"`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/seq.html" class="external-link"><code class="sourceCode r R">seq_len</code></a></span><span class="op">`(`</span><span class="va">`n`</span><span class="op">`)`</span><span class="op">`)`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/paste.html" class="external-link"><code class="sourceCode r R">paste0</code></a></span><span class="op">`(`</span><span class="st">`"g"`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/seq.html" class="external-link"><code class="sourceCode r R">seq_len</code></a></span><span class="op">`(`</span><span class="va">`g`</span><span class="op">`)`</span><span class="op">`)`</span><span class="op">`)`</span>\
\
<span class="va">`x`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`as_cell_program_data`](../reference/as_cell_program_data.html)</span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/list.html" class="external-link"><code class="sourceCode r R">list</code></a></span><span class="op">`(`</span>`CT1 `<span class="op">`=`</span>` `<span class="va">`Y`</span><span class="op">`)`</span><span class="op">`)`</span>\
<span class="fu">[`validate_cell_program_data`](../reference/validate_cell_program_data.html)</span><span class="op">`(`</span><span class="va">`x`</span><span class="op">`)`</span>

</div>

</div>

<div class="section level2">

## Fit with the torch backend<a href="#fit-with-the-torch-backend" class="anchor" aria-label="anchor"></a>

<div id="cb3" class="sourceCode">

\
<span class="va">`fit`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`fit_celltype_programs`](../reference/fit_celltype_programs.html)</span><span class="op">`(`</span><span class="va">`x`</span>`, backend `<span class="op">`=`</span>` `<span class="st">`"cebmf_torch"`</span>`,`\
`                             max_factors `<span class="op">`=`</span>` `<span class="fl">`6`</span>`, seed `<span class="op">`=`</span>` `<span class="fl">`1`</span><span class="op">`)`</span>\
<span class="va">`fl`</span>` `<span class="op">`<-`</span>` `<span class="va">`fit`</span><span class="op">`$`</span><span class="va">`fits`</span><span class="op">`$`</span><span class="va">`CT1`</span><span class="op">`$`</span><span class="va">`flash`</span>\
<span class="va">`fl`</span><span class="op">`$`</span><span class="va">`backend`</span>\
<span class="va">`fl`</span><span class="op">`$`</span><span class="va">`n_factors`</span>\
<span class="fu"><a href="https://rdrr.io/r/utils/head.html" class="external-link"><code class="sourceCode r R">head</code></a></span><span class="op">`(`</span><span class="va">`fl`</span><span class="op">`$`</span><span class="va">`objective_trace`</span><span class="op">`)`</span>`   `<span class="co">`# convergence record (the torch backend's ELBO analogue)`</span>

</div>

The registered mapping contract: the torch backend does not produce an ELBO, so [`program_summary()`](../reference/program_summary.html) reports `NA` there and nothing is manufactured to fill the gap; `objective_trace` is the honest convergence record. Loadings and scores arrive in the same shapes as the flashier backend:

<div id="cb4" class="sourceCode">

\
<span class="fu"><a href="https://rdrr.io/r/base/dim.html" class="external-link"><code class="sourceCode r R">dim</code></a></span><span class="op">`(`</span><span class="va">`fl`</span><span class="op">`$`</span><span class="va">`L_pm`</span><span class="op">`)`</span>`   `<span class="co">`# samples x factors`</span>\
<span class="fu"><a href="https://rdrr.io/r/base/dim.html" class="external-link"><code class="sourceCode r R">dim</code></a></span><span class="op">`(`</span><span class="va">`fl`</span><span class="op">`$`</span><span class="va">`F_pm`</span><span class="op">`)`</span>`   `<span class="co">`# genes x factors`</span>

</div>

</div>

<div class="section level2">

## Determinism at a fixed seed<a href="#determinism-at-a-fixed-seed" class="anchor" aria-label="anchor"></a>

<div id="cb5" class="sourceCode">

\
<span class="va">`fit2`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`fit_celltype_programs`](../reference/fit_celltype_programs.html)</span><span class="op">`(`</span><span class="va">`x`</span>`, backend `<span class="op">`=`</span>` `<span class="st">`"cebmf_torch"`</span>`,`\
`                              max_factors `<span class="op">`=`</span>` `<span class="fl">`6`</span>`, seed `<span class="op">`=`</span>` `<span class="fl">`1`</span><span class="op">`)`</span>\
<span class="fu"><a href="https://rdrr.io/r/base/identical.html" class="external-link"><code class="sourceCode r R">identical</code></a></span><span class="op">`(`</span><span class="va">`fit`</span><span class="op">`$`</span><span class="va">`fits`</span><span class="op">`$`</span><span class="va">`CT1`</span><span class="op">`$`</span><span class="va">`flash`</span><span class="op">`$`</span><span class="va">`L_pm`</span>`, `<span class="va">`fit2`</span><span class="op">`$`</span><span class="va">`fits`</span><span class="op">`$`</span><span class="va">`CT1`</span><span class="op">`$`</span><span class="va">`flash`</span><span class="op">`$`</span><span class="va">`L_pm`</span><span class="op">`)`</span>

</div>

</div>

<div class="section level2">

## var_floor passes through<a href="#var_floor-passes-through" class="anchor" aria-label="anchor"></a>

`var_floor` (a positive scalar) is forwarded to the torch fit as a known noise standard deviation, replacing the column-wise noise estimate:

<div id="cb6" class="sourceCode">

\
<span class="va">`fit3`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`fit_celltype_programs`](../reference/fit_celltype_programs.html)</span><span class="op">`(`</span><span class="va">`x`</span>`, backend `<span class="op">`=`</span>` `<span class="st">`"cebmf_torch"`</span>`,`\
`                              max_factors `<span class="op">`=`</span>` `<span class="fl">`6`</span>`, seed `<span class="op">`=`</span>` `<span class="fl">`1`</span>`, var_floor `<span class="op">`=`</span>` `<span class="fl">`0.5`</span><span class="op">`)`</span>\
<span class="fu"><a href="https://rdrr.io/r/base/all.html" class="external-link"><code class="sourceCode r R">all</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/is.finite.html" class="external-link"><code class="sourceCode r R">is.finite</code></a></span><span class="op">`(`</span><span class="va">`fit3`</span><span class="op">`$`</span><span class="va">`fits`</span><span class="op">`$`</span><span class="va">`CT1`</span><span class="op">`$`</span><span class="va">`flash`</span><span class="op">`$`</span><span class="va">`objective_trace`</span><span class="op">`)`</span><span class="op">`)`</span>

</div>

</div>

<div class="section level2">

## Covariate adjustment (residualize mode)<a href="#covariate-adjustment-residualize-mode" class="anchor" aria-label="anchor"></a>

Sample-level covariates are regressed out of the pseudobulk before factorization — pass column names of `sample_metadata` (or a one-sided formula). To show the adjustment doing real work, plant a known batch effect (a global shift on a nuisance gene set for the `b2` samples), fit with and without adjustment, and compare how much program activity aligns with batch:

<div id="cb7" class="sourceCode">

\
<span class="va">`batch_genes`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/sort.html" class="external-link"><code class="sourceCode r R">sort</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/sample.html" class="external-link"><code class="sourceCode r R">sample.int</code></a></span><span class="op">`(`</span><span class="va">`g`</span>`, `<span class="fl">`60`</span><span class="op">`)`</span><span class="op">`)`</span>\
<span class="va">`batch`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/rep.html" class="external-link"><code class="sourceCode r R">rep</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link"><code class="sourceCode r R">c</code></a></span><span class="op">`(`</span><span class="st">`"b1"`</span>`, `<span class="st">`"b2"`</span><span class="op">`)`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link"><code class="sourceCode r R">c</code></a></span><span class="op">`(`</span><span class="fl">`40`</span>`, `<span class="fl">`20`</span><span class="op">`)`</span><span class="op">`)`</span>\
<span class="va">`Yb`</span>` `<span class="op">`<-`</span>` `<span class="va">`Y`</span>\
<span class="va">`Yb`</span><span class="op">`[`</span><span class="va">`batch`</span>` `<span class="op">`==`</span>` `<span class="st">`"b2"`</span>`, `<span class="va">`batch_genes`</span><span class="op">`]`</span>` `<span class="op">`<-`</span>` `<span class="va">`Yb`</span><span class="op">`[`</span><span class="va">`batch`</span>` `<span class="op">`==`</span>` `<span class="st">`"b2"`</span>`, `<span class="va">`batch_genes`</span><span class="op">`]`</span>` `<span class="op">`+`</span>` `<span class="fl">`1.2`</span>\
<span class="va">`meta`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/data.frame.html" class="external-link"><code class="sourceCode r R">data.frame</code></a></span><span class="op">`(`</span>`observation_id `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">rownames</code></a></span><span class="op">`(`</span><span class="va">`Yb`</span><span class="op">`)`</span>`, batch `<span class="op">`=`</span>` `<span class="va">`batch`</span><span class="op">`)`</span>\
<span class="va">`x_batch`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`as_cell_program_data`](../reference/as_cell_program_data.html)</span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/list.html" class="external-link"><code class="sourceCode r R">list</code></a></span><span class="op">`(`</span>`CT1 `<span class="op">`=`</span>` `<span class="va">`Yb`</span><span class="op">`)`</span>`, sample_metadata `<span class="op">`=`</span>` `<span class="va">`meta`</span><span class="op">`)`</span>\
\
<span class="va">`fit_raw`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`fit_celltype_programs`](../reference/fit_celltype_programs.html)</span><span class="op">`(`</span><span class="va">`x_batch`</span>`, backend `<span class="op">`=`</span>` `<span class="st">`"cebmf_torch"`</span>`,`\
`                                 max_factors `<span class="op">`=`</span>` `<span class="fl">`6`</span>`, seed `<span class="op">`=`</span>` `<span class="fl">`1`</span><span class="op">`)`</span>\
<span class="va">`fit_adj`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`fit_celltype_programs`](../reference/fit_celltype_programs.html)</span><span class="op">`(`</span><span class="va">`x_batch`</span>`, backend `<span class="op">`=`</span>` `<span class="st">`"cebmf_torch"`</span>`,`\
`                                 covariates `<span class="op">`=`</span>` `<span class="st">`"batch"`</span>`,`\
`                                 max_factors `<span class="op">`=`</span>` `<span class="fl">`6`</span>`, seed `<span class="op">`=`</span>` `<span class="fl">`1`</span><span class="op">`)`</span>\
\
<span class="va">`batch_num`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/numeric.html" class="external-link"><code class="sourceCode r R">as.numeric</code></a></span><span class="op">`(`</span><span class="va">`batch`</span>` `<span class="op">`==`</span>` `<span class="st">`"b2"`</span><span class="op">`)`</span>\
<span class="va">`max_batch_corr`</span>` `<span class="op">`<-`</span>` `<span class="kw">`function`</span><span class="op">`(`</span><span class="va">`fit`</span><span class="op">`)`</span>` `<span class="op">`{`</span>\
`  `<span class="va">`Z`</span>` `<span class="op">`<-`</span>` `<span class="va">`fit`</span><span class="op">`$`</span><span class="va">`scores`</span><span class="op">`$`</span><span class="va">`CT1`</span>\
`  `<span class="kw">`if`</span>` `<span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/NULL.html" class="external-link"><code class="sourceCode r R">is.null</code></a></span><span class="op">`(`</span><span class="va">`Z`</span><span class="op">`)`</span>` `<span class="op">`||`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/nrow.html" class="external-link"><code class="sourceCode r R">ncol</code></a></span><span class="op">`(`</span><span class="va">`Z`</span><span class="op">`)`</span>` `<span class="op">`==`</span>` `<span class="fl">`0`</span><span class="op">`)`</span>` `<span class="kw"><a href="https://rdrr.io/r/base/function.html" class="external-link"><code class="sourceCode r R">return</code></a></span><span class="op">`(`</span><span class="cn">`NA_real_`</span><span class="op">`)`</span>\
`  `<span class="fu"><a href="https://rdrr.io/r/base/Extremes.html" class="external-link"><code class="sourceCode r R">max</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/MathFun.html" class="external-link"><code class="sourceCode r R">abs</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/apply.html" class="external-link"><code class="sourceCode r R">apply</code></a></span><span class="op">`(`</span><span class="va">`Z`</span>`, `<span class="fl">`2`</span>`, `<span class="va">`cor`</span>`, y `<span class="op">`=`</span>` `<span class="va">`batch_num`</span><span class="op">`)`</span><span class="op">`)`</span><span class="op">`)`</span>\
<span class="op">`}`</span>\
<span class="fu"><a href="https://rdrr.io/r/base/Round.html" class="external-link"><code class="sourceCode r R">round</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link"><code class="sourceCode r R">c</code></a></span><span class="op">`(`</span>`unadjusted `<span class="op">`=`</span>` `<span class="fu">`max_batch_corr`</span><span class="op">`(`</span><span class="va">`fit_raw`</span><span class="op">`)`</span>`,`\
`        adjusted   `<span class="op">`=`</span>` `<span class="fu">`max_batch_corr`</span><span class="op">`(`</span><span class="va">`fit_adj`</span><span class="op">`)`</span><span class="op">`)`</span>`, `<span class="fl">`3`</span><span class="op">`)`</span>

</div>

The unadjusted fit devotes its strongest component to the batch shift; after `covariates = "batch"` the alignment drops and the planted programs — which are independent of batch — carry the factorization instead.

Two contract notes: `covariate_mode = "fixed"` (covariates as fixed loadings inside the factorization) is **flashier-only** — the torch backend aborts with guidance rather than silently residualizing. And the covariate design must leave at least 10 residual degrees of freedom in every cell type; the fit fails closed otherwise.

</div>

<div class="section level2">

## Downstream is backend-agnostic<a href="#downstream-is-backend-agnostic" class="anchor" aria-label="anchor"></a>

<div id="cb8" class="sourceCode">

\
<span class="va">`fit`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`canonicalize_programs`](../reference/canonicalize_programs.html)</span><span class="op">`(`</span><span class="va">`fit`</span><span class="op">`)`</span>\
<span class="fu">[`program_summary`](../reference/program_summary.html)</span><span class="op">`(`</span><span class="va">`fit`</span><span class="op">`)`</span><span class="op">`[`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link"><code class="sourceCode r R">c</code></a></span><span class="op">`(`</span><span class="st">`"cell_type"`</span>`, `<span class="st">`"n_obs"`</span>`, `<span class="st">`"n_programs"`</span>`, `<span class="st">`"R2"`</span>`, `<span class="st">`"elbo"`</span><span class="op">`)`</span><span class="op">`]`</span>\
\
<span class="va">`pid`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/utils/head.html" class="external-link"><code class="sourceCode r R">head</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">colnames</code></a></span><span class="op">`(`</span><span class="fu">[`program_scores`](../reference/program_scores.html)</span><span class="op">`(`</span><span class="va">`fit`</span><span class="op">`)`</span><span class="op">`)`</span>`, `<span class="fl">`1`</span><span class="op">`)`</span>\
<span class="fu">[`top_genes`](../reference/top_genes.html)</span><span class="op">`(`</span><span class="va">`fit`</span>`, `<span class="va">`pid`</span>`, n `<span class="op">`=`</span>` `<span class="fl">`8`</span><span class="op">`)`</span>

</div>

<div id="cb9" class="sourceCode">

\
<span class="fu">[`plot_top_genes`](../reference/plot_top_genes.html)</span><span class="op">`(`</span><span class="va">`fit`</span>`, programs `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/utils/head.html" class="external-link"><code class="sourceCode r R">head</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">colnames</code></a></span><span class="op">`(`</span><span class="fu">[`program_scores`](../reference/program_scores.html)</span><span class="op">`(`</span><span class="va">`fit`</span><span class="op">`)`</span><span class="op">`)`</span>`, `<span class="fl">`2`</span><span class="op">`)`</span>`, n `<span class="op">`=`</span>` `<span class="fl">`8`</span><span class="op">`)`</span>

</div>

</div>

<div class="section level2">

## When to prefer which backend<a href="#when-to-prefer-which-backend" class="anchor" aria-label="anchor"></a>

- **flashier (default)**: pure R, ELBO-based null checks with a fallback ladder, mature on small pseudobulk problems, and the only backend with a `covariate_mode = "fixed"` path.
- **cebmf_torch**: torch acceleration, covariate-EBMF machinery, and loud failure when the fit is unusable (no silent fallback). Covariates are handled identically for both backends — residualized in the R preprocessing stage (`covariate_mode = "residualize"`, the default).

<!-- -->

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
    #>  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0<br>
    #>  [5] xfun_0.61         cachem_1.1.0      knitr_1.52        htmltools_0.5.9<br>
    #>  [9] rmarkdown_2.32    lifecycle_1.0.5   cli_3.6.6         sass_0.4.10<br>
    #> [13] pkgdown_2.2.1     textshaping_0.4.0 jquerylib_0.1.4   systemfonts_1.3.2
    #> [17] compiler_4.4.2    tools_4.4.2       ragg_1.5.2        bslib_0.12.0<br>
    #> [21] evaluate_1.0.5    yaml_2.3.12       otel_0.2.0        jsonlite_2.0.0<br>
    #> [25] rlang_1.3.0       fs_2.1.0          htmlwidgets_1.6.4

</div>
