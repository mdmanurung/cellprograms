<div class="section level1">

<div class="page-header">

# cellprograms<a href="#cellprograms" class="anchor" aria-label="anchor"></a>

</div>

<a href="https://github.com/mdmanurung/cellprograms/actions/workflows/R-CMD-check.yaml" class="external-link"><img src="https://github.com/mdmanurung/cellprograms/actions/workflows/R-CMD-check.yaml/badge.svg" alt="R-CMD-check" /></a>

Discover donor-level, cell-type-specific transcriptional programs from pseudobulk expression using empirical-Bayes matrix factorization. Fit locally within each cell type; investigate cross-type coordination afterwards rather than assuming it.

**Development software, not a scientific release certificate.** The supported primary route is ordinary local EBMF with explicit preprocessing, named axes, effective estimator settings and failure accounting. Experimental integrations and historical benchmark summaries are not established by ordinary-core tests. See [support and migration](IMPLEMENTATION_SUPPORT.html) and the installed `validation/CORE_SUPPORT.md` record.

<div class="section level2">

## Installation<a href="#installation" class="anchor" aria-label="anchor"></a>

<div id="cb1" class="sourceCode">

\
<span class="fu">`remotes`</span><span class="fu">`::`</span><span class="fu"><a href="https://remotes.r-lib.org/reference/install_github.html" class="external-link"><code class="sourceCode r R">install_github</code></a></span><span class="op">`(`</span><span class="st">`"mdmanurung/cellprograms"`</span><span class="op">`)`</span>

</div>

Requires R \>= 4.3. GitHub dependency revisions are pinned in `DESCRIPTION`. For reproducibility, install a **full reviewed commit SHA**, not moving `main`, and use the exact-source build/install receipt tooling in <a href="https://github.com/mdmanurung/cellprograms-reproducibility" class="external-link">cellprograms-reproducibility</a>. A package version string alone does not identify installed source.

</div>

<div class="section level2">

## Ordinary local fit<a href="#ordinary-local-fit" class="anchor" aria-label="anchor"></a>

<div id="cb2" class="sourceCode">

\
<span class="kw"><a href="https://rdrr.io/r/base/library.html" class="external-link"><code class="sourceCode r R">library</code></a></span><span class="op">`(`</span><span class="va">`cellprograms`</span><span class="op">`)`</span>\
\
<span class="co">`# One observations x genes matrix per cell type, with unique named axes.`</span>\
<span class="co">`# metadata contains observation_id and any explicitly selected covariates.`</span>\
<span class="va">`x`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`as_cell_program_data`](reference/as_cell_program_data.html)</span><span class="op">`(`</span><span class="va">`pseudobulk`</span>`, sample_metadata `<span class="op">`=`</span>` `<span class="va">`metadata`</span><span class="op">`)`</span>\
<span class="fu">[`validate_cell_program_data`](reference/validate_cell_program_data.html)</span><span class="op">`(`</span><span class="va">`x`</span><span class="op">`)`</span>\
\
<span class="va">`fit`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`fit_celltype_programs`](reference/fit_celltype_programs.html)</span><span class="op">`(`</span>\
`  `<span class="va">`x`</span>`, max_factors `<span class="op">`=`</span>` `<span class="fl">`10`</span>`, feature_scale `<span class="op">`=`</span>` `<span class="st">`"logcounts"`</span>`, seed `<span class="op">`=`</span>` `<span class="fl">`1`</span>\
<span class="op">`)`</span>\
<span class="va">`fit`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`canonicalize_programs`](reference/canonicalize_programs.html)</span><span class="op">`(`</span><span class="va">`fit`</span><span class="op">`)`</span>\
\
<span class="va">`Z`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`program_scores`](reference/program_scores.html)</span><span class="op">`(`</span><span class="va">`fit`</span><span class="op">`)`</span>\
<span class="va">`W`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`program_loadings`](reference/program_loadings.html)</span><span class="op">`(`</span><span class="va">`fit`</span><span class="op">`)`</span>\
<span class="fu">[`program_summary`](reference/program_summary.html)</span><span class="op">`(`</span><span class="va">`fit`</span><span class="op">`)`</span>\
<span class="va">`fit`</span><span class="op">`$`</span><span class="va">`fits`</span><span class="op">`[[`</span><span class="fl">`1`</span><span class="op">`]`</span><span class="op">`]`</span><span class="op">`$`</span><span class="va">`effective_config`</span>

</div>

Declare units truthfully; the software cannot infer that a matrix represents counts, log counts, normalized expression or some other scale. Choice of rank, preprocessing and covariates is analysis-specific, not a universal recommendation established by the default arguments.

The normal activity / point-Laplace loading estimator inherits its native flashier/ebnm behavior. Fallback attempts are warned and recorded; inspect `effective_config` and `fallback`. Numeric `var_floor` is a **fixed known residual standard deviation**, not a lower bound on an estimated variance. Known-noise fits use a binding native path; they must not silently ignore it.

</div>

<div class="section level2">

## Frozen preprocessing and projection<a href="#frozen-preprocessing-and-projection" class="anchor" aria-label="anchor"></a>

<div id="cb3" class="sourceCode">

\
<span class="va">`spec`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`freeze_preprocessing`](reference/freeze_preprocessing.html)</span><span class="op">`(`</span>\
`  `<span class="va">`x`</span>`, features `<span class="op">`=`</span>` `<span class="st">`"variable"`</span>`, n_variable_genes `<span class="op">`=`</span>` `<span class="fl">`2000`</span>`,`\
`  feature_scale `<span class="op">`=`</span>` `<span class="st">`"logcounts"`</span>\
<span class="op">`)`</span>\
<span class="va">`fit`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`canonicalize_programs`](reference/canonicalize_programs.html)</span><span class="op">`(`</span><span class="fu">[`fit_celltype_programs`](reference/fit_celltype_programs.html)</span><span class="op">`(`</span>\
`  `<span class="va">`x`</span>`, preprocessing_spec `<span class="op">`=`</span>` `<span class="va">`spec`</span>`, max_factors `<span class="op">`=`</span>` `<span class="fl">`10`</span>`, seed `<span class="op">`=`</span>` `<span class="fl">`1`</span>\
<span class="op">`)`</span><span class="op">`)`</span>\
<span class="va">`projected`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`project_programs_to_cells`](reference/project_programs_to_cells.html)</span><span class="op">`(`</span>\
`  `<span class="va">`fit`</span>`, expression `<span class="op">`=`</span>` `<span class="va">`query_matrix`</span>`, cell_types `<span class="op">`=`</span>` `<span class="va">`query_cell_types`</span>`,`\
`  feature_scale `<span class="op">`=`</span>` `<span class="st">`"logcounts"`</span>\
<span class="op">`)`</span>

</div>

`query_matrix` has query observations in rows and genes in columns; `query_cell_types` is named by query observation IDs. Supply donor/sample metadata explicitly when the trained design requires it.

Frozen specifications verify their digest, feature identity, transformations and units. Trained covariate transformations are replayed, not learned again from the query composition. Query-time `feature_scale` describes **query** units; it cannot relabel a reference with unknown units. Unknown-unit references, legacy residualized fits without trained replay, or missing native manifests require a truthful refit. Current round-trip tests do not certify every historical package-version pair.

</div>

<div class="section level2">

## Stability and honest denominators<a href="#stability-and-honest-denominators" class="anchor" aria-label="anchor"></a>

<div id="cb4" class="sourceCode">

\
<span class="va">`fit`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`assess_program_stability`](reference/assess_program_stability.html)</span><span class="op">`(`</span><span class="va">`fit`</span>`, n_boot `<span class="op">`=`</span>` `<span class="fl">`20`</span>`, seed `<span class="op">`=`</span>` `<span class="fl">`1`</span>`, progress `<span class="op">`=`</span>` `<span class="cn">`FALSE`</span><span class="op">`)`</span>\
<span class="va">`fit`</span><span class="op">`$`</span><span class="va">`stability`</span><span class="op">`$`</span><span class="va">`replicate_log`</span>\
<span class="va">`fit`</span><span class="op">`$`</span><span class="va">`stability`</span><span class="op">`$`</span><span class="va">`per_program`</span>

</div>

Inspect requested, successful and failed refits separately. The legacy `recovery_freq` is `n_hit / n_boot`; `conditional_recovery_freq` uses successful refits, including successful zero-rank fits. For a rate conditional on a matchable cell type, use `n_hit / n_eligible` when positive. All-failure rates are undefined (`NA`), not zero. Stability evaluates the full fitting procedure, including its documented fallback ladder, not a guarantee of the same selected prior on every resample.

</div>

<div class="section level2">

## Coordination and sharing: interpretation limits<a href="#coordination-and-sharing-interpretation-limits" class="anchor" aria-label="anchor"></a>

[`find_coordinated_programs()`](reference/find_coordinated_programs.html) computes fixed-input permutation-max-adjusted **node** p-values. Programs are not refitted under permutation. Validity depends on the declared null, exchangeability, missingness and selection design; this does not establish general learned-factor, edge-family or merged-CMP-set error control. Ordinary row shuffling is not patient-cluster permutation for repeated visits. Preserve uncallable nodes and reasons rather than interpreting them as negative evidence.

[`program_sharing()`](reference/program_sharing.html) and block memberships are **descriptive**. Activity in both conditions does not prove equal effects, common structure, private biology or calibrated inferential sharing. Prediction and invalid patient-pairing claims in historical case studies remain withdrawn.

</div>

<div class="section level2">

## Optional and historical surfaces<a href="#optional-and-historical-surfaces" class="anchor" aria-label="anchor"></a>

- `cebmf_torch` is a distinct experimental estimator, not numerical equivalence to ordinary flashier. Its supported controls and native environment must be validated independently; unsupported combinations fail early.
- SuperCell, SEACells/mcRigor, optional covariate integrations and metacell projection have separate dependency and aggregation-scale limitations. Ordinary checks do not certify them.
- Historical COVID, DIVAS, immune-aging and comparator material is retained for provenance. Maintained reconstruction and study status belong in the reproducibility repository. Neither a historical report nor an executed tiny fixture is new biological, calibration or independent-reproduction evidence.
- The [previous README](reports/historical/README_before_completion.html) is preserved as an explicitly superseded historical record, not current claims.

</div>

<div class="section level2">

## Validation and documentation<a href="#validation-and-documentation" class="anchor" aria-label="anchor"></a>

``` sh
Rscript -e 'testthat::test_local(".", stop_on_failure = TRUE)'
R CMD build .
# Use the actual tarball produced above, not an assumed release filename.
R CMD check --no-manual cellprograms_<version>.tar.gz
```

The installed test entrypoint executes ordinary native-fit and actual-refit regressions as well as algebra, design replay, units, support and failure tests. Expected warnings are asserted; unexpected or missing messages fail their tests. Optional tests must identify their applicability instead of becoming empty mandatory checks. CI defines Linux current/minimum-R, macOS and Windows lanes; only actual run results establish coverage.

- [Quickstart](vignettes/quickstart.Rmd)
- [NEWS](NEWS.html)
- [Support and migration](IMPLEMENTATION_SUPPORT.html)
- <a href="https://github.com/mdmanurung/cellprograms-reproducibility" class="external-link">Reproducibility repository</a>

</div>

<div class="section level2">

## License and citation<a href="#license-and-citation" class="anchor" aria-label="anchor"></a>

MIT for the package’s licensed code; preserve third-party notices. Dataset, comparator and archived-result rights are separate. Use `citation("cellprograms")` for package attribution and cite inherited methods separately.

</div>

</div>
