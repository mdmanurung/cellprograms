# Quick start: from pseudobulk to shared programs

`cellprograms` discovers transcriptional programs **de novo within each
cell type** from pseudobulk single-cell RNA-seq, then investigates
cross-cell-type coordination and sharing afterwards — never assuming it.
This vignette runs the whole pipeline on a small simulated dataset.

## Simulate pseudobulk data

[`simulate_pseudobulk()`](https://mdmanurung.github.io/cellprograms/reference/simulate_pseudobulk.md)
plants private and shared programs with known ground truth. The result
is a named list of observations x genes matrices, one per cell type.

\
[`library`](https://rdrr.io/r/base/library.html)`(``cellprograms``)`\
\
`sim`` ``<-`` `[`simulate_pseudobulk`](https://mdmanurung.github.io/cellprograms/reference/simulate_pseudobulk.md)`(``scenario ``=`` ``"mixed"``, n_samples ``=`` ``60``,`\
`                           n_genes ``=`` ``800``, program_size ``=`` ``40``, seed ``=`` ``1``)`\
[`str`](https://rdrr.io/r/utils/str.html)`(``sim``$``pseudobulk``, max.level ``=`` ``1``)`\
`#> List of 5`\
`#>  $ B   : num [1:60, 1:800] 1.146 -2.403 0.573 0.375 -0.425 ...`\
`#>   ..- attr(*, "dimnames")=List of 2`\
`#>  $ CD4 : num [1:60, 1:800] -1.474 0.59 -0.13 0.372 1.452 ...`\
`#>   ..- attr(*, "dimnames")=List of 2`\
`#>  $ CD8 : num [1:60, 1:800] 0.0679 0.9149 0.8625 -1.1996 -2.2594 ...`\
`#>   ..- attr(*, "dimnames")=List of 2`\
`#>  $ NK  : num [1:60, 1:800] -0.0478 1.5893 1.175 0.0746 -0.9374 ...`\
`#>   ..- attr(*, "dimnames")=List of 2`\
`#>  $ Mono: num [1:60, 1:800] -0.7613 -0.4367 -1.4659 0.0574 -1.4384 ...`\
`#>   ..- attr(*, "dimnames")=List of 2`

## Build the data object and fit programs

\
`x`` ``<-`` `[`as_cell_program_data`](https://mdmanurung.github.io/cellprograms/reference/as_cell_program_data.md)`(``sim``$``pseudobulk``)`\
[`validate_cell_program_data`](https://mdmanurung.github.io/cellprograms/reference/validate_cell_program_data.md)`(``x``)`\
\
`fit`` ``<-`` `[`fit_celltype_programs`](https://mdmanurung.github.io/cellprograms/reference/fit_celltype_programs.md)`(``x``, loading_prior ``=`` ``"point_laplace"``,`\
`                             max_factors ``=`` ``6``, features ``=`` ``"variable"``,`\
`                             n_variable_genes ``=`` ``400``)`\
`fit`` ``<-`` `[`canonicalize_programs`](https://mdmanurung.github.io/cellprograms/reference/canonicalize_programs.md)`(``fit``)`\
`fit`\
`#> ``cell_program_data`\
`#> Cell types: 5 ``"B"``, ``"CD4"``, ``"CD8"``, ``"NK"``, and ``"Mono"`\
`#> B: 60 observations x 800 genes`\
`#> CD4: 60 observations x 800 genes`\
`#> CD8: 60 observations x 800 genes`\
`#> NK: 60 observations x 800 genes`\
`#> Mono: 60 observations x 800 genes`\
`#> Shared genes across all cell types: 800`

[`fit_celltype_programs()`](https://mdmanurung.github.io/cellprograms/reference/fit_celltype_programs.md)
fits one empirical-Bayes matrix factorization per cell type
(‘flashier’);
[`canonicalize_programs()`](https://mdmanurung.github.io/cellprograms/reference/canonicalize_programs.md)
resolves sign and scale ambiguity so programs are comparable downstream.
Real-data fits should additionally pass `covariates`
(e.g. `c("site", "sex")`) and use `max_factors = 30` when the goal is
CMP detection (see the CMP benchmark note in `README`).

## Inspect programs

\
[`program_summary`](https://mdmanurung.github.io/cellprograms/reference/program_summary.md)`(``fit``)``          ``# one row per cell type: K, R2, ELBO`\
`#>   cell_type n_obs n_genes n_programs        R2      elbo`\
`#> 1         B    60     400          2 0.2439779 -36065.09`\
`#> 2       CD4    60     400          1 0.1707640 -35709.99`\
`#> 3       CD8    60     400          1 0.2010540 -35859.83`\
`#> 4        NK    60     400          2 0.2948849 -36081.50`\
`#> 5      Mono    60     400          2 0.2997384 -35859.12`\
[`top_genes`](https://mdmanurung.github.io/cellprograms/reference/top_genes.md)`(``fit``, ``"B_1"``, n ``=`` ``8``)``  ``# top loadings of program B_1`\
`#>    gene   loading          lfsr`\
`#> 1  G474  4.582835 2.630908e-244`\
`#> 3  G733  3.164218 1.075143e-117`\
`#> 2  G422  2.724563  7.611642e-88`\
`#> 4  G197 -2.688759  0.000000e+00`\
`#> 6  G462  2.555132  1.603963e-77`\
`#> 8  G349 -2.213170  0.000000e+00`\
`#> 9  G514  2.067245  1.752200e-51`\
`#> 10 G720  2.031627  8.535883e-50`

Scores (donor x program activities) and loadings (gene x program) are
one call away:

\
`Z`` ``<-`` `[`program_scores`](https://mdmanurung.github.io/cellprograms/reference/program_scores.md)`(``fit``)``    ``# donors x all programs (wide matrix)`\
`W`` ``<-`` `[`program_loadings`](https://mdmanurung.github.io/cellprograms/reference/program_loadings.md)`(``fit``)``  ``# genes x all programs`\
[`dim`](https://rdrr.io/r/base/dim.html)`(``Z``)``; `[`dim`](https://rdrr.io/r/base/dim.html)`(``W``)`\
`#> [1] 60  8`\
`#> NULL`

## Stability

[`assess_program_stability()`](https://mdmanurung.github.io/cellprograms/reference/assess_program_stability.md)
refits on donor subsamples and reports how often each program is
recovered.

\
`fit`` ``<-`` `[`assess_program_stability`](https://mdmanurung.github.io/cellprograms/reference/assess_program_stability.md)`(``fit``, n_boot ``=`` ``4``, resample_frac ``=`` ``0.7``,`\
`                                seed ``=`` ``1``, progress ``=`` ``FALSE``)`\
`fit``$``stability``$``per_program`\
`#>   program_id cell_type recovery_freq median_score_cor median_loading_cos`\
`#> 1        B_1         B             1        0.9999571          0.9974227`\
`#> 2        B_2         B             1        0.9997754          0.9939955`\
`#> 3      CD4_1       CD4             1        0.9999864          0.9983564`\
`#> 4      CD8_1       CD8             1        0.9999776          0.9980614`\
`#> 5       NK_1        NK             1        0.9999427          0.9970874`\
`#> 6       NK_2        NK             1        0.9999653          0.9954930`\
`#> 7     Mono_1      Mono             1        0.9999674          0.9958291`\
`#> 8     Mono_2      Mono             1        0.9999071          0.9953091`\
`#>   sign_consistency top_gene_jaccard`\
`#> 1                1        1.0000000`\
`#> 2                1        0.9047619`\
`#> 3                1        0.9047619`\
`#> 4                1        0.9047619`\
`#> 5                1        0.9047619`\
`#> 6                1        1.0000000`\
`#> 7                1        0.9047619`\
`#> 8                1        1.0000000`

## Coordinated multicellular programs (CMPs)

[`find_coordinated_programs()`](https://mdmanurung.github.io/cellprograms/reference/find_coordinated_programs.md)
tests whether programs from *different* cell types co-vary across
donors, with family-wise error control by permutation. Gene overlap
plays no role in the call.

\
`cmp`` ``<-`` `[`find_coordinated_programs`](https://mdmanurung.github.io/cellprograms/reference/find_coordinated_programs.md)`(``fit``, n_perm ``=`` ``199``, seed ``=`` ``1``)`\
`cmp`\
`#> <cmp_result> 1 coordinated multicellular program(s) at FWER 0.05 (mode: marginal, ambient: drop)`\
`#>   CMP 1: B_1, CD4_1 (-), NK_2 (-), Mono_1 (-)`

## Shared vs specific programs

[`program_sharing()`](https://mdmanurung.github.io/cellprograms/reference/program_sharing.md)
goes the other way: it matches programs across cell types by loading
correlation (permutation-calibrated), then asks with ‘mashr’ whether the
matched programs share their gene *effects* — correlated loadings alone
do not imply shared effects.

\
`if`` ``(`[`requireNamespace`](https://rdrr.io/r/base/ns-load.html)`(``"mashr"``, quietly ``=`` ``TRUE``)``)`` ``{`\
`  ``shr`` ``<-`` `[`program_sharing`](https://mdmanurung.github.io/cellprograms/reference/program_sharing.md)`(``fit``, n_perm ``=`` ``99``, seed ``=`` ``1``)`\
`  ``shr`\
`}`\
`#> program_sharing: 0 matched group(s), 8 program(s) assessed`\
`#>   8 unmatched program(s): private (no significant match)`

## Artifact genes

Immunoglobulin/TCR V(D)J segments, known donor-genotype markers, and
sex-chromosome genes can dominate programs for technical or genetic
rather than biological reasons. Flag them in existing fits, or exclude
them at fit time (exclusion happens before variable-gene selection):

\
[`flag_artifact_genes`](https://mdmanurung.github.io/cellprograms/reference/flag_artifact_genes.md)`(``fit``)`\
`#>   cell_type program_id vdj_frac vdj_flag genotype_marker_frac`\
`#> 1         B        B_1        0    FALSE                    0`\
`#> 2         B        B_2        0    FALSE                    0`\
`#> 3       CD4      CD4_1        0    FALSE                    0`\
`#> 4       CD8      CD8_1        0    FALSE                    0`\
`#> 5        NK       NK_1        0    FALSE                    0`\
`#> 6        NK       NK_2        0    FALSE                    0`\
`#> 7      Mono     Mono_1        0    FALSE                    0`\
`#> 8      Mono     Mono_2        0    FALSE                    0`\
`#>   genotype_marker_flag`\
`#> 1                FALSE`\
`#> 2                FALSE`\
`#> 3                FALSE`\
`#> 4                FALSE`\
`#> 5                FALSE`\
`#> 6                FALSE`\
`#> 7                FALSE`\
`#> 8                FALSE`\
\
`# fit-time masking:`\
`# excl <- artifact_genes(colnames(sim$pseudobulk$B))`\
`# fit_masked <- fit_celltype_programs(x, exclude_genes = excl, ...)`

## Next steps

- The real-data case-study article (pkgdown site) walks through the
  Stephenson and COMBAT COVID cohorts, including the artifact-masking
  sensitivity analysis.
- `benchmarks/` contains the pre-registered simulation benchmarks behind
  the default settings and the CMP caller.
