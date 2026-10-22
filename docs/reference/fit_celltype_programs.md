# Fit cell-type-specific programs with flashier

Fit cell-type-specific programs with flashier

## Usage

``` r
fit_celltype_programs(
  x,
  loading_prior = "point_laplace",
  max_factors = 10,
  center = TRUE,
  scale = FALSE,
  features = "all",
  n_variable_genes = 2000,
  exclude_genes = NULL,
  exclude_artifacts = FALSE,
  var_type = 1,
  backfit = TRUE,
  nullcheck = TRUE,
  covariates = NULL,
  covariate_mode = c("residualize", "fixed"),
  adjust_batches = NULL,
  seed = 1,
  flash_control = list(),
  feature_scale = "unknown",
  preprocessing_spec = NULL,
  backend = c("flashier", "cebmf_torch"),
  var_floor = "none"
)
```

## Arguments

- x:

  A \`cell_program_data\` object.

- loading_prior:

  Gene-side EBNM prior family: "point_laplace" (default),
  "point_normal", or "unimodal".

- max_factors:

  Maximum candidate factors per cell type (greedy_Kmax). Must leave at
  least one residual degree of freedom per cell type: the fit aborts if
  \`max_factors \> n - q - 1\`, where n is the number of donors for the
  cell type and q the number of covariate columns adjusted for (safety
  guard; never binds at n \>= 50 with \`max_factors \<= 10\`). Default
  10: benchmarked on COMBAT pseudobulk (32-config factorial), K=10
  dominated all prior/variance/backfit combinations; K \>= 30 overfits
  and collapses disease-signal recovery. For coordinated-program
  detection (\`find_coordinated_programs()\`), use \`max_factors = 30\`:
  it raised CMP membership F1 in the cmp-bench-v1.2 held-out benchmark.

- center:

  Center genes within each cell type (default TRUE).

- scale:

  Unit-variance scale genes (default FALSE; not recommended).

- features:

  Feature selection mode: "all", "variable", or a named list of gene ID
  vectors per cell type.

- n_variable_genes:

  Number of most-variable genes when \`features = "variable"\`.

- exclude_genes:

  Optional character vector of gene IDs removed from every cell type's
  feature universe BEFORE feature selection, so excluded genes do not
  consume \`n_variable_genes\` slots (with \`features = "variable"\`)
  and never enter the fit. Exclusion changes the feature panel and is
  NOT the default; it is intended for artifact-gene classes such as
  V(D)J or common-germline-variation genes – see \`artifact_gene_sets\`
  and \`artifact_genes()\`. With a named-list \`features\`, excluded
  genes are dropped from the requested vectors after validation.

- exclude_artifacts:

  Convenience preset: when TRUE, union \`artifact_genes()\` (V(D)J,
  genotype-marker, sex and canonical-marker sets) over the pseudobulk
  gene universe into \`exclude_genes\`, so the preset blacklist never
  enters the fit and does not consume \`n_variable_genes\` slots.
  Default FALSE: exclusion changes the feature panel, so it is opt-in;
  explicit \`exclude_genes\` entries are kept.

- var_type:

  flashier residual-variance type: 0 = constant, 1 = one variance per
  row/donor (default), 2 = per-gene (unstable on log-normalized
  pseudobulk), or \`c(1, 2)\` = Kronecker rank-one structure s_ij = a_i
  \* b_j (most flexible; substantially slower).

- backfit:

  Run flashier backfitting (default TRUE).

- nullcheck:

  Run flashier nullcheck (default TRUE).

- covariates:

  Optional sample-level covariates to regress out of each cell type's
  pseudobulk (genes x samples fit separately per cell type, over the
  donors observed in that cell type) before factorization: a character
  vector of \`sample_metadata\` column names or a one-sided formula,
  e.g. \`~ age + sex + batch\`. Programs are then fit to the residuals.
  Do not include the outcome you plan to test. Variable-gene selection
  is done before adjustment. The covariate design must leave \>= 10
  residual degrees of freedom in every cell type (fail-closed abort
  otherwise); levels carried by a single donor trigger a leverage
  warning.

- covariate_mode:

  How \`covariates\` are handled: \`"residualize"\` (default; hard
  adjustment: regress out of \`Y\` before factorization) or \`"fixed"\`
  (soft, SOFA-style: the covariate design columns enter flashier as
  fixed sample-side factors with unshrunk gene effects, and the K free
  programs are fit jointly with them; free programs can still share
  variance with a covariate. Requires \`center = TRUE\`, ignores
  \`flash_control\`). Fixed columns are not returned as programs;
  \`Y_used\` is \`Y\` minus the fitted covariate part. Note: \`"fixed"\`
  fits cannot be projected to single cells —
  \`project_programs_to_cells()\` refuses them because the fixed
  covariate adjustment cannot be replayed on a query.

- adjust_batches:

  Optional name (character scalar) of a batch column in
  \`sample_metadata\`. When given, \`sva::ComBat\` removes the batch
  location/ scale from each cell type's panel before covariate
  residualization. Intended for small-batch designs where OLS
  residualization is underpowered (the batch benchmark found ComBat
  better at ~5 donors/batch but with residual batch structure). Not the
  default; not supported with a frozen \`preprocessing_spec\`;
  \`project_programs_to_cells()\` refuses batch-adjusted fits (the
  adjustment cannot yet be replayed on a query).

- seed:

  Random seed for reproducibility.

- flash_control:

  Optional list overriding flashier::flash arguments.

- feature_scale:

  Declared measurement units of the pseudobulk expression values,
  recorded in the fit manifest and required for
  \`project_programs_to_cells()\` (the projection oracle refuses queries
  whose declared units differ from the source fit's). "unknown"
  (default) keeps the fit unprojectable; pass e.g. "logcounts" for
  log-normalized pseudobulk. Only used on the flashier.utils backend.

- preprocessing_spec:

  Optional frozen preprocessing specification from
  \`freeze_preprocessing()\` (MR-0). When supplied, the spec governs
  feature selection and adjustment: \`features\`, \`n_variable_genes\`,
  \`exclude_genes\`, \`exclude_artifacts\`, \`covariates\`,
  \`covariate_mode\`, \`center\`, and \`scale\` must all be left at
  their defaults (any explicit conflict is an error, so a refit can
  never silently drift from the frozen panel), and the spec's digest is
  recorded in the fit's preprocessing provenance. The spec's cell types
  must cover \`x\` and all frozen genes must be present in the dataset.

- backend:

  EBMF implementation: \`"flashier"\` (default; the flashier.utils R
  backend, with a nullcheck fallback ladder) or \`"cebmf_torch"\` (the
  covariate-EBMF torch backend, which fails loudly instead of falling
  back). Both backends share preprocessing, program extraction, and the
  output contract, with two exceptions: torch fits carry no posterior
  SDs (so \`program_sharing()\` is unavailable) and cannot be projected
  to single cells with \`project_programs_to_cells()\` (no
  flashier.utils manifest). The torch backend runs in a pinned Python
  environment: set options(cellprograms.cebmf_python_env = ...) to point
  at an existing environment, otherwise a basilisk environment is used.

- var_floor:

  Residual SD handling: \`"none"\` (default; residual SDs are estimated
  as usual), \`"vignette"\` (fixed at the flashier single-cell
  vignette's S = sd(log(Poisson(1/n) + 1)), drawn once per cell type
  with a fixed local seed, so the value is deterministic and independent
  of \`seed\`), or a positive numeric scalar S (residual SDs fixed at
  S). Fixed S binds on the \`cebmf_torch\` backend. On the \`flashier\`
  backend the value is passed through but does not currently bind
  (CAND-13: fits with and without a fixed scalar S are bit-identical at
  the pinned flashier version); a fail-closed remediation is pending.

## Value

The input object augmented with class \`cell_program_fit\` and
components \`fits\` (per-cell-type backend fits with preprocessing
replay records), \`preprocessing\`, and \`provenance\`. Program
\`scores\` and \`loadings\` are added by \`canonicalize_programs()\` —
run it before downstream analysis. (There is no \`programs\` component;
an earlier version of this documentation listed one in error.)

## Examples

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- fit_celltype_programs(x, max_factors = 3, features = "variable",
                             n_variable_genes = 150)
fit <- canonicalize_programs(fit)
program_summary(fit)
#>   cell_type n_obs n_genes n_programs        R2      elbo
#> 1         B    24     150          2 0.4107044 -5719.888
#> 2       CD4    24     150          1 0.3463505 -5544.727
#> 3       CD8    24     150          1 0.2399331 -5506.694
#> 4        NK    24     150          2 0.3943826 -5716.026
#> 5      Mono    24     150          2 0.4113416 -5653.769
```
