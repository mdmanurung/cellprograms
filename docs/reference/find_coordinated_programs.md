<div class="page-header">

# Find coordinated multicellular programs

<div class="d-none name">

`find_coordinated_programs.Rd`

</div>

</div>

<div class="ref-description section level2">

Detects groups of programs from different cell types whose donor-level scores co-vary (coordinated multicellular programs, CMPs), with permutation-max-adjusted node p-values from within-strata donor permutation. General learned-factor, edge-family and CMP-set error control are not established by this calculation. Gene overlap plays no role in the call; programs in a CMP may be thematically unrelated.

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
find_coordinated_programs(
  fit,
  strata = NULL,
  context = NULL,
  mode = c("marginal", "within"),
  ambient = c("drop", "flag", "none"),
  ambient_threshold = 0.3,
  profiles = NULL,
  alpha = 0.05,
  n_perm = 199L,
  tau = 2,
  link_threshold = 3,
  min_shared = 20L,
  min_donors = 20L,
  seed = 1L
)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

fit<a href="#arg-fit" class="anchor" aria-label="anchor"></a><br>
A canonicalized \`cell_program_fit\` (see \`canonicalize_programs()\`).

strata<a href="#arg-strata" class="anchor" aria-label="anchor"></a><br>
Optional blocking variable (e.g. site or batch): a \`sample_metadata\` column name or a vector named by observation ID. Scores are residualized on it and donors are permuted only within it.

context<a href="#arg-context" class="anchor" aria-label="anchor"></a><br>
Optional donor-level variable (e.g. disease status) used by \`mode = "within"\`; numeric variables with more than 5 distinct values are binned into quintiles.

mode<a href="#arg-mode" class="anchor" aria-label="anchor"></a><br>
\`"marginal"\` (default): any shared donor variation counts, including variation driven by \`context\`. \`"within"\`: only coordination that persists within levels of \`context\` (residualized on, and permuted within, strata x context).

ambient<a href="#arg-ambient" class="anchor" aria-label="anchor"></a><br>
Treatment of ambient-like programs (see \`ambient_program_scores()\`): \`"drop"\` (default; excluded before testing), \`"flag"\` (tested, flagged in the output) or \`"none"\`.

ambient_threshold<a href="#arg-ambient-threshold" class="anchor" aria-label="anchor"></a><br>
Ambient score above which a program is ambient-like (default 0.3).

profiles<a href="#arg-profiles" class="anchor" aria-label="anchor"></a><br>
Optional per-cell-type mean log-expression profiles for the ambient score; defaults to column means of \`fit\$pseudobulk\`.

alpha<a href="#arg-alpha" class="anchor" aria-label="anchor"></a><br>
Nominal node rejection threshold (default 0.05); rejection uses strict \`node_p \< alpha\`. This is not a certified CMP-set error rate.

n_perm<a href="#arg-n-perm" class="anchor" aria-label="anchor"></a><br>
Number of permutations (default 199).

tau<a href="#arg-tau" class="anchor" aria-label="anchor"></a><br>
Node-statistic offset; only edges with \|T\| \> \`tau\` contribute (default 2).

link_threshold<a href="#arg-link-threshold" class="anchor" aria-label="anchor"></a><br>
Minimum \|T\| for an edge to join two significant programs into one CMP (default 3).

min_shared<a href="#arg-min-shared" class="anchor" aria-label="anchor"></a><br>
Minimum donors shared by two cell types for their programs to be compared (default 20).

min_donors<a href="#arg-min-donors" class="anchor" aria-label="anchor"></a><br>
Minimum donors for a cell type to be used (default 20).

seed<a href="#arg-seed" class="anchor" aria-label="anchor"></a><br>
Random seed.

</div>

<div class="section level2">

## Value<a href="#value" class="anchor" aria-label="anchor"></a>

An object of class \`cmp_result\`, a list with

cmps<br>
one row per (CMP, program): \`cmp\`, \`cell_type\`, \`program_id\`, \`node_p\` (permutation-max-adjusted), \`sign\` (orientation relative to the CMP's first program, propagated along linking edges as the product of edge signs along the path; CAND-05) and \`ambient_flag\`.

edges<br>
edges among members of the same CMP (\`stat\` = \|T\|).

nodes<br>
all tested programs with node statistic, adjusted p-value and ambient score.

ambient<br>
ambient scores of all programs (or NULL).

settings<br>
the arguments used.

</div>

<div class="section level2">

## Observed support<a href="#observed-support" class="anchor" aria-label="anchor"></a>

Each program is rank-transformed and projected only on its observed donors; unused blocking levels are removed before a rank-aware QR projection. Program callability requires at least \`min_donors\` observations, positive residual degrees of freedom and nonconstant residuals. Edges additionally require at least \`max(4, min_shared)\` pairwise-complete donors and a defined correlation. Unavailable edges and nodes without any callable edge have NA statistics/p-values, not scientific zeros. \`program_support\` records counts, design rank, residual df and callability; \`edge_support\` records every candidate edge's overlap and callability, including unavailable edges. These guards are eligibility checks, not a new calibration guarantee under arbitrary missingness or disconnected donor support.

</div>

<div class="section level2">

## Inferential scope<a href="#inferential-scope" class="anchor" aria-label="anchor"></a>

The input score matrices, selected programs and ambient filter are held fixed during permutation; the full discovery pipeline is not refitted. The proposed fixed-score node null is independence from the joint score collection of all other cell types, conditional on the declared blocks. Its validity requires an appropriate exchangeability argument; residualizing on blocks is not a blanket guarantee. Row permutation does not implement patient-cluster permutation for repeated visits. A future known-score calibration cannot be relabeled as strong end-to-end learned-factor FWER. Fixed-threshold links and greedy grouping add no edge- or CMP-set error guarantee; two genuine groups can be falsely bridged without rejecting a true-null node. Calibration campaigns and scientific signoff remain pending.

</div>

<div class="section level2">

## Program cap<a href="#program-cap" class="anchor" aria-label="anchor"></a>

CMP detection depends on each member program existing upstream. In the historical cmp-bench-v1.2 aggregate reports, cap 30 performed better than cap 10 on the reported membership endpoint. Per-replicate raw outputs have not been recovered in the bounded implementation audit, so these reports are historical summary-only evidence, not an independently reproduced general cap recommendation or error-control certificate. Select the cap by a prespecified training-only rule for the intended analysis. The fitting default remains 10. A one-time message is shown when a fit with \`max_factors \< 30\` has a cell type at its cap.

</div>

<div class="section level2">

## Missing donor scores<a href="#missing-donor-scores" class="anchor" aria-label="anchor"></a>

Programs with missing donor scores (NA) are handled by pairwise-complete donors: rank transformation and strata residualization keep NA donors as NA, cross-cell-type correlations use the donors observed in both programs, and edge \`n\` reports the pairwise-complete donor count for that program pair. A warning lists the affected programs. Missingness reduces power (fewer effective donors); it does not manufacture high-scoring donors.

</div>

<div class="section level2">

## Examples<a href="#ref-examples" class="anchor" aria-label="anchor"></a>

``` r
# Toy canonicalized fit: program 1 of cell types A and B shares donor
# activity z (a planted coordinated program) and top genes G1-G15.
set.seed(1)
ids <- paste0("d", seq_len(40))
genes <- paste0("G", seq_len(100))
z <- rnorm(40)
mk <- function(ct) {
  S <- matrix(rnorm(40 * 2), 40, 2, dimnames = list(ids, paste0(ct, "_", 1:2)))
  S[, 1] <- z + rnorm(40, sd = 0.3)
  W <- matrix(rnorm(length(genes) * 2, sd = 0.1), length(genes), 2,
              dimnames = list(genes, colnames(S)))
  W[1:15, 1] <- W[1:15, 1] + 1
  list(S = S, W = W)
}
ta <- mk("A"); tb <- mk("B")
fit <- structure(list(
  scores = list(A = ta$S, B = tb$S), loadings = list(A = ta$W, B = tb$W),
  cell_types = c("A", "B"), observation_ids = ids,
  sample_metadata = data.frame(observation_id = ids,
                               group = rep(c("ctrl", "stim"), each = 20))),
  class = c("cell_program_fit", "cell_program_data", "list"))
res <- find_coordinated_programs(fit, n_perm = 99, ambient = "none")
res$cmps
```

</div>
