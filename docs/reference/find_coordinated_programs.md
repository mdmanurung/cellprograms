# Find coordinated multicellular programs

Detects groups of programs from different cell types whose donor-level
scores co-vary (coordinated multicellular programs, CMPs), with
family-wise error control by within-strata donor permutation. Gene
overlap plays no role in the call; programs in a CMP may be thematically
unrelated.

## Usage

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

## Arguments

- fit:

  A canonicalized \`cell_program_fit\` (see
  \`canonicalize_programs()\`).

- strata:

  Optional blocking variable (e.g. site or batch): a \`sample_metadata\`
  column name or a vector named by observation ID. Scores are
  residualized on it and donors are permuted only within it.

- context:

  Optional donor-level variable (e.g. disease status) used by \`mode =
  "within"\`; numeric variables with more than 5 distinct values are
  binned into quintiles.

- mode:

  \`"marginal"\` (default): any shared donor variation counts, including
  variation driven by \`context\`. \`"within"\`: only coordination that
  persists within levels of \`context\` (residualized on, and permuted
  within, strata x context).

- ambient:

  Treatment of ambient-like programs (see \`ambient_program_scores()\`):
  \`"drop"\` (default; excluded before testing), \`"flag"\` (tested,
  flagged in the output) or \`"none"\`.

- ambient_threshold:

  Ambient score above which a program is ambient-like (default 0.3).

- profiles:

  Optional per-cell-type mean log-expression profiles for the ambient
  score; defaults to column means of \`fit\$pseudobulk\`.

- alpha:

  Family-wise error level (default 0.05).

- n_perm:

  Number of permutations (default 199).

- tau:

  Node-statistic offset; only edges with \|T\| \> \`tau\` contribute
  (default 2).

- link_threshold:

  Minimum \|T\| for an edge to join two significant programs into one
  CMP (default 3).

- min_shared:

  Minimum donors shared by two cell types for their programs to be
  compared (default 20).

- min_donors:

  Minimum donors for a cell type to be used (default 20).

- seed:

  Random seed.

## Value

An object of class \`cmp_result\`, a list with

- cmps:

  one row per (CMP, program): \`cmp\`, \`cell_type\`, \`program_id\`,
  \`node_p\` (FWER-adjusted), \`sign\` (direction relative to the CMP's
  first program) and \`ambient_flag\`.

- edges:

  edges among members of the same CMP (\`stat\` = \|T\|).

- nodes:

  all tested programs with node statistic, adjusted p-value and ambient
  score.

- ambient:

  ambient scores of all programs (or NULL).

- settings:

  the arguments used.

## Program cap

CMP detection depends on each member program existing upstream. In the
cmp-bench-v1.2 held-out benchmark, fits with \`max_factors = 30\` raised
pooled membership F1 from 0.571 to 0.705 over the default \`max_factors
= 10\` (paired difference 0.127, 95 guardrails \<= 0.06. The
\`fit_celltype_programs()\` default of 10 was tuned for a different goal
(donor-level disease signal in the concatenated representation), so it
is not changed; use \`max_factors = 30\` when the fit is intended for
\`find_coordinated_programs()\`. A one-time message is shown when a fit
with \`max_factors \< 30\` has a cell type at its cap.

## Examples

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
#>   cmp cell_type program_id node_p sign ambient_flag
#> 1   1         A        A_1   0.01    1           NA
#> 2   1         B        B_1   0.01    1           NA
```
