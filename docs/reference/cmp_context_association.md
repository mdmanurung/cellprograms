# Association of coordinated programs with a donor variable

Regresses each CMP's consensus donor score (\`cmp_scores()\`) on a
donor-level variable, optionally adjusting for covariates.

## Usage

``` r
cmp_context_association(fit, cmp, variable, covariates = NULL)
```

## Arguments

- fit:

  A canonicalized \`cell_program_fit\`.

- cmp:

  A \`cmp_result\` or its \`cmps\` table.

- variable:

  A \`sample_metadata\` column name or a vector named by observation ID.

- covariates:

  Optional character vector of \`sample_metadata\` columns to adjust
  for.

## Value

Data frame: \`cmp\`, \`term\`, \`estimate\`, \`t\`, \`p\`, \`n\`.

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
cmp_context_association(fit, res, "group")
#> Error in cmp_context_association(fit, res, "group"): could not find function "cmp_context_association"
```
