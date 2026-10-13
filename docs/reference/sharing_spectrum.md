# Pairwise sharing spectra across all cell types

Pairwise sharing spectra across all cell types

## Usage

``` r
sharing_spectrum(
  fit,
  pairs = NULL,
  space = c("scores", "loadings"),
  n_perm = NULL,
  ...
)
```

## Arguments

- fit:

  A \`cell_program_fit\` object.

- pairs:

  Optional list of length-2 cell-type vectors; default all pairs.

- space:

  Passed to \`principal_angles\` (default \`"scores"\`).

- n_perm:

  Permutations per pair. Default (NULL) scales with the number of pairs,
  \`max(999, 40 \* \#pairs)\`, so that the permutation floor
  \`1/(n_perm+1)\` is below the BH level of a single true pair (\`0.05 /
  \#pairs\`); a smaller value warns. Pairs get distinct seeds derived
  from \`seed\`.

- ...:

  Passed to \`principal_angles\` (e.g. \`seed\`, \`strata\`).

## Value

Per-pair \`principal_angles\` objects plus a summary table. \`p_pair\`
is a single permutation p-value per pair (statistic sum of squared
cosines); \`z_pair\` and \`excess_frac\` are graded effect sizes against
the permutation null (z-score of the statistic; fraction of the smaller
subspace's dimensions shared beyond chance); \`q_pair\` is BH across
callable pairs; \`underpowered\` marks pairs whose chance-level cosine
\`(sqrt(ka)+sqrt(kb))/sqrt(n)\` is \>= 1 (excluded from BH, never
called); \`shared\` = \`q_pair \< 0.05\`. Needs \`n_perm \> 0\`.

## Examples

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- fit_celltype_programs(x, max_factors = 3, features = "variable", n_variable_genes = 150)
fit <- canonicalize_programs(fit)
ss <- sharing_spectrum(fit, n_perm = 199, seed = 1)
ss$summary
#>           pair k_a k_b effective_dim min_angle_deg n_shared_05 p_pair
#> 1     B vs CD4   2   1            24          11.4           1  0.005
#> 2     B vs CD8   2   1            24          88.5           0  1.000
#> 3      B vs NK   2   2            24          12.0           1  0.005
#> 4    B vs Mono   2   2            24          10.9           1  0.005
#> 5   CD4 vs CD8   1   1            24          88.1           0  0.905
#> 6    CD4 vs NK   1   2            24           9.0           1  0.005
#> 7  CD4 vs Mono   1   2            24          10.6           1  0.005
#> 8    CD8 vs NK   1   2            24          10.9           1  0.005
#> 9  CD8 vs Mono   1   2            24          88.2           0  1.000
#> 10  NK vs Mono   2   2            24           8.6           1  0.005
#>        z_pair excess_frac underpowered      q_pair shared
#> 1  11.0233569  0.95722189        FALSE 0.007142857   TRUE
#> 2  -1.0524466 -0.09714698        FALSE 1.000000000  FALSE
#> 3   6.8369595  0.42814977        FALSE 0.007142857   TRUE
#> 4   7.3284086  0.44892048        FALSE 0.007142857   TRUE
#> 5  -0.7759779 -0.04408835        FALSE 1.000000000  FALSE
#> 6  12.8766472  0.97320021        FALSE 0.007142857   TRUE
#> 7  11.1441564  0.96258391        FALSE 0.007142857   TRUE
#> 8   9.5325177  0.96016475        FALSE 0.007142857   TRUE
#> 9  -1.1733655 -0.11093361        FALSE 1.000000000  FALSE
#> 10  7.6289088  0.43957052        FALSE 0.007142857   TRUE
```
