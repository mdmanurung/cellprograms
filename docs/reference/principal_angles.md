# Principal angles between the program subspaces of two cell types

Principal angles between the program subspaces of two cell types

## Usage

``` r
principal_angles(
  fit,
  ct_a,
  ct_b,
  space = c("scores", "loadings"),
  weight = c("none", "residual_sd"),
  n_perm = 200L,
  seed = NULL,
  strata = NULL
)
```

## Arguments

- fit:

  A \`cell_program_fit\` object.

- ct_a, ct_b:

  Cell-type labels.

- space:

  \`"scores"\` (activity coordination; default) or \`"loadings"\`
  (shared gene usage).

- weight:

  Loadings space only: \`"none"\` or \`"residual_sd"\` whitening.

- n_perm:

  Permutations for the label-permutation null (0 skips).

- seed:

  Random seed for permutations.

- strata:

  Scores space only: named vector (names = observation ids) of stratum
  labels, e.g. batch or institute. Observations of \`ct_b\` are permuted
  only within strata, so the null preserves stratum-level structure
  shared by both cell types (blocked permutation).

## Value

An object of class \`principal_angles\`.

## Examples

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- fit_celltype_programs(x, max_factors = 3, features = "variable", n_variable_genes = 150)
fit <- canonicalize_programs(fit)
# CD8 and NK share a planted donor-activity program in this simulation
pa <- principal_angles(fit, "CD8", "NK", n_perm = 99, seed = 1)
#> Error in principal_angles(fit, "CD8", "NK", n_perm = 99, seed = 1): could not find function "principal_angles"
pa
#> Error: object 'pa' not found
```
