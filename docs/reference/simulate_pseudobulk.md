# Simulate cell-type pseudobulk matrices with planted programs

Simulate cell-type pseudobulk matrices with planted programs

## Usage

``` r
simulate_pseudobulk(
  scenario = "mixed",
  n_samples = 150,
  n_genes = 1500,
  program_size = 50,
  signal_sd = 1.5,
  noise_sd = 1,
  cell_types = c("B", "CD4", "CD8", "NK", "Mono"),
  seed = 1
)
```

## Arguments

- scenario:

  One of "private" (Scenario 1), "mixed" (Scenario 4, principal),
  "same_genes_independent" (Scenario 7 negative control), or "global"
  (Scenario 2).

- n_samples:

  Number of biological observations.

- n_genes:

  Total gene pool size.

- program_size:

  Number of active genes per program.

- signal_sd:

  Program activity standard deviation.

- noise_sd:

  Residual noise standard deviation.

- cell_types:

  Character vector of cell-type names.

- seed:

  Random seed.

## Value

A list with \`pseudobulk\` (named list of observations x genes
matrices), \`truth\` (per-cell-type list of true score vectors, loading
indices, and program labels), and \`scenario\`.

## Examples

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
names(sim$pseudobulk)
#> [1] "B"    "CD4"  "CD8"  "NK"   "Mono"
dim(sim$pseudobulk$B)
#> [1]  24 300
```
