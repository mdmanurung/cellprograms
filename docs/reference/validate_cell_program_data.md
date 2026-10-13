# Validate a cell program data object

Runs the full data-contract check and returns the object invisibly,
emitting warnings for soft-failure conditions (small N, near-zero
variance genes, extreme sparsity).

## Usage

``` r
validate_cell_program_data(x)
```

## Arguments

- x:

  A \`cell_program_data\` object.

## Examples

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
validate_cell_program_data(x)
```
