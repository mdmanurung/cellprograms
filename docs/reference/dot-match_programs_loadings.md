# Candidate cross-cell-type program matches from loading correlation

For every pair of cell types, correlates (absolute Pearson) each
program's canonicalized loading vector over the shared gene universe
(intersection of the two feature sets) and calibrates the statistic
against a permutation null: gene labels are permuted within each cell
type's loading matrix (\`n_perm\` permutations; per-pair p = (1 +
exceedances) / (n_perm + 1)). \`q\` is BH across all tested program
pairs.

## Usage

``` r
.match_programs_loadings(fit, n_perm = 199L, seed = 1L)
```

## Arguments

- fit:

  A canonicalized \`cell_program_fit\`.

- n_perm:

  Permutations for the null (default 199; \>= 1 required).

- seed:

  Master seed; each cell-type pair gets a deterministic sub-seed.

## Value

A data.frame with columns ct_a, prog_a, ct_b, prog_b, cor, p, q, sorted
by q then -\|cor\|. Rows with q \< alpha are candidate matches.
