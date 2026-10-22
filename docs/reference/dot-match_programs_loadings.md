<div class="page-header">

# Candidate cross-cell-type program matches from loading correlation

<div class="d-none name">

`dot-match_programs_loadings.Rd`

</div>

</div>

<div class="ref-description section level2">

For every pair of cell types, correlates (absolute Pearson) each program's canonicalized loading vector over the shared gene universe (intersection of the two feature sets) and calibrates the statistic against a permutation null: gene labels are permuted within each cell type's loading matrix (\`n_perm\` permutations; per-pair p = (1 + exceedances) / (n_perm + 1)). \`q\` is BH across all tested program pairs.

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
.match_programs_loadings(fit, n_perm = 199L, seed = 1L)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

fit<a href="#arg-fit" class="anchor" aria-label="anchor"></a><br>
A canonicalized \`cell_program_fit\`.

n_perm<a href="#arg-n-perm" class="anchor" aria-label="anchor"></a><br>
Permutations for the null (default 199; \>= 1 required).

seed<a href="#arg-seed" class="anchor" aria-label="anchor"></a><br>
Master seed; each cell-type pair gets a deterministic sub-seed.

</div>

<div class="section level2">

## Value<a href="#value" class="anchor" aria-label="anchor"></a>

A data.frame with columns ct_a, prog_a, ct_b, prog_b, cor, p, q, sorted by q then -\|cor\|. Rows with q \< alpha are candidate matches.

</div>
