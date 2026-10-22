<div class="page-header">

# Pairwise sharing spectra across all cell types

<div class="d-none name">

`sharing_spectrum.Rd`

</div>

</div>

<div class="ref-description section level2">

Pairwise sharing spectra across all cell types

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
sharing_spectrum(
  fit,
  pairs = NULL,
  space = c("scores", "loadings"),
  n_perm = NULL,
  seed = 1L,
  ...
)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

fit<a href="#arg-fit" class="anchor" aria-label="anchor"></a><br>
A \`cell_program_fit\` object.

pairs<a href="#arg-pairs" class="anchor" aria-label="anchor"></a><br>
Optional list of length-2 cell-type vectors; default all pairs.

space<a href="#arg-space" class="anchor" aria-label="anchor"></a><br>
Passed to \`principal_angles\` (default \`"scores"\`).

n_perm<a href="#arg-n-perm" class="anchor" aria-label="anchor"></a><br>
Permutations per pair. Default (NULL) scales with the number of pairs, \`max(999, 40 \* \#pairs)\`, so that the permutation floor \`1/(n_perm+1)\` is below the BH level of a single true pair (\`0.05 / \#pairs\`); a smaller value warns. Pairs get distinct seeds derived from \`seed\`.

seed<a href="#arg-seed" class="anchor" aria-label="anchor"></a><br>
Master seed (default 1, so default runs are reproducible); each pair gets a deterministic sub-seed. \`NULL\` leaves the RNG unseeded.

...<a href="#arg--" class="anchor" aria-label="anchor"></a><br>
Passed to \`principal_angles\` (e.g. \`strata\`).

</div>

<div class="section level2">

## Value<a href="#value" class="anchor" aria-label="anchor"></a>

Per-pair \`principal_angles\` objects plus a summary table. \`p_pair\` is a single permutation p-value per pair (statistic sum of squared cosines); \`z_pair\` and \`excess_frac\` are graded effect sizes against the permutation null (z-score of the statistic; fraction of the smaller subspace's dimensions shared beyond chance); \`q_pair\` is BH across callable pairs; \`underpowered\` marks pairs whose chance-level cosine \`(sqrt(ka)+sqrt(kb))/sqrt(n)\` is \>= 1 (excluded from BH, never called); \`shared\` = \`q_pair \< 0.05\`. Needs \`n_perm \> 0\`.

</div>

<div class="section level2">

## Examples<a href="#ref-examples" class="anchor" aria-label="anchor"></a>

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- fit_celltype_programs(x, max_factors = 3, features = "variable", n_variable_genes = 150)
fit <- canonicalize_programs(fit)
ss <- sharing_spectrum(fit, n_perm = 199, seed = 1)
ss$summary
```

</div>
