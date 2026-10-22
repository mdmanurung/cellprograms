<div class="page-header">

# Principal angles between the program subspaces of two cell types

<div class="d-none name">

`principal_angles.Rd`

</div>

</div>

<div class="ref-description section level2">

Principal angles between the program subspaces of two cell types

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

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

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

fit<a href="#arg-fit" class="anchor" aria-label="anchor"></a><br>
A \`cell_program_fit\` object.

ct_a, ct_b<a href="#arg-ct-a-ct-b" class="anchor" aria-label="anchor"></a><br>
Cell-type labels.

space<a href="#arg-space" class="anchor" aria-label="anchor"></a><br>
\`"scores"\` (activity coordination; default) or \`"loadings"\` (shared gene usage).

weight<a href="#arg-weight" class="anchor" aria-label="anchor"></a><br>
Loadings space only: \`"none"\` or \`"residual_sd"\` whitening.

n_perm<a href="#arg-n-perm" class="anchor" aria-label="anchor"></a><br>
Permutations for the label-permutation null (0 skips).

seed<a href="#arg-seed" class="anchor" aria-label="anchor"></a><br>
Random seed for permutations.

strata<a href="#arg-strata" class="anchor" aria-label="anchor"></a><br>
Scores space only: named vector (names = observation ids) of stratum labels, e.g. batch or institute. Observations of \`ct_b\` are permuted only within strata, so the null preserves stratum-level structure shared by both cell types (blocked permutation).

</div>

<div class="section level2">

## Value<a href="#value" class="anchor" aria-label="anchor"></a>

An object of class \`principal_angles\`.

</div>
