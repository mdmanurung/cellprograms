<div class="page-header">

# CorShrink-style null-correlation estimate for a matched program set

<div class="d-none name">

`dot-estimate_V.Rd`

</div>

</div>

<div class="ref-description section level2">

Estimates the correlation of loading estimation noise across cell types (shared donors, ambient leakage) from the z-scores (loading / posterior SD) of null genes — genes with \|z\| \< \`z_thresh\` in EVERY condition of the matched set — and shrinks the correlation matrix toward the identity with a Ledoit-Wolf-style analytic intensity: lambda = sum of plug-in sampling variances (1 - r^2)^2 / (n - 1) over the sum of squared off-diagonal correlations, clipped to \[0, 1\]. With fewer than \`min_null\` null genes the identity is returned (conservative; warns).

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
.estimate_V(fit, matches, z_thresh = 1, min_null = 20L, mats = NULL)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

fit<a href="#arg-fit" class="anchor" aria-label="anchor"></a><br>
A canonicalized \`cell_program_fit\`.

matches<a href="#arg-matches" class="anchor" aria-label="anchor"></a><br>
One matched set: a data.frame with columns \`cell_type\` and \`program_id\` (one row per cell type), or a single-row pairwise match with columns ct_a, prog_a, ct_b, prog_b.

z_thresh<a href="#arg-z-thresh" class="anchor" aria-label="anchor"></a><br>
\|z\| below this in every condition marks a null gene.

min_null<a href="#arg-min-null" class="anchor" aria-label="anchor"></a><br>
Minimum number of null genes for a usable estimate.

</div>

<div class="section level2">

## Value<a href="#value" class="anchor" aria-label="anchor"></a>

The shrunk C x C correlation matrix (dimnames = cell types), with attributes \`n_null\` and \`lambda\`.

</div>
