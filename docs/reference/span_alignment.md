<div class="page-header">

# Principal-angle alignment between two loading subspaces

<div class="d-none name">

`span_alignment.Rd`

</div>

</div>

<div class="ref-description section level2">

Computes the principal angles between the column spans of two loading matrices (genes x programs) on the same gene panel — e.g. a discovery fit vs a stability refit, or the same cohort fit by two methods. Each span is orthonormalized (rank-aware via SVD), and the singular values of the cross-product of the bases are the cosines of the principal angles between the spans.

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
span_alignment(W1, W2, k = NULL)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

W1<a href="#arg-w-" class="anchor" aria-label="anchor"></a><br>
Loading matrix (genes x programs) on the common gene panel.

W2<a href="#arg-w-" class="anchor" aria-label="anchor"></a><br>
Second loading matrix (genes x programs) on that same panel. Row names should match; rows are matched by name when both have row names, otherwise by position.

k<a href="#arg-k" class="anchor" aria-label="anchor"></a><br>
Optional number of leading principal cosines to average for \`mean_overlap_topk\` (default: all).

</div>

<div class="section level2">

## Value<a href="#value" class="anchor" aria-label="anchor"></a>

A list with \`principal_cosines\` (descending cosines of the principal angles, length = min(rank1, rank2)), \`mean_overlap\`, \`mean_overlap_topk\`, and \`dim\` (the two span ranks).

</div>

<div class="section level2">

## Details<a href="#details" class="anchor" aria-label="anchor"></a>

The registered tier-2 use (immune-aging application, 2026-10-16) is the \*\*stability-normalized conservation ratio\*\*: the between-cohort span overlap (\`mean_overlap\`) divided by the within-cohort bootstrap ceiling — the same quantity computed between discovery and bootstrap refits of the SAME cohort (from \[assess_program_stability()\] refits). A ratio near 1 means cross-cohort alignment is as strong as the fit's own sampling noise ceiling; near 0 means the shared subspace is no larger than noise.

</div>

<div class="section level2">

## Examples<a href="#ref-examples" class="anchor" aria-label="anchor"></a>

``` r
W <- matrix(rnorm(100 * 3), 100, 3)
span_alignment(W, W)$mean_overlap          # 1: identical span
W2 <- W[, 1:2, drop = FALSE]
span_alignment(W, W2)$mean_overlap_topk    # subset span fully contained
```

</div>
