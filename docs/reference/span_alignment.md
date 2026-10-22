# Principal-angle alignment between two loading subspaces

Computes the principal angles between the column spans of two loading
matrices (genes x programs) on the same gene panel — e.g. a discovery
fit vs a stability refit, or the same cohort fit by two methods. Each
span is orthonormalized (rank-aware via SVD), and the singular values of
the cross-product of the bases are the cosines of the principal angles
between the spans.

## Usage

``` r
span_alignment(W1, W2, k = NULL)
```

## Arguments

- W1, :

  W2 Loading matrices (genes x programs) on the same gene panel. Row
  names should match; rows are matched by name when both have row names,
  otherwise by position.

- k:

  Optional number of leading principal cosines to average for
  \`mean_overlap_topk\` (default: all).

## Value

A list with \`principal_cosines\` (descending cosines of the principal
angles, length = min(rank1, rank2)), \`mean_overlap\`,
\`mean_overlap_topk\`, and \`dim\` (the two span ranks).

## Details

The registered tier-2 use (immune-aging application, 2026-10-16) is the
\*\*stability-normalized conservation ratio\*\*: the between-cohort span
overlap (\`mean_overlap\`) divided by the within-cohort bootstrap
ceiling — the same quantity computed between discovery and bootstrap
refits of the SAME cohort (from \[assess_program_stability()\] refits).
A ratio near 1 means cross-cohort alignment is as strong as the fit's
own sampling noise ceiling; near 0 means the shared subspace is no
larger than noise.

## Examples

``` r
W <- matrix(rnorm(100 * 3), 100, 3)
span_alignment(W, W)$mean_overlap          # 1: identical span
#> [1] 1
W2 <- W[, 1:2, drop = FALSE]
span_alignment(W, W2)$mean_overlap_topk    # subset span fully contained
#> [1] 1
```
