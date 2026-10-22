<div class="page-header">

# Numerically degenerate programs, two kinds: 1. flashier can retain factors whose loadings underflow to numerical noise (median \|loading\| many orders of magnitude below the cell type's strongest program; e.g. sparse factors whose support is a few genes and whose remaining entries sit at ~1e-16). Their posterior SDs underflow as well, producing garbage z-scores, and correlation-based matching is scale-invariant, so they spuriously match every strong program. 2. Fixed, never-updated factors: real-scale loadings whose posterior SD column collapsed to machine precision (observed for the smallest cell type in both real cohorts). Point estimates without usable uncertainty: mashr cannot assess them. A program is degenerate when median \|loading\| \< \`rel_tol\` x the cell type's largest, or when its median posterior SD (scaled to the loadings) is \<= \`sd_floor\`.

<div class="d-none name">

`dot-degenerate_programs.Rd`

</div>

</div>

<div class="ref-description section level2">

Numerically degenerate programs, two kinds: 1. flashier can retain factors whose loadings underflow to numerical noise (median \|loading\| many orders of magnitude below the cell type's strongest program; e.g. sparse factors whose support is a few genes and whose remaining entries sit at ~1e-16). Their posterior SDs underflow as well, producing garbage z-scores, and correlation-based matching is scale-invariant, so they spuriously match every strong program. 2. Fixed, never-updated factors: real-scale loadings whose posterior SD column collapsed to machine precision (observed for the smallest cell type in both real cohorts). Point estimates without usable uncertainty: mashr cannot assess them. A program is degenerate when median \|loading\| \< \`rel_tol\` x the cell type's largest, or when its median posterior SD (scaled to the loadings) is \<= \`sd_floor\`.

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
.degenerate_programs(fit, rel_tol = 1e-06, sd_floor = 1e-12)
```

</div>
