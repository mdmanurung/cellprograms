<div class="page-header">

# Coordination calls from a block_joint fit (benchmark output contract)

<div class="d-none name">

`block_clusters.Rd`

</div>

</div>

<div class="ref-description section level2">

A shared factor (class "global" or "partial") with \>= 2 member cell types is one coordination call; membership = block mass \> \`mass_thresh\` (0.1 matches the classifier's partial rule). Private factors are single-cell-type and never become calls.

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
block_clusters(bfit, mass_thresh = 0.1)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

bfit<a href="#arg-bfit" class="anchor" aria-label="anchor"></a><br>
A \`block_joint\` object from \[fit_joint_block()\].

mass_thresh<a href="#arg-mass-thresh" class="anchor" aria-label="anchor"></a><br>
Block-mass threshold for member cell types.

</div>

<div class="section level2">

## Value<a href="#value" class="anchor" aria-label="anchor"></a>

data.frame(cluster, cell_type, program_id, stat, fwer_p) with one row per member cell type, plus \`activity\`: donor x cluster state matrix.

</div>
