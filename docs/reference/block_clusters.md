# Coordination calls from a block_joint fit (benchmark output contract)

A shared factor (class "global" or "partial") with \>= 2 member cell
types is one coordination call; membership = block mass \>
\`mass_thresh\` (0.1 matches the classifier's partial rule). Private
factors are single-cell-type and never become calls.

## Usage

``` r
block_clusters(bfit, mass_thresh = 0.1)
```

## Arguments

- bfit:

  A \`block_joint\` object from \[fit_joint_block()\].

- mass_thresh:

  Block-mass threshold for member cell types.

## Value

data.frame(cluster, cell_type, program_id, stat, fwer_p) with one row
per member cell type, plus \`activity\`: donor x cluster state matrix.
