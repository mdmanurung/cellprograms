<div class="page-header">

# Block-structured joint EBMF (Level 3 candidate A)

<div class="d-none name">

`fit_joint_block.Rd`

</div>

</div>

<div class="ref-description section level2">

Block-structured joint EBMF (Level 3 candidate A)

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
fit_joint_block(
  fit,
  level = c("scores", "expression"),
  order = c("shared_first", "private_first"),
  max_factors = 20L,
  max_private = 10L,
  private_threshold = 0.8,
  global_entropy = 0.8,
  seed = NULL,
  verbose = 0L
)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

fit<a href="#arg-fit" class="anchor" aria-label="anchor"></a><br>
A \`cell_program_fit\` (uses its processed matrices / scores).

level<a href="#arg-level" class="anchor" aria-label="anchor"></a><br>
\`"expression"\` or \`"scores"\`.

order<a href="#arg-order" class="anchor" aria-label="anchor"></a><br>
\`"shared_first"\` (default) or \`"private_first"\` (expression level only). Scores/private-first is refused before fitting: subtracting the stage-1 fit from its own scores has no defined residual semantics.

max_factors<a href="#arg-max-factors" class="anchor" aria-label="anchor"></a><br>
Candidate factors for the joint pass.

max_private<a href="#arg-max-private" class="anchor" aria-label="anchor"></a><br>
Candidate factors per cell type in the private pass.

private_threshold<a href="#arg-private-threshold" class="anchor" aria-label="anchor"></a><br>
Block-mass threshold for private classification.

global_entropy<a href="#arg-global-entropy" class="anchor" aria-label="anchor"></a><br>
Normalized-entropy threshold for global classification.

seed<a href="#arg-seed" class="anchor" aria-label="anchor"></a><br>
Random seed.

verbose<a href="#arg-verbose" class="anchor" aria-label="anchor"></a><br>
Verbosity.

</div>

<div class="section level2">

## Value<a href="#value" class="anchor" aria-label="anchor"></a>

An integration-style list with method \`"block"\`.

</div>
