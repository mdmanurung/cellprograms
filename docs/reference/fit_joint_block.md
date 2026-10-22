# Block-structured joint EBMF (Level 3 candidate A)

Block-structured joint EBMF (Level 3 candidate A)

## Usage

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

## Arguments

- fit:

  A \`cell_program_fit\` (uses its processed matrices / scores).

- level:

  \`"expression"\` or \`"scores"\`.

- order:

  \`"shared_first"\` (default) or \`"private_first"\` (ablation).

- max_factors:

  Candidate factors for the joint pass.

- max_private:

  Candidate factors per cell type in the private pass.

- private_threshold:

  Block-mass threshold for private classification.

- global_entropy:

  Normalized-entropy threshold for global classification.

- seed:

  Random seed.

- verbose:

  Verbosity.

## Value

An integration-style list with method \`"block"\`.
