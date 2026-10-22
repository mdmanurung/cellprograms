# Association of coordinated programs with a donor variable

Regresses each CMP's consensus donor score (\`cmp_scores()\`) on a
donor-level variable, optionally adjusting for covariates.

## Usage

``` r
cmp_context_association(fit, cmp, variable, covariates = NULL)
```

## Arguments

- fit:

  A canonicalized \`cell_program_fit\`.

- cmp:

  A \`cmp_result\` or its \`cmps\` table.

- variable:

  A \`sample_metadata\` column name or a vector named by observation ID.

- covariates:

  Optional character vector of \`sample_metadata\` columns to adjust
  for.

## Value

Data frame: \`cmp\`, \`term\`, \`estimate\`, \`t\`, \`p\`, \`n\`.
