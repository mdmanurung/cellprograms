# EB-CP coordination programs from a fitted \`cell_program_fit\`

Runs \[fit_eb_cp()\] on the adjusted expression matrices (\`Y_used\`) of
a \`cell_program_fit\`. If \`strata\` is missing, it is derived from the
covariates used at fit time (interaction of all covariate columns), so
the held-out split respects the adjustment covariates by default.

## Usage

``` r
eb_cp_from_fit(fit, strata, ...)
```

## Arguments

- fit:

  A \`cell_program_fit\` object (see \[fit_celltype_programs()\]).

- strata:

  Optional named donor strata overriding the default derived from the
  fit's covariates.

- ...:

  Forwarded to \[fit_eb_cp()\] (\`Kmax\`, \`heldout_z\`, \`seed\`, ...).
