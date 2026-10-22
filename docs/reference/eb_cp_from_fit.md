<div class="page-header">

# EB-CP coordination programs from a fitted \`cell_program_fit\`

<div class="d-none name">

`eb_cp_from_fit.Rd`

</div>

</div>

<div class="ref-description section level2">

Runs \[fit_eb_cp()\] on the adjusted expression matrices (\`Y_used\`) of a \`cell_program_fit\`. If \`strata\` is missing, it is derived from the categorical covariates used at fit time. Automatic strata from numeric covariates are refused: explicitly choose a categorical split strategy rather than silently generating singleton age strata.

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
eb_cp_from_fit(fit, strata, ...)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

fit<a href="#arg-fit" class="anchor" aria-label="anchor"></a><br>
A \`cell_program_fit\` object (see \[fit_celltype_programs()\]).

strata<a href="#arg-strata" class="anchor" aria-label="anchor"></a><br>
Optional named donor strata overriding the default derived from the fit's covariates.

...<a href="#arg--" class="anchor" aria-label="anchor"></a><br>
Forwarded to \[fit_eb_cp()\] (\`Kmax\`, \`heldout_z\`, \`seed\`, ...).

</div>
