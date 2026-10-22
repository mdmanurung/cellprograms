<div class="page-header">

# Association of coordinated programs with a donor variable

<div class="d-none name">

`cmp_context_association.Rd`

</div>

</div>

<div class="ref-description section level2">

Regresses each CMP's consensus donor score (\`cmp_scores()\`) on a donor-level variable, optionally adjusting for covariates.

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
cmp_context_association(fit, cmp, variable, covariates = NULL)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

fit<a href="#arg-fit" class="anchor" aria-label="anchor"></a><br>
A canonicalized \`cell_program_fit\`.

cmp<a href="#arg-cmp" class="anchor" aria-label="anchor"></a><br>
A \`cmp_result\` or its \`cmps\` table.

variable<a href="#arg-variable" class="anchor" aria-label="anchor"></a><br>
A \`sample_metadata\` column name or a vector named by observation ID.

covariates<a href="#arg-covariates" class="anchor" aria-label="anchor"></a><br>
Optional character vector of \`sample_metadata\` columns to adjust for.

</div>

<div class="section level2">

## Value<a href="#value" class="anchor" aria-label="anchor"></a>

Data frame: \`cmp\`, \`term\`, \`estimate\`, \`t\`, \`p\`, \`n\`.

</div>
