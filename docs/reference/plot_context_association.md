<div class="page-header">

# Plot CMP context associations

<div class="d-none name">

`plot_context_association.Rd`

</div>

</div>

<div class="ref-description section level2">

Dot-and-whisker plot of the context-variable effect on each CMP's donor scores (linear model from \`cmp_context_association\`), with 95 the model standard errors.

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
plot_context_association(assoc, p_alpha = 0.05)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

assoc<a href="#arg-assoc" class="anchor" aria-label="anchor"></a><br>
Output of \`cmp_context_association\` (data.frame with \`cmp\`, \`term\`, \`estimate\`, \`t\`, \`p\`, \`n\`).

p_alpha<a href="#arg-p-alpha" class="anchor" aria-label="anchor"></a><br>
Significance level for solid vs open points (default 0.05).

</div>

<div class="section level2">

## Value<a href="#value" class="anchor" aria-label="anchor"></a>

A ggplot object.

</div>

<div class="section level2">

## Examples<a href="#ref-examples" class="anchor" aria-label="anchor"></a>

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 3, seed = 1))
cmp <- find_coordinated_programs(fit, n_perm = 200, seed = 1)
# context variable: high/low donor score of the first planted program
z <- sim$truth$B[[1]]$z > median(sim$truth$B[[1]]$z)
ctx <- setNames(factor(z, c(FALSE, TRUE), c("low", "high")), paste0("obs", seq_along(z)))
plot_context_association(cmp_context_association(fit, cmp, ctx))
```

</div>
