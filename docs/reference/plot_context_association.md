# Plot CMP context associations

Dot-and-whisker plot of the context-variable effect on each CMP's donor
scores (linear model from \`cmp_context_association\`), with 95 the
model standard errors.

## Usage

``` r
plot_context_association(assoc, p_alpha = 0.05)
```

## Arguments

- assoc:

  Output of \`cmp_context_association\` (data.frame with \`cmp\`,
  \`term\`, \`estimate\`, \`t\`, \`p\`, \`n\`).

- p_alpha:

  Significance level for solid vs open points (default 0.05).

## Value

A ggplot object.

## Examples

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
