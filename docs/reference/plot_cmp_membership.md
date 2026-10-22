<div class="page-header">

# Plot CMP membership across settings

<div class="d-none name">

`plot_cmp_membership.Rd`

</div>

</div>

<div class="ref-description section level2">

Alluvial view of how each program's coordination membership changes across a named list of \`cmp_result\` objects — typically the same fit run with different \`find_coordinated_programs()\` settings, e.g. \`ambient = "drop"\` versus \`ambient = "flag"\`. Strata are membership classes ("CMP k", "significant, no CMP", "not significant", "excluded (ambient)"); ribbons are programs, colored by their destination class, so a setting change that calls or drops programs shows up as a colored ribbon leaving the grey baseline.

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
plot_cmp_membership(cmp_list, alpha = 0.05, label_programs = FALSE)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

cmp_list<a href="#arg-cmp-list" class="anchor" aria-label="anchor"></a><br>
A named list of \`cmp_result\` objects; the names become the x-axis labels, in order.

alpha<a href="#arg-alpha" class="anchor" aria-label="anchor"></a><br>
Membership significance level used to classify programs (default 0.05; each result's own \`settings\$alpha\` takes precedence).

label_programs<a href="#arg-label-programs" class="anchor" aria-label="anchor"></a><br>
When TRUE and at most 12 programs are present, annotate ribbons whose membership changed with the program IDs carried.

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
cmps <- list(
  drop = find_coordinated_programs(fit, n_perm = 200, seed = 1),
  flag = find_coordinated_programs(fit, n_perm = 200, seed = 1, ambient = "flag")
)
plot_cmp_membership(cmps)
```

</div>
