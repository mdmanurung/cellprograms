# Plot CMP membership across settings

Alluvial view of how each program's coordination membership changes
across a named list of \`cmp_result\` objects — typically the same fit
run with different \`find_coordinated_programs()\` settings, e.g.
\`ambient = "drop"\` versus \`ambient = "flag"\`. Strata are membership
classes ("CMP k", "significant, no CMP", "not significant", "excluded
(ambient)"); ribbons are programs, colored by their destination class,
so a setting change that calls or drops programs shows up as a colored
ribbon leaving the grey baseline.

## Usage

``` r
plot_cmp_membership(cmp_list, alpha = 0.05, label_programs = FALSE)
```

## Arguments

- cmp_list:

  A named list of \`cmp_result\` objects; the names become the x-axis
  labels, in order.

- alpha:

  Membership significance level used to classify programs (default 0.05;
  each result's own \`settings\$alpha\` takes precedence).

- label_programs:

  When TRUE and at most 12 programs are present, annotate ribbons whose
  membership changed with the program IDs carried.

## Value

A ggplot object.

## Examples

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
