<div class="page-header">

# Plot program scores across donors

<div class="d-none name">

`plot_program_scores.Rd`

</div>

</div>

<div class="ref-description section level2">

Tile heatmap of the canonicalized program score matrix (donors x programs), optionally with donor-annotation strips (e.g. severity, batch) aligned underneath.

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
plot_program_scores(fit, cell_types = NULL, annotation = NULL, zscore = TRUE)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

fit<a href="#arg-fit" class="anchor" aria-label="anchor"></a><br>
A canonicalized \`cell_program_fit\` object.

cell_types<a href="#arg-cell-types" class="anchor" aria-label="anchor"></a><br>
Cell types to include (default: all with \>= 1 program).

annotation<a href="#arg-annotation" class="anchor" aria-label="anchor"></a><br>
Optional character vector of column names in \`fit\$sample_metadata\` to display as annotation strips, or a data.frame with an \`observation_id\` column plus the columns to display.

zscore<a href="#arg-zscore" class="anchor" aria-label="anchor"></a><br>
Rescale each program's scores to mean 0 / SD 1 before plotting (default TRUE; keeps the color scale comparable across programs).

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
plot_program_scores(fit)
# with an annotation strip built from a planted program's donor scores
sev <- sim$truth$B[[1]]$z > median(sim$truth$B[[1]]$z)
md <- data.frame(observation_id = paste0("obs", seq_along(sev)),
                 severity = factor(sev, c(FALSE, TRUE), c("low", "high")))
fit2 <- as_cell_program_data(sim$pseudobulk, sample_metadata = md)
plot_program_scores(canonicalize_programs(fit_celltype_programs(fit2, max_factors = 3, seed = 1)),
                    annotation = "severity")
```

</div>
