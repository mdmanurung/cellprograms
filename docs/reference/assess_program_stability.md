# Assess program stability by donor subsampling

Assess program stability by donor subsampling

## Usage

``` r
assess_program_stability(
  fit,
  n_boot = 100L,
  match_method = c("subspace", "correlation"),
  sim_threshold = 0.5,
  top_n = 20L,
  resample_frac = 0.8,
  seed = NULL,
  progress = TRUE
)
```

## Arguments

- fit:

  A canonicalized \`cell_program_fit\` (refit uses its stored
  preprocessing settings, gene sets and covariates).

- n_boot:

  Number of subsample refits.

- match_method:

  Matching method for \`match_programs\`.

- sim_threshold:

  Minimum \|cosine\| for a factor-level match.

- top_n:

  Top-gene set size for Jaccard recovery.

- resample_frac:

  Fraction of donors kept per cell type per refit.

- seed:

  Random seed.

- progress:

  Print progress messages.

## Value

The fit with \`fit\$stability\` populated: per-program recovery
frequency, median score correlation, median loading cosine, sign
consistency, and top-gene Jaccard.

## Examples

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- fit_celltype_programs(x, max_factors = 3, features = "variable", n_variable_genes = 150)
fit <- canonicalize_programs(fit)
# few subsamples for a quick example; use more (e.g. 100) in practice
fit <- assess_program_stability(fit, n_boot = 2, seed = 1, progress = FALSE)
fit$stability$per_program
#>   program_id cell_type recovery_freq median_score_cor median_loading_cos
#> 1        B_1         B             1        0.9999650          0.9971983
#> 2        B_2         B             1        0.9997030          0.9870352
#> 3      CD4_1       CD4             1        0.9998824          0.9921341
#> 4      CD8_1       CD8             1        0.9993324          0.9828631
#> 5       NK_1        NK             1        0.9999204          0.9946079
#> 6       NK_2        NK             1        0.9999242          0.9946153
#> 7     Mono_1      Mono             1        0.9999864          0.9985128
#> 8     Mono_2      Mono             1        0.9998608          0.9920544
#>   sign_consistency top_gene_jaccard
#> 1                1        0.9047619
#> 2                1        0.8614719
#> 3                1        0.8614719
#> 4                1        0.9090909
#> 5                1        0.9523810
#> 6                1        0.9090909
#> 7                1        0.9523810
#> 8                1        0.7424242
```
