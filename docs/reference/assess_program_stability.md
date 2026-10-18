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

The fit with \`fit\$stability\` populated:

- per_program:

  One row per program: \`n_hit\` (refits whose match recovered the
  program), \`n_eligible\` (successful refits in which the program's
  cell type produced at least one factor, i.e. matching was possible),
  \`recovery_freq\` (\*\*legacy\*\*: unconditional hit fraction over all
  REQUESTED iterations, \`n_hit / n_boot\` — not an estimate of
  biological recovery conditional on successful refitting),
  \`conditional_recovery_freq\` (\`n_hit / valid_refits\`, NA when no
  refit succeeded), median score correlation, median loading cosine,
  sign consistency, and top-gene Jaccard.

- valid_refits, failed_refits:

  Refit denominator accounting. Failed refits are never silently
  dropped: any failure warns, and an all-fail run yields NA rates
  (undefined), never an artificial 0.

- replicate_log:

  One row per requested refit: \`replicate\`, \`seed\`, \`status\`
  ("ok"/"error"), \`error_class\`/\`error_message\`,
  \`elapsed_seconds\`, \`n_obs\` (sampled observations per cell type),
  \`n_factors\` (refit factor count per cell type), and \`eligible\`
  (whether matching could be attempted). A successful refit with no
  matched factor is distinguishable from a failed refit.

## Examples

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- fit_celltype_programs(x, max_factors = 3, features = "variable", n_variable_genes = 150)
fit <- canonicalize_programs(fit)
# few subsamples for a quick example; use more (e.g. 100) in practice
fit <- assess_program_stability(fit, n_boot = 2, seed = 1, progress = FALSE)
fit$stability$per_program
#>   program_id cell_type n_hit n_eligible recovery_freq conditional_recovery_freq
#> 1        B_1         B     2          2             1                         1
#> 2        B_2         B     2          2             1                         1
#> 3      CD4_1       CD4     2          2             1                         1
#> 4      CD8_1       CD8     2          2             1                         1
#> 5       NK_1        NK     2          2             1                         1
#> 6       NK_2        NK     2          2             1                         1
#> 7     Mono_1      Mono     2          2             1                         1
#> 8     Mono_2      Mono     2          2             1                         1
#>   median_score_cor median_loading_cos sign_consistency top_gene_jaccard
#> 1        0.9999650          0.9971983                1        0.9047619
#> 2        0.9997030          0.9870352                1        0.8614719
#> 3        0.9998824          0.9921341                1        0.8614719
#> 4        0.9993324          0.9828631                1        0.9090909
#> 5        0.9999204          0.9946079                1        0.9523810
#> 6        0.9999242          0.9946153                1        0.9090909
#> 7        0.9999864          0.9985128                1        0.9523810
#> 8        0.9998608          0.9920544                1        0.7424242
```
