<div class="page-header">

# Assess program stability by donor subsampling

<div class="d-none name">

`assess_program_stability.Rd`

</div>

</div>

<div class="ref-description section level2">

Assess program stability by donor subsampling

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

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

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

fit<a href="#arg-fit" class="anchor" aria-label="anchor"></a><br>
A canonicalized \`cell_program_fit\` (refit uses its stored preprocessing settings, gene sets and covariates).

n_boot<a href="#arg-n-boot" class="anchor" aria-label="anchor"></a><br>
Number of subsample refits.

match_method<a href="#arg-match-method" class="anchor" aria-label="anchor"></a><br>
Matching method for \`match_programs\`.

sim_threshold<a href="#arg-sim-threshold" class="anchor" aria-label="anchor"></a><br>
Minimum \|cosine\| for a factor-level match.

top_n<a href="#arg-top-n" class="anchor" aria-label="anchor"></a><br>
Top-gene set size for Jaccard recovery.

resample_frac<a href="#arg-resample-frac" class="anchor" aria-label="anchor"></a><br>
Fraction of donors kept per cell type per refit.

seed<a href="#arg-seed" class="anchor" aria-label="anchor"></a><br>
Random seed.

progress<a href="#arg-progress" class="anchor" aria-label="anchor"></a><br>
Print progress messages.

</div>

<div class="section level2">

## Value<a href="#value" class="anchor" aria-label="anchor"></a>

The fit with \`fit\$stability\` populated:

per_program<br>
One row per program: \`n_hit\` (refits whose match recovered the program), \`n_eligible\` (successful refits in which the program's cell type produced at least one factor, i.e. matching was possible), \`recovery_freq\` (\*\*legacy\*\*: unconditional hit fraction over all REQUESTED iterations, \`n_hit / n_boot\` — not an estimate of biological recovery conditional on successful refitting), \`conditional_recovery_freq\` (\`n_hit / valid_refits\`, NA when no refit succeeded; the denominator counts successful refits, not the replicates where the program's cell type retained a factor — use \`n_hit / n_eligible\` for the per-program conditional rate), median score correlation, median loading cosine, sign consistency, and top-gene Jaccard.

valid_refits, failed_refits<br>
Refit denominator accounting. Failed refits are never silently dropped: any failure warns, and an all-fail run yields NA rates (undefined), never an artificial 0.

replicate_log<br>
One row per requested refit: \`replicate\`, \`seed\`, \`status\` ("ok"/"error"), \`error_class\`/\`error_message\`, \`elapsed_seconds\`, \`n_obs\` (sampled observations per cell type), \`n_factors\` (refit factor count per cell type), and \`eligible\` (whether matching could be attempted). A successful refit with no matched factor is distinguishable from a failed refit. \`refit_settings\` records the intended estimator arguments, including backend/noise controls, for every success and failure. Frozen specs govern preprocessing and units; legacy fits use their stored units. \`effective_refit\` retains actual per-type effective settings, fallback, known noise and backend convergence when available; NULL after failure. Stability targets the full fitting procedure including documented fallback behavior, not a fixed selected estimator.

</div>

<div class="section level2">

## Examples<a href="#ref-examples" class="anchor" aria-label="anchor"></a>

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- fit_celltype_programs(x, max_factors = 3, features = "variable", n_variable_genes = 150)
fit <- canonicalize_programs(fit)
# few subsamples for a quick example; use more (e.g. 100) in practice
fit <- assess_program_stability(fit, n_boot = 2, seed = 1, progress = FALSE)
fit$stability$per_program
```

</div>
