<div class="page-header">

# Assess CMP stability under synchronized donor resampling

<div class="d-none name">

`assess_cmp_stability.Rd`

</div>

</div>

<div class="ref-description section level2">

Draws a donor subset once per replicate (optionally within donor-level strata), applies the identical donor set to all cell types, and reruns \[find_coordinated_programs()\] with the caller's configuration. Reports per-CMP membership recovery (Jaccard), assignment ambiguity, sign consistency, donor-score concordance, cell-type survival, and honest failure denominators. Replicates that error are ledgered and excluded from conditional rates; rates are \`NA\` when no valid replicate exists rather than an artificial zero.

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
assess_cmp_stability(
  fit,
  donor_col = NULL,
  visit_col = NULL,
  donor_strata = NULL,
  n_boot = 100L,
  subsample_frac = 0.8,
  jaccard_threshold = 0.5,
  cmp_args = list(),
  seed = NULL,
  progress = TRUE
)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

fit<a href="#arg-fit" class="anchor" aria-label="anchor"></a><br>
A canonicalized fit (see \[find_coordinated_programs()\]).

donor_col<a href="#arg-donor-col" class="anchor" aria-label="anchor"></a><br>
Optional \`sample_metadata\` column giving the biological donor when it differs from \`observation_ids\`. All observations of a drawn donor are kept together in every cell type where present (cluster subsampling); missing cell types/visits are never synthesized.

visit_col<a href="#arg-visit-col" class="anchor" aria-label="anchor"></a><br>
Optional \`sample_metadata\` column identifying the visit/timepoint within donor. Required when \`donor_col\` is given and any donor contributes more than one observation: the estimand is then the donor-synchronized draw retaining all of a donor's visits, and the per-visit index is preserved in replicate metadata.

donor_strata<a href="#arg-donor-strata" class="anchor" aria-label="anchor"></a><br>
Optional \`sample_metadata\` column with a donor-level stratification (must be constant within donor). The draw is taken within each stratum separately.

n_boot<a href="#arg-n-boot" class="anchor" aria-label="anchor"></a><br>
Number of resampling replicates (default 100).

subsample_frac<a href="#arg-subsample-frac" class="anchor" aria-label="anchor"></a><br>
Fraction of donors drawn per replicate, without replacement (default 0.8).

jaccard_threshold<a href="#arg-jaccard-threshold" class="anchor" aria-label="anchor"></a><br>
Minimum membership Jaccard between a reference CMP and its best-matching replicate CMP to count as recovered (default 0.5).

cmp_args<a href="#arg-cmp-args" class="anchor" aria-label="anchor"></a><br>
Named list of arguments passed to \[find_coordinated_programs()\] unchanged (the predeclared caller configuration), e.g. \`list(strata = "site", n_perm = 199, tau = 2, link_threshold = 3)\`. \`seed\` is always overridden by the per-replicate seed recorded in the replicate log.

seed<a href="#arg-seed" class="anchor" aria-label="anchor"></a><br>
Base seed for reproducibility. Replicate \`b\` uses seed \`seed + b - 1\` for both the donor draw and the CMP permutation null.

progress<a href="#arg-progress" class="anchor" aria-label="anchor"></a><br>
Print progress every 10 replicates.

</div>

<div class="section level2">

## Value<a href="#value" class="anchor" aria-label="anchor"></a>

An object of class \`cmp_stability\`, a list with

summary<br>
one row per reference CMP: \`cmp\`, \`n_members\`, \`members\`, \`valid_refits\`, \`recovered_n\`, \`recovery_freq\` (conditional on valid replicates; \`NA\` when none), \`mean_jaccard\`, \`mean_ambiguity_margin\` (best minus runner-up Jaccard), \`sign_concordance\` (fraction of recovered replicates whose member signs match the reference orientation up to one global flip), \`donor_score_concordance\` (mean correlation of CMP donor scores on shared donors among recovered replicates), \`many_to_one_n\` (replicates where this reference CMP shared its best-matching replicate CMP with another reference CMP).

members<br>
one row per reference CMP member: recovery frequency in the matched replicate CMP and sign concordance given presence.

cell_type_survival<br>
per cell type, fraction of valid replicates passing the \`min_donors\` filter.

replicates<br>
per-replicate ledger: \`replicate\`, \`seed\`, \`status\`, \`error_class\`, \`error_message\`, \`elapsed_seconds\`, \`n_donors_drawn\`, \`n_obs\`, \`n_cmps\`, \`eligible_rows\` (per cell type), \`surviving_cell_types\`.

reference<br>
the reference \[find_coordinated_programs()\] result.

settings<br>
the arguments used.

When the reference fit yields no CMPs, stability is undefined: the function returns an empty \`summary\`, zero-row \`replicates\`, and \`status = "no_reference_cmps"\` without spending replicates.

</div>

<div class="section level2">

## Identity levels<a href="#identity-levels" class="anchor" aria-label="anchor"></a>

\`observation_ids\` are the score rows (observation level). \`donor_col\` promotes resampling to the donor level; \`visit_col\` disambiguates repeated observations per donor. Cohort/stratification effects on the null remain the responsibility of \`cmp_args\$strata\`.

</div>

<div class="section level2">

## Examples<a href="#ref-examples" class="anchor" aria-label="anchor"></a>

``` r
set.seed(1)
ids <- paste0("d", seq_len(40))
z <- rnorm(40)
mk <- function(ct) {
  S <- matrix(rnorm(40 * 2), 40, 2, dimnames = list(ids, paste0(ct, "_", 1:2)))
  S[, 1] <- z + rnorm(40, sd = 0.3)
  S
}
fit <- structure(list(
  scores = list(A = mk("A"), B = mk("B")),
  loadings = list(A = matrix(0, 5, 2), B = matrix(0, 5, 2)),
  cell_types = c("A", "B"), observation_ids = ids,
  sample_metadata = data.frame(observation_id = ids)),
  class = c("cell_program_fit", "cell_program_data", "list"))
stab <- assess_cmp_stability(fit, n_boot = 5, subsample_frac = 0.8,
                             cmp_args = list(n_perm = 99, ambient = "none",
                                             min_shared = 10, min_donors = 10))
stab$summary
```

</div>
