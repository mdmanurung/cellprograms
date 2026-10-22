# Post-hoc shared-vs-specific assessment of fitted programs (mashr)

Assesses which fitted programs are the SAME biological program across
cell types (shared gene loadings) versus cell-type-private, using
mashr-style adaptive shrinkage: for each matched program group (one
program per cell type) the G x C matrix of canonicalized loadings and
their posterior SDs go to \[mashr::mash()\] with the canonical plus
data-driven (PCA/ED) covariance set and the CorShrink-style
null-correlation matrix \[\`.estimate_V\`\]. The estimated mixture
weights over covariance patterns ARE the empirical sharing patterns
(e.g. \`equal_effects\` = shared across all conditions,
\`singletons\_\*\` = one cell type only).

## Usage

``` r
program_sharing(
  fit,
  matches = NULL,
  alpha = 0.05,
  n_top = 50L,
  lfsr_thresh = 0.05,
  n_perm = 199L,
  seed = 1L,
  min_genes = 50L,
  min_null = 20L,
  z_thresh = 1,
  data_driven = TRUE,
  npc = 3L,
  verbose = FALSE,
  ...
)
```

## Arguments

- fit:

  A canonicalized \`cell_program_fit\` with stored loading posterior SDs
  (fits from older versions must be refit).

- matches:

  \`NULL\` (default): run \[\`.match_programs_loadings\`\] and group
  matches significant at FDR \`alpha\`. Alternatively a pairwise match
  table (columns ct_a/prog_a/ct_b/prog_b, optionally q), one matched
  set, or a list of matched sets (data.frames with columns cell_type,
  program_id).

- alpha:

  FDR level for candidate matches and the strong-gene subset for
  data-driven covariances (default 0.05).

- n_top:

  Number of top-\|loading\| genes per group for the summary fractions
  (default 50). Choose relative to the smallest shared support you want
  to detect: a group sharing S genes can reach \`frac_top_shared\` at
  most S/\`n_top\`, so with the 0.5 "shared" threshold any support
  smaller than \`n_top\`/2 (25 at the default) classifies at best
  "partial" — never "shared" — no matter how clean the signal
  (confirm-scale evidence: planted 10-gene sex/genotype-shared programs
  were never promoted in any arm; CP-4 Amendment 1, limitation register
  L3). For narrow shared modules (VDJ/sex/genotype scale, ~10-30 genes),
  lower \`n_top\` accordingly, and report \`frac_top_shared\`,
  \`frac_top_private\` and \`p_all_shared\` directly rather than the
  ternary label alone.

- lfsr_thresh:

  lfsr threshold calling a gene active in a condition (default 0.05).

- n_perm, seed:

  Permutation calibration for the matching step. CAND-07 caveat: the
  gene-permutation null treats genes as exchangeable; block-correlated
  (co-regulated) genes violate that assumption, so matching q-values are
  anticonservative when a large correlated gene block drives both
  programs (audit D-1: 67.0 equicorrelation 0.55 vs 2.0 registered
  calibration benchmark; interpret marginal matches accordingly.

- min_genes:

  Minimum shared-panel size per group (default 50).

- min_null:

  Minimum null genes for V estimation (else identity).

- z_thresh:

  \|z\| below this in every condition marks a null gene for V estimation
  (default 1). With many conditions the intersection at 1 can be empty;
  raise to ~1.5 for large matched groups.

- data_driven:

  Add PCA/ED data-driven covariances to the canonical set (default TRUE;
  skipped with a warning when too few strong genes).

- npc:

  Number of PCA covariance components (capped at \#conditions).

- verbose:

  Let mashr print progress (default FALSE).

- ...:

  Passed to \[mashr::mash()\].

## Value

An object of class \`program_sharing\`: \* \`matches\`: the pairwise
candidate-match table (NULL for explicit groups); \* \`groups\`: list of
matched-set data.frames (cell_type, program_id); \* \`results\`:
per-group list with \`B_hat\`/\`S_hat\` (mash inputs), \`V\`,
\`posterior_mean\` and \`lfsr\` (G x C), \`mixture_weights\` (mashr's pi
over the scaled covariance grid), \`pattern_weights\` (pi aggregated
over base covariance patterns), \`p_all_shared\` (mean posterior weight
of the equal_effects pattern over the top genes), \`frac_top_shared\` /
\`frac_top_private\` (fraction of the top \`n_top\` genes with lfsr \<
\`lfsr_thresh\` in \>= 2 / exactly 1 conditions), \`classification\`
("shared" if frac_top_shared \>= 0.5, "private" if frac_top_private \>=
0.5, else "partial"). D-2a caveat: "shared" means active-in-both, not
proven shared structure — at high loading SNR, independent loadings also
reach frac_top_shared \>= 0.5 (forced-pair diagnostic: 30/30 false
"shared" at loading SE = 0.125 while mashr's mixture put 96 the identity
pattern, p_all_shared ~ 0.001). Check \`p_all_shared\`. Each per-group
result also carries the full \`mash\` fit; \* \`summary\`: one row per
program. Programs in no matched group are classified "private" with
evidence "no_significant_match".

## Details

This concerns shared GENE LOADINGS only. Coordinated DONOR ACTIVITY (CMP
detection, \[find_coordinated_programs()\]) is a separate axis and a
separate workstream; a shared program need not have correlated donor
scores.

## Examples

``` r
sim <- simulate_pseudobulk(n_samples = 20, n_genes = 150, program_size = 20, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- fit_celltype_programs(x, max_factors = 3, features = "variable", n_variable_genes = 100)
fit <- canonicalize_programs(fit)
# Assess one matched set explicitly (skips the permutation matching step);
# with matches = NULL, candidate matches are data-driven.
if (requireNamespace("mashr", quietly = TRUE)) {
  grp <- data.frame(cell_type = c("CD8", "NK"), program_id = c("CD8_1", "NK_1"))
  res <- program_sharing(fit, matches = grp, n_top = 20)
  res$summary
}
#> Warning: `matches` supplied: matched sets are user-forced, not data-driven. The 'shared' classification counts top genes active (lfsr < 0.05) in >= 2 conditions and can label independent loadings 'shared' at high loading SNR; check `p_all_shared` (mixture weight of the equal-effects pattern) before interpretation.
#>  - Computing 73 x 241 likelihood matrix.
#>  - Likelihood calculations took 0.00 seconds.
#>  - Fitting model with 241 mixture components.
#>  - Model fitting took 0.06 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.00 seconds.
#>   program_id cell_type group_id n_conditions classification frac_top_shared
#> 1        B_1         B       NA            1        private              NA
#> 2        B_2         B       NA            1        private              NA
#> 3      CD4_1       CD4       NA            1        private              NA
#> 4      CD8_1       CD8        1            2        private            0.05
#> 5       NK_1        NK        1            2        private            0.05
#> 6       NK_2        NK       NA            1        private              NA
#> 7     Mono_1      Mono       NA            1        private              NA
#> 8     Mono_2      Mono       NA            1        private              NA
#>   frac_top_private p_all_shared             evidence
#> 1               NA           NA no_significant_match
#> 2               NA           NA no_significant_match
#> 3               NA           NA no_significant_match
#> 4              0.6            0                mashr
#> 5              0.6            0                mashr
#> 6               NA           NA no_significant_match
#> 7               NA           NA no_significant_match
#> 8               NA           NA no_significant_match
```
