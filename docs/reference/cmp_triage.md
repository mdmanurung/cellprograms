# One-row-per-CMP triage table

Collects the checks needed before interpreting a CMP: the recurrent top
genes of its members and how many of them fall in nuisance classes
(\`cmp_nuisance_genes()\`), the fraction of ambient-flagged members, an
optional association with a donor variable, optional sign replication of
member-pair co-variation within levels of a grouping variable (e.g.
site), and optional TF footprints from \`annotate_cmp_tf()\`.

## Usage

``` r
cmp_triage(
  fit,
  cmp,
  context = NULL,
  covariates = NULL,
  replicate_by = NULL,
  tf = NULL,
  nuisance = cmp_nuisance_genes(),
  top_n = 15L,
  n_recurrent = 12L,
  max_ambient_frac = 0.5,
  max_vdj_frac = 0.5,
  min_marker_genes = 2L,
  min_concordance = 0.75,
  min_level_donors = 15L,
  min_ieg = 3L
)
```

## Arguments

- fit:

  A canonicalized \`cell_program_fit\`.

- cmp:

  A \`cmp_result\` or its \`cmps\` table (needs \`ambient_flag\` for the
  ambient check; use \`ambient = "flag"\` in
  \`find_coordinated_programs()\`).

- context:

  Optional donor variable (column name or named vector) for
  \`cmp_context_association()\`.

- covariates:

  Optional \`sample_metadata\` columns to adjust for.

- replicate_by:

  Optional grouping variable (column name or named vector), e.g. site.
  Within each level with \>= \`min_level_donors\` donors, member-pair
  Pearson correlations of CMP-oriented scores are computed; concordance
  is the fraction that are positive.

- tf:

  Optional \`cmp_tf\` result.

- nuisance:

  Named regex vector (default \`cmp_nuisance_genes()\`).

- top_n, n_recurrent:

  Top genes per member (15) and recurrent genes kept (12).

- max_ambient_frac, max_vdj_frac, min_marker_genes, min_concordance:

  Flag thresholds.

- min_level_donors:

  Minimum donors per \`replicate_by\` level (default 15).

- min_ieg:

  Minimum number of immediate-early genes among the recurrent top genes
  for the \`stress_ieg\` flag (default 3).

## Value

Data frame, one row per CMP.

## Details

Recurrent top genes: the \`top_n\` genes with the largest CMP-oriented
loadings in each member, ranked by how many members share them; the
first \`n_recurrent\` are kept. Flags (all thresholds are arguments):
\`ambient\` (ambient fraction \>= \`max_ambient_frac\`),
\`platelet_rbc\` (\>= \`min_marker_genes\` contamination markers),
\`vdj\` (V(D)J fraction of recurrent genes \>= \`max_vdj_frac\`),
\`sex_chromosome\` and \`genotype_candidate\` (\>= \`min_marker_genes\`
marker genes), \`stress_ieg\` (\>= \`min_ieg\` IEGs), plus soft notes
\`genotype_marker_single\` (one marker gene) and \`low_replication\`
(minimum within-level concordance \< \`min_concordance\`). \`label\` is
the first hard flag in that order, else \`"candidate"\`, or
\`"candidate_associated"\` when the context association has p \< 0.05.
Labels are heuristics for prioritization, not classifications.

## Examples

``` r
# Toy canonicalized fit: program 1 of cell types A and B shares donor
# activity z (a planted coordinated program) and top genes G1-G15.
set.seed(1)
ids <- paste0("d", seq_len(40))
genes <- paste0("G", seq_len(100))
z <- rnorm(40)
mk <- function(ct) {
  S <- matrix(rnorm(40 * 2), 40, 2, dimnames = list(ids, paste0(ct, "_", 1:2)))
  S[, 1] <- z + rnorm(40, sd = 0.3)
  W <- matrix(rnorm(length(genes) * 2, sd = 0.1), length(genes), 2,
              dimnames = list(genes, colnames(S)))
  W[1:15, 1] <- W[1:15, 1] + 1
  list(S = S, W = W)
}
ta <- mk("A"); tb <- mk("B")
fit <- structure(list(
  scores = list(A = ta$S, B = tb$S), loadings = list(A = ta$W, B = tb$W),
  cell_types = c("A", "B"), observation_ids = ids,
  sample_metadata = data.frame(observation_id = ids,
                               group = rep(c("ctrl", "stim"), each = 20))),
  class = c("cell_program_fit", "cell_program_data", "list"))
res <- find_coordinated_programs(fit, n_perm = 99, ambient = "none")
cmp_triage(fit, res, context = "group", replicate_by = "group")
#>   cmp n_members cell_types members                       recurrent_top_genes
#> 1   1         2        A,B A_1,B_1 G1,G10,G11,G12,G13,G14,G15,G2,G3,G4,G5,G6
#>   ambient_frac vdj_frac n_genotype_markers genotype_markers n_sex_markers
#> 1           NA        0                  0                              0
#>   sex_markers n_ieg n_platelet_rbc context_p context_estimate min_concordance
#> 1                 0              0  0.726823       -0.1092466               1
#>   concordance_by_level tf_specific tf_shared flags     label
#> 1  ctrl:1.00;stim:1.00        <NA>      <NA>       candidate
```
