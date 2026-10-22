<div class="page-header">

# One-row-per-CMP triage table

<div class="d-none name">

`cmp_triage.Rd`

</div>

</div>

<div class="ref-description section level2">

Collects the checks needed before interpreting a CMP: the recurrent top genes of its members and how many of them fall in nuisance classes (\`cmp_nuisance_genes()\`), the fraction of ambient-flagged members, an optional association with a donor variable, optional sign replication of member-pair co-variation within levels of a grouping variable (e.g. site), and optional TF footprints from \`annotate_cmp_tf()\`.

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

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

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

fit<a href="#arg-fit" class="anchor" aria-label="anchor"></a><br>
A canonicalized \`cell_program_fit\`.

cmp<a href="#arg-cmp" class="anchor" aria-label="anchor"></a><br>
A \`cmp_result\` or its \`cmps\` table (needs \`ambient_flag\` for the ambient check; use \`ambient = "flag"\` in \`find_coordinated_programs()\`).

context<a href="#arg-context" class="anchor" aria-label="anchor"></a><br>
Optional donor variable (column name or named vector) for \`cmp_context_association()\`.

covariates<a href="#arg-covariates" class="anchor" aria-label="anchor"></a><br>
Optional \`sample_metadata\` columns to adjust for.

replicate_by<a href="#arg-replicate-by" class="anchor" aria-label="anchor"></a><br>
Optional grouping variable (column name or named vector), e.g. site. Within each level with \>= \`min_level_donors\` donors, member-pair Pearson correlations of CMP-oriented scores are computed; concordance is the fraction that are positive.

tf<a href="#arg-tf" class="anchor" aria-label="anchor"></a><br>
Optional \`cmp_tf\` result.

nuisance<a href="#arg-nuisance" class="anchor" aria-label="anchor"></a><br>
Named regex vector (default \`cmp_nuisance_genes()\`).

top_n, n_recurrent<a href="#arg-top-n-n-recurrent" class="anchor" aria-label="anchor"></a><br>
Top genes per member (15) and recurrent genes kept (12).

max_ambient_frac, max_vdj_frac, min_marker_genes, min_concordance<a href="#arg-max-ambient-frac-max-vdj-frac-min-marker-genes-min-concordance" class="anchor" aria-label="anchor"></a><br>
Flag thresholds.

min_level_donors<a href="#arg-min-level-donors" class="anchor" aria-label="anchor"></a><br>
Minimum donors per \`replicate_by\` level (default 15).

min_ieg<a href="#arg-min-ieg" class="anchor" aria-label="anchor"></a><br>
Minimum number of immediate-early genes among the recurrent top genes for the \`stress_ieg\` flag (default 3).

</div>

<div class="section level2">

## Value<a href="#value" class="anchor" aria-label="anchor"></a>

Data frame, one row per CMP.

</div>

<div class="section level2">

## Details<a href="#details" class="anchor" aria-label="anchor"></a>

Recurrent top genes: the \`top_n\` genes with the largest CMP-oriented loadings in each member, ranked by how many members share them; the first \`n_recurrent\` are kept. Flags (all thresholds are arguments): \`ambient\` (ambient fraction \>= \`max_ambient_frac\`), \`platelet_rbc\` (\>= \`min_marker_genes\` contamination markers), \`vdj\` (V(D)J fraction of recurrent genes \>= \`max_vdj_frac\`), \`sex_chromosome\` and \`genotype_candidate\` (\>= \`min_marker_genes\` marker genes), \`stress_ieg\` (\>= \`min_ieg\` IEGs), plus soft notes \`genotype_marker_single\` (one marker gene) and \`low_replication\` (minimum within-level concordance \< \`min_concordance\`). \`label\` is the first hard flag in that order, else \`"candidate"\`, or \`"candidate_associated"\` when the context association has p \< 0.05. Labels are heuristics for prioritization, not classifications.

</div>

<div class="section level2">

## Examples<a href="#ref-examples" class="anchor" aria-label="anchor"></a>

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
```

</div>
