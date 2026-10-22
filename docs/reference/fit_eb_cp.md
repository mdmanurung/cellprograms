<div class="page-header">

# EB-CP tensor factorization across cell types

<div class="d-none name">

`fit_eb_cp.Rd`

</div>

</div>

<div class="ref-description section level2">

Fits a greedy rank-greedy three-mode tensor factorization (genes x donors x cell types) to per-cell-type pseudobulk matrices with EBNM normal-means updates per mode (point-laplace on gene and cell-type modes, normal on the donor mode) and per-(cell-type, donor) residual variance. Components are selected by an experimental held-out screening statistic on a donor-synchronous 50/50 split (stratified by \`strata\` when supplied), refit on the full data, and backfit. This is not a calibrated likelihood test or a guarantee of CP uniqueness. Memberships are descriptive, not error-controlled CMP sets. See the file header for the model and selection details.

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
fit_eb_cp(
  mats,
  strata = NULL,
  Kmax = 15L,
  ebnm_gene = ebnm::ebnm_point_laplace,
  ebnm_donor = ebnm::ebnm_normal,
  ebnm_ct = ebnm::ebnm_point_laplace,
  member_thresh = 0.2,
  heldout_z = 3,
  max_iter = 100L,
  tol = 1e-06,
  backfit_iter = 2L,
  seed = 1L
)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

mats<a href="#arg-mats" class="anchor" aria-label="anchor"></a><br>
Named list of (donors_ct x genes) matrices, one per cell type; common gene columns are intersected, rownames are donor ids. Inputs must already be covariate-adjusted and gene-centered (e.g. the \`Y_used\` matrices of a \`cell_program_fit\`).

strata<a href="#arg-strata" class="anchor" aria-label="anchor"></a><br>
Optional named factor/vector of donor strata (e.g. the batching covariate) used to stratify the held-out split. Strongly recommended whenever the input was residualized on a covariate.

Kmax<a href="#arg-kmax" class="anchor" aria-label="anchor"></a><br>
Maximum number of greedy components.

ebnm_gene<a href="#arg-ebnm-gene" class="anchor" aria-label="anchor"></a><br>
EBNM function for the gene mode.

ebnm_donor<a href="#arg-ebnm-donor" class="anchor" aria-label="anchor"></a><br>
EBNM function for the donor mode.

ebnm_ct<a href="#arg-ebnm-ct" class="anchor" aria-label="anchor"></a><br>
EBNM function for the cell-type mode.

member_thresh<a href="#arg-member-thresh" class="anchor" aria-label="anchor"></a><br>
Membership threshold on max-normalized \|W\| (stored on the result; used by \[eb_cp_programs()\]).

heldout_z<a href="#arg-heldout-z" class="anchor" aria-label="anchor"></a><br>
Held-out z threshold for keeping a component.

max_iter<a href="#arg-max-iter" class="anchor" aria-label="anchor"></a><br>
Maximum per-component coordinate iterations.

tol<a href="#arg-tol" class="anchor" aria-label="anchor"></a><br>
Per-component convergence tolerance.

backfit_iter<a href="#arg-backfit-iter" class="anchor" aria-label="anchor"></a><br>
Backfitting sweeps after the greedy pass.

seed<a href="#arg-seed" class="anchor" aria-label="anchor"></a><br>
Random seed (split and EBNM reproducibility).

</div>

<div class="section level2">

## Value<a href="#value" class="anchor" aria-label="anchor"></a>

A list of class \`eb_cp_fit\` with gene loadings \`U\` (genes x components), donor activities \`V\` (donors x components), cell-type loadings \`W\` (cell types x components), the retained \`genes\`/\`donors\`/ \`cts\`, selection diagnostics (\`kept_heldout_z\`, \`trace_z\`, \`frac_unexplained\`, \`n_components\`), and the frozen rules.

</div>
