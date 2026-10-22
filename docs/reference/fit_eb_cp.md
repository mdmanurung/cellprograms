# EB-CP tensor factorization across cell types

Fits a greedy rank-greedy three-mode tensor factorization (genes x
donors x cell types) to per-cell-type pseudobulk matrices with EBNM
normal-means updates per mode (point-laplace on gene and cell-type
modes, normal on the donor mode) and per-(cell-type, donor) residual
variance. Components are selected by a held-out likelihood test on a
donor-synchronous 50/50 split (stratified by \`strata\` when supplied),
refit on the full data, and backfit. See the file header for the model
and selection details.

## Usage

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

## Arguments

- mats:

  Named list of (donors_ct x genes) matrices, one per cell type; common
  gene columns are intersected, rownames are donor ids. Inputs must
  already be covariate-adjusted and gene-centered (e.g. the \`Y_used\`
  matrices of a \`cell_program_fit\`).

- strata:

  Optional named factor/vector of donor strata (e.g. the batching
  covariate) used to stratify the held-out split. Strongly recommended
  whenever the input was residualized on a covariate.

- Kmax:

  Maximum number of greedy components.

- ebnm_gene, ebnm_donor, ebnm_ct:

  EBNM functions for the gene, donor, and cell-type modes.

- member_thresh:

  Membership threshold on max-normalized \|W\| (stored on the result;
  used by \[eb_cp_programs()\]).

- heldout_z:

  Held-out z threshold for keeping a component.

- max_iter, tol:

  Per-component coordinate ascent controls.

- backfit_iter:

  Backfitting sweeps after the greedy pass.

- seed:

  Random seed (split and EBNM reproducibility).

## Value

A list of class \`eb_cp_fit\` with gene loadings \`U\` (genes x
components), donor activities \`V\` (donors x components), cell-type
loadings \`W\` (cell types x components), the retained
\`genes\`/\`donors\`/ \`cts\`, selection diagnostics
(\`kept_heldout_z\`, \`trace_z\`, \`frac_unexplained\`,
\`n_components\`), and the frozen rules.
