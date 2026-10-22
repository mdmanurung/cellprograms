# Refit program loadings on a common (union) gene panel

When each cell type selected its own features (e.g. \`features =
"variable"\`), the per-cell-type loading matrices live on different gene
sets and cannot be compared across cell types. This helper takes the
union of the fitted gene panels and, for each cell type, regresses the
union-panel expression on that cell type's program scores by least
squares. The fits and scores are not modified: this is route B of
ROADMAP item C5 (the \`genes_after\` arm of the model-selection
benchmark showed that fitting on the union panel from the start, route
A, does not improve recovery). Genes a cell type did not select get a
least-squares loading on its programs instead of being missing; for the
genes a cell type did fit, the refit loadings are close to (but not
identical to) the flashier posterior means, which are shrunk.

## Usage

``` r
union_panel_loadings(fit, genes = NULL)
```

## Arguments

- fit:

  A canonicalized \`cell_program_fit\`.

- genes:

  Optional character vector restricting the panel; default is the union
  of the per-cell-type fitted genes. Genes absent from a cell type's
  pseudobulk get \`NA\` loadings there, as do genes that are degenerate
  (zero-variance after centering) in that cell type when scaling.

## Value

A list with \`loadings\`, a named list of gene x program matrices on the
common panel (one per cell type, 0 columns for skipped cell types), and
\`genes\`, the panel used.

## Details

The fit's preprocessing is replayed on the panel before the regression:
covariate residualization (residualize mode) or removal of the fitted
covariate part (fixed mode), then the recorded centering/scaling.
