The default `cellprograms` pipeline fits each cell type separately and
only then asks whether programs are shared. Two experimental
(v0.2.0.9001) extensions invert the ordering and factorize across cell
types directly:

- **EB-CP** (`fit_eb_cp()` / `eb_cp_from_fit()`): a tensor-style greedy
  EBMF over the stacked cell-type x donor matrices. Each rank-1
  component is a gene loading (`U`), a donor activity vector (`V`)
  shared across the cell types the component spans, and a cell-type
  loading vector (`W`). Components are accepted by a held-out z gate on
  a donor-synchronous split, so the number of components is selected by
  out-of-sample predictive evidence, not by a variance cutoff.
- **Block EBMF** (`fit_joint_block()` / `block_clusters()`): a
  shared-first flashier pass over the stacked expression matrix,
  followed by per-cell-type private passes on the residuals.
  `block_clusters()` collapses shared factors into gene-set clusters by
  loading overlap and returns donor activity per cluster.

This vignette runs both on a small simulated dataset and states the
scope caveats measured in the coordination benchmark (see
`benchmarks/coordination/` in the repository).

    library(cellprograms)

    sim <- simulate_pseudobulk(scenario = "mixed", n_samples = 60,
                               n_genes = 800, program_size = 40, seed = 1)
    x <- as_cell_program_data(sim$pseudobulk)
    fit <- fit_celltype_programs(x, loading_prior = "point_laplace",
                                 max_factors = 10L, seed = 1)
    fit <- canonicalize_programs(fit)

## EB-CP: greedy components with held-out selection

`eb_cp_from_fit()` runs EB-CP on the adjusted per-cell-type matrices of
a fitted `cell_program_fit`. If `strata` is not supplied, the
donor-synchronous train/held-out split is derived from the covariates
used at fit time, so the split respects the adjustment; pass `strata`
explicitly when you want a different stratification.

    cp <- eb_cp_from_fit(fit)
    cp$n_components
    #> [1] 8
    cp$kept_heldout_z          # held-out z of each kept component
    #> [1] 141.06089 154.97144 103.76617 107.39676  94.15278  63.45945  93.65694
    #> [8]  51.93119
    cp$trace_z                 # z trace including the first failed component
    #> [1] 141.06089 154.97144 103.76617 107.39676  94.15278  63.45945  93.65694
    #> [8]  51.93119

`eb_cp_programs()` turns components into program calls. Components
spanning &gt;= `min_members` cell types at the membership threshold are
flagged as coordination programs; the rest are kept as private programs,
flagged rather than dropped.

    pr <- eb_cp_programs(cp)
    pr$programs
    #>   program cell_type n_members is_coordination heldout_z
    #> 1    CP01      Mono         1           FALSE 141.06089
    #> 2    CP02    CD8+NK         2            TRUE 154.97144
    #> 3    CP03         B         1           FALSE 103.76617
    #> 4    CP04       CD4         1           FALSE 107.39676
    #> 5    CP05        NK         1           FALSE  94.15278
    #> 6    CP06      Mono         1           FALSE  63.45945
    #> 7    CP07    CD8+NK         2            TRUE  93.65694
    #> 8    CP08         B         1           FALSE  51.93119
    head(pr$activity[, 1:2, drop = FALSE])   # donor activity per coordination component
    #>        cluster_02  cluster_07
    #> obs1   2.53382698 -2.76698019
    #> obs10  0.49417548 -0.54670036
    #> obs11 -0.36878957  0.34052851
    #> obs12 -1.12515663  1.19204711
    #> obs13 -0.07259629  0.05518669
    #> obs14 -0.71861307  0.58362302

## Block EBMF: shared-first pass plus private residuals

`fit_joint_block()` stacks the per-cell-type matrices and fits shared
factors first, then per-cell-type private factors on the residuals.

    bfit <- fit_joint_block(fit, level = "expression", order = "shared_first",
                            max_factors = 10L, seed = 1, verbose = 0L)
    bfit$shared_fit$n_factors
    #> [1] 4
    bc <- block_clusters(bfit, mass_thresh = 0.1)
    head(bc$clusters, 10)
    #>   cluster cell_type program_id      stat fwer_p
    #> 1       1         B       SF01 0.3283760     NA
    #> 2       1       CD4       SF01 0.3283760     NA
    #> 3       1        NK       SF01 0.3283760     NA
    #> 4       1      Mono       SF01 0.3283760     NA
    #> 5       2       CD8       SF02 0.5278789     NA
    #> 6       2        NK       SF02 0.5278789     NA

Each row of `bc$clusters` is a gene-set cluster (a group of shared
factors with overlapping loadings); `bc$activity` holds donor activity
per cluster. The block-mass rule (`mass_thresh`) keeps clusters whose
membership spreads across at least two cell types.

## Scope caveats from the coordination benchmark

Both extensions trade robustness for directness, and the preregistered
coordination benchmark (`benchmarks/coordination/`, tag
`prereg-coord-v1`) quantifies where:

- **Block EBMF (`block_E`) is powerful but over-calls without
  calibration.** On null data the shared-first pass still returns
  factors whose loadings spread across blocks, and the block-mass rule
  alone called spurious shared programs in 16/16 null replicates. A
  permutation-calibrated significance layer (within-strata donor
  permutation, max-T family-wise control) is available in the benchmark
  harness and is being productized; until then, treat uncalibrated
  `block_clusters()` output as a candidate list, not as evidence of
  sharing.
- **EB-CP is conservative.** The held-out z gate keeps specificity high
  (near-zero false positives on null data) but loses power when shared
  programs have cell-type-specific gene loadings, because the model’s
  shared donor-activity assumption is only approximately met.
- **Low-K regime.** Both methods were validated with few recurring
  programs (1-3) against many private ones. With dozens of true
  per-cell-type programs, set `max_factors` above the default and expect
  the shared-first ordering to matter more.

For the default-by-construction path (per-cell-type fits, then
`sharing_spectrum()` / `plot_sharing_spectrum()`), see
`vignette("quickstart")`.
