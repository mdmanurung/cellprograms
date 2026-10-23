# Real-cohort runbook: applying the cellprograms default path

Status: **no cohort data was available in this session**, so this runbook
documents the exact data contract and the validated default path
(`pa_cc`: per-cell-type EBMF, then sharing assessment) so the application run
can be executed as soon as a cohort is supplied. Everything below reflects the
prereg-coord-v1 benchmark verdicts (see `DECISION_prereg_coord_v1.md`).

## 1. Data requirements

| Item | Requirement | Benchmark reference |
|---|---|---|
| Unit of analysis | Donor-level pseudobulk per cell type (sum/mean of counts per donor x cell type), one numeric matrix per cell type, rows = donors, cols = genes | `gen_coord` / `simulate_pseudobulk` |
| Cell types | Named list; names must be unique, non-empty | `as_cell_program_data()` |
| Donors per cell type | >= 50 recommended; the benchmark used 50-100 donors per cell type with a shared 30-donor core across small types | sizes 50/75/100 |
| Genes | >= 2,000 after filtering; columns must be shared (identical gene IDs) across cell-type matrices | G = 2000 in the benchmark |
| Missing values | None allowed: NA/NaN/Inf are rejected at ingest | `as_cell_program_data()` CAND-03 check |
| Donor metadata | One row per donor: donor ID plus any covariates to adjust (batch, institute, sequencing round). Used for `covariates` residualization and, for the joint methods, as permutation strata | `institute` (3 levels) |
| Value scale | Normalized expression (CPM/log1p or similar); the fit residualizes covariates and gene-centers internally | preregistered pre-stage-1 adjustment |

## 2. Ingest and validation

```r
library(cellprograms)

# pseudobulk: named list of donors x genes matrices, one per cell type
x <- as_cell_program_data(pseudobulk, sample_metadata = meta)
validate_cell_program_data(x)

# optional but recommended: flag artifact genes before fitting
flagged <- flag_artifact_genes(x)
```

## 3. Default path (adopted by the benchmark decision)

The preregistered decision kept `pa_cc` (per-cell-type fits + sharing
spectrum) as the default; no benchmarked arm displaced it.

```r
fit <- fit_celltype_programs(x, max_factors = 10L,
                             covariates = "institute",      # your covariate column(s)
                             covariate_mode = "residualize",
                             seed = 1)
fit <- canonicalize_programs(fit)                            # required before sharing_spectrum()

sp  <- sharing_spectrum(fit, n_perm = 999L, seed = 1)        # pairwise sharing evidence
plot_sharing_spectrum(sp)

scores  <- program_scores(fit, cell_type = "CD4")            # donor activities
genes   <- top_genes(fit, cell_type = "CD4", program = 1, n = 20)
```

Interpretation guardrails from the benchmark:

- **Low-K regime.** Validated with 1-3 recurring programs per cohort against
  many private ones. If you expect dozens of shared programs, raise
  `max_factors` and treat conclusions as outside the validated envelope.
- **Activity-correlated private programs are the known trap.** Two cell types
  with genuinely distinct programs but correlated donor activity were called
  shared by the ungated default in 46% of trap replicates. Any multi-cell-type
  sharing call should be checked for gene-set overlap, not activity
  correlation alone (`match_programs()` / `align_programs()` help).
- **Do not use the experimental joint methods as defaults.** `block_E`
  (block EBMF) had strong power but called spurious shared programs in 16/16
  null replicates without the (not yet productized) permutation calibration.
  `eb_cp` was specific but underpowered on cell-type-specific-gene scenarios.
  Both are documented in `vignette("coordination_joint")` as experimental.

## 4. What is deliberately out of scope here

- Confirmatory-scale replication of the block_E calibration result: gated on
  the bounded rerun outcome (see `DECISION_prereg_coord_v1.md` follow-ups).
- Any cohort-specific preprocessing (cell-type annotation, doublet removal,
  pseudobulk aggregation choices): use your pipeline of choice; the package
  contract starts at the pseudobulk matrices.
