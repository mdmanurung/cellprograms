# Cross-cohort immune-aging analysis (Terekhova et al., Immunity 2026)

Application of the cellprograms pipeline (EBMF donor-level programs per cell
type, MR-0 freeze / MR-1 projection) to the 8-cohort PBMC aging resource
(2,609 healthy donors; ABF300, OneK1K, AIDA, CIMA, KCL, UCSF, SoundLife,
SPAC). This directory contains everything needed to reproduce the analysis
locally: scripts, small result tables, and exact commands.

## Data access

All eight cohorts are served by the authors' portal (open access; no Synapse
or EGA account needed for the curated h5ads):

1. Cohort config (donor counts, cohort list):
   `https://artyomovlab.wustl.edu/human-immune-aging-trajectories/cohorts.config.json`
2. File listing: `https://artyomovlab.wustl.edu/scn-m/scn/getFiles?token=<TOKEN>`
   with one token per cohort directory, e.g. `OneK1K_new_annotation_clean`,
   `ABF300_new_annotation_clean`, `SPAC_new_annotation_clean`, `KCL_annotation`,
   `CIMA_new_annotation_clean`, `AIDA_new_annotation_clean`,
   `SoundLife_new_annotation_clean`, `UCSF_new_annotation_clean`.
3. Download: `https://artyomovlab.wustl.edu/scn-m` + the listed path
   (`.../files/<Cohort>_reannotated.h5ad`), ~1.9-41 GB per cohort.

Each h5ad carries `.X` (log-normalized) and `.raw/X` (integer counts) with a
13-column harmonized `obs` schema (donor_id, age, sex, ethnicity, CMV, BMI,
Dataset, barcode, annotation, author_celltype, major/minor_celltypes,
total_counts). `annotation` holds the 59 harmonized cell types; drop the
"Bad quality" group at pseudobulking. CMV is informative only in KCL
(Neg/Pos/No_info) and SoundLife (Neg/Pos).

## Pipeline (run in this order)

Environment: R >= 4.4 with `cellprograms` (this repo), `rhdf5`, `clue`,
`sva` (benchmark only); Python >= 3.10 with `h5py`, `numpy`, `scipy`.
All scripts take explicit paths; `<dir>` below is a per-cohort output
directory.

```bash
# 1. Audit schema/units of each cohort (optional but recommended)
python scripts/audit_h5ad.py <cohort>.h5ad

# 2. Pseudobulk: sum .raw/X counts over all cells per donor x cell type
#    (min 10 cells/donor-ct; drops "Bad quality"); writes pseudobulk.h5,
#    sample_metadata.csv, ct_summary.csv
python scripts/pseudobulk_h5ad.py <cohort>.h5ad <dir> <COHORT> --min-cells 10

# 3. Fit EBMF programs per cell type (top-2000 variable genes, log1p-CP10K,
#    max_factors=10, point_laplace, var_type=1, seed=20261016, NO covariate
#    residualization — see amendment note below)
Rscript scripts/fit_pilot.R <dir> <COHORT> <dir>/<cohort>_full_fit.rds ""

# 4. Associations: lm(score ~ age + sex + CMV) per program, BH-FDR within ct
Rscript scripts/assoc_programs.R <dir>/<cohort>_full_fit.rds <dir> \
        <dir>/<cohort>_full_assoc.csv

# 5. Projection (per query cohort vs the discovery fit): MR-1 projection of
#    discovery programs onto query pseudobulks + matched score concordance
Rscript scripts/project_cohort.R <discovery_dir>/<disc>_full_fit.rds <dir> \
        <QueryName> <dir>/<query>_vs_<disc>

# 6. Association-profile concordance (matched refit programs)
Rscript scripts/assoc_concordance.R <disc>_full_assoc.csv <q>_full_assoc.csv \
        <dir>/<query>_vs_<disc>_concordance.csv <dir>/<query>_assoc_concordance.csv

# 7. Projected-score age replication (discovery programs tested on projected
#    scores in each query; sidesteps decomposition mismatch)
Rscript scripts/projected_age_replication.R

# 8. Assemble the 3-tier conservation panel + gene-set overlap + age summaries
Rscript scripts/assemble_conservation.R
Rscript scripts/gene_overlap.R
Rscript scripts/age_summary.R
Rscript scripts/age_topgenes.R   # top genes of key age programs (OneK1K)
```

**Registered operating point (do not silently change):** max_factors=10,
point_laplace, var_type=1, features="variable" (top 2000 variable genes per
cell type), feature_scale="logcounts", seed=20261016, canonicalize after fit.
**Amendment (pilot-falsified):** fits use NO covariate residualization —
residualizing age+sex+CMV inside the fit makes program scores orthogonal to
age by construction and the association test vacuous. Associations are tested
post hoc on unresidualized scores.

## Result tables (committed, small)

| Table | Content |
|---|---|
| `tables/conservation_panel.csv` | 3-tier panel per cohort pair: matched score concordance (median \|cor\|, frac >= 0.5, sign flips) + association-profile Spearman/sign agreement |
| `tables/age_by_cohort.csv` | programs tested and fraction age-associated (FDR<0.05) per cohort |
| `tables/age_by_celltype_consistency.csv` | strongest-age program per cell type: cross-cohort sign consistency |
| `tables/projected_age_replication.csv` | per-program age betas of discovery programs on projected scores in each query |
| `tables/gene_overlap_top20.csv` | top-20 gene Jaccard for matched program pairs |

## Headline results (as run 2026-10-16)

- Matched score concordance: median \|cor\| 0.49-0.63 across the five
  ABF300-discovery pairs and the KCL->SoundLife pair.
- Age-associated programs (FDR<0.05): 44-56% in the 400-1000 donor cohorts
  vs 3-19% below 200 donors (power scaling).
- Thymic-involution axis: CD8+ RTE strongest-age program negative in 7/7
  cohorts (median -0.29 SD/10y); CD4+ RTE 6/7 (-0.32); gd T naive 6/7 (-0.32);
  gd T Vd2 GZMK+ 7/7 (-0.28).
- Matched age-profile Spearman is null-to-weak everywhere (-0.07..0.17)
  including at 981 donors — decomposition mismatch, not noise (see vignette).
- Projected-score replication of discovery age hits: 36-68% per query cohort
  (null ~4%); CIMA 68%, OneK1K 52%.
- Top-20 gene Jaccard of matched programs: median 0.03-0.13 — conservation
  lives at the score/subspace level, not in top-N gene lists.

## Large artifacts (not committed)

Per-cohort fits (`<cohort>_full_fit.rds`, 0.2-1.6 GB), pseudobulks
(`pseudobulk.h5`, 0.2-7 GB), and raw h5ads are reproducible from the portal
via the commands above. The full fits and projection outputs are also
distributed alongside this repository (see the session results archive,
`aging_pbmc_repro/`).
