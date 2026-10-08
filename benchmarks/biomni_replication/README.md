# Biomni pseudobulk replication

`01_pseudobulk.py --cohort {combat,stephenson,hlca,onek1k} [--source {X,raw}]` (run via
`sbatch -J pb_<cohort> 01_pseudobulk.sbatch <cohort> [--source X|raw]`; `--partition=all`, 4 cpu, 128G, ~2-3 min per cohort).
Census pinned to 2025-11-08; dataset resolved from the Census datasets table; source h5ad fetched with
`cellxgene_census.download_source_h5ad` (author labels are not in Census obs), streamed in row chunks with a sparse group-sum.
Outputs: `mat/ct_*.csv`, `mat/global.csv`, `mat/group_*.csv` (COMBAT only), `donor_meta.csv`, `composition.csv`, `pb_var.csv`,
`pb_counts.parquet`, `pb_coldata.csv`, `qc.json` (records `source`).

## Two arms
- `--source X` (default) -> `data/<cohort>/`. Sums the log-normalized `X` (log1p CP10k) per donor x cell type, truncates the sums to
  integers (`np.floor`), then applies Biomni's `prep_matrix` (CPM + log2 on those sums, top-3000 variance genes). This is what
  Biomni's `pb_counts.parquet` actually contained (the name is misleading: it is not raw counts). Validated against
  `ref/execution_trace/pseudobulk.ipynb` cell 4: all 120 CT x marker cells reproduced, max abs diff 0.00 (rounding instead of
  truncating: max 0.76; no integer cast: 0.84; mean of X used as-is: values 0-4.5, nowhere near). Top-3000 panels vs Biomni's
  `program_annotations.csv` top genes: 100% overlap in all 10 COMBAT cell types. Check scripts: `data/_validation/val.py`, `val2.py`.
- `--source raw` -> `data/<cohort>_raw/`. Sums raw counts (`raw.X`; `X` when there is no `raw`). Kept for comparison: it gives GNLY/NK 13.57
  (trace 10.94) and only 59-95% of Biomni's top program genes inside our panels (4-42% missing per cell type, PB worst).
- Note: other scripts that read `data/<cohort>/` now get the X arm. OneK1K has no `raw` and its `X` is integer counts, so its two arms are
  numerically identical (same `pb_counts`, same `mat/`); only donor_meta/composition differ (see RISKS).
- `data/_h5ad/` is shared by both arms.

## Choices
- COMBAT: label `major_subset`, dataset ebc2e1ff-c8f9-466a-acf4-9d291afaf8b3. Exact recipe; string "nan" labels dropped from
  pseudobulk but kept as a `nan` column in composition (18 columns, as in Biomni's file).
- Stephenson: `author_cell_type`, CTs with >=30 donors at >=50 cells (19). Age is not an obs column of the source h5ad;
  `Age` is INFERRED from `development_stage` ("N-year-old stage"; NaN for decade-only stages). donor_meta = first row per donor.
- HLCA (INFERRED): core dataset 066943a2 (collection 6f6d381a, title "(core)"); label = Census `cell_type` (matches the 9 CTs and all
  n_obs exactly); `is_primary_data` is False for every core cell so the primary filter is skipped for this cohort;
  donor_meta = per-donor mode (tissue_level_2 varies within donor); `Age` = `age_or_mean_of_age_range`; technical: `study`, `dataset`,
  `sequencing_platform`.
- OneK1K (INFERRED): dataset 3faad104, all donors "normal"; label = Census `cell_type` mapped to CD4_CM etc.; cap = 300 donors
  chosen by `default_rng(42).permutation(sorted donors)[:300]`. The exact Biomni draw is unknown: B_Naive/B_Mem/CD4_CM
  match, the other CTs are off by 1-9 donors (best of ~60 tried seeded-draw variants). `pool_number` kept as technical covariate.
- The earlier claim that the COMBAT gene-order difference was due to ties was wrong. Cause: the raw-count arm used a different input
  (raw counts instead of Biomni's truncated sums of log-normalized X), which changes both the gene ranking and which genes enter the panel.

## QC (X arm; donor counts per CT vs Biomni `program_summary.csv`)
COMBAT, Stephenson, HLCA: all CTs match (COMBAT: CD4 121, CD8 119, B 109, NK 118, cMono 120, ncMono 105, DC 67, PB 52, GDT 50, DP 57).
OneK1K: 3 of 9 CTs match (B_Naive, B_Mem, CD4_CM); CD14_Mono 72/68, CD4_EM 52/54, CD4_Naive 292/290, CD8_EM 265/259, CD8_Naive 114/105,
NK 280/281 (inferred 300-donor draw). X and raw arms give identical donor counts for every cohort.

## RISKS
- `donor_meta` with `agg="first"` (COMBAT, Stephenson, OneK1K) takes the first row per donor; arbitrary for variables that vary within a
  donor (e.g. Stephenson `Days_from_onset`, `Collection_Day`, `Status`, `sample_id` for resampled donors; COMBAT `scRNASeq_sample_ID`, `GEX_region`).
- `composition.csv`/`donor_meta.csv` use all primary cells of the cohort. For OneK1K the X arm now restricts both to the 300 capped donors;
  `data/onek1k_raw/` was written before this fix and still has 981 donor rows.
- Stephenson `Age` is NaN for 62/120 donors (decade-only `development_stage`); any model using Age drops or imputes half the donors.
- HLCA `donor_meta` uses the per-donor mode, whereas Biomni's trace used the first row per donor; categorical covariates can differ for
  donors with several sites/samples.
- OneK1K donor cap is an inferred seeded draw; 6 of 9 CTs differ from Biomni's n_obs by 1-9 donors.

