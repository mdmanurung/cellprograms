# As-run scripts (exact provenance)

These are the scripts **exactly as executed** for the real-data MOFAcellulaR
comparison, preserved verbatim (including the analysis machine's `/workspace`
paths), so both committed vignette caches
(`vignettes/articles/mofacellular_stephenson_cache{,_excl}.rds`) have an exact
provenance record.

## Default arm (274-gene vdj + genotype_marker blacklist)

| Script | Produced |
|---|---|
| `run_masked_refit.R` | `stephenson_K30_masked_fit.rds` (B1 masked refit; also ran the combat cohort and the mashr sharing rerun in the same loop) |
| `run_cmp_stephenson.R` | `stephenson_K30_cmp.rds` |
| `run_mofa_stephenson.R` | `stephenson_mofa.hdf5` |
| `build_cache.R` | `mofacellular_stephenson_cache.rds` |

## Sensitivity arm (exclude_artifacts = TRUE, all 4 artifact sets)

| Script | Produced |
|---|---|
| `run_excl_refit.R` | `stephenson_K30_masked_fit_excl.rds` + `stephenson_K30_cmp_excl.rds` (refit and its CMP run) |
| `run_mofa_stephenson_excl.R` | `stephenson_mofa_excl.hdf5` |
| `build_cache_excl.R` | `mofacellular_stephenson_cache_excl.rds` |

The maintained, portable equivalents are the parameterized scripts in
`scripts/` (`01`–`05`, with `VIEW_DIMS`); use those to reproduce from a fresh
checkout. Known deltas between the two:

- `build_cache.R` is the executed version, which defines `view_dims` inline.
  Two earlier drafts were superseded in place and are not reproduced here:
  a `build_cache.R` without the `view_dims` definition (not runnable
  standalone), and an inline cache builder inside `run_mofa_stephenson.R`
  (field set `r2_df`, no `qc_note`/`n_donors_input`) whose output was
  overwritten by the final cache.
- The as-run default cache carries a hardcoded `qc_note` string; the
  parameterized `05_build_cache.R` derives the donor-retention numbers from
  the data instead. The as-run sensitivity cache has no `qc_note`; the
  parameterized script emits it for both arms. The vignette reads neither.
- `build_cache_excl.R` additionally records `blacklist_n_default = 274L`
  (unused by the vignette, which hardcodes the 274/313 labels).
- `05_build_cache.R` originally hardcoded the `prior_hvg` settings label;
  it now derives the label from the fit's `exclude_artifacts` flag so one
  script serves both arms.
