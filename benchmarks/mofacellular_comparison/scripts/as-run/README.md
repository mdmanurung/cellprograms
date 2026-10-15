# As-run sensitivity scripts (exact provenance)

These are the three scripts **exactly as executed** for the
`exclude_artifacts = TRUE` sensitivity arm, preserved verbatim (including the
analysis machine's `/workspace` paths) so the committed sensitivity cache
(`vignettes/articles/mofacellular_stephenson_cache_excl.rds`) has an exact
provenance record.

| Script | Produced |
|---|---|
| `run_excl_refit.R` | `stephenson_K30_masked_fit_excl.rds` + `stephenson_K30_cmp_excl.rds` (refit and its CMP run) |
| `run_mofa_stephenson_excl.R` | `stephenson_mofa_excl.hdf5` |
| `build_cache_excl.R` | `mofacellular_stephenson_cache_excl.rds` |

The maintained, portable equivalents are the parameterized scripts in
`scripts/` (`04_fit_excl.R`, `02_cmp.R`, `03_mofa.R`, `05_build_cache.R` with
`VIEW_DIMS`); use those to reproduce from a fresh checkout. Known trivial
deltas between the two:

- `build_cache_excl.R` additionally records `blacklist_n_default = 274L`
  (unused by the vignette, which hardcodes the 274/313 labels).
- `05_build_cache.R` originally hardcoded the `prior_hvg` settings label;
  it now derives the label from the fit's `exclude_artifacts` flag so one
  script serves both arms.
