# CMP benchmark harness (`cmp_bench`)

The harness used for the cellprograms-vs-MOFAcellulaR comparison and the
sharing/subspace calibration studies (preregistered in `PREREG.md`). Copied
verbatim from the session workspace so it is preserved alongside the package.
The numbered pipeline scripts that produced the committed comparison results
live in `../mofacellular_comparison/scripts/`, with the exact executed
versions under `../mofacellular_comparison/scripts/as-run/`.

## Layout

- **Harness modules (this folder)**: shared infrastructure (`common.R`,
  `config.R`, `registry.R`, `fit_cache*.R`), candidate methods and evaluation
  (`cand_*.R`, `candidates.R`, `run_cand*.R`, `eval_cmp.R`, `v2_score.R`),
  real-data drivers (`real_*.R`, `cross_cohort*.R`, `joint_fit.R`,
  `paired.R`), figures (`figs*.R`, `annot_*.R`), checks (`parity.R`,
  `rngcheck.R`, `prereg_check*.R`, `ctx_test.R`, `split_half.R`,
  `sim_cmp.R`, `ambient.R`, `lr_eval.R`), and `mofa_helper.py` (MOFA2 helper).
- **`drivers/`**: top-level scripts run from the workspace root that source
  the harness or operate on its cached fits.

## Drivers

| Script | What it does |
|---|---|
| `run_sharing.R` | `program_sharing()` on the real Stephenson/COMBAT K30 fits (permutation-calibrated matching on stored posterior SDs; no refit needed). |
| `run_subspace.R` | Multi-candidate evaluation driver reading fits from a local mirror; same semantics as `run_cand_v2.R` but reads each fit once and evaluates all candidates on it (avoids re-reading the slow S3-FUSE store). |
| `run_b2_b3.R` | B2: V-estimation `z_thresh` sweep on the unmasked K30 fits; B3: matching calibration. Outputs `v_sensitivity.csv`, `matching_calibration.csv`. |
| `smoke_sim.R` | Simulation smoke/timing loop over bases (`param`, `semi`) x CMP scenarios via the harness. |
| `smoke_subspace.R` | Subspace-candidate smoke run through the harness. |
| `prof_eval.R` / `prof_eval2.R` | Profiling wrappers around the harness candidate evaluation (`eval_cmp.R`). |
| `prof_subspace.R` | Profiling wrapper around the subspace candidate path. |
| `summarize_sharing.R` | Empirical sharing patterns from the real K30 `program_sharing()` runs (reads `{co}_K30_sharing.rds` + artifact flags, writes summary figures). |
| `breadth_winner.R` | Signal-breadth summary across candidate methods from the simulation scorecards. |

## Notes

- Everything here is **as-run provenance**, not a portable pipeline: scripts
  contain absolute paths to the original session layout
  (`/mnt/shared-workspace/cmp/bench/`, `/workspace/.Rlib`, cached `.rds`
  files). The harness originally lived at `/mnt/shared-workspace/cmp/bench/`;
  drivers `source()` that path, so rerunning requires restoring the harness
  at that location or editing the paths.
- `run_sharing.R` and `run_b2_b3.R` load cellprograms directly via
  `pkgload::load_all()` and require `mashr` for `program_sharing()`.
