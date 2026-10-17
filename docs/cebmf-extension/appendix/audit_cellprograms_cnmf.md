# Audit: cellprograms (R + benchmarks) and cNMF — integration map for a multiresolution EBMF workflow

Audited at commit `b191a5a` (default branch `main`, shallow clone) of `github.com/mdmanurung/cellprograms`, and cNMF at `/workspace/.venv/lib/python3.11/site-packages/cnmf/cnmf.py` (1,298 lines). All paths relative to the repo root unless prefixed `cnmf.py:`.

## (a) Repo layout map

```
R/                      package core (17 files, see inventory below)
benchmarks/
  biomni_replication/   Census pseudobulk builder + replication of the Biomni trace
                        (01_pseudobulk.py, 02*_metrics.R, 03_fit_cp.R, ...);
                        fork/ = vendored copy of the ebmf package (drift risk, see (e));
                        ref/ = reference cellprograms_benchmark outputs (combat/hlca/onek1k/stephenson)
  cmp_bench/            CMP detection benchmark (candidates, joint_fit, eval)
  cmp_detection/        same scripts, as-run copy under scripts/
  batch_adjustment/     ComBat-style batch-adjustment benchmark
  coordination/         CMP coordination pilot (methods.R, scoring.R)
  divas_comparison/     DIVAS comparison
  model_selection/      factorial model-selection benchmark (arms.R, sim_semi.R)
  mofacellular_comparison/, c2_matching/
experiments/            (README only)
reports/                ROADMAP.md, HANDOFF.md, STAGE_B_ROADMAP.md (MR-1..MR-4 staged plan), audit/
python/                 benchmark python drivers (run_combat_*.py, run_stephenson_*.py, mofa_fit.R, ...)
biomni/                 REVIEW.md, REVIEW_round2.md (review docs only)
analyses/, docs/, inst/, man/, tests/, vignettes/
```

## (b) Existing capabilities inventory (file:line)

**Data contract — pseudobulk is an INPUT, not built by the R package**
- `R/data-model.R:23` `as_cell_program_data(pseudobulk, sample_metadata, ...)` — named list of genes x observations numeric matrices, one per cell type; validation at `R/data-model.R:113` (`validate_cell_program_data`); NA-free requirement at `R/data-model.R:41`.
- `R/simulation.R:29` `simulate_pseudobulk()` — synthetic only.
- **No R function constructs pseudobulk from single-cell counts.** All construction lives in the benchmark layer.

**Pseudobulk construction (benchmark layer) — `benchmarks/biomni_replication/01_pseudobulk.py`**
- Aggregation unit: **donor x cell type** (sample key `"donor||ct"`, line 121); condition/disease are NOT part of the unit — they are carried as donor-level covariates in `donor_meta.csv` (`donor_meta`, line 156).
- `select_samples` (line 93): keep (donor, ct) groups with >= `MIN_CELLS=50` cells (line 23); per-cohort config `CFG` (lines 39-67): COMBAT fixed 10 `major_subset` types; Stephenson `author_cell_type` with >= 30 donors/CT; HLCA `cell_type` (9 types, `agg="mode"`); OneK1K `cell_type` capped at 300 healthy donors (seeded, line 47).
- `group_sum` (line 119): streamed 100k-row chunks, sparse indicator matmul sums X (log-normalized, floored, line 145) or raw.X per (donor, ct).
- `prep_matrix` (line 149): verbatim COMBAT recipe — genes expressed (>=1) in >=10% of samples, CPM, `log2(CPM+1)`, top 3000 HVG by variance.
- Outputs (main, line 180): `pb_counts.parquet`, `pb_coldata.csv`, per-cell-type `mat/ct_<CT>.csv`, plus `mat/global.csv` and `mat/group_<G>.csv` (group-level pseudobulks, lines 208-212) — the group level is a crude second resolution that already exists.
- Census pinned: `CENSUS = "2025-11-08"` (line 22).

**Metacell: zero code, placeholder only**
- No hits for metacell/SEACell in R/, benchmarks/, python/ (only reference PDFs/citations).
- `reports/STAGE_B_ROADMAP.md:51-58` — "MR-3 — Metacells design doc (STAGED)": metacells named as the middle ground between donor pseudobulk and raw cells; SEACell/Metacell-2 as grouting candidates; "no code until the doc is reviewed".
- Related staged items in the same file: MR-2 cNMF baseline (lines 42-49, registered comparison rule), MR-4 mcRigor/anchored refinement (lines 60-70, to be built on flashier.utils `fit_initialized_ebmf()`).

**Program projection onto new cells — exists and is mature**
- `R/cell-projection.R:80` `project_programs_to_cells(fit, expression, cell_types, ...)` — replays the fit's stored preprocessing (per-gene OLS residualization coefficients, then center/scale; lines 190-206), calls `flashier.utils::project_factor_activity` (line ~216), re-applies the canonicalization RMS-scale/sign multipliers (lines 219-221), returns cells x programs activity + coverage + provenance.
- `R/cell-projection.R:230` `.replay_design` — rebuilds the covariate design for query donors, fails closed on level mismatch.
- Limits: dense matrices only ("sparse not supported in v1", line 87), per-cell-type loop with no cell batching, requires the flashier.utils manifest (line 178), `min_coverage=0.9` intersection gene matching.

**Downstream interpretation**
- `R/cmp-interpret.R:28` `prepare_collectri(raw)` — OmniPath CollecTRI -> signed `source/target/mor` network; NFKB/AP1 COMPLEX collapsing (lines 34-38).
- `R/cmp-interpret.R:144` `annotate_cmp_tf(fit, cmp, net, ...)` — hand-rolled decoupler-style ULM of signed program loadings on TF regulons, with gene-permutation `padj` and a self-contained sign-flip test (`p_self`); `cmp_triage` at line 321; `cmp_nuisance_genes` at line 250. (No decoupler package dependency — ULM is implemented in-base.)
- CMP detection: `R/coordination.R:269` `find_coordinated_programs`, `cmp_scores` at 470, `ambient_program_scores` at 407.
- Sharing/stability/matching: `R/program-sharing.R:399` `program_sharing` (mashr-based); `R/stability.R:31` `assess_program_stability`; `R/principal-angles.R:153` `principal_angles`, `:294` `sharing_spectrum`, `:354` `align_programs`, `:392` `match_programs` (cross-cohort program matching).

**Fitting backend (the seam cebmf_torch would plug into)**
- `R/ebmf.R:94` `fit_celltype_programs(x, loading_prior="point_laplace", ...)` — per-cell-type flashier fit `Y_c = Z_c W_c' + E_c`.
- `R/ebmf.R:467` `.fit_via_flashier_utils(args, seed, feature_scale)` — ordinary fits route through `flashier.utils::fit_ebmf()`; legacy `do.call(flashier::flash, args)` at line 293; escape hatch `options(cellprograms.use_flashier_utils = FALSE)` (lines 125-128).
- `R/ebmf.R:545` `canonicalize_programs(fit)` — RMS/sign canonicalization that projection depends on.
- `DESCRIPTION:20,22,34-35` — Imports `flashier` (pinned commit `fd44811`) and `flashier.utils` (`mdmanurung/flashier-utils`).

## (c) cNMF consensus internals (`cnmf.py`)

- `consensus(k, density_threshold=0.5, local_neighborhood_size=0.30, ...)` — line 823.
  1. Load merged spectra (all iterations of all runs, n_iter*k rows); L2-normalize each spectrum to unit length (line ~872).
  2. Local density: `n_neighbors = int(0.30 * n_spectra / k)`; Euclidean distances among L2 spectra (`euclidean_distances`); local density = mean distance to the n_neighbors nearest neighbors via `argpartition` (lines 885-899, cached in `local_density_cache`).
  3. Density filter: keep spectra with `local_density < density_threshold` (line 903); `>= 2.0` disables filtering; aborts if nothing remains (line 906).
  4. `KMeans(n_clusters=k, n_init=10, random_state=1)` on the filtered L2 spectra (line 910); consensus spectrum per cluster = per-gene **median** of member spectra, then row-normalized to sum 1 (lines 914-917).
  5. Usages: `refit_usage` (line 776) = NNLS via sklearn NMF with `H` fixed (`update_H=False`) on the normalized counts; a final usage refit on std-scaled HVG TPM runs when `refit_usage=True` (lines ~985-995).
  6. `spectra_tpm` via `refit_spectra` (line 805, NNLS with usage fixed on TPM); `gene_spectra_score` = OLS coefficients of Z-scored TPM on usages (`efficient_ols_all_cols`, line 55).
- **K-selection statistic** — `k_selection_plot` (line 1119) calls `consensus(k, skip_density_and_return_after_stats=True)` per k; in that mode the density filter is skipped (`density_threshold_str='2'`, line 877) and the statistic is `silhouette_score(l2_spectra.values, kmeans_cluster_labels, metric='euclidean')` (line 923): the **mean Euclidean silhouette of the KMeans(k) clustering of the L2-normalized (unfiltered) merged spectra** — i.e., stability of replicate-iteration clustering, borrowed from Alexandrov et al. 2013 (line 1122). Companion statistic: `prediction_error` = squared Frobenius norm of `norm_counts - usage.dot(median_spectra)` (lines 925-930). Plotted as stability (blue) vs error (red) vs K (lines 1142-1150).
- **Reusable objects for a fair multi-method benchmark**:
  - `load_results(K, density_threshold, n_top_genes=100, norm_usage=True)` — line 1161: returns `(usage [cells x K, optionally sum-to-1], spectra_scores [K x genes, Z-score coefficients], spectra_tpm [K x genes], top_genes [genes x K])` read from the saved txt artifacts. This is the clean handoff surface for benchmarking cNMF programs against EBMF programs on identical inputs.
  - `refit_usage(X, spectra)` / `refit_spectra(X, usage)` (lines 776/805) — generic fixed-factor NNLS; the natural cNMF-side "projection onto new cells" (NNLS refit of usages given frozen spectra), directly comparable to `project_programs_to_cells`.
  - Consensus artifacts saved at paths defined in `_initialize_dirs` (lines 312-315): `*.spectra.k_%d.dt_%s.consensus`, `*.usages.k_%d.dt_%s.consensus` (npz + txt).
  - Caveat: cNMF expects non-negative input (counts/TPM pipeline in `prepare`/`get_norm_counts`, lines 333/487), whereas cellprograms fits log-normalized pseudobulk — the benchmark harness must feed each method its native scale rather than mutating either package.

## (d) Net-new vs build-on for the four capabilities

| Capability | Status | Evidence / seam |
|---|---|---|
| **Pseudobulk discovery** | **BUILD-ON** (contract exists; constructor is benchmark-only) | Input contract `R/data-model.R:23`; fitting seam `R/ebmf.R:94` + backend routing `R/ebmf.R:467` (swap flashier.utils backend for cebmf_torch). Constructor `01_pseudobulk.py:119` (`group_sum`) and `:149` (`prep_matrix`) work but are cohort-hardcoded and Census-pinned — promote to a reusable module; donor x cell-type unit and group-level pseudobulks (`main`, lines 208-212) are reusable patterns. |
| **Metacell discovery** | **NET-NEW** | Zero code anywhere; only the staged MR-3 design doc (`reports/STAGE_B_ROADMAP.md:51-58`) naming SEACell/Metacell-2. The existing list-of-matrices data contract absorbs metacell pseudobulk unchanged (metacells become the "observations" of a cell-type block). |
| **Single-cell batched projection** | **BUILD-ON for the wrapper, NET-NEW for the numerics** | Wrapper template `R/cell-projection.R:80` (preprocessing replay + canonicalization + coverage/provenance) is the contract to copy; but it is dense-only, unbatched, and bound to `flashier.utils::project_factor_activity`. A cebmf_torch projection oracle and cell batching are net-new. For the cNMF arm, `cnmf.py:776` `refit_usage` is the ready NNLS projection analog. |
| **Refinement (anchored refit at single-cell/metacell resolution)** | **MOSTLY NET-NEW** (registered plan + hook exist, no code) | MR-4 (`reports/STAGE_B_ROADMAP.md:60-70`) specifies anchored refinement via flashier.utils `fit_initialized_ebmf()` with a registered adoption rule; no implementation in R/. cebmf_torch-side initialization from Stage-A programs is net-new; `R/stability.R:31` and the model_selection benchmark provide the drift/ablation harness. |

**Natural module boundaries** (per the stated responsibility split):
- `cebmf_torch` — numerics only: EBMF fit, canonicalization math, projection oracle, anchored refit. No I/O formats, no cell-type semantics.
- `flashier-utils` (R-facing backend interface) — the only layer cellprograms talks to: `fit_ebmf`, `project_factor_activity`, `fit_initialized_ebmf`; manifests/replay records live here. cebmf_torch becomes a backend behind this interface (or its successor) so `R/ebmf.R:467` routing is the single change point.
- `cellprograms` — data prep (pseudobulk/metacell construction promoted out of `benchmarks/`), workflow orchestration (discovery -> projection -> refinement), and biological validation (`R/cmp-interpret.R`, `R/coordination.R`, `R/program-sharing.R`).

## (e) Circular-dependency risks

1. **cellprograms <-> flashier-utils is already a hard bidirectional coupling in spirit**: cellprograms Imports `flashier.utils` (DESCRIPTION:22,35) and its projection path fails closed without it (`R/cell-projection.R:169-171`). If the new flashier-utils wraps cebmf_torch, the layering cebmf_torch -> flashier-utils -> cellprograms must stay strictly one-way; cellprograms must never import cebmf_torch directly or the backend abstraction inverts.
2. **Vendored fork drift**: `benchmarks/biomni_replication/fork/R/*.R` is a full copy of the ebmf package inside the repo. New workflow code must import the package, not the fork, or fixes diverge silently.
3. **Version pin**: `DESCRIPTION:34` pins flashier to commit `fd44811`; a cebmf_torch-backed backend bump must be a deliberate pin change, not an implicit install.
4. **Benchmark-scale mismatch**: cNMF's non-negative counts/TPM contract vs cellprograms' log-normalized pseudobulk (see (c) caveat) — if the multiresolution benchmark feeds one matrix to both, one method is silently advantaged; the MR-2 registration rule (`reports/STAGE_B_ROADMAP.md:42-49`) should state the per-method input scale explicitly.
5. **Metacell grouting before projection**: MR-3/MR-4 assume programs are fit before projecting remaining cells; if metacell assignment is computed from a projection of the same programs, a fit/project/regrout loop can become self-referential — keep grouting independent of the programs being evaluated (the roadmap's "projection of remaining cells" wording already implies this).
