<div class="page-header">

# Changelog

</div>

<div class="section level2">

## cellprograms 0.2.0.9001 (development)<a href="#cellprograms-0209001-development" class="anchor" aria-label="anchor"></a>

<div class="section level3">

### Repository completion<a href="#repository-completion-0-2-0-9001" class="anchor" aria-label="anchor"></a>

- Promoted the real ordinary/native-fit, serialization and actual stability acceptance fixtures into the maintained installed test suite.
- Replaced placeholder maintainer metadata with the existing public owner contact. The development README now separates supported core behavior, experimental integrations, historical evidence and open scientific gates.
- Reconciled concurrent upstream cNMF/semi-base fixes and preserved its negative coordination decision; these are not new validation claims by this change.
- Ordinary CI installs hard dependencies explicitly, limits threads, and permits only documented optional-Suggests notes. Cross-platform coverage requires actual CI outcomes; development checks do not constitute CRAN acceptance.

</div>

<div class="section level3">

### Bounded ordinary-core acceptance<a href="#bounded-ordinary-core-acceptance-0-2-0-9001" class="anchor" aria-label="anchor"></a>

- Projection now rejects unknown or missing reference measurement units before native inference and requires refitting with truthful units declared. Projection-time `feature_scale` declares query units only; it cannot relabel reference provenance (delegated implementation decision TR03-D001).
- The installed suite includes the byte-identical historical node-count accounting oracle. Its fixture no longer becomes an empty installed test. This is accounting validation, not a calibration campaign.
- Known warning fixtures assert their complete expected messages and recorded fallback settings rather than leaving extra warnings unexplained.

</div>

<div class="section level3">

### CPUB correctness and support corrections (local candidate, 2026-10-22)<a href="#cpub-correctness-and-support-corrections-local-candidate-0-2-0-9001" class="anchor" aria-label="anchor"></a>

- EB-CP now orients two modes, preserving each component’s tensor product. Reconstruction diagnostics are computed from the returned factors. Held-out selection is experimental screening, not a winner’s-curse or set-error guarantee; this supersedes stronger wording in the historical entry below. Unusable continuous/singleton strata and insufficient per-type split halves are refused before solving.
- Covariate residualization records trained terms, prediction variables, factor levels, contrasts and numeric scaling constants. Linear, factor, polynomial and spline query designs use training state. Unsupported custom transforms and legacy residualized models without state require explicit preprocessing/refitting, rather than query-dependent recomputation.
- Frozen preprocessing schema 1.1 records feature units and validates its content digest on consumption. Intact schema 1.0 remains readable with unknown units.
- Known residual noise now uses `S` together with `var_type = NULL` internally; passing `S` alone did not disable noise estimation in the pinned engine. Effective engine, priors, variance mode and fallback controls are recorded.
- Experimental `cebmf_torch` is explicitly a distinct Laplace/Laplace estimator: use `loading_prior="laplace"`, `var_type=2`, `backfit=TRUE`, `nullcheck=FALSE`. Incompatible controls fail before Python. Loss increases cannot establish plateau convergence. Only `flash_control$verbose` remains accepted; estimator controls must use their explicit arguments.
- Score-level block fits retain block masses; score/private-first is refused. Observed-support CMP residualization uses rank-aware QR and exposes program and edge callability. These guards do not establish new calibration.
- Stability refits forward backend/noise controls in frozen and unfrozen paths and retain intended settings for failed as well as successful replicates.
- Optional metacell reduction preparation returns the updated Seurat object with a consistent PCA key. Full SuperCell/mcRigor integration remains optional.
- Added the installed-package test entrypoint. CI now configures minimum R, current Linux/macOS/Windows and an upstream-dependency lane; configuration is not evidence that these platforms have executed successfully.

</div>

<div class="section level3">

### Batch adjustment, subspace metrics, and covariate-audit tests (2026-10-22)<a href="#batch-adjustment-subspace-metrics-and-covariate-audit-tests-0-2-0-9001" class="anchor" aria-label="anchor"></a>

- New `adjust_batches=` argument in [`fit_celltype_programs()`](../reference/fit_celltype_programs.html): `sva::ComBat` removes batch location/scale from each cell type’s panel before covariate residualization. Intended for small-batch designs where OLS residualization is underpowered (~5 donors/batch). Not the default; unsupported with a frozen `preprocessing_spec`; [`project_programs_to_cells()`](../reference/project_programs_to_cells.html) refuses batch-adjusted fits (the adjustment cannot yet be replayed on a query).
- New [`span_alignment()`](../reference/span_alignment.html): principal-angle subspace metrics between two loading sets (mean overlap, top-k overlap) — the discovery-vs-refit span-conservation diagnostic.
- [`fit_joint_block()`](../reference/fit_joint_block.html) no longer crashes when the shared-first pass finds zero factors: the zero reconstruction now carries dimnames.
- New Suggests: sva (ComBat reference implementation, used by the covariate-audit regression tests and cross-checked against `adjust_batches`).
- Ported the covariate-audit regression suite (`test-covariates.R`) and a new batch-subspace suite (`test-batch-subspace.R`, 15 checks) covering the replay center/scale, NA-coefficient, fixed-mode centering, and interaction-drop fixes already in main since 2026-10-19.
- Benchmark-only: `block_E_cal` arm in `benchmarks/coordination/` — a no-refit FWER permutation layer on block_E’s calls (within-strata donor-row permutation, max-T across clusters). Reduces but does not eliminate the raw block-mass rule’s null false positives (double-dipping: member loadings are fit on member expression). Confirmatory run preregistered as `benchmarks/coordination/PREREG_confirm_v1.md`.
- Registrations filed: `PREREG_confirm_v1.md` (block_E_cal confirmatory) and `benchmarks/cmp_detection/PREREG_mr2-cnmf-v1.md` (cNMF baseline vs the Stage A operating point); `res/hallmark.csv` and `res/lr_consensus.csv` reconstructed (msigdbr 26.1.1; CellChatDB, 1020 ligand-receptor pairs).
- Corrected the Stephenson cohort accession in the mofacellular benchmark README: E-MTAB-8731 (wrong, unrelated mouse study) -\> E-MTAB-10026.

</div>

<div class="section level3">

### Metacells-to-programs vignette (2026-10-22)<a href="#metacells-to-programs-vignette-0-2-0-9001" class="anchor" aria-label="anchor"></a>

- New article `vignettes/articles/metacells_to_programs.Rmd`: the full downstream workflow after SuperCell 2.0 + mcRigor metacell construction — metacell counts -\> donor pseudobulk (sum mode is lossless at the donor level) -\> `fit_celltype_programs(feature_scale = "logcounts")` -\> canonicalization -\> stability -\> CMP detection (planted coordinated pair recovered, private programs not called) -\> [`program_sharing()`](../reference/program_sharing.html) -\> [`project_programs_to_cells()`](../reference/project_programs_to_cells.html) back onto the original single cells (per-cell activities recover the planted donor gradient).

</div>

<div class="section level3">

### EB-CP tensor factorization and block EBMF (2026-10-21)<a href="#eb-cp-tensor-factorization-and-block-ebmf-0-2-0-9001" class="anchor" aria-label="anchor"></a>

- New [`fit_eb_cp()`](../reference/fit_eb_cp.html): greedy rank-1 empirical-Bayes CANDECOMP/PARAFAC factorization across cell types (genes x donors x cell types) with EBNM normal-means updates per mode (point-Laplace on gene and cell-type modes, normal on donors) and per-(cell-type, donor) residual variance — the tensor analogue of flashier `var_type = 1`. Components are selected by a winner’s-curse-proof held-out likelihood test on a donor-synchronous 50/50 split and refit on the full data with backfitting.
- New [`eb_cp_from_fit()`](../reference/eb_cp_from_fit.html): runs the tensor factorization on the adjusted expression matrices (`Y_used`) of a `cell_program_fit`; the held-out split is stratified on the fit’s covariates by default (stratification is load-bearing: an unstratified split leaves residual covariate imbalance that the greedy pass spends components on).
- New [`eb_cp_programs()`](../reference/eb_cp_programs.html): per-component cell-type memberships (member iff \|W\| \>= 0.2 of the component maximum) and donor activities; components with \>= 2 member cell types are coordination calls.
- New [`fit_joint_block()`](../reference/fit_joint_block.html) / [`block_clusters()`](../reference/block_clusters.html): joint block EBMF at expression or scores level (shared-first greedy + private pass), with block-mass-based membership extraction.
- New benchmark suite `benchmarks/coordination/`: simulator with correlated-private decoys (`priv_corr`) and each-alone breadth families, six methods (pa_cc, node_cor, subspace_sum, gene_match, block_E, eb_cp) + oracle floor, scorer with Hungarian matching, resumable replicate harness, runtime/amplitude pilots, and the filled preregistration (`benchmarks/coordination/PREREG.md`, prereg-coord-v1).

</div>

<div class="section level3">

### SuperCell backend (2026-10-22)<a href="#supercell-backend-0-2-0-9001" class="anchor" aria-label="anchor"></a>

- New default R-native metacell path implementing the preregistered benchmark verdict (REPLACE): [`MetacellsByGroupsSuperCell()`](../reference/MetacellsByGroupsSuperCell.html), [`OptimizeGammaSuperCell()`](../reference/OptimizeGammaSuperCell.html) and [`GetMetacellSuperCell()`](../reference/GetMetacellSuperCell.html) wrap SuperCell 2.0 (`SCimplify_for_Seurat_v5`, pinned `GfellerLab/SuperCell@89b34c0`) with the same grouping semantics, mcRigor QC and return structure as the SEACells trio. No Python/basilisk required, and the end-to-end tests run in the plain R test suite (the SEACells e2e stays `CPMC_INTEGRATION`-gated). `target_metacells` is reached by log-gamma bisection; [`OptimizeGammaSuperCell()`](../reference/OptimizeGammaSuperCell.html) optimizes native cells-per-metacell by the mcRigor score. The user-facing surface grows from 42 to 45 verbs.
- [`MetacellsByGroupsSEACells()`](../reference/MetacellsByGroupsSEACells.html) and the SEACells backend remain exported and unchanged.
- mcRigor wrapper fix (both backends): `.mcrigor_detect_group()` now passes MC-prefix string labels named by cell id to `mcRigor_DETECT`. Integer codes are coerced to NA by DETECT’s named-membership filter at mcRigor c70cd6f and aborted with “Cannot find cells provided”; the string convention is what the benchmark validated. \# cellprograms 0.2.0.9000 (development)

</div>

<div class="section level3">

### Metacell benchmark (2026-10-22)<a href="#metacell-benchmark-0-2-0-9001" class="anchor" aria-label="anchor"></a>

- Preregistered benchmark of SEACells vs SuperCell 2.0 (pure R) with mcRigor QC across 2 scales x 3 seeds x 4 granularities (`benchmarks/metacell_comparison/`). Verdict: **REPLACE** - SuperCell 2.0 + mcRigor is the default R-native metacell path (mean \|ARI gap\| 0.019 \<= 0.05; mean relative fidelity gap 5.4% \<= 10%; membership construction ~200x faster; downstream program fit is N- and method-insensitive). SEACells remains available and unchanged.

</div>

<div class="section level3">

### API streamlining (2026-10-21)<a href="#api-streamlining-0-2-0-9001" class="anchor" aria-label="anchor"></a>

- Eight internal helpers are no longer exported (still callable via `cellprograms:::`; nothing is deleted): `program_metadata` (derivable from [`program_loadings()`](../reference/program_loadings.html) column names), `principal_angles` (pairwise primitive behind [`align_programs()`](../reference/align_programs.html) / [`match_programs()`](../reference/match_programs.html)), `cmp_scores`, `cmp_gene_overlap`, `cmp_context_association` and `cmp_nuisance_genes` (components orchestrated by [`cmp_triage()`](../reference/cmp_triage.html)), `union_panel_loadings` and `ambient_program_scores` (coordination machinery). The user-facing surface drops from 50 to 42 verbs.
- The pkgdown reference index is reorganized from a flat 50-function list into 8 curated family sections (Data & simulation / Fit & canonicalize / Extract / Compare, align & share / Coordination & annotation / Metacells / Program plots / Coordination & sharing plots), following the mixOmics reference layout.

</div>

<div class="section level3">

### Bug fixes and hardening (2026-10-21 all-exports smoke sweep)<a href="#bug-fixes-and-hardening-all-exports-smoke-sweep-0-2-0-9001" class="anchor" aria-label="anchor"></a>

- [`plot_cmp_network()`](../reference/plot_cmp_network.html): no longer crashes on CMP results computed with `ambient = "none"` — those results carry no `ambient_score` column (nothing was excluded), and the plot assumed it existed (`setNames` length mismatch). Nodes are now unflagged for ambient in that case.
- Fit consumers fail closed on a non-canonicalized fit: [`program_summary()`](../reference/program_summary.html), [`program_metadata()`](../reference/program_metadata.html), [`program_scores()`](../reference/program_scores.html), [`program_loadings()`](../reference/program_loadings.html), [`top_genes()`](../reference/top_genes.html), [`assess_program_stability()`](../reference/assess_program_stability.html), [`principal_angles()`](../reference/principal_angles.html), [`sharing_spectrum()`](../reference/sharing_spectrum.html), [`align_programs()`](../reference/align_programs.html) and [`match_programs()`](../reference/match_programs.html) now abort with “run canonicalize_programs() first” instead of the opaque “argument is of length zero” (`ncol(NULL)` on the unpopulated scores list) or silently empty tables. Three entry points already had this guard ([`find_coordinated_programs()`](../reference/find_coordinated_programs.html), [`project_programs_to_cells()`](../reference/project_programs_to_cells.html), [`plot_variance_explained()`](../reference/plot_variance_explained.html)); the message is now uniform across the API.
- Regression test: every guarded consumer rejects a raw [`fit_celltype_programs()`](../reference/fit_celltype_programs.html) output with the fix named, and passes a canonicalized fit through unchanged.

</div>

<div class="section level3">

### Bug fixes and hardening (2026-10-20 edge-case segment)<a href="#bug-fixes-and-hardening-edge-case-segment-0-2-0-9001" class="anchor" aria-label="anchor"></a>

- [`OptimizeGammaSEACells()`](../reference/OptimizeGammaSEACells.html): `gammas` is now validated (positive integers) and deduplicated/sorted at entry — duplicate candidates previously produced identically named membership columns for `mcRigor_OPTIMIZE`.
- `.mc_get_group_ids()` (both metacell entry points): fails closed when a `group.by` column is missing or contains NA — <a href="https://rdrr.io/r/base/paste.html" class="external-link"><code>paste()</code></a> previously rendered NA as the literal string `"NA"`, silently lumping unassigned cells into a spurious group.
- [`OptimizeGammaSEACells()`](../reference/OptimizeGammaSEACells.html) membership coding: cells with NA membership labels in any candidate gamma are excluded with a warning before scoring; <a href="https://rdrr.io/r/base/unique.html" class="external-link"><code>unique()</code></a> previously counted NA as a level, inflating the degeneracy check and leaking NA codes into mcRigor.
- `.cebmf_use_python_env()`: an existing directory without `bin/python` now aborts with the registered message instead of falling through to reticulate (<a href="https://rdrr.io/r/base/files.html" class="external-link"><code>file.exists()</code></a> is TRUE for directories).
- Regression test: the cebmf basilisk pip spec is checked for the invariants the CI shakedown established (no flag entries, no whitespace, `==` pins or commit-pinned direct references only).
- Packaging: `hdWGCNA` removed from Suggests (never imported; design attribution lives in comments). Its presence broke every CI dependency solve — pak cannot resolve GitHub-only packages from CRAN/Bioconductor.
- CI: workflow actions bumped to Node 24 majors (checkout v7, upload-artifact v7, github-script v9, upload-pages-artifact v5, deploy-pages v5); Node 20 was removed from runners on 2026-09-23. The drift issue now also auto-closes on a green `workflow_dispatch` run, not just scheduled runs.

</div>

<div class="section level3">

### Bug fixes and hardening (2026-10-20 cebmf basilisk segment)<a href="#bug-fixes-and-hardening-cebmf-basilisk-segment-0-2-0-9001" class="anchor" aria-label="anchor"></a>

- The cebmf_torch backend’s basilisk environment was never constructible and the defect was masked by every audit run setting the venv option. Three fixes: the pip spec now pins `numpy`/`scipy` (which `cebmf_ext` imports directly and torch does not provide), the CPU torch wheel is selected via `--extra-index-url` as its own pip token (index flags are global pip options, not requirement syntax), and the commit-pinned `cebmf_torch` direct reference is written unspaced so the shell cannot split it when basilisk passes pip entries through `system2` unquoted. The basilisk version-spec check is relaxed only while this one environment is created (it cannot parse PEP 508 direct references; the commit pin is fully reproducible). The env now builds and the backend fits through it with no option set.

</div>

<div class="section level3">

### Bug fixes and hardening (2026-10-19 implementation-review segment)<a href="#bug-fixes-and-hardening-implementation-review-segment-0-2-0-9001" class="anchor" aria-label="anchor"></a>

Full implementation review follow-through: 2 critical correctness fixes, 12 major fixes, and 50 minor fixes. Regression coverage for every fixed finding lives in `tests/testthat/test-review-fixes.R` and `tests/testthat/test-artifact-genes.R`.

<div class="section level4">

#### Critical<a href="#critical-0-2-0-9001" class="anchor" aria-label="anchor"></a>

- [`GetMetacellSEACells()`](../reference/GetMetacellSEACells.html): multi-group soft-assignment matrices are now combined block-diagonally (one block per group) instead of row-binding group-local column indices, which crashed or silently misassigned cells to metacells when more than one group was processed.
- [`cmp_context_association()`](../reference/cmp_context_association.html): the context term is selected exactly through the model-matrix term assignment instead of a `^x` prefix grep that also captured adjustment covariates whose names start with “x” (e.g. `xsite`). Covariates named `y` or `.ctx` are refused (reserved internal names).

</div>

<div class="section level4">

#### Major<a href="#major-0-2-0-9001" class="anchor" aria-label="anchor"></a>

- [`fit_celltype_programs()`](../reference/fit_celltype_programs.html): within-cell-type constant genes no longer crash the gene filter under `scale = TRUE`; their non-finite column mass is treated as degenerate and the genes are dropped.
- `fit_celltype_programs(backend = "cebmf_torch")` warns when a cell type’s fit exhausts `maxit` without converging (the flag was stored but never surfaced).
- [`assess_program_stability()`](../reference/assess_program_stability.html) warns and skips cell types with fewer than 5 donors instead of crashing during resampling.
- [`plot_program_scores()`](../reference/plot_program_scores.html) handles per-cell-type observation sets (union of donors, NA where a cell type was not observed) instead of requiring equal row counts; z-scoring is NA-aware.
- [`plot_sharing_spectrum()`](../reference/plot_sharing_spectrum.html) rewritten: significance class (shared / n.s. / underpowered) drives a fixed colorblind-safe palette, q-values are carried in the pair labels, and \>9 pairs no longer break the color scale.
- [`plot_artifact_flags()`](../reference/plot_artifact_flags.html) detects any `*_frac` column set (not just the built-in screens) and accepts CMP-shaped flag tables.
- [`union_panel_loadings()`](../reference/union_panel_loadings.html) leaves zero-overlap or constant panel genes NA (with warnings) instead of aborting, and replays rank-deficient fixed covariate designs via least squares instead of crashing in <a href="https://rdrr.io/r/base/solve.html" class="external-link"><code>solve()</code></a>.
- [`project_programs_to_cells()`](../reference/project_programs_to_cells.html) refuses a `feature_scale` declaration that contradicts the fit’s recorded units (previously silently accepted).
- [`flag_artifact_genes()`](../reference/flag_artifact_genes.html): the V(D)J screen now also covers J/D segments, IGLL-family genes, and JCHAIN.
- Metacells: `max_iter` documentation corrected (wrapped SEACells default is 20); `use_sparse` warns when set (ignored by the wrapped SEACells version); per-metacell QC columns are joined into the result metadata.
- Metacell IO coercions use the non-deprecated `CsparseMatrix` form throughout: the pinned flashier’s `.onLoad` sets `Matrix.warnDeprecatedCoerce = 2`, under which the previous `dgCMatrix` coercions of `readMM()` output error out (exposed by the first `CPMC_INTEGRATION=1` end-to-end run, which now passes).
- `withr` moved from Suggests to Imports: <a href="https://withr.r-lib.org/reference/with_seed.html" class="external-link"><code>withr::local_seed()</code></a> is called unconditionally by the fitting code since the scoped-RNG fix.
- [`OptimizeGammaSEACells()`](../reference/OptimizeGammaSEACells.html) now works under Seurat v5: mcRigor’s `buildmc` aligns per-metacell sizes through `AggregateExpression` colnames, which mangle membership labels (digit-leading values gain a “g” prefix; “\_” becomes “-”) and only undo the digit case — the `"mc_<k>"` labels used here never realigned, so every group failed with “No cell overlap between new meta data and Seurat object”. Membership is now integer-coded per gamma column (the same contract as the `mcRigor_DETECT` path) and gamma candidates whose membership collapses to a single metacell are dropped with a warning (`NA` in `scores`) instead of failing the whole group’s optimization.
- [`MetacellsByGroupsSEACells()`](../reference/MetacellsByGroupsSEACells.html) no longer crashes when SEACells collapses a group to fewer metacells than archetypes fitted (degenerate fits on weak structure): the soft-assignment matrix has one column per archetype and cannot be renamed to the realized metacell names, so it is dropped for that group with a warning instead of erroring with “length of Dimnames\[\[2\]\] is not equal to Dim\[2\]”.

</div>

<div class="section level4">

#### Minor (selected)<a href="#minor-selected-0-2-0-9001" class="anchor" aria-label="anchor"></a>

- All internal RNG is now scoped with <a href="https://withr.r-lib.org/reference/with_seed.html" class="external-link"><code>withr::local_seed()</code></a>; no <a href="https://rdrr.io/r/base/Random.html" class="external-link"><code>set.seed()</code></a> leaks into the user session.
- [`plot_variance_explained()`](../reference/plot_variance_explained.html) maps PVEs to programs through the flashbridge native-factor mapping and skips (with a warning) cell types whose PVE vector cannot be reconciled, instead of positionally truncating.
- [`plot_context_association()`](../reference/plot_context_association.html) omits confidence intervals for associations with a zero t statistic (previously infinite whiskers) and labels the context term plainly after the `.ctx` internal rename.
- [`plot_program_scores()`](../reference/plot_program_scores.html) fails closed when an annotation data.frame lacks an `observation_id` column (previously a silent all-NA strip).
- [`match_programs()`](../reference/match_programs.html) documents that `space` only affects the subspace angles (factor-level matches are always loading-cosine based) and fails closed when `space = "scores"` has fewer than 2 shared observations.
- [`as_mofacellular_input()`](../reference/as_mofacellular_input.html) aborts on composite-name collisions under the `<cell_type>_<observation_id>` scheme, on duplicate `cell_counts` keys, and on zero-observation cell types.
- Metacell IO: missing counts-like layers warn or abort, `dims` are validated, and numeric-like cell names abort (they would be mangled by the Python bridge).
- CMP machinery: sign-concordance ties and all-NA sign vectors no longer score deterministically; [`prepare_collectri()`](../reference/prepare_collectri.html) parses integer 0/1 `mor` correctly; NA `mor` edges are dropped with a warning; [`cmp_gene_overlap()`](../reference/cmp_gene_overlap.html) is schema-stable on empty and single-member input.
- Sharing/alignment: degenerate-program NA propagation fixed; [`sharing_spectrum()`](../reference/sharing_spectrum.html) gains an explicit `seed` argument; the variance estimate reuses the sharing matrices (single computation, consistent panels); principal-angle p-value NA handling no longer biases the permutation null.

</div>

</div>

<div class="section level3">

### New features (2026-10-18 feature-completion segment)<a href="#new-features-feature-completion-segment-0-2-0-9001" class="anchor" aria-label="anchor"></a>

- Metacell construction is folded into `cellprograms` from the retired `r_packages/cellprogramsMC` sub-package: [`MetacellsByGroupsSEACells()`](../reference/MetacellsByGroupsSEACells.html) (group-stratified SEACells metacells with mcRigor QC), [`OptimizeGammaSEACells()`](../reference/OptimizeGammaSEACells.html) (mcRigor granularity selection), and [`GetMetacellSEACells()`](../reference/GetMetacellSEACells.html) (result accessor). The basilisk Python environment is now created lazily on first use (internal `.seacells_env()`, same pattern as the cebmf backend), so `basilisk` stays in Suggests and installing `cellprograms` never requires basilisk, Python, Seurat, or mcRigor; the exported `seacells_env` object from the sub-package is retired (nothing in the project referenced it). `Seurat`, `SeuratObject`, and `mcRigor` are Suggests with explicit fail-closed guards; `Matrix` and `methods` move to Imports. Note: the basilisk environment is re-registered under `cellprograms`, so the first metacell run after this change rebuilds the SEACells conda environment (one-time, needs network).
- [`fit_celltype_programs()`](../reference/fit_celltype_programs.html) gains `backend = c("flashier", "cebmf_torch")`. The `cebmf_torch` backend runs the pure-PyTorch cEBMF implementation (pinned upstream `william-denault/cebmf_torch@34abac5`) through the `cebmf_ext` adapter, with a registered numerical guard set (sebetahat floor 1e-4, tau precision cap 1e8, explicit seeding, seed-family fallback with diagnostics) and a strict mapping contract: the bridge never manufactures flashier-only quantities, so `elbo`, `L_psd`, `F_psd`, and `F_lfsr` are NULL on this backend and [`program_summary()`](../reference/program_summary.html) reports `elbo = NA`. Fixed-mode covariates are refused on this backend (no fixed-loadings path); use the default `covariate_mode = "residualize"`. Python resolution: `options(cellprograms.cebmf_python_env = ...)` (python binary or venv with `cebmf_ext` importable), else basilisk (pinned `torch==2.7.1` CPU — the pin is load-bearing, 2.14.1+cpu diverges on inputs where 2.7.1+cpu is green), else a fail-closed error. No flashier-style fallback ladder: backend failures abort loudly per cell type.
- [`fit_celltype_programs()`](../reference/fit_celltype_programs.html) gains `var_floor = "none" | "vignette" | <positive scalar>`: fix the residual SD instead of estimating it. `"vignette"` uses the flashier single-cell vignette’s floor, `S = sd(log(Poisson(1/n) + 1))` with `n` the cell type’s donor count, drawn once per cell type at a fixed local seed (deterministic, fit-seed-independent). A `var_floor`-derived S is not mappable through <a href="https://mdmanurung.github.io/flashier-utils/reference/fit_ebmf.html" class="external-link"><code>flashier.utils::fit_ebmf()</code></a>, so on the fast path the fit reroutes through <a href="https://rdrr.io/pkg/flashier/man/flash.html" class="external-link"><code>flashier::flash()</code></a> with a warning; the cebmf_torch backend passes S through as KNOWN noise. The default `"none"` preserves existing behavior exactly. Calibration evidence: CAND-13 demonstrator (`reports/audit/E_programs/`).
- Covariate-audit patches (R2/F1-F4): gene-drop blocks now subset the projection-replay `center`/`scale` vectors to retained genes (F1); `covariate_mode = "fixed"` with `center = FALSE` is refused even when the combination arrives through a frozen preprocessing spec (F3); covariate terms referencing within-cell-type-constant variables — including interaction terms — are dropped instead of crashing <a href="https://rdrr.io/r/stats/model.matrix.html" class="external-link"><code>model.matrix()</code></a> (F4); aliased rank-deficient replay designs zero their NA coefficients at projection instead of propagating NaN activities (F2).

</div>

<div class="section level3">

### New features (spec v2 execution segment)<a href="#new-features-spec-v2-execution-segment-0-2-0-9001" class="anchor" aria-label="anchor"></a>

- [`assess_cmp_stability()`](../reference/assess_cmp_stability.html) (CP-2): stability of coordinated multicellular programs under SYNCHRONIZED donor resampling. Donors are drawn once per replicate (optionally within donor-level `donor_strata`) and the identical donor set is applied to every cell type — independent per-type subsampling ([`assess_program_stability()`](../reference/assess_program_stability.html)) changes donor overlap arbitrarily and cannot answer whether the same donors support a CMP. Resampling is subsampling without replacement (`subsample_frac`, default 0.8): duplicated donor rows would pseudo-replicate the donor-level evidence the permutation null evaluates. `donor_col` promotes resampling to the biological donor when it differs from `observation_ids` (all of a donor’s observations are kept together; missing cell types/visits are never synthesized); `visit_col` is required when donors contribute \> 1 observation, making the paired-observation estimand explicit. Each replicate reruns [`find_coordinated_programs()`](../reference/find_coordinated_programs.html) with the caller’s predeclared `cmp_args` and reports per-CMP membership recovery (Jaccard), assignment-ambiguity margins, sign concordance up to one global flip, donor-score concordance, cell-type survival, and CP-1-style failure accounting (`valid_refits` / `failed_refits`, conditional rates, NA on all-fail). When the reference fit yields no CMPs, stability is reported as ineligible rather than zero.

</div>

<div class="section level3">

### Audit remediations<a href="#audit-remediations-0-2-0-9001" class="anchor" aria-label="anchor"></a>

- [`find_coordinated_programs()`](../reference/find_coordinated_programs.html): missing donor scores (NA) are no longer silently converted into top ranks. Rank transformation and strata residualization keep NA donors as NA, cross-cell-type correlations use pairwise-complete donors, edge `n` reports the pairwise-complete donor count for each program pair, and a warning lists programs with missing scores. Previously, <a href="https://rdrr.io/r/base/rank.html" class="external-link"><code>rank()</code></a>’s default `na.last=TRUE` manufactured fake top-scoring donors, which either destroyed detection power (donor-specific missingness) or inflated significance (missingness shared across programs); a planted pair with 20% missingness could be missed entirely with no warning. Regression tests: `tests/testthat/test-cmp-missingness.R`.

</div>

<div class="section level3">

### Remediation queue (2026-10-18): CAND-01, CAND-05, CAND-07, D-2a, CAND-03+D-3b, CAND-12<a href="#remediation-queue-cand-cand-cand-da-canddb-cand-0-2-0-9001" class="anchor" aria-label="anchor"></a>

- [`assess_program_stability()`](../reference/assess_program_stability.html) (CAND-01/CP-1, P1): subsample refit failures are now counted and reported (`fit$stability$valid_refits`/`failed_refits`, warning on any failure) with a per-replicate `replicate_log` (status, error class/message, elapsed time, sampled observations and factor counts per cell type, matching eligibility). `recovery_freq` keeps its LEGACY meaning (unconditional hit fraction over requested iterations, `n_hit / n_boot`); the new `conditional_recovery_freq` is defined over SUCCESSFUL refits (`n_hit / valid_refits`) — previously the silent denominator was `n_boot`, so a 30% refit-failure rate read as a 30% recovery loss with zero user-facing signal, and an all-fail run was indistinguishable from zero stability. All-fail now yields NA rates with a distinct warning; per-program `n_hit` and `n_eligible` added. Subsampled refit objects now carry consistent `observation_ids`/`sample_metadata`/`dims` (previously stale full-data values).
- [`find_coordinated_programs()`](../reference/find_coordinated_programs.html) (CAND-05, P2): CMP member orientation now propagates along the linking edges that formed the CMP (product of edge signs along the path from the reference member). Chain-joined members previously inherited the sign of the direct reference–member edge from the full edge table — which can be a non-linking noise edge, i.e. a random +/- that propagated into [`cmp_scores()`](../reference/cmp_scores.html), [`annotate_cmp_tf()`](../reference/annotate_cmp_tf.html) and [`cmp_triage()`](../reference/cmp_triage.html) direction calls.
- [`program_sharing()`](../reference/program_sharing.html) (D-2a/D-2c, P1): `matches=` references are validated (unknown programs abort informatively instead of crashing in the mashr layer) and every forced-matches call warns that the “shared” label counts top genes ACTIVE (lfsr \< thresh) in \>= 2 conditions — which independent loadings can also reach at high loading SNR (audit D-2a: 30/30 false “shared” for forced independent pairs at loading SE 0.125, while mashr’s own mixture put 96% weight on the identity pattern, p_all_shared ~ 0.001); `p_all_shared` is the diagnostic. The classification rule itself is unchanged pending a registered benchmark of any replacement.
- Matching-step documentation (CAND-07, P1): the gene-permutation null’s exchangeability caveat is now stated where matching q-values are produced (mirroring `cmp-interpret.R`): block-correlated genes inflate false matches (audit D-1: 67.0% type-I under block equicorrelation 0.55 vs 2.0% iid). The null itself is UNCHANGED — a registered mini-benchmark (null calibration at block equicorrelation + planted-truth power) precedes any null change.
- [`as_cell_program_data()`](../reference/as_cell_program_data.html) (CAND-03, P2): rejects ALL non-finite entries (NA/NaN/Inf) with a count-bearing abort; Inf previously passed the NA-only check and crashed validate/fit with an opaque “missing value where TRUE/FALSE needed”.
- `.greedy_match()` (D-3b, P2): non-finite similarities are demoted to “unmatched” with a warning instead of crashing on <a href="https://rdrr.io/r/base/Extremes.html" class="external-link"><code>max()</code></a> of an NA-bearing matrix.
- [`fit_celltype_programs()`](../reference/fit_celltype_programs.html) + [`assess_program_stability()`](../reference/assess_program_stability.html) (CAND-12, P2): the full frozen `preprocessing_spec` is stored on the fit (not just its digest) and stability refits REPLAY it through the MR-0 fail-closed checks; `feature_scale` joined the MR-0 conflict list and is passed through in legacy-path refits (previously silently “unknown” in refit manifests).

</div>

</div>

<div class="section level2">

## cellprograms 0.2.0<a href="#cellprograms-020" class="anchor" aria-label="anchor"></a>

<div class="section level3">

### Stage B: frozen preprocessing and cell-level projection<a href="#stage-b-frozen-preprocessing-and-cell-level-projection-0-2-0" class="anchor" aria-label="anchor"></a>

- [`freeze_preprocessing()`](../reference/freeze_preprocessing.html) (MR-0): resolve the feature panel and adjustment decisions once into a content-addressed specification (`cell_program_preproc_spec`, sha256 digest); refit any cohort carrying the same genes with `fit_celltype_programs(preprocessing_spec = ...)`, which replays the spec fail-closed (conflicting explicit arguments are refused, and the spec’s digest is recorded in the fit’s preprocessing provenance).
- [`project_programs_to_cells()`](../reference/project_programs_to_cells.html) (MR-1): project fitted programs onto individual cells. Query expression is adjusted with the fit’s recorded preprocessing (residualization coefficients, centering, scale), activities are inferred against the frozen gene loadings via <a href="https://mdmanurung.github.io/flashier-utils/reference/project_factor_activity.html" class="external-link"><code>flashier.utils::project_factor_activity()</code></a>, and raw activities are canonicalized with the same RMS/sign multipliers as the donor scores. Requires a canonicalized fit with declared `feature_scale`; fits with residualized covariates require `donors` + `sample_metadata`, and the rebuilt covariate design must match the fit’s design columns exactly.
- [`fit_celltype_programs()`](../reference/fit_celltype_programs.html) gains `feature_scale` (declared measurement units, recorded in the fit manifest and required for projection) and `preprocessing_spec` (MR-0 replay); fits now record a `preproc_replay` block per cell type so the exact adjustment can be replayed on query data.
- New roadmap: `reports/STAGE_B_ROADMAP.md` (MR-0..MR-4 spec, staged main-run costs, registration drafts, gates).

</div>

<div class="section level3">

### Backend synchronization (Stage A)<a href="#backend-synchronization-stage-a-0-2-0" class="anchor" aria-label="anchor"></a>

- Ordinary EBMF fits now delegate to <a href="https://mdmanurung.github.io/flashier-utils/reference/fit_ebmf.html" class="external-link"><code>flashier.utils::fit_ebmf()</code></a> and canonicalization uses `flashier.utils::standardize_factors(scaling = "rms")`, with numerical parity verified to machine precision against the previous direct <a href="https://rdrr.io/pkg/flashier/man/flash.html" class="external-link"><code>flashier::flash()</code></a> backend. The public workflow, output fields, and scientific results are unchanged. `flashier`, `ebnm`, and `flashier.utils` move to `Imports` (R \>= 4.3 required); the fixed-covariate path (`covariate_mode = "fixed"`) and degenerate fits remain on the legacy backend. Set `options(cellprograms.use_flashier_utils = FALSE)` to force the legacy backend. See `reports/SYNC_ROADMAP.md` for the version pair, parity evidence, and rollback plan.

</div>

<div class="section level3">

### New features<a href="#new-features-0-2-0" class="anchor" aria-label="anchor"></a>

- [`find_coordinated_programs()`](../reference/find_coordinated_programs.html) and the CMP toolkit: coordinated multicellular program detection from donor-level program scores with within-strata permutation FWER control, plus [`cmp_scores()`](../reference/cmp_scores.html), [`ambient_program_scores()`](../reference/ambient_program_scores.html), [`annotate_cmp_lr()`](../reference/annotate_cmp_lr.html), [`annotate_cmp_themes()`](../reference/annotate_cmp_themes.html), [`cmp_gene_overlap()`](../reference/cmp_gene_overlap.html), [`cmp_context_association()`](../reference/cmp_context_association.html), [`annotate_cmp_tf()`](../reference/annotate_cmp_tf.html) (CollecTRI footprints via [`prepare_collectri()`](../reference/prepare_collectri.html)), [`cmp_triage()`](../reference/cmp_triage.html), and [`cmp_nuisance_genes()`](../reference/cmp_nuisance_genes.html).
- [`program_sharing()`](../reference/program_sharing.html): post-hoc shared-vs-specific assessment of matched programs across cell types with ‘mashr’ — permutation-calibrated loading matching, greedy agglomeration into disjoint matched groups, CorShrink-style null-correlation estimation (`z_thresh` configurable), and per-gene lfsr calls. Fits now store loading posterior SDs (`F_psd`/`L_psd`).
- Artifact-gene masking: `artifact_gene_sets` (V(D)J, genotype-marker, sex sets), [`artifact_genes()`](../reference/artifact_genes.html), [`flag_artifact_genes()`](../reference/flag_artifact_genes.html), and fit-time `exclude_genes` in [`fit_celltype_programs()`](../reference/fit_celltype_programs.html) (applied before variable-gene selection).

</div>

<div class="section level3">

### Bug fixes and robustness<a href="#bug-fixes-and-robustness-0-2-0" class="anchor" aria-label="anchor"></a>

- [`program_sharing()`](../reference/program_sharing.html) handles degenerate real-data fits: numerically degenerate programs (noise-floor loadings; fixed factors with collapsed posterior SDs) are excluded from matching, and zero-SD point-laplace posteriors are floored relative to each condition’s loading scale instead of tripping mashr’s zero check.
- [`simulate_pseudobulk()`](../reference/simulate_pseudobulk.html) no longer fails with custom `cell_types` under the “mixed” and “same_genes_independent” scenarios (shared-program objects are created by the first present member instead of a hardcoded cell type).

</div>

<div class="section level3">

### Documentation and infrastructure<a href="#documentation-and-infrastructure-0-2-0" class="anchor" aria-label="anchor"></a>

- Full roxygen2 documentation: `man/` pages for all exported functions with runnable examples; NAMESPACE is now roxygen-managed.
- R-CMD-check CI (ubuntu + macOS) and a pkgdown site skeleton.
- Quickstart vignette and a precomputed real-data case-study article.

</div>

<div class="section level3">

### Benchmarks<a href="#benchmarks-0-2-0" class="anchor" aria-label="anchor"></a>

- `benchmarks/cmp_detection/`: pre-registered v1.1/v1.2 confirmations of the frozen CMP caller (`max_factors = 30` upstream; held-out seeds 1051-1100, all 7 criteria passed) and dev verdicts for four challenger arms (joint-subspace rejected; two-stage identical to winner; Cauchy closest miss; EB shrinkage no-op).

</div>

<div class="section level3">

### Cross-cohort replication case study (CP-5, 2026-10-19)<a href="#cross-cohort-replication-case-study-cp-0-2-0" class="anchor" aria-label="anchor"></a>

- New evaluated article `vignettes/articles/cross_cohort_replication.Rmd`: an end-to-end replication workflow on a three-cohort simulated fixture (conserved outcome-linked program, cohort-specific program, nuisance-only batch effect) orchestrating existing exports only — no new API. Shows the fixed-basis pathway ([`project_programs_to_cells()`](../reference/project_programs_to_cells.html), learned transform applied) versus the de-novo pathway ([`freeze_preprocessing()`](../reference/freeze_preprocessing.html) + independent refit + [`match_programs()`](../reference/match_programs.html)), with predeclared same-sign-as-discovery replication rules, abstention accounting (attempted/projected/matched/ineligible/replicated denominators), three negative controls (200-fold permuted outcome at null ~0.05; scrambled gene IDs abstain loudly; batch-only effect does not replicate as biology), and a leave-one-cohort-out check in which the cohort-specific program correctly abstains. `confirm_programs_across_cohorts()` is deliberately NOT exported: no multi-user need is proven.

</div>

</div>
