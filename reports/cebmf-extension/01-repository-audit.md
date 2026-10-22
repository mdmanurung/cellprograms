# 01 — Repository Audit: cebmf_torch and the cellprograms Seam

Workstream G / Milestone M0. Date: 2026-10-17. Auditor: coordinating agent (G), with three
parallel inspection passes. Full source-grounded reports are preserved verbatim in
`appendix/audit_cebmf_torch.md`, `appendix/audit_flashier_utils.md`,
`appendix/audit_cellprograms_cnmf.md`. This document consolidates their findings into a
capability matrix and flags what is verified vs untested. All file:line references were
verified against the pinned revisions listed below.

Pinned revisions:
- `william-denault/cebmf_torch` — cloned 2026-10-17 (SHA recorded in appendix report header)
- `william-denault/flashier` @ fd44811 (load-bearing pin; CRAN 1.0.7 breaks flashier-utils)
- `william-denault/flashier-utils` @ 271f583
- `mdmanurung/cellprograms` @ e00ca52 (local; 5 commits ahead of origin, awaiting push)
- cNMF v2.x (Kotliar et al. 2019, eLife)

Publication verification: the package implements the method of Denault, Tayeb, Carbonetto,
Willwerscheid, Stephens (2025), "Covariate-moderated Empirical Bayes Matrix Factorization,"
NeurIPS 38, doi:10.52202/085713-1573. The covariate-assisted prior families (cash, lcash,
po_lcash, cgb) are the paper's core contribution. Foundation: Wang & Stephens (2021), JMLR
(doi:10.48550/arxiv.1802.06931) for standard EBMF.

## 1. cebmf_torch capability matrix

| Capability | Status | Evidence |
|---|---|---|
| Variational EBMF objective (negative ELBO, not Frobenius) | VERIFIED | `cebmf/cebmf.py:418-442` |
| Coordinate-wise EM: per-factor L then F update against residual-with-k-added-back | VERIFIED | `cebmf/cebmf.py:306-326` |
| Per-iteration cost O(NP·K²) (3 full residual matmuls per factor) + O(K·(N+P)·G) prior fits | VERIFIED (code reading) | `cebmf/cebmf.py:258-268` |
| Mini-batching affects ONLY ASH mixture-weight EM; factor/residual/tau updates full-batch; statistical objective unchanged | VERIFIED | `utils/mixture.py:7-120`, `ebnm/ash.py:38` |
| README claim "Adam on logits" | INACCURATE (documented limitation) | mixture weights use EM or L-BFGS; Adam only an off-by-default warm start in point-Laplace (`point_laplace.py:73-77`) |
| 11 registered priors incl. covariate-assisted (cash, cgb, cgb_sharp, lcash, po_lcash, emdn, spiked_emdn) | VERIFIED | `priors/registry.py:72-81` |
| "exp" prior name collision: point-exponential overwrites ASH exponential mixture | CONFIRMED DEFECT (P2, upstream) | runtime-verified; see appendix 1 §priors |
| Posterior 2nd moments computed internally (L2/F2) but dropped at Prior boundary; CEBMFResult = {L, F, tau, history_obj} only | CONFIRMED GAP | `cebmf.py:21-39, 502-504`; `priors/base.py:19-33` |
| No LFSR anywhere | CONFIRMED GAP | repo-wide search, appendix 1 |
| NaN-native missing-data mask | VERIFIED | `cebmf.py:499-500` |
| Noise models: CONSTANT / ROW_WISE / COLUMN_WISE / KNOWN(S) | VERIFIED | `cebmf.py:523-608` |
| Stopping: fixed maxit only; no ELBO-based convergence check | CONFIRMED GAP | `cebmf.py:210` |
| Sparse input / chunking / streaming / AnnData: none (dense float32 only) | CONFIRMED LIMITATION | `cebmf.py:185` |
| Determinism (CPU, fixed seed): bit-identical refits | VERIFIED (runtime smoke test, 200×50, K=3, maxit=20, 5.1 s) | appendix 1 §determinism |
| License MIT; deps torch>=2.0 + tqdm only | VERIFIED | `LICENCE.md`, `pyproject.toml` |

## 2. The cellprograms seam (what a Python backend must satisfy)

Exact drop-in point: the `run_flash` closure at `R/ebmf.R:290-293` (`do.call(flashier::flash, args)`
at `:293`); flashier-utils route `:467-495`; fixed-L/F route `:429-455`; fallback ladder `:300-338`
gates on finite ELBO (`.fit_ok :374`); `canonicalize_programs :545-627` consumes `L_pm/F_pm`.

Required-fields contract for a backend fit (from audit 2, §contract):
- Per cell type: `fits[[ct]]$flash` list with `L_pm, F_pm, elbo` (required), `L_psd, F_psd,
  F_lfsr` (optional, tolerated as NULL), `fixed_L/fixed_F` when applicable; plus `$Y_used,
  $genes, $preproc_replay` (center, scale, residualize_coefficients, design_columns),
  `$fallback`, `$skipped`.
- Top level: `cell_types`, `observation_ids`, `sample_metadata` (keyed by `observation_id`),
  `preprocessing`, `provenance`; post-canonicalization scores/loadings with colnames `<ct>_<k>`.
- `flashbridge_metadata` optional for fitting, **mandatory for projection**
  (`R/cell-projection.R:141-142`).

Top 3 translation gaps (audit 2):
1. **Uncertainty**: no lfsr/posterior-sd equivalent in cebmf_torch → `top_genes(by="lfsr")`
   errors (tolerated), CMP loading-shrinkage silently degrades. Mitigation: expose internal
   2nd moments upstream (PR candidate) or ship NULL and document.
2. **Projection oracle**: `project_programs_to_cells` needs a live flashier object + manifest →
   Python fits need a Python-side projection implementation or lose cell-level projection.
3. **ELBO scale**: not comparable across backends; `.fit_ok` finiteness is portable; fallback
   triggers and multistart ranking are not.

## 3. cellprograms / cNMF integration facts (audit 3)

- Pseudobulk construction exists only in `benchmarks/biomni_replication/01_pseudobulk.py`
  (group_sum `:119`, prep_matrix `:149`; "donor||ct" key, ≥50 cells) — needs promotion to a
  reusable module before the extension can depend on it.
- Metacell discovery is net-new repo-wide (MR-3 staged at `reports/STAGE_B_ROADMAP.md:51-58`,
  design-doc-only gate).
- Single-cell projection wrapper exists (`R/cell-projection.R:80` →
  `flashier.utils::project_factor_activity`): dense-only, unbatched, backend-bound.
- Refinement is net-new (MR-4 `:60-70`; `fit_initialized_ebmf()` registered, no implementation;
  drift harness `R/stability.R:31`).
- cNMF consensus internals pinned for fair benchmarking: L2-normalize merged replicate spectra →
  local-density filter (density_threshold=0.5) → KMeans(k, n_init=10, seed=1) → per-cluster
  median spectrum → NNLS usage refit (`cnmf.py:823`); K statistic = mean silhouette of KMeans on
  unfiltered merged spectra + Frobenius prediction error (`cnmf.py:923`, `:1119`).
- **Layering rule (hard constraint)**: cellprograms → flashier-utils → cebmf_torch, strictly
  one-way. No direct cellprograms→cebmf_torch imports. Benchmark arms must not feed different
  input scales to count-based (cNMF) vs log-based (EBMF/flashier) methods.

## 4. Untested claims carried into the design docs

- GPU execution path of cebmf_torch: code-reviewed only; no CUDA available in this sandbox
  (verified: no device nodes, torch CPU-only). All performance claims in this workstream are
  CPU-only.
- Chunked-EM equivalence (the core M2 design claim): derived from the objective's sufficient
  statistics, to be validated numerically (test 5 in the plan), not yet demonstrated.
- Covariate-assisted priors (cash/lcash/po_lcash): API-verified, not yet exercised on real
  data in this session; behavior under donor-level covariate dependence is a design question
  (02-statistical-design.md §B), not a verified property.
