# 03 — Computational Design: Scaling cEBMF to Large Cohorts (CPU-First)

Workstream G / Milestone M0. Companion to `01-repository-audit.md` and
`02-statistical-design.md`. All performance statements are CPU-only (no CUDA in the sandbox;
verified). The GPU path of cebmf_torch remains code-reviewed-but-unexecuted.

## 1. Baseline cost model (from the audit, verified against code)

Per EM iteration of the full fit (N obs × P genes × K factors):
- 3 full residual matmuls per factor (coordinate-wise L then F update against
  residual-with-k-added-back): O(N·P·K) per factor pass, O(N·P·K²) per iteration
  (`cebmf/cebmf.py:306-326, 258-268`).
- Prior fits: O(K·(N+P)·G) with G mixture components (ASH/Laplace EM per factor column).
- ELBO: O(N·P·K).

Memory (dense float32, per matrix entry): Y (4 B) + residual R (4) + mask (1) + L2/F2 (4) +
L/F (8) ≈ 21 B/entry persistent; temporaries during updates push peak to ~45-55 B/entry
(audit 1 §memory). At 5M cells × 3000 genes (1.5e10 entries): ~315 GB persistent, ~700+ GB
peak — impossible on a 64 GB machine, and the motivation for the chunked design.

## 2. Key structural fact: the EM updates need only sufficient statistics

The likelihood factorizes over rows given (F, tau) and over columns given (L, tau). All
quantities the coordinate updates consume are accumulable across row-chunks:

- L-update given F: per-row weighted regression against F with per-gene precision tau_j.
  Sufficient per-chunk statistics: Y_chunk @ (F * tau) (n×K) and (F * tau) @ F^T (K×K,
  shared across chunks).
- F-update given L: per-gene regression against L. Sufficient per-chunk statistics:
  Y_chunk^T @ (L * tau_i) (P×K) and (L * tau_i) @ L^T (K×K).
- tau (COLUMN_WISE): masked mean squared residual per gene — sum of squared masked residuals
  and masked counts per column, both accumulable per chunk.
- ELBO: loglik term from the same masked residual sums; KL terms depend only on the
  (small) prior parameters and L/F.

Consequence: an EXACT chunked EM is realizable with working memory O(chunk·P + K·P + K²)
plus the sparse count matrix, instead of O(N·P). The per-factor coordinate structure
(residual-with-k-added-back) is preserved by reconstructing R_chunk = Y_chunk - L F^T +
l_k f_k^T per chunk — no full-size residual is ever materialized. This is the core design
claim of M2 and is validated numerically (test 5: chunked vs dense L agreement < 1e-4
relative on 10k cells; validation runs in float64 to separate chunking error from float32
roundoff).

What is NOT chunked in M2: the prior-fit step (EBNM per factor column) operates on the
K·(N+P) parameter vectors — small, full-batch, unchanged. Mini-batching inside ASH mixture
EM (upstream) is untouched: it does not change the objective (audit 1).

## 3. Equivalence-vs-approximation map (binding for claims)

| Operation | Status | Condition |
|---|---|---|
| Chunked L-update, F and tau fixed | EXACT in the likelihood; APPROXIMATE in the prior | rows are independent given (F, tau), but the package refits the L-prior per call on the current statistics (cebmf.py:347-360), so chunked priors see chunk-level values. MEASURED (M2, synthetic 1000x60): signal-factor shape error vs dense is converged (not a convergence artifact) and scales with chunk count: 3.0-10.9% at 10 chunks, 1.9-6.4% at 4, 1.3-4.7% at 2, 0% at 1 chunk. At 10k-row benchmark chunks the effect is smaller. Immaterial for usage scoring; material for exact posterior reproduction. Upstream feature request: posterior-only updates with frozen prior state. |
| Chunked tau update from accumulated masked residual sums | EXACT | sufficient statistics |
| Chunked F-update from accumulated statistics | EXACT | same argument over columns |
| Full chunked EM (L, F, tau all updated chunk-wise) | EXACT in likelihood; prior-refit caveat as above | per-chunk residual reconstruction preserves coordinate structure |
| Projection with F fixed, few passes | APPROXIMATE as inference | it conditions on F (ignores F uncertainty); exact for the conditional L-update only |
| Refinement from init | FULL fit | initialization-dependent; compared against from-scratch at matched maxit (02 §C) |

## 4. Adapter architecture (M1-M2 code layout, cellprograms repo `python/`)

- `cebmf_ext/io.py` — sparse loaders (mtx/h5ad → CSR), gene/obs id handling, no dense
  materialization at load.
- `cebmf_ext/preproc.py` — replayable shifted-log recipe (CP10K → log1p → HVG subset →
  center/scale), serialized into the result contract; chunk-wise replay for projection.
- `cebmf_ext/fit.py` — thin wrapper over `cEBMF` (pinned API: `cEBMF(data, K, prior_L,
  prior_F, internal_epoch=10, allow_backfitting=True, prune_thresh=1-1e-3, noise_type,
  S=None, device=None)`; `fit(maxit)` → `CEBMFResult(L, F, tau, history_obj)`); adds
  ELBO-plateau stopping (adapter-level; the objective is untouched), effective-K detection
  from pruned factors, runtime/peak-memory capture (psutil), and the result contract.
- `cebmf_ext/project.py` — chunked projection (M2): fixed F, sparse CSR input, per-chunk
  preprocessing replay, tau from accumulated statistics; returns full-N L without ever
  holding dense N×P.
- `cebmf_ext/contract.py` — result contract schema (factor ids, gene/cell ids, prior info,
  config, objective trace, runtime, memory, backend versions, preprocessing provenance;
  unavailable fields serialized as null — no manufactured uncertainties).
- Tests: `tests/test_determinism.py`, `test_nan_masking.py`, `test_pruning.py`,
  `test_elbo_monotone.py`, `test_chunk_equivalence.py`, `test_contract_schema.py`.

Boundary rule: this package imports cebmf_torch, never the reverse; cellprograms (R) will
consume the adapter's outputs through the documented contract (audit 2 §2), not through
direct imports.

## 5. Benchmark plan (Phase F scope, after the M3 gate)

- Small 10k × 3000: full dense fit and chunked path; correctness cross-check; ~10 min.
- Medium 100k × 3000: chunked path; memory profile; ~1-2 h.
- Large 1M × 3000: chunked path, background with checkpoints to shared workspace every
  completed stage; ~8-24 h estimated from measured Medium throughput.
- 5M × 3000: NOT run. Design-validated on a 100k shard; reported as measured-throughput
  extrapolation with the memory model above (the 5M dense baseline is included in the
  extrapolation table as infeasible-on-64GB, which is itself a result).
- All arms: fixed seeds, pinned revisions, recorded versions; repeated-run variability
  reported for timing claims (≥ 3 repeats on Small/Medium).

## 6. Checkpointing and failure policy

Long runs (Medium/Large) checkpoint L, F, tau, ELBO trace, config JSON to the shared
workspace after each completed stage; a resumed run verifies checkpoint integrity (shapes +
config hash) before continuing. Any checkpoint claimed as preserved is verified by
enumeration and size check before being reported.
