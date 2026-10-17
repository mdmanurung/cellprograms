# 02 — Statistical Design: Three Approaches with Pre-Registered Acceptance Criteria

Workstream G / Milestone M0. Companion to `01-repository-audit.md` (evidence base) and
`03-computational-design.md` (computational realization). Acceptance criteria below are
pre-registered: they are fixed before any M1-M3 implementation is benchmarked and must not be
relaxed post hoc (standing rule of the audit program).

## Common preprocessing contract (all approaches)

Shifted-log normalization, fixed recipe (matches the registered operating point of
workstreams E/F): counts → CP10K → log1p → subset to top-2000 HVGs → center genes, scale by
gene sd. The recipe is recorded in the result contract (`preprocessing` block) so any backend
can replay it exactly. Rationale: cebmf_torch and flashier are both Gaussian-likelihood
matrix factorizations; feeding either raw counts or mismatched scales would bias any
comparison against cNMF (count model) or EBMF (log model). Benchmark arms therefore share one
preprocessed matrix per granularity.

## Approach A — Standard EBMF (cEBMF core, covariate-free)

**Model.** Y (N obs × P genes) ≈ L F^T + E, E_ij ~ N(0, tau_i or tau_j), with point-Laplace
sparsity priors on L and F columns (cebmf_torch `laplace/laplace`), COLUMN_WISE noise
(per-gene variance), K=10 with backfitting and pruning (prune_thresh = 1-1e-3). This mirrors
the flashier reference fit (point-Laplace priors, var_type=1) used by cellprograms at the
registered operating point (max_factors=10, seed=20261016).

**Non-identifiability handling (binding).** EBMF is invariant to factor rotation, sign, and
column scaling. Exact numerical equality of L/F across backends is NOT expected and is NOT an
acceptance criterion. All cross-backend comparisons use matched quantities: Hungarian matching
on |Pearson correlation| of factor scores, top-gene set Jaccard, and held-out prediction
error. Within-backend determinism (bit-identical refits) IS an acceptance criterion because
cebmf_torch is deterministic given fixed init.

**Known objective difference (documented, not a defect).** flashier places priors on
standardized loadings with per-row/per-column variance options and estimates uncertainty
(lfsr, psd); cebmf_torch's public result carries point estimates only. Downstream consumers
that require lfsr degrade gracefully to NULL (verified in audit 2). This is a capability
difference, recorded in the comparison, not a defect to fix in M1-M3.

**Pre-registered acceptance criteria (A).**
- A1: On clean synthetic data (no dropout, K_true ≤ 5 planted programs), cebmf_torch recovers
  all planted programs with Hungarian-matched |cor| ≥ 0.9 between true and estimated scores.
- A2: On the same data, flashier and cebmf_torch agree: matched |cor| ≥ 0.8 on shared
  programs, and top-50 gene Jaccard ≥ 0.6 per matched program.
- A3: Refit determinism: bit-identical L/F across two runs with identical inputs and seed.
- A4: ELBO finite and non-decreasing (within 1e-6 relative tolerance per iteration) on all
  benchmark arms; the adapter records the full objective trace.
- A5: Effective K detected by pruning within ±1 of the planted K on synthetic data with 2
  null factors added.

## Approach B — Covariate-moderated EBMF (the cEBMF paper's contribution)

**Model.** Prior parameters on L (or F) depend on covariates X via the cash/lcash/po_lcash
families: e.g., program usage may scale with donor age, or program loadings may be
moderated by stimulation status. This is prior moderation (the covariate changes the
sparsity/scale of what the factor CAN express), not post-hoc regression on scores.

**Leakage rules (binding).**
- B-L1: Covariates must be external to the factorized matrix (donor age, sex, condition,
  batch) or derived by a leakage-safe pipeline (e.g., reference-trained cell-type scores).
  No covariate may be computed from the same matrix being factorized (e.g., library size,
  program scores themselves).
- B-L2: Covariates enter at the observation level of the factorized matrix. For pseudobulk
  (donor × cell type) matrices the observation is the donor-pseudobulk sample; cell-level
  covariates must never be treated as independent replicates of donor-level covariates
  (donor is the unit of dependence; consistent with the CMP design and mcRigor's donor-blind
  caveat noted in workstream E).
- B-L3: Evaluation of moderation claims uses held-out donors: fit with covariates on a donor
  training split, evaluate moderation (prior-parameter vs covariate association) on held-out
  donors. No moderation claim is accepted from in-sample fit statistics alone.

**Pre-registered acceptance criteria (B).**
- B1: On synthetic data with a planted covariate-moderated program (usage amplitude a
  function of a continuous covariate), the cash/lcash fit recovers the moderation (Spearman
  of estimated usage vs covariate ≥ 0.8 on held-out donors) while a covariate-free fit
  (approach A) does not systematically (|Spearman| ≤ 0.5 on the same program).
- B2: Moderation does not destroy recovery: planted non-moderated programs still meet A1
  thresholds under the covariate-moderated fit.
- B3: Any claim that cEBMF moderation outperforms "A + post-hoc regression of scores on
  covariates" requires the held-out comparison B1/B3; if post-hoc regression matches the
  moderation result, the report must say so (the paper's value proposition is tested, not
  assumed).

## Approach C — Multiresolution discovery (pseudobulk/metacell init → single-cell refine)

**Stage specification (per-stage objective, prior, normalization, convergence).**

| Stage | Data | Objective | Prior | Convergence |
|---|---|---|---|---|
| C1 discovery | pseudobulk (donor×ct) or metacell matrix | full cEBMF EM (ELBO) | laplace/laplace, K=10, backfitting | ELBO plateau: relative improvement < 1e-5 for 3 consecutive checks, maxit cap 100 |
| C2 projection | single-cell counts (sparse, chunked) | L-update with F fixed from C1 (exact per-chunk; see 03 §3) | prior from C1 (fixed) | single pass (or ≤3 passes); tau from accumulated sufficient statistics |
| C3 refinement | single-cell matrix | full cEBMF EM initialized from C1/C2 (user-provided factors) | laplace/laplace re-fitted | ELBO plateau as C1, maxit cap 50 |
| C4 rare-state scan | single-cell residuals | residual analysis + factor addition (backfitting on) | as C3 | as C3 |

**Approximation map (what is exact vs approximate).**
- C2 projection with F and tau fixed is EXACT for the L-update: the Gaussian likelihood
  factorizes over rows, so chunked inference equals full-batch inference (validated by test 5
  in the plan, < 1e-4 relative agreement).
- C2 with tau re-estimated from accumulated per-column sufficient statistics (masked residual
  sums) is also exact; tau updated once per full pass, not per chunk.
- C3 refinement is a full posterior fit (not an approximation) but is initialization-dependent;
  the multiresolution claim is therefore about cost and rare-state recovery, not about
  reaching a better optimum than a from-scratch single-cell fit. The comparison is explicit:
  C3-from-init vs A-from-scratch at matched maxit.
- Known risk from workstream E: ebmf K-collapse under high dropout (ZeroRate). C1 on
  pseudobulk/metacell data (lower ZeroRate than single-cell) is the designed mitigation;
  the ZeroRate budget of the input granularity is recorded per arm (mcRigor grounding,
  Liu & Li 2025, doi:10.1038/s41467-025-63626-5).

**Pre-registered acceptance criteria (C).**
- C1: C3-from-init reaches ELBO within 1% of A-from-scratch ELBO on the same single-cell
  matrix at matched maxit, at ≥ 5× lower wall-clock (the cost claim).
- C2: Rare-state recovery: on synthetic data with a rare program (≤ 5% of cells), the
  C1→C2→C3 pipeline recovers it (matched |cor| ≥ 0.7) in ≥ 8/10 seeds; from-scratch A fit is
  the reference arm (its own recovery rate reported either way).
- C3: Program concordance across granularities (pseudobulk vs metacell vs single-cell
  pipelines) meets or exceeds the workstream F pilot concordance of the flashier backend
  (matched |cor| and top-50 Jaccard, same metrics, same data).
- C4: No granularity arm may claim "discovery" without reporting its ZeroRate and effective K
  after pruning (dropout-collapse transparency, from the workstream E finding).

## Relationship to the registered Stage B gates

MR-3 (metacell design doc) and MR-4 (anchored refinement + ablation) overlap with approach C.
This workstream does NOT supersede those gates: the C-stage specs above are the design
content MR-3 requires, and the C3 refinement arm is the fit_initialized_ebmf() realization
MR-4 registered. Reconciliation: when Stage B opens, MR-3/MR-4 adopt or cite these specs
rather than re-deriving them; until then, workstream G results are pilot evidence, labeled as
such, distinct from the registered MR-2 cNMF baseline benchmark.
