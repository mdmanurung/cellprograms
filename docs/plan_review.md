# CellPrograms — Plan Review

Review of `cellprograms_package_scope_implementation_plan.md` and
`cellprograms_benchmarking_plan.md`, prioritized by risk. Overall verdict: both
plans are unusually rigorous and buildable as written; the findings below are
corrections and clarifications, not structural objections.

## A. Statistical / method soundness (fix before coding)

1. **Stability score is ill-defined under bootstrap resampling (Section 13 of the
   implementation plan) — correctness bug.** Bootstrap refits produce score
   matrices on *resampled* rows, so `cor(Z_k, Z_l)` between reference and
   bootstrap factors is undefined across different row sets. Spec should state:
   match factors primarily by gene-loading cosine (always comparable), and
   compute score correlation only on rows shared between reference and resample
   (with duplicate weights), or drop the score term when row sets differ.
2. **flashier API facts (verified in this session).** `flashier::flash()` takes
   no `seed` argument — use `set.seed()` before the call. Fitted factors are
   `fit$L_pm` / `fit$F_pm` (posterior means), not `fit$loadings[[...]]`. The
   plan's `ebnm_fn = list(ebnm_normal, ebnm_point_normal)` correctly maps
   Normal→sample side, point-Normal→gene side for an observations × genes
   matrix.
3. **Canonicalization (Section 11) is correct** — RMS-scale + sign preserves
   reconstruction exactly (verified to 1.8e-15). Add a deterministic tie-break
   for the largest-|loading| gene (magnitude → gene name).
4. **Missingness vs. benchmark metrics.** The package correctly refuses to
   zero-fill missing program scores, but the benchmark's KNN-based
   biology-retention metrics need a defined available-case rule (e.g.,
   masked-distance KNN over observed program scores). Specify it, or metric
   implementations will each invent their own.
5. **Optimized defaults (this session's evidence).** Per-gene variance
   (`var_type = 1`) with a point-Laplace gene prior beat the plan's Section 9
   defaults consistently across scenarios, held-out data, and sample sizes,
   with the largest edge at small N (+0.0064 at N=50). The plan's `var_type = 2`
   should be the fallback, not the default.

## B. Benchmark design gaps

1. **Fairness regime A is not literally matched for distributional methods.**
   GloScope/PILOT consume cell-level latent representations, not pseudobulk
   matrices. State explicitly that matched-input applies to pseudobulk methods;
   distributional methods get the same sample filters and a standard cell-level
   PCA embedding.
2. **Comparator availability risk.** mc-ASTRA/MOFA-FLEX and EV-BIDIFAC may not
   be publicly/legally usable. Define the fallback now: MOFA2 with explicit
   grouped views as the joint-factorization comparator; keep EV-BIDIFAC
   optional (as the plan already does).
3. **No uncertainty on benchmark conclusions.** Add bootstrap CIs (over
   samples) or per-dataset sign tests for the headline local-EBMF vs
   mc-ASTRA comparisons, so the Section 20 regime map is statistically
   supported.
4. **Strengths worth keeping:** separate biology/technical retention panels,
   no single aggregate score, the 7-scenario simulation suite, and the explicit
   ablation ladder (Section 16) are all sound and above field standard.

## C. Engineering feasibility

1. **flashier is GitHub-only** (stephenslab/flashier); ebnm is on CRAN. Pin
   versions in provenance (this session: flashier 1.1.42, ebnm 1.0.59, R 4.4.2).
2. **Zero-factor edge case (found by this session's noise control).** On pure
   noise, flashier retains 0 factors and `L_pm` is NULL — all extraction code
   must guard K = 0. This is desirable behavior (the pipeline does not
   hallucinate structure) but will crash naive implementations.
3. API, object structure, and the Phase 0–10 roadmap are buildable and
   correctly sequenced; the MVP (Section 28) is well chosen. Minor: document
   the `gene_variance` flag in `fit_celltype_programs()`.

## D. What was actually done this session

- Phase 0 package skeleton (`cellprograms` 0.1.0): data model, EBMF backend,
  canonicalization, simulation module; 17/17 unit tests pass.
- Benchmark-first optimization loop (tusoskill, 50 substantive iterations):
  evaluator frozen at eval-v3 after two bug fixes (gene-universe alignment;
  R-1-based vs Python-0-based indexing), then a 15+13+10-candidate config
  search with guardrails (noise control, seed sensitivity, held-out scenario,
  N=50/300, missing cell type).
- Locked default config: `center=TRUE, features=all, prior=point_laplace,
  var_type=1, greedy_Kmax=30, backfit=TRUE, nullcheck=TRUE, no PVE threshold,
  no merging`. Frozen-benchmark aggregate 0.9212; held-out 0.9468.
