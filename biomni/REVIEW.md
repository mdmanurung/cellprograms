# Review of Biomni deliverables (2026-10-07)

Scope: `Biomni_lab_downloads_20261007_115149.zip` (PMD integration on COMBAT) and
`..._115157.zip` (4-cohort benchmark, Biomni's `cellprograms/` fork, trace notebooks).
Read-only review of the zip contents; no repo code was changed.

## Verdict

The work is careful in places: the PMD permutation-null bug was found and fixed, and
the final COMBAT run (exec 21) postdates the fix (exec 20). But the benchmark headlines
are **not yet trustworthy** (dimension-biased metric, zero-filled missingness, superseded
defaults), and **4 of the 5 PMD "multicellular programs" look like donor-level
confounders**, not coordination.

## Benchmark findings

1. **`max |Spearman|` favours high-dimensional representations.**
   Null simulation, n = 122, mean max |rho| over d random columns:
   d = 10 → 0.167, d = 30 → 0.209, d = 266 → 0.274 (p95 0.25 / 0.30 / 0.32).
   The random floor was built at 10 dims (`rep_random(donors, 10L)`; observed 0.200).
   COMBAT Age: cp_L1groups 0.464 vs perct_pca 0.422 (30 dims) is a tie or reversal
   against dimension-matched floors (~0.27 vs ~0.21). Same issue for HLCA Age.
   Stephenson Age (0.653 vs <=0.48) likely survives. The "cp wins are conservative"
   caveat holds only for kNN macro-F1 (sex), not for max-Spearman.
2. **Missing donors are zero-filled** (`align_donors`). Only 9/124 COMBAT donors have all
   10 cell types, so kNN neighbours partly encode which cell types a donor lacks (tracks
   abundance, hence disease; PB and DC). Larger effect for cp (more blocks). Contradicts
   `docs/plan_review.md` A.4 (available-case rule). Applies to every Panel A/B number.
3. **Superseded defaults.** Benchmark used K = 30 (COMBAT) / 25 (rollout) with the fork's
   `point_normal`, `var_type = 2` defaults. Repo commit `0c3baef` replaced these (K = 10,
   `point_laplace`, `var_type = 1`); the repo's own result is Stephenson disease kNN
   0.468 -> 0.648. "PCA wins Source/Status" may be a configuration artifact. The PMD
   run used K = 10, so the two Biomni deliverables are mutually inconsistent.
4. **Panel B reading needs care.** "cp sits at the technical floor" is only informative
   alongside biology. Spot-check of Stephenson: cp_block Status 0.833 with Site 0.558
   (floor 0.444) is a genuinely good result; cp_L0 Status 0.460 (floor 0.265) with Site
   0.684 is weak. COMBAT cp Source ~0.28 (floor 0.099) at Institute floor is mixed.
   Per-representation, not blanket, conclusions.
5. **Provenance gap.** The COMBAT trace stops after baselines/comparators (exec 67);
   cp evaluation, the other 3 cohorts, stability, cross-cohort matching and
   `benchmarks/R/metrics.R` are not in the zip. Most numbers cannot be audited.

## PMD integration findings

6. **MCP biology (COMBAT, K = 10).** All five have perm p = 1/101 (floor), all 10 cell types
   active, view contributions ~8 of max 9 (mean pairwise r ~0.9).
   - MCP01: Ig V-genes lead in CD4/CD8/NK; ambient/plasmablast contamination likely.
   - MCP02: sex (p = 5e-18); donor-level by construction.
   - MCP03: interferon; the plausible real MCP (Source p = 1.5e-6, HV lowest).
   - MCP04: JUN/ATF3/NR4A2/DUSP; Institute p = 2e-4 (St George's +1.0); processing batch.
   - MCP05: GSTM1/CHCHD10/HOXB2 lead; possible germline axis. The state histogram is
     unimodal, so this is unconfirmed.
   The value-permutation null tests "any shared donor-level variable", not coordination.
   **Fix:** residualize per-CT scores on sex + Institute (+ optional Ig score) before PMD.
7. Smaller: p-value saturates at 1/101 so rho tuning degenerates to max objective (use
   a null z-score); `n_active_celltypes` is always "all" under SUMCOR; p ignores rho
   selection. Positives: simulation validation (S1/S2/S4/S7, 30% missing) is solid, and
   the contribution-weighted state fixes private leakage (S4: 0.82 -> 0.95).
8. **Not drop-in:** depends on fork-only `integrate_programs`, `program_scores(.., "list")`,
   `trace_state`, `.stopf`/`.warnf`, `%||%`.

## Minor

- 6 of the 13 "strong" COMBAT<->Stephenson matches are sex programs; no null for the
  |cos| threshold; one subtype per cell type (e.g. CD4 -> CD4.CM only). Counts verified:
  28 pairs, 13 with |cos| >= 0.8.
- OneK1K sex 0.978 is a Y-gene positive control.
- L2 "negative result" (PB R2 0.57 -> 0.13, ELBO -256.7k -> -264.1k) looks like an
  implementation defect: a borrowing step that lowers ELBO should fall back to L0.
- Held up: L1-global ~ L1-groups, scores-space principal-angle degeneracy
  (k_a + k_b > n), L1-global runtime blow-up.

## Port decision: fork -> repo

Repo fits already expose `fit$scores` / `fit$loadings` after `canonicalize_programs()`
(`R/ebmf.R:237-256`). Gaps: program IDs (`ct_k` vs `ct__F001`), `fit$programs`,
`fit$input$observation_ids` (~20-line shim). Do **not** replace the repo's `fit_celltype_programs`
with the fork's: the repo's is ahead (tuned defaults, point_laplace fallback, K=0 handling).

| Fork piece | Decision | Why |
|---|---|---|
| `principal-angles.R` (`match_programs`, `principal_angles`, loadings-space `sharing_spectrum`) | Port now | Needed for stability and cross-cohort matching |
| `stability.R` (`assess_program_stability`) | Port now | Promised in repo DESCRIPTION; COMBAT 73% stable core |
| `extract.R` (`program_loadings`, `top_genes`, `program_summary`) | Port now | Thin accessors the above use |
| `block-joint.R` (`fit_joint_block`) | After K = 10 re-run | Clean sims, but validated only at K = 25-30 |
| `integrate-pmd.R` + `integrate.R` scaffolding | After residualization fix | Solver sound; 4/5 real MCPs confounded |
| `simulate.R` (S1-S7) | With block/PMD | Test harness for integration |
| L1 shared priors | Don't port | No rep-level effect; 3-18x runtime |
| L2 init borrowing | Don't port | Collapses rare CTs; likely bug |
| `evbidifac.R` | Don't port | Batch-contaminated, slow, convergence-sensitive |
| ICA / PVA | Don't port | No evidence shown |

## Suggested next step

Phase 1 port: principal-angles + stability + extract + shim, with tests (~1-2 h).
R is not on PATH in this shell, so tests need an R environment.
