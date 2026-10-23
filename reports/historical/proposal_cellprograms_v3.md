# cellprograms v3 — roadmap and proposal

**Date:** 2026-10-12
**Status:** proposal for user decision; nothing here is preregistered or authorized for implementation yet
**Baseline:** `feat/coordinated-programs` (PR #3 open): cap-30 CMP detection, frozen `n_agg_noanc` caller, TF footprints + triage; 46 tests passing. cmp-bench-v1.2 held-out: pooled F1 0.705 (param 0.475, semi 0.934), all 7 prereg criteria passed.

---

## 0. Where v1.2 left us

**Validated.** CMP caller with exact FWER control (pooled null FWER 0.033 param / 0.000 semi); near-ceiling caller performance on realistic data (semi oracle ~1.0, caller F1 0.93–0.99 at every breadth); cap-30 fits double real-data CMP yield (13→27 Stephenson, 11→25 COMBAT) and improve cross-cohort matching (4→7).

**Known gaps, ranked by scientific cost:**

| # | Gap | Evidence |
|---|---|---|
| G1 | Narrow CMPs (2–3 of 12 cell types) missed ~71% on weak-signal data | param R20 miss 0.712 at both caps; oracle ceiling 0.57–0.59 — the bottleneck is **upstream fitting**, not the caller |
| G2 | ~Half of new cap-30 real-data CMPs carry V(D)J/genotype-marker genes | 7–9 of 13 new per cohort; unvalidated heuristic; genotype origin is a hypothesis, not a finding |
| G3 | No principled shared-vs-cell-type-specific program assessment | current `cmp_gene_overlap()` is a pairwise top-100-gene heuristic; the coupled-EBMF gene-sharing plan is blocked on GS-00 (ELBO derivation, structural criterion) |
| G4 | Downstream association is ad hoc | severity p-values computed per-CMP without mixed models, composition covariates, or multiple-testing discipline |
| G5 | Cap tension | K ≥ 30 helps CMP detection but collapses disease-signal recovery in the concatenated representation (COMBAT model-selection benchmark); cap-30 recommendation is scoped to CMP detection only |
| G6 | Fitting costs ~20 min/dataset; no reuse across cohorts | every new cohort refits from scratch |
| G7 | Ambient decoys still called in 4%/2% of ambient datasets | profile-correlation filter at 0.3 |

**Rejected, do not revisit without new evidence:** MOFA-style joint fitting (F1 0.047, broken guardrails), C1/C2 caller add-ons (≤ +0.004), all-genes/4000-HVG fitting (−0.21), point-normal prior (−0.091), L1/L2 fit borrowing, EV-BIDIFAC, ICA/PVA, replacing `fit_celltype_programs`.

**In flight (this session):** `subspace_sum` joint-subspace arm on dev (360 datasets); param nullx fit regeneration (seeds 31–130) for pooled null calibration.

---

## 1. v3 objectives and non-goals

**Objectives:**
1. Close the narrow-CMP power gap (G1) or prove it is irreducible at realistic donor counts.
2. Make the real-data CMP catalog credible (G2): genotype/VDJ artifacts characterized, validated, or filtered.
3. Ship a principled shared-vs-specific program module (G3) that also unblocks the coupled-EBMF workstream's required comparator.
4. Ship a defensible CMP–phenotype association module (G4).
5. Cut the cost of applying the pipeline to a new cohort from ~20 min × donors to seconds (G6).

**Non-goals for v3:** joint latent-factor fitting as a *detection* method; changing the frozen v1.2 caller without a new preregistered benchmark; single-cell-level (non-pseudobulk) coordination; the coupled-EBMF optimizer itself (GS-00 derivation is a separate workstream — v3 supplies its comparator, not its optimizer).

---

## 2. Workstreams

### WS0 — Close out v1.2 (precondition for everything)

| Item | Action | Exit criterion |
|---|---|---|
| subspace_sum arm | Finish dev → gate selection → heldout2 + nullx scoring vs `n_agg_noanc` | Adopt/reject recorded in `tradeoffs.csv` with paired CI; adopt only if held-out F1 improves **and** all guardrails ≤ 0.10 |
| Param nullx fits | Regeneration running (8/100 done) | Pooled param null FWER reportable for any new arm |
| PR #3 | Merge after user review | Branch on `main` |
| `strata` verdict (M2 leftover) | User decision: rejected-as-registered vs adopted-with-deviation | Recorded in `benchmarks/model_selection/DEVIATIONS.md` |
| `feat/sofa-comparator` | PR or archive | No unpushed branches |

### WS1 — Caller statistics (cmp-bench v1.3)

**Rationale.** The caller is at ceiling on semi but the param R20 regime may admit a better statistic. Two untried ideas remain from the prereg review, plus one new one:

1. **Two-stage pair-gate → node test** (review recommendation, not yet tried): screen edges at liberal threshold, test nodes only within surviving neighborhoods — fewer tests, more power for narrow CMPs.
2. **Cauchy/omnibus combination** of per-partner p-values at each node, replacing the hard tau threshold: powerful when 2–3 edges sit just below threshold. Same permutation null, exact FWER preserved.
3. **subspace_sum** verdict from WS0 feeds here.

**Gate (cmp-bench v1.3 prereg, same skeleton as v1.2):** dev seeds 1–30 + nullx for selection; heldout2 seeds 1051–1100 for confirmation; adopt only if paired F1 CI excludes 0 against `n_agg_noanc` **and** max guardrail ≤ 0.10 **and** runtime ≤ 3×. Explicitly powered for the R20 stratum (report per-breadth F1, not just pooled).

**Cost:** ~1 day of HPC per arm (fits are cached; caller evals only).
**Risk:** low — worst case is a clean rejection recorded in `tradeoffs.csv`.

### WS2 — Upstream power for narrow CMPs (the real G1 bottleneck)

**Rationale.** Param oracle is 0.57–0.59: the member program is never fitted in ~40% of narrow-CMP datasets. No caller beats that. Three candidate fixes, in increasing order of ambition:

1. **Score-uncertainty propagation.** flashier already computes posterior second moments; the caller currently treats scores as point estimates. Down-weight high-variance donor scores in the edge/node statistics. *Cheap; no refitting; testable as a caller arm on existing fits.*
2. **Adaptive per-cell-type caps.** Replace global `max_factors` with a per-cell-type K chosen by stability-penalized PVE elbow. May simultaneously help narrow-CMP recovery (more factors where identifiable) and the cap tension G5 (fewer factors where not). **Dual gate required:** must clear cmp-bench v1.3 *and* not degrade the model-selection benchmark's disease-recovery endpoint beyond a preregistered tolerance (this is the lesson of the cap-30 deviation: fit changes need both benchmarks' verdicts, scoped explicitly).
3. **Guided initialization.** Initialize factors from Hallmark/Reactome gene sets, then let EBMF rotate freely. Not the rejected joint fit — no coordination assumed, only within-cell-type gene-set priors. *Highest risk; only if (1) and (2) leave the gap.*

**Gate:** each variant is an upstream arm in cmp-bench v1.3; adoption requires the same criteria as v1.2 (7 criteria) plus the dual-gate rule for fit changes.

**Cost:** (1) ~0.5 day; (2) refits needed — ~2 days HPC for dev+nullx+heldout2 at K-adaptive; (3) ~1 week including a new prereg.

### WS3 — Shared-vs-specific program module (mashr-based, post-hoc)

**Rationale.** G3. The coupled-EBMF plan is blocked on GS-00; meanwhile the package's own doctrine — *fit independently, assess sharing subsequently* — supports a post-hoc answer now. mashr (Urbut et al.) is built for "which conditions share this effect": units = genes, conditions = cell types, effects = sign-aligned loadings, SEs = loading posterior SDs.

**Design:**
- `program_sharing()`: candidate program matching across cell types via loading correlation / principal angles (exists) with permutation-calibrated match statistic → mashr on aligned loadings → per-gene lfsr + data-driven sharing patterns ("all-shared", "lymphoid-only", …) with mixture weights.
- **V estimation** (null correlation of loading z-scores across cell types, induced by shared donors/ambient): CorShrink-style estimator implemented in-package (~50 lines; CorShrink itself is GitHub-only and lightly maintained — do not take the dependency).
- **Storage change:** retain `F_psd`/`L_psd` in `fit$fits` (one line in `ebmf.R`; flashier computes them already).
- **Coordinate-system contract** from the gene-sharing plan §3.1 applies: loadings mapped to a common gene coordinate system before comparison; view-specific scaling is not absorbable.

**Hard caveats to encode in the API:** mashr does not solve alignment (rotation ambiguity) — false matches become false sharing claims, so the match calibration is part of the shipped function, not a user responsibility. Shared loadings ≠ coordinated donor activity; the two axes stay separate in outputs.

**Gate:** new gene-sharing benchmark (the planned `benchmarks/gene_sharing/`): planted fully-shared / partially-shared / private programs; metrics = false-sharing rate, false-membership rate, private-program loss. This module also enters that benchmark as the **strong independent-EBMF comparator** the coupled-EBMF plan's primary hypothesis is measured against — implementing it unblocks GS-00-adjacent evaluation without deriving the coupled ELBO.

**Cost:** ~2–3 days implementation + benchmark prereg; benchmark run ~1 day HPC.

### WS4 — Real-data artifact validity (G2)

1. **`genotype_aware` fitting mode:** option to exclude IG/TR variable genes and a curated CNV/deletion/HLA-haplotype marker list at fit time (or mask from CMP calling), plus regression of donor ancestry PCs when available. Ship the marker list as documented data with provenance.
2. **Validation study:** test whether genotype-marked CMP donor scores track tagging SNPs/eQTLs in a cohort with genotypes (e.g., OneK1K). Converts the biggest real-data caveat into either a finding or a cleanup rule. *Requires external data access; can start as an eQTL lookup (public summary stats) before any new data download.*
3. **SoupX-style per-donor ambient estimates** replacing the profile-correlation heuristic (G7); gate: ambient decoy rate ≤ 0.02 on cmp-bench without F1 loss.

**Gate:** (1) and (3) are cmp-bench v1.3 arms (ambient scenario is the relevant guardrail); (2) is an analysis study with its own pre-analysis plan (outcomes: fraction of genotype-marked CMPs with a significant cis-eQTL/tagging-SNP association vs matched null CMPs).

### WS5 — CMP–phenotype association module (G4)

- `cmp_associate()`: per-CMP donor-score association with a phenotype via mixed models (site/batch random effects, cell-type-composition covariates), BH across CMPs.
- `cmp_mash()`: per-cell-type or per-cohort effect estimates + SEs → mashr → posterior effects, lfsr, sharing-by-magnitude/sign. Replaces the ad-hoc severity p-values and formalizes cross-cohort conservation (Stephenson × COMBAT as two conditions) in place of the |r| ≥ 0.5 heuristic.
- Error notions stay complementary: FWER for calling CMPs (caller), lfsr for sign/sharing claims (associations).

**Gate:** simulation check on null phenotypes (calibration of p-values/lfsr under site confounding and composition shifts); then reanalysis of the Stephenson/COMBAT severity associations as the worked example.

**Cost:** ~2 days + validation runs.

### WS6 — Engineering and adoption (G6)

1. **Project-and-call mode:** fit programs once on a reference atlas; project new cohorts (loadings fixed, scores estimated) → caller. Seconds per new cohort. Gate: projected-score CMP calls on a held-out cohort match from-scratch calls (F1 ≥ 0.9 against the from-scratch labels on that cohort).
2. **CMP atlas:** run the pipeline across public PBMC cohorts (CELLxGENE), cluster recurring CMPs into a named reference taxonomy, ship as a data package. Depends on WS4 (artifact filtering) and WS6.1.
3. **CI + pkgdown + vignette** ("CMP detection in 10 minutes"); the 46-test suite exists but has no visible CI.
4. **Nextflow export of cmp-bench** for end-to-end reproducibility.

---

## 3. Sequencing

```
WS0 (close v1.2) ──► WS1 (caller arms) ─┐
                    WS2.1 (score uncertainty, caller-level) ─┤
                    WS3 (sharing module) ────────────────────┤──► cmp-bench v1.3 prereg
                    WS4.1/4.3 (artifact arms) ───────────────┘      + gene-sharing prereg
WS2.2 (adaptive caps; needs dual-gate prereg) ──► after v1.3 dev results
WS2.3 (guided init) ──► only if WS2.1+2.2 leave the R20 gap
WS4.2 (genotype validation study) ──► independent, can start immediately
WS5 (association module) ──► after WS3 (shares mashr machinery)
WS6.1 (project-and-call) ──► after caller frozen post-v1.3
WS6.2 (atlas) ──► after WS4 + WS6.1
```

**Milestone exit criteria:**

| Milestone | Exit criterion |
|---|---|
| V3-M1 (WS0) | subspace verdict recorded; PR #3 merged; no unpushed branches; param nullx complete |
| V3-M2 (WS1+WS2.1) | cmp-bench v1.3 prereg frozen; caller arms decided on dev; held-out confirmation done |
| V3-M3 (WS3) | `program_sharing()` shipped with calibrated matching; gene-sharing benchmark has its comparator arm |
| V3-M4 (WS2.2±2.3) | narrow-CMP miss rate on param R20 reduced from 0.71 with dual gate passed, **or** gap documented as irreducible at n≈100 donors with a power analysis to prove it |
| V3-M5 (WS4) | genotype hypothesis resolved (validated or filtered); ambient decoy ≤ 0.02 |
| V3-M6 (WS5+WS6) | association module validated on null phenotypes; project-and-call shipped; atlas v0 |

---

## 4. Decision framework

Everything enters `R/` only through a preregistered rule — the discipline that made v1.2 defensible stays:

1. **cmp-bench v1.3 prereg** (caller + upstream + artifact arms): same 8 scenarios, same guardrails (≤ 0.10 each), paired bootstrap CIs, per-breadth reporting, dev → heldout2 flow, `tradeoffs.csv` logging of every tried variant including rejections.
2. **Gene-sharing prereg** (WS3): new benchmark; false-sharing, false-membership, private-loss metrics; the mashr comparator is itself an arm.
3. **Dual-gate rule for fit changes** (WS2.2+): cmp-bench improvement AND model-selection-benchmark non-degradation, with the scope of any recommendation stated explicitly (the cap-30 deviation, generalized into a rule).
4. **Deviation discipline:** any post-hoc rule change is logged in DEVIATIONS.md *before* confirmation runs.

## 5. Top risks

| Risk | Mitigation |
|---|---|
| Narrow-CMP gap is irreducible at ~100 donors (power, not method) | WS2 power analysis is an explicit acceptable outcome (V3-M4) — a negative result with a number is a deliverable |
| mashr alignment step imports rotation-ambiguity errors into sharing claims | Match calibration is inside `program_sharing()`, permutation-tested; conservative default |
| Adaptive caps re-open the model-selection can of worms | Dual-gate rule; scoped recommendations; no default flips without both benchmarks |
| Genotype validation finds most new CMPs are artifacts | That is a *good* outcome for credibility — the filter then ships with evidence |
| Scope creep (6 workstreams) | WS0–WS3 are the core; WS5–WS6 can slip without blocking the science |

## 6. What v3 success looks like

1. Param R20 miss ≤ 0.50 (from 0.71) with all guardrails intact — or a power analysis showing the floor.
2. A shipped, calibrated `program_sharing()` answering shared-vs-specific with lfsr — and the coupled-EBMF workstream unblocked.
3. Every real-data CMP carries an artifact annotation (V(D)J/genotype/ambient) with a validated interpretation.
4. `cmp_associate()`/`cmp_mash()` replace all ad-hoc phenotype p-values.
5. A new cohort goes from pseudobulk to CMP calls in under a minute (project-and-call).
