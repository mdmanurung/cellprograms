# CP-4 preregistration — gene-loading sharing and artifact sensitivity fixtures

Frozen 2026-10-18 before any confirmatory run. Mirrors the CP-3A protocol style:
fixed scenarios, fixed seeds, fixed decision rules, Monte Carlo uncertainty reported.

## Layer under test

`program_sharing()` (matching: `.match_programs_loadings()`; assessment: mashr on
canonicalized loading posteriors with `.estimate_V` null correlation). Simulations are
at the **loading level**: synthetic `cell_program_fit` objects carry known loading
matrices and posterior SDs directly. Justification (same argument as CP-3A): the
sharing layer consumes (B, S) matrices; EBMF fitting is upstream and separately tested.
Nothing here re-implements `program_sharing()`; every scenario drives the production
entry point.

## Fixed universe and design

- 440 genes: 300 background (`G001`-`G300`), 40 IG/ambient (matching the `vdj` regex),
  20 TCR (`TRBV*`/`TRAV*`, matching `vdj`), 10 sex (`artifact_gene_sets$sex` members),
  10 genotype (`artifact_gene_sets$genotype_marker` members), 60 state (`ST01`-`ST60`,
  in no artifact set).

**Amendment 1 (2026-10-18, shakeout phase, before any confirmatory run):** the state
set was enlarged from 20 to 60 genes. Shakeout with 20 state genes showed the positive
control systematically classified "partial" because `frac_top_shared` is computed over
`n_top = 50` genes: a 20-gene shared support can contribute at most 20/50 = 0.4. That
n_top/support-size interaction is itself a logged sensitivity finding (a program whose
support is smaller than `n_top` cannot reach the "shared" threshold on the top-fraction
metric); the positive control is instead given a support larger than `n_top` so the
end-to-end pipeline has an unambiguous shared case to recover.
- 4 cell types (CD8T, NK, B, Mono), 3 programs each. Program support genes load at
  N(1.2, 0.1); all other entries N(0, 0.02). Posterior SDs: 0.08 on support genes,
  0.15 elsewhere. Mono P3 is constructed degenerate (loadings ~1e-18) to exercise
  degeneracy logging.
- Caller configuration: `alpha = 0.05`, `n_perm = 199`, `n_top = 50`, `min_genes = 50`,
  `data_driven = TRUE`, defaults otherwise. Replicates: 20 per scenario
  (smoke mode runs reps 1-3 only). Seeds: `31000 + scenario_index * 100 + replicate`.

## Scenarios

| name | planted structure | expectation |
|---|---|---|
| P1_shared_biology | P1 of all 4 types shares the 20 state genes | matched group of 4, "shared", high `p_all_shared` |
| N0_unrelated | all programs independent | no group promoted to "shared" |
| N1_ambient_ig | P1 of all types shares the 40 IG genes (ambient contamination) | full arm fires AND is artifact-flagged; prespecified arm does not promote |
| N2_vdj_lineage | P1 of CD8T only carries the 20 TCR genes (real lineage biology) | no cross-type promotion; prespecified vdj exclusion removes it (documented cost) |
| N3_sex_shared | P1 of all types shares the 10 sex genes | as N1: fires with flag; removed by prespecified exclusion |
| N4_genotype_shared | P1 of all types shares the 10 genotype genes | as N1 |

## Arms (analysed as separate analyses, never pooled)

1. **full**: the 400-gene universe.
2. **prespecified**: universe minus `artifact_genes(universe, c("vdj","genotype_marker","sex"))`
   (320 genes), exclusion decided before seeing results.
3. **posthoc**: full arm first; genes flagged in any matched group's top-12 loadings
   (`flag_artifact_genes`, sets vdj/genotype_marker/sex) are then removed and the
   analysis rerun. The selection is data-dependent by construction and is logged as such.

## Logged per replicate x arm

feature-panel size and hash; matched-pool size and min q; per-group size,
classification, `p_all_shared`, `frac_top_shared`, `frac_top_private`; artifact-overlap
fraction of group top genes; degenerate programs; V identity vs estimated with null-gene
count; data-driven covariance use vs fallback; errors.

## Sub-fixture Q1 (two-axis distinction)

P1-style shared loadings with (a) independent donor scores across types, (b) correlated
donor scores (shared latent driver, rho = 0.6, N = 60 donors). CMP detection
(`find_coordinated_programs`, n_perm = 199, min_shared = min_donors = 10,
ambient = "none") must be negative in (a) and positive in (b) while sharing fires in
both: loading sharing and donor-activity coordination are distinct axes.

## Decision rules

- **N0**: promotion rate (any group classified "shared") with binomial 95% CI; a defect
  investigation is triggered only if the CI excludes 5% on the high side.
- **P1**: recovery (4-type group classified "shared" with `p_all_shared` > 0.5) should be
  >= 80%; otherwise a sensitivity defect is investigated.
- **N1/N3/N4**: full-arm "shared" promotion is EXPECTED (loading similarity is real);
  the gate is that every promoted group carries a nonzero artifact-overlap flag and the
  prespecified arm does not promote the nuisance program. A missing flag is a
  documentation/flagging gap, not a mashr calibration defect.
- **N2**: the lineage program must not be promoted cross-type in the full arm; its
  removal under prespecified vdj exclusion is reported as the documented cost of
  fit-time exclusion (spec: avoid automatic removal of real lineage signals — the
  fixture quantifies the tradeoff; no code change).
- **Q1**: any disagreement with the expected 2x2 pattern triggers investigation.
- All findings are reported with Monte Carlo CIs. Code corrections (CP-4 has none
  preplanned) only if a defect is demonstrated; otherwise results become documented
  limitations/flags.
