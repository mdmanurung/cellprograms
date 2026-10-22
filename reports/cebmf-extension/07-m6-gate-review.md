# 07 — M6 Gate Review: adopt / revise / abandon recommendation

Status: FINAL (2026-10-18). Scope: the cebmf_torch extension (M0-M3 evidence),
the metacell plan (W1-W7), and the real-data case studies. Baseline preserved
throughout (413e290); every code change ledgered in reports/audit/TASKS.md.

## Recommendation: ADOPT WITH REVISIONS

Adopt the cebmf_torch backend + the cebmf_ext wrapper as the program-discovery
engine for cellprograms, conditional on the revision set below. Abandon
nothing outright; retire the unsupervised k-selection path from
recommendation status (evidence below).

## Evidence table

| # | Evidence | Result | Verdict |
|---|---|---|---|
| E1 | M3 cebmf battery (270 fits: 3 granularities x 3 dropout rates x 10 seeds, Poisson UMI truth) | gates + divergence burden + heldout recorded (battery_cebmf_aggregated.csv) | engine recovers planted structure at registered operating point |
| E2 | DIVAS D1 from-scratch (30 fits) + D2 anchored refinement + SVD-30 ceilings + A4 normalization extension (50 fits pflog/clr) | tables reproduce ledgered means; ceilings = max-per-column convention | parity with published DIVAS workflow |
| E3 | S5 cNMF anchor battery (30 fits) | cNMF density-silhouette picks k_hat=2 vs 4 planted at dropout 0/0.1 (10/10 seeds), rare 5% program never recovered there; 4/10 seeds at 0.3 | unsupervised k-selection underselects -> fixed-K + pruning design supported |
| E4 | Granularity workstream F (12 arms) | cebmf route concordance 0.82 (metacell vs single-cell); hdwgna arms 0.96 internal, 0.72-0.76 vs cebmf; cnmf ~0.00 vs everything; ebmf route-robustness 0.19 | cEBMF is the most route-robust engine tested |
| E5 | Metacell builder (cellprogramsMC, SEACells 0.3.3 + mcRigor) | smoke PASSED end-to-end after 8 ledgered bridge fixes; builder health clean (convergence, sizes, alignment) | builder usable as companion package (basilisk isolation; GPLv2/GPLv3 deps never source-copied) |
| E6 | W5 transfer gates (30 combos) | G1 FAIL (rare-lift 3/10 at rate 0.0 near-miss, 0/10 under dropout), G2 FAIL (all-pairs mean diluted by noise columns; real programs 0.67-0.86), G4 FAIL (soft-projection 0.365 mean; Jensen compression between sum-aggregated and per-cell log scales), G3 not directly measured (smoke covers it) | metacells trade rare-state fidelity + projection exactness for denoising/scale; transparency requirements mandatory |
| E7 | CD8T aging case study (8 subsets x 2 routes, 99 donors) | 19/160 FDR<0.05; strongest: effector program (CST7/PRF1/KLRB1) declining with age in naive CD8T (beta=-6.88, FDR=7e-5); GZMA/GZMH increasing in GZMB+ Tem; ROR1/CTSV in RTE; honest caveats (ribosomal/MT content blocks; weak cross-subset sharing: 0 edges at >=0.8) | real-data utility demonstrated; interpretation requires the content-block caveat |
| E8 | fit_celltype_programs covariate audit | 46 PASS / 0 FAIL; 3 P2 + 1 P3 defects with proposed patches (NOT applied; user decision pending); T11b resolved FALSE-POSITIVE | covariate path scientifically sound; patches recommended before release |
| E9 | Numerical robustness | two NaN modes root-caused and guarded: tau -> 1e12 cap (backend), sebetahat -> 0 -> laplace-EBNM NaN (se-floor, deviation 7); seed-family fallbacks | guards required for production use (revision R1) |

## Required revisions (R1-R4)

- **R1 (adopted, in the wrapper)**: numerical guard set — tau precision cap
  1e8, sebetahat floor 1e-4 before laplace EBNM fits, explicit seed application
  (np.random.seed + torch.manual_seed; FitConfig.seed is a no-op), seed-family
  fallbacks. All registered as deviations with diagnostics.
- **R2 (proposed, pending user approval)**: covariate-audit patches F1-F4
  (replay center/scale subsetting; aliased-column fail-closed handling;
  frozen-spec center guard; interaction-term validation). Recommend applying
  before release.
- **R3 (design requirement)**: metacell-route transparency — fixed-K + pruning
  with ZeroRate/K_eff reporting per fit (E3, E6); never present metacell-route
  rare-program results without the G1-style transfer caveat.
- **R4 (documentation)**: aggregation-scale mismatch (E6 G4) documented in the
  metacell vignette: soft-projection scores are not interchangeable with
  exact per-cell projections when metacells aggregate by sum.

## Retired / not recommended

- Unsupervised k-selection (cNMF density-silhouette) as the primary K choice
  for rare-state discovery (E3).
- Plain ebmf as the route-robustness benchmark (E4: 0.19 pseudobulk-vs-
  single-cell vs cebmf 0.82).
- Source-copying SEACells/hdWGCNA into the MIT package (license rule;
  dependency/wrapper only).

## Open queue (not blocking adoption)

M4 1M-cell chunked run; M5 500k shard scaling; E-M3.2/E-M3.4; Zenodo 403 on
the DIVAS RDS (D2(i) blocked); upstream frozen-prior request (weak-factor
projection fidelity, M2); covariate patch application (R2).
