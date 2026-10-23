# Spec-v2 completion evidence pack (REL-1, 2026-10-19)

Curated evidence for the completed spec-v2 milestone queue. The canonical
code state is on GitHub at immutable SHAs; this pack carries the evidence
that justifies the milestone claims, not a code copy.

## Final state

| Repo | SHA | GitHub |
|---|---|---|
| cellprograms | `95ead12` | https://github.com/mdmanurung/cellprograms |
| flashier-utils | `738eba1` | https://github.com/mdmanurung/flashier-utils |
| cellprograms-reproducibility | `c03ea3c` | https://github.com/mdmanurung/cellprograms-reproducibility |

(cellprograms code/docs are identical to `cce06d9`; `95ead12` adds only the
audit-ledger completion entries. The reproducibility bundle pins `cce06d9`.)

Pinned backend: `flashier@fd44811` (1.0.60); `cebmf_torch@34abac5`;
`torch==2.7.1+cpu`; `SEACells==0.3.3`; mcRigor = JSB-UCLA/mcRigor.

## Final suites (REL-1, at the release SHAs)

| Package | Result | Log |
|---|---|---|
| flashier-utils | FAIL 0, WARN 11 (pre-existing max-iteration notices in projection tests), SKIP 6 | `suite_rel1_flashier_utils.log` |
| cellprograms | FAIL 0, WARN 3 (documented var_floor legacy-backend routing), SKIP 1 (opt-in SEACells integration) | `suite_rel1_cellprograms.log` |

## Milestone outcomes

- **FU-0/FU-1/FU-2B** (earlier): backend preflight, FND oracles, grouped
  specificity validation — see NEWS_flashier_utils.md.
- **FU-2A/FU-2C**: opt-in membership burden; `center_factor_activity()`;
  `factor_features_global()` + `plot_factor_features(select="global")`.
  53 oracle-pinned expectations.
- **FU-3A**: RNA→ADT fixed-activity projection validated against a native
  oracle (per-protein effect cor > 0.98; held-out-cell reconstruction;
  explicit rank-deficiency/absent-cell/undeclared-unit refusals). Vignette
  section added. No ImmGenT numeric parity claimed.
- **FU-3B/CP-6**: exact chunked projection is a NO-GO (query-set-dependent
  empirical-Bayes priors; pooled variance at var_type 0/2). Decision record:
  flashier-utils `reports/FU3B_CP6_CHUNKING_DECISION.md`. Documented in both
  projection functions.
- **CP-3A/CP-3B**: CMP calibration confirm, 5500 replicates — A0/A1
  calibrated, B1 strong FWER control supported, C-family boundary FWER
  documented (docs-only), D1 power 1.000 / D2 0.825. Analysis:
  `../spec_v2_implementation/CONFIRM_ANALYSIS.md`.
- **CP-4**: sharing sensitivity confirm, 20 reps/arm — all preregistered
  gates PASS; N3/N4 small-support dilution documented (Amendment 1).
  Analysis: `../spec_v2_implementation/CP4_CONFIRM_ANALYSIS.md`.
- **CP-5**: cross-cohort replication case-study article (evaluated;
  fixed-basis vs de-novo pathways, abstention accounting, three negative
  controls, leave-one-cohort-out). No new export — existing API suffices.
- **Metacells fold-in**: cellprogramsMC merged into cellprograms
  (SEACells + mcRigor QC, basilisk Suggests).
- **REL-1**: NEWS merged/finalized, limitation registers
  (`LIMITATIONS_*.md`), confirm-analysis label fixes, pkgdown site rebuilt
  and committed, repro-repo pins updated.

## Contents

- `NEWS_cellprograms.md`, `NEWS_flashier_utils.md` — release notes at REL-1
- `LIMITATIONS_cellprograms.md`, `LIMITATIONS_flashier_utils.md` —
  normative limitation registers
- `suite_rel1_*.log` — final full-suite logs
- Confirm-scale benchmark evidence: `../spec_v2_implementation/`
  (preregistrations live in the package repo under `benchmarks/`)

## Post-release supplement (2026-10-20)

The REL-1 snapshot above is unchanged. The post-release hardening round
(64 review findings, SEACells integration, exports sweep, mcRigor
single-metacell root cause, OptimizeGammaSEACells Seurat-v5 fix,
archetype-collapse guard) is documented in `POST_RELEASE_hardening.md`
with its own suite logs (`suite_post_hardening_cellprograms.log`,
`suite_integration_cellprograms.log`,
`suite_integration_metacells_planted.log`). Current `main` HEAD
supersedes the SHAs above.
