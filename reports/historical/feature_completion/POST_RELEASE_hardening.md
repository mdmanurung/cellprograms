# Post-release hardening round (2026-10-19/20)

Supplement to the REL-1 spec-v2 pack: everything below happened AFTER the
release SHAs pinned in MANIFEST.md and is tracked on GitHub `main`.

## Code fixes on cellprograms (in order)

| SHA | Fix |
|---|---|
| `6d39811` | All 64 code-review findings (report: `../report_cellprograms_code_review.md`) |
| `5143beb` | First CPMC_INTEGRATION=1 end-to-end run; 4 bugs (CsparseMatrix coercion, withr import, Seurat v5 variable features, skipped_groups row count) |
| `47841bb` | All-exports execution sweep: 50/50 exports, 48 PASS / 2 by-design WARN / 0 FAIL |
| `dc82460` | mcRigor fallback root cause: single-metacell effect (AggregateExpression ignores 1-value groupings) |
| `239abc3` | OptimizeGammaSEACells: integer-coded membership for mcRigor (Seurat v5 label mangling killed EVERY group) + single-metacell gamma drop with warning |
| `c1df192` | man page regen for the dropped-gamma return contract |
| (this commit) | MetacellsByGroupsSEACells archetype-collapse guard (soft assignments dropped with warning instead of dimnames crash) |

## Suite baselines after the hardening round

| Mode | Result | Log |
|---|---|---|
| Baseline | FAIL 0, WARN 11 (REL-1 warning set unchanged), SKIP 3 (opt-in integration) | `suite_post_hardening_cellprograms.log` |
| CPMC_INTEGRATION=1 | FAIL 0, SKIP 0, WARN 31 (historical integration baseline held) | `suite_integration_cellprograms.log` |
| metacells file only, integration | FAIL 0 (22 tests incl. planted multi-metacell QC + gamma optimization + collapse guard) | `suite_integration_metacells_planted.log` |

## New verification assets

- Planted-archetype integration fixture: 5 subpopulations per 50-cell
  stratum (private 6-gene marker blocks) so SEACells produces genuine
  multi-metacell output and the mcRigor QC path runs end-to-end; random
  Poisson data collapses to 1 metacell and cannot exercise it.
- OPTIMIZE full-path smoke test: mcRigor_OPTIMIZE on integer-coded
  membership (the post-fix buildmc contract) through
  `.mcrigor_optimize_pick`.
- Collapse-guard test: degenerate random-data fit drops misaligned soft
  assignments gracefully.

## Session infrastructure (not repo changes)

- Basilisk SEACells env persisted under `/workspace/.basilisk` via
  `BASILISK_EXTERNAL_DIR` (set by `/workspace/bin/R` wrappers on the test
  machine); fresh-session pickup 5.1 s vs 10-30 min rebuild. `$HOME`
  caches do not survive machine hibernation.
- worker-2 retired; its unique 1.9 GB `UCSF_reannotated.h5ad` checkpointed
  to session-shared storage (md5-verified) before deletion.
