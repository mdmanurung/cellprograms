# C1 implementation checks

Checkout: isolated feat/c1-complete, baseline 584c502, API commit 0ffa453.
Runtime: R4_51, R 4.5.1, flashier 1.0.60, ebnm 1.1.42. Exact versions are recorded in benchmarks/coordination/runtime.csv.

- R CMD build and R CMD check --no-manual: Status OK, including installed tests and runnable coordination example. See c1_package_build.log and c1_package_check.log.
- testthat::test_local: 180 passing expectations, zero failures or skips. Four expected warnings cover the existing nullcheck fallback, small donor count, zero gene variance and small gene intersection. See c1_source_tests.log. The installed testthat was built under R 4.5.2 and emits a version warning.
- Direct-source test_dir was unsuitable for the existing namespace-mocking test: one harness error, not a passing suite. Preserved in c1_direct_source_inapplicable.log. The corrected package-loaded run and installed-package run pass.
- Evaluator oracle, maximum-weight counterexample, fragment handling and nuisance/activity identity-conflict checks pass: c1_evaluator_tests.log.
- Qualification pipeline self-check: 1,400 synthetic records qualify; altered false controls reject; missing output rejects. c1_summary_selfcheck.log. This tests report plumbing, not scientific qualification.
- Existing EBMF scores, loadings and sharing summaries match baseline on the strongest small fixture at tolerance 1e-12. Detailed receipt: c1_compatibility.rds. This is fixture invariance, not real-cohort invariance; no numerical changes to existing fitting/sharing functions were made.
- Initial discarded pilots: scheduler receipt c1_initial_pilot_scheduler.txt preserves six failed tasks with trailing parser errors after their RDS finalization. Structurally valid exploratory records may guide pilot settings, never main qualification. Extension 25739463 uses the pilot3_sources.csv manifest and unchanged sources during execution.

Held-out benchmark, scientific qualification and final independent review remain pending. A green package check does not label C1 done. PDF manual compilation was not run.
