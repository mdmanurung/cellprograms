# Independent C1 plan and pre-freeze review

Reviewer: `c1_plan_review`, a separate read-only agent. The source was inspected independently; no reviewer production edits.

## Plan review
Initial verdict: revise before preregistration. Accepted findings:
1. Complete-linkage and signed consistency are needed before calling a cluster a coordinated program.
2. Full private/nuisance truth and deterministic one-to-one matching need explicit definitions; adjusted nuisance truth is zero.
3. True-score reference power is independent of expression amplitude. Pilot amplitudes using fitted detection; ideal references are algorithm-specific, not universal optimal bounds.
4. Candidate selection belongs in the pilot; held-out main data confirm one selected candidate. Null and independent-activity/same-gene controls gate separately.
5. Categorical within-block adjustment permits exchangeability assumptions; arbitrary continuous adjustment and heterogeneous residuals need separate validation.
6. Absolute node precision/recovery, empty-call handling, concrete paired margins, and R100 non-loss are required.

After the explicit resolutions in C1_PLAN.md: **implementation may proceed**. Precision is pooled correct/called nodes per cell with dataset-bootstrap intervals; both R50/R20 gains must pass; signed cycles must be flagged or rejected. Those clarifications are adopted.

## Pre-freeze source review
Reviewed R/coordination.R, tests, simulator, evaluator and runner. The reviewer independently reproduced a duplicate score-list-name defect. Public validation now rejects duplicate names and cardinality mismatch; its regression fixture passes. Main records now retain their verified source manifest, and summarization requires an identical manifest. Macro uncertainty now resamples within fixed equal-weight cells.

The detector separates edge inference from cluster summaries and enforces complete, distinct-cell-type, sign-consistent memberships. The evaluator uses deterministic assignments and separate pa_cc membership-only scoring; its independent-oracle checks passed under R4_51. Root added fragmentation and nuisance-loading/activity-conflict fixtures. Constant-on-intersection score columns are excluded individually.

This review does not establish scientific qualification. Held-out evaluation and independent final review remain required before C1 can be labelled done.

## Prospective pilot extension
The scalar 1/2 same-gene/small-R20 pilots failed selection. Reviewer approved ten fresh datasets per gene mode for ordered profiles (1,1,4), then (2,2,4), keeping the original per-mode thresholds, backend, detector, matching and main qualification gates. Retain rejected scalar results and their descriptive stress cells. Explicitly limit any eventual support to the selected signal regime. The reviewer independently confirmed column-level intersection-constant exclusion and retracted the earlier pair-level wording suggestion. Parser failures after exploratory output finalization must remain visible as task failures; main sources remain frozen throughout execution.

## Conditional rank refinement and pilot futility
Independent discarded-fit diagnostics found concentrated unmatched score factors (effective donor support about 1.6–2.1) inflating Pearson max-T tails. On the same failed fit, a scratch per-intersection Spearman calculation recovered all three planted programs without fragments; this is a mechanism diagnostic, not validation.

If the balanced profile also fails, reviewer approved a prospective rank candidate: use Spearman for both observed and permuted score correlations; retain complete/sign-consistent cliques, native Pearson truth identity matching, backend and every main gate. The Fisher-like scaling orders permutation statistics and supplies no analytic p-value. Label monotonic rank activity coordination explicitly; correlation/min_abs_cor become Spearman values, ties use R's average ranks, provenance becomes node_rank_complete. Require fresh scalar 1/2 activity pilots (ten per mode) and fresh null/genes_only/raw nuisance pilots before selection. Reference: R stats cor documentation, https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cor.html.

Reviewer also approved deterministic pilot futility stopping: with ten planned replicates and at most one allowed spurious outcome, two validated spurious records make that profile ineligible regardless of remaining results. Once both profiles are irreversibly rejected, cancel only remaining owned pilot tasks, preserve observed/planned counts and cancellation receipts, and report incomplete sampling and rejection by futility. Sourced files cannot change until every pilot task is terminal. No main stopping rule or qualification gate changes are authorized by this review.

## Conservative operating threshold before main
Before any main output, reviewer approved locking candidate/default alpha=.025 with every qualification gate unchanged. At actual false-call probability .05, exact probability of passing <=10/200 in both independent controls is .339967; at .025 it is .975013. This motivates a conservative operating point below the empirical .05 ceiling, at a potential power cost governed by the unchanged recovery gates. Reuse the discarded rank pilot's stored probabilities/null draws to recluster candidate and ideal reference at .025, preserving immutable .05 records and labelling derived results. Historical pa_cc, subspace_sum and gene_match remain at their existing .05 levels (different inferential units/error corrections); disclose that distinction. No alpha selection from main outcomes. Public default/docs will be .025 before final checks and registration; sourced pilot4 files remain unchanged until terminal.
