# C1 side clone: audit trail

These files come from a separate clone (`.c1-work/`, branch `feat/c1-complete`) that tried a
different coordinated-program detector (`node_rank`: Spearman edges, one joint donor
permutation, complete sign-consistent cliques, alpha .025) with a qualification harness. It
never reached a held-out run. The package default caller on `release/0.3.0` is `n_agg_noanc`
(`benchmarks/coordination/results/DECISION_default_caller.md`); C1 code conflicts with it and
is not part of the package.

What is here, and only as a record:

| Path | What it holds |
|---|---|
| `C1_PLAN.md`, `C1_PLAN_REVIEW.md` | The C1 plan, the review, and the alpha = .025 lock with its exact-binomial rationale (`C1_PLAN_REVIEW.md`, section "Conservative operating threshold before main"). |
| `pearson_pilot/` | The discarded Pearson-statistic pilot: the frozen spec, per-cell metrics and summary. Scalar amplitudes 1 and 2 failed the same-gene and small-R20 selection thresholds; the profiles (1,1,4) and (2,2,4) were then piloted. Read `pilot_pearson_summary.csv` for their results; I did not summarize them here. |
| `validation/` | Test, build and check logs of the C1 implementation and scheduler receipts of the pilots, including the retained failures. |
| `code/evaluate.R`, `code/test_evaluate.R` | A node-to-truth scorer (exact one-to-one assignment of fitted nodes to true score columns, `|r| >= .5`). **Not runnable on release/0.3.0**: it reads `d$truth$scores`, `d$truth$loadings` and a `nuisance_group` covariate that only the C1 simulator provides. Porting it is an open row in the prereg-coord-v2 draft. |

Not copied (132 MB): the raw pilot records under `.c1-work/benchmarks/coordination/results/`.
They are in `/exports/para-lipg-hpc/mdmanurung/cellprograms_backup_2026-11-17/untracked_local_material.tar.gz`.
The two C1 commits are kept as git tags `archive/c1-0ffa453` and `archive/c1-e95f1d7`.

## One finding that applies to release

`benchmarks/coordination/sim_coord.R:78-79` sets each member's activity to
`act_cor * z_shared + sqrt(1 - act_cor^2) * noise`. `act_cor` is therefore each member's
correlation with the shared driver, and the pairwise correlation between two members is
`act_cor^2` (0.49 at `act_cor = 0.7`). The primary endpoint of prereg-coord-v1 uses
`act_cor = 1`, where both agree. The 0.9 and 0.7 robustness axis means pairwise 0.81 and 0.49.
