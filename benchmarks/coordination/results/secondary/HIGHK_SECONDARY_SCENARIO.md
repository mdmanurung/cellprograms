# P3c — High-K secondary scenario (exploratory, outside prereg grid)

**Question.** Does the prereg conclusion (pa_cc default) hold when the private-program
burden is raised from 1-2 to 25 private programs per cell type (G=4000, stage-1
`max_factors=30`)?

**Design.** 24 datasets: 16 signal (`act_only`, breadth=`all`, amp=0.5, seeds
9000+38*(rep-1)) + 8 null (seeds 9500+38*(rep-1)). Arms: pa_cc, node_cor,
subspace_sum, gene_match, oracle. block_E and eb_cp were skipped: their frozen
caps (20/15 factors) are structurally misspecified at this burden and both are
under active diagnosis (P3a/P3b).

## Results

**P(found) on signal** (recall>=0.5 & precision>=0.5):

| method | P(found) | mean recall \| found | mean calls | ds with spurious |
|---|---|---|---|---|
| oracle | 1.000 | 1.00 | 3.00 | 0.00 |
| subspace_sum | 0.250 | 0.63 | 1.50 | 0.75 |
| pa_cc | 0.188 | 0.79 | 0.94 | 0.375 |
| node_cor | 0.021 | 0.50 | 0.13 | 0.06 |
| gene_match | 0.000 | — | 0.00 | 0.00 |

**Null (8 reps):** pa_cc, node_cor, gene_match, oracle all 0 spurious;
subspace_sum 1/8 datasets with 1 spurious call. No traps for any arm.

**Realized K (stage-1 flash, cap 30):** median 9 (IQR 5-13), max 21 — the cap
never binds even at 25 privates/ct. Flash's greedy pass self-terminates far
below the cap, independently corroborating the P3b eb_cp diagnosis
(self-termination is a property of the greedy pass, not a cap artifact).

**Runtime:** stage-1 dominates (~870 s median); ~2560 s total per dataset.

## Interpretation

1. The scenario is recoverable in principle (oracle 100%), so the power loss is
   in private/shared disentangling, not signal existence.
2. High private burden degrades every experimental method; pa_cc degrades most
   gracefully (0 null spurious retained, 19% strict-found on signal). Its calls
   are often partial (81% of signal datasets get >=1 call but only 19% meet the
   strict found bar).
3. subspace_sum trades precision for power (25% found, 75% of signal datasets
   with >=1 spurious call) — consistent with its main-run profile.
4. The prereg decision (pa_cc default) is unchanged; the REAL_COHORT_RUNBOOK
   guardrail "low-K regime" is reinforced: at high private burden, expect
   substantially reduced precision from all joint methods.

Artifacts: `highk_dataset_level.csv`, `highk_program_level.csv`,
`highk_realized_k.csv` in this folder; raw per-replicate CSVs in
`raw_replicates.tar.gz` under `out/highk/`.
