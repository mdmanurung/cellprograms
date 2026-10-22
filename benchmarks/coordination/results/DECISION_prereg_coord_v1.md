# prereg-coord-v1 — Main Run Decision Table

**Status:** main run complete (272/272 datasets, exit 0) · summarized with the frozen
summarizer (`summarize_main.R`, `OUT_SUBDIR=main`) · decision by the frozen rule in
`benchmarks/coordination/PREREG.md` (prereg-coord-v1, FILLED 2026-10-21).

**Correction log:** the summarizer's paired-primary section originally pooled the
act_cor 1/0.9/0.7 sweep into the primary contrast (the filter omitted `act_cor == 1` and
the pairing key omitted `act_cor`). Fixed in `summarize_main.R` and rerun; the table below
uses the corrected, frozen endpoint (act_cor = 1 only). Per-family tables were never affected.

**Scale caveat (preregistered):** development scale — 32 replicates per primary family,
16 null, 24 priv_corr|R20. The prereg specifies 96 reps/primary cell for a confirmatory run.

## Frozen decision rule

- Hierarchy: `eb_cp` → `block_E` → attribution arms (`subspace_sum`, `gene_match`, `node_cor`) vs `pa_cc` baseline.
- Adopt an arm iff paired recall gain ≥ **+0.10** AND paired precision not worse by > 0.05
  AND FPR gate: P(any spurious call | null) ≤ **0.10 absolute** AND ≤ pa_cc + 0.03.
- Primary endpoint: membership recall/precision/breadth, R50 + R20 each-alone, `act_only`,
  act_cor = 1 (n = 64 paired). pa_cc baseline recall on this endpoint: **0.592**.
- Oracle floor (node-based arms only); act_cor sweep {1, 0.9, 0.7} as robustness axis.

## Decision table (primary endpoint, n = 64 paired replicates)

| Arm | Paired recall vs pa_cc (margin +0.10) | Paired precision (not worse by >0.05) | FPR null (gate ≤0.10, ≤pa_cc+0.03) | Trap priv_corr | Verdict |
|---|---|---|---|---|---|
| `eb_cp` | **−0.584** [−0.692, −0.477] | −0.641 [−0.755, −0.526] | 0.0625 (1/16) ✓ | 0.000 ✓ | **Not adopted** — calibrated but no power |
| `block_E` | **+0.130** [+0.004, +0.255] | +0.068 [−0.030, +0.167] | **1.000** (16/16, mean 15.9 spurious calls) ✗ | **1.000** (24/24) ✗ | **Not adopted** — passes both margins, fails specificity |
| `subspace_sum` | −0.297 [−0.394, −0.200] | — | 0.000 ✓ | 0.000 ✓ | Not adopted |
| `gene_match` | −0.592 [−0.696, −0.488] | — | 0.000 ✓ | 0.000 ✓ | Not adopted |
| `node_cor` | −0.389 [−0.501, −0.277] | — | 0.0625 ✓ | 0.042 ✓ | Not adopted |
| `oracle` (floor) | +0.408 [+0.304, +0.512] | — | 0.0625 | 1.000 (by construction) | Floor: signal detectable |
| `pa_cc` (baseline) | — (0.592 mean recall) | — | 0.0625 | 0.458 | Remains the default |

**Preregistered conclusion: no arm meets the adoption rule at this operating point.
`pa_cc` remains the package default for coordination calls.**

## Key secondary results

- **block_E is the power direction, but its margin is modest on the frozen endpoint.**
  +0.130 paired recall [+0.004, +0.255] — the point estimate clears the +0.10 adoption
  margin, with acceptable precision (+0.068). It dominates every family with shared
  activity: found-rate 0.97/1.00 on R20/R20small act_only (pa_cc 0.41/0.00), 1.00 at
  act_cor 0.7 R20 (pa_cc 0.00), 0.88 on R100. It is disqualified *purely* by specificity:
  on null data the shared-first pass returns ~16 factors and the block-mass > 0.1
  membership rule (no significance layer) calls all of them; on priv_corr it traps the
  correlated-private pair in 24/24 replicates. This is a calibration defect, not a
  factorization defect — a permutation-calibrated block-mass threshold (FWER layer) is
  the indicated fix before any re-adoption attempt.
- **eb_cp's failure mode (preregistered interpretation axis): the gene-sharing assumption,
  and beyond it, greedy selection under private-factor competition.** On `act_only`
  (disjoint member genes) the shared gene mode is structurally diluted — expected. But even
  on `act+genes` R50 (identical shared genes, where `gene_match` scores 1.00) `eb_cp`
  reaches only 0.125, and its held-out trace stops early once private structure dominates
  half A. The EB-CP implementation is well-calibrated (null 1/16, trap 0/24) — the
  bottleneck is component selection, not calibration.
- **pa_cc profile:** strong on R50 (0.94) and act+genes R50 (0.97), weak on R20 (0.41),
  zero on R20small and at act_cor 0.7; traps the correlated-private pair in 46% of
  priv_corr replicates (baseline is not gated).
- **Oracle floor** is 1.00 on all primary and act+genes families and ~0 on `genes_only`
  (correct: no shared activity exists there) — the simulator, scorer, and operating point
  are sound; shortfalls are method-side.
- **gene_match is the gene-reuse specialist:** 1.00 on act+genes R50, 0.94 on genes_only
  R50, 0.00 on all act_only families — exactly the secondary-endpoint complement to the
  activity-based arms.
- **Runtime (median s/dataset):** block_E 908 (dominates; 985 on priv_corr), gene_match 77,
  pa_cc 13, eb_cp 12, subspace_sum 8, node_cor 5, oracle 6. Total run ≈ 6.5 h for 272
  datasets on 12 shards.

## Files

- `summary_rates.csv` — per (family, method): spurious/decoy/trap rates, runtime.
- `summary_programs.csv` — per (family, method, program): found rate, recall, precision, breadth error.
- `summary_paired_primary.csv` — paired recall vs pa_cc with 95% CIs (primary endpoint, act_cor = 1).
