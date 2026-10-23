# P3a — block_E permutation calibration: negative result (4 probe rounds)

**Goal.** Give block_E a calibrated FWER layer so its raw recall advantage
(+0.130 vs pa_cc) could be retained while fixing its disqualifying false
positive rates (null FPR 1.000, trap rate 1.000 in the prereg main run).

**Outcome.** No calibration tested achieves the prereg gates. The raw block_E
recall advantage is attributable to winner's-curse selection; removing the
curse removes the power. The prereg decision (pa_cc default) stands, now with a
mechanistic explanation. The 104-dataset rerun was **not** launched (gated on
probe pass), and the P3d confirmatory-scale run is **unmotivated** for a
negative result.

## Probe history (all on the prereg probe seeds: null R50 r1=5774,
## act_only R50 r1=1038, priv_corr R20 r1=6366)

| round | design | null spurious | signal found | privcorr traps | verdict |
|---|---|---|---|---|---|
| 1 | mass statistic, colnames bug | crash | crash | crash | fixed (flashier F_pm has no colnames) |
| 2 | name-based donor matching | 13/18 | TRUE | 5 | NO-OP: permuted rows carried rownames, re-sorted back |
| 3 | cross-block alignment stat, labels fixed | 7/18 | FALSE | 1 | winner's curse: full-data fit's noise factors inflate observed alignment |
| 4 | split-half (fit A, test B) | 0/18 | FALSE | 0 | sound but zero power for the true effect |

Round-4 detail (split-half, the principled design mirroring eb_cp's held-out
selection): refit `fit_joint_block` on stratified half-A donors, test
cross-block alignment of half-A factors on held-out half-B donors with max-T
permutation (199 perms, alpha 0.05).

- **null**: 18 raw spurious calls -> 0 calibrated; all factor p >= 0.805. Gate PASS.
- **priv_corr**: 17 raw calls / 5 traps -> 0 calls, 0 traps; every decoy-pair
  factor p >= 0.255. Gate PASS (decoy-ambiguity limitation moot at 0 calls).
- **signal (act_only R50, amp 0.5)**: raw found=TRUE (17 spurious) ->
  calibrated 0 calls, found=FALSE; best factor p = 0.555. Gate **FAIL**.

## Mechanistic diagnosis (diag_splitpower.R, act_only R50 r1 seed 1038)

The split-half design works as intended — the power loss is the family-wise
correction itself, not a residual curse:

- The half-A fit fragments the true program (5 member cell types -> best
  cluster holds 3, Jaccard 0.6 across 16 half-A clusters), but its factor is
  genuinely detectable: T_obs = 0.571 on held-out half-B donors, **uncorrected
  permutation p = 0.010**.
- The max-T correction over the ~15 callable factors is what kills it: the
  permutation null's max-T reaches 0.613 (95th pct 0.326), above the true
  factor's 0.571, so its family-wise corrected p is 0.555.
- Ceiling check: the true program's own z vectors have cross-block alignment
  1.0 on half-B — the effect is present; the FWER control over the noise-factor
  family consumes the power.
- For contrast, the same factor's full-data alignment is 0.984 — the raw
  path's in-sample selection value, i.e. the winner's curse the calibration
  removes.

At this effect size (amp 0.5, breadth R50, 100 donors), the true effect's
out-of-sample alignment sits inside the max-T null range: FWER control and
useful power are not jointly achievable for block_E here.

## Decision

- Keep block_E **out of the package** in calibrated form; the raw variant stays
  documented as exploratory with its FPR failure explicitly recorded.
- pa_cc remains the default joint method (prereg decision table unchanged).
- The activity-ambiguity limitation of decoy pairs (no activity-based test can
  reject cor-0.7 activity on disjoint genes) is superseded: the binding
  constraint is power, not decoy ambiguity.
- P3d confirmatory-scale replication (96 reps/primary cell): **not launched** —
  there is no positive calibrated result to replicate.

Artifacts: probe logs (rounds 1-4) in the session record;
`diag_splitpower.R`, `probe_blockfix.R`, `methods.R` (split-half
implementation), `blockfix_run.R`, `analyze_blockfix.R` committed under
`benchmarks/coordination/` in the repo; `blockfix_grid.csv` (104 rows, unused)
alongside.
