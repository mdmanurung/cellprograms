# Confirmatory decision — coordination benchmark (PREREG_confirm_v1.md)

**Date**: 2026-10-22 · **Run**: 288/288 jobs, exit 0 (12 shards, serial BLAS, 64 GB machine)
**Grid**: 96 act_only R50 + 96 act_only R20 + 48 null + 48 priv_corr, amp 0.5, act_cor 1, seeds from 1038.

## §E gate evaluation for `block_E_cal`

| Gate | Criterion | Result | Verdict |
|---|---|---|---|
| 1 | recall ≥ +0.10 vs pa_cc | diff = **−0.040** (95% CI [−0.113, +0.033]) | **FAIL** |
| 2 | precision not worse by > 0.05 | R20: +0.063 (pass); R50: **−0.203** | **FAIL** |
| 3 | null FPR ≤ 0.10 abs and ≤ pa_cc + 0.03 | block_E_cal made **no calls** on null (P_any_spurious NA); pa_cc = 0 | PASS (vacuous) |
| 4 | trap ≤ pa_cc + 0.03 | block_E_cal = 0.646 vs pa_cc = 0.500 (**+0.146**) | **FAIL** |

## Decision

**`block_E_cal` is REJECTED.** It fails three of four gates. The within-strata permutation
max-T calibration does not rescue the block_E family: it reduces raw block_E's catastrophic
spurious rate (mean 15.6 → 5.0 spurious calls) but at the cost of recall *below* the pa_cc
default and a trap rate 0.146 above pa_cc on priv_corr.

This corroborates the independent split-half probe (other session): max-T calibration on the
block_E statistic kills the true factor. `pa_cc` remains the default coordination caller.

## Reference points (same run)

- `oracle`: recall 1.000, precision ≥ 0.997, null FPR 0 — ceiling confirmed.
- `block_E` (raw): recall 0.734 but P_any_spurious = 1.0 on every family incl. null — disqualified.
- `pa_cc`: recall 0.579, null FPR 0, trap 0.500 — retained as default.

## Caveats

- block_E_cal null row is NA (zero calls) — the FPR gate passes vacuously; the method is
  over-conservative on null, not well-calibrated.
- The double-dipping caveat in PREREG §E stands: the same replicates inform both the
  calibration diagnostic and the gate evaluation.
