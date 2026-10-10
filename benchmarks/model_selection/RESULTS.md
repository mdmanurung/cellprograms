# M2 main run: results (2026-10-10)

Harness at commit `584c502` (tag `prereg-v1` plus DEVIATIONS D3-D7). 900 datasets (450 `param`, 450 `semi`; S0 200, S1 50, S2 50, S3 100, S7 50), 8 arms, all paired. No failed tasks. `summarize.R` output: `results/decisions.csv`, `results/decisions_base_calibration.csv`.

## Decisions as registered (guardrails include D4)
**No arm adopted.** M3 (needs >= 2 adopted arms) is not triggered.

| Arm | Endpoint | `param` est [CI] | `semi` est [CI] | Outcome |
|---|---|---|---|---|
| `k_n` | dP on S3 (>= +0.10) | +0.20 [0.09, 0.31] | +0.50 [0.39, 0.61] | effect passes; guards fail: `semi` dF on S0 = +0.565 (F 0.14 -> 0.705), `semi` dR = -0.075, `param` untested 0.223 |
| `k_20` | dP on S3 | 0.00 | -0.05 [-0.10, 0.00] | no effect; `semi` time 1.8x |
| `genes_after` | dR on S1+S2 (>= +0.02) | -0.003 [-0.034, 0.028] | 0.000 [-0.007, 0.006] | no effect |
| `vt_12` | dR on S1+S2 | -0.408 [-0.464, -0.353] | -0.541 [-0.586, -0.493] | harmful |
| `strata` | dF on S0 | -0.240 [-0.305, -0.175] | -0.075 [-0.115, -0.035] | effect passes; rejected only by the D4 `frac_untested` guard on `param` (0.223 > 0.05) |
| `k_stab`, `w_S` | n/a | NA | NA | dropped from the run (DEVIATIONS D6) |

## Calibration and confounder check (S0, share of datasets with any shared call)
| | `base` | `none` | `none_strata` | `strata` |
|---|---|---|---|---|
| `param` | 0.295 (CP 0.23-0.36) | **1.000** | 0.045 | 0.055 |
| `semi` | 0.140 (CP 0.10-0.20) | 0.090 | 0.060 | 0.065 |

- `param`: the confounder check passes (`none` is fully inflated) and the strata null restores calibration.
- `semi`: the check **fails** (`none` is not inflated). A weak institute effect and an unaligned one cannot be told apart from this data. `summarize.R` does not compute this check; the numbers above were computed separately from the `rep*_pair.csv` files.

## Untested pairs (D3/D4)
`param`: 12-27% of (arm, scenario) pairs are untested for every arm, `base` included (a cell type kept no program). `semi`: 0%. The registered absolute guard (<= 0.05) therefore fails every `param` arm, including ones that are otherwise fine. A guard relative to `base` would have been the better design.

## Open decision
Treat `strata` as rejected (D4 as registered), or report it as adopted under a post-hoc `base`-relative guard, labelled as a deviation made after seeing results. Not decided at the time of writing; this file records the registered outcome.
