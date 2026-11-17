# Decision: default caller in `find_coordinated_programs()` (0.3.0)

**Status:** decided by the user on 2026-11-17. Supersedes the line "`pa_cc` remains the
package default for coordination calls" in `DECISION_prereg_coord_v1.md` as far as the
package default goes. `prereg-coord-v1` itself stands as recorded.

**Decision:** `find_coordinated_programs()` keeps its implemented caller `n_agg_noanc`
(rank-Gauss scores residualized on strata, ambient filter, node statistic
`A_v = sum_c' (max_u |T_vu| - tau)_+^2`, max over nodes under within-strata permutation).
`pa_cc` stays the baseline of `prereg-coord-v1` and stays reachable through
`principal_angles()`, `sharing_spectrum()` and `align_programs()`.

## Why the two records differ

| | `prereg-coord-v1` (`benchmarks/coordination/`) | `cmp-bench` v1.1 / v1.2 (`benchmarks/cmp_detection/`) |
|---|---|---|
| Arms | `eb_cp`, `block_E`, `subspace_sum`, `gene_match`, `node_cor` vs `pa_cc` | `n_agg_noanc` vs `b_node_cor`, `b_pa_cc` and cap variants |
| `n_agg_noanc` an arm? | no | yes (the caller under test) |
| Endpoint | membership recall, `act_cor = 1`, R50 and R20 | pooled membership F1 over 8 scenarios |
| Data | generated, 32 replicates per family (development scale) | generated parametric base plus semi-synthetic base from Stephenson pseudobulk, 600 positive datasets |

Neither run compares the two callers on the same simulator, so the user decision rests on the
`cmp-bench` held-out run and not on a head-to-head under `prereg-coord-v1`.

## Evidence used (held-out seeds 1051-1100, `PREREG_cmp-bench-v1.2.md`)

table fallback: the source figures (`tusoskill/cmp_v2/figures/`) are not in this repository,
and making a new plot needs an approved grill plan. The numbers below come from
`benchmarks/cmp_detection/results/heldout_v1.2_scorecard.csv`, `heldout_v1.2_criteria.csv`
and `heldout_v1.2_paired.csv`.

| Run (cap 10 fits unless noted) | Pooled F1 | Parametric F1 | Semi F1 | Param ambient decoy |
|---|---|---|---|---|
| `n_agg_noanc`, cap 30 | 0.7045 | 0.4747 | 0.9343 | 0.04 |
| `n_agg_noanc`, cap 10 | 0.5712 | 0.4718 | 0.6706 | 0.04 |
| `b_pa_cc` | 0.3554 | 0.3104 | 0.4004 | 1.00 |
| `b_node_cor` | 0.3487 | 0.0182 | 0.6791 | 0.78 |

- At cap 10 the caller beats `b_pa_cc` on pooled F1 (0.5712 vs 0.3554) on the same fits.
- Paired difference, cap-30 `n_agg_noanc` minus `b_pa_cc`: +0.351 (95% CI 0.328 to 0.373, n = 600).
- Guardrails of `n_agg_noanc` (cap 30), parametric / semi: null FWER 0.04 / 0; ambient decoy 0.04 / 0.02;
  within-context false 0 / 0; genes-only false 0.06 / 0. All 7 preregistered criteria passed.

## What this does not show

- All evidence is simulated or semi-synthetic. No real-data validity claim follows.
- The gain from cap 30 is on the semi-synthetic base (+0.251). The parametric base is neutral (+0.003).
  The package fitting default stays `max_factors = 10`.
- Per-replicate raw outputs of the cmp-bench runs were not recovered. The reports are summary-only.
- Edge |T| under cluster permutation is not a calibrated statistic. Node permutation p-values are.
- `n_agg_noanc` was never run as an arm in `prereg-coord-v1`. If you want the same-simulator
  head-to-head, port it as arm `node_agg` in `benchmarks/coordination/methods.R` (a 0.4 item).
