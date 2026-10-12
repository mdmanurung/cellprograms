# Caller-arm dev verdicts (cmp-bench-v1.2, dev seeds 1-30, 180 datasets/arm)

All arms run on the winner's exact prep (`n_agg_noanc`: rank-first + ambient drop 0.3) so the
comparison isolates the caller statistic. Adoption rule (frozen): paired CI excludes 0 in the
positive direction AND guardrails <= 0.10. Scorecards: `dev_subspace_vs_winner.csv`,
`dev_ws1_ws21_vs_winner.csv`. Ledger: tusoskill iterations 61-65.

## Verdicts

| arm | semi F1 | paired vs winner [CI] | R20 | R50 | R100 | CTX | verdict |
|---|---|---|---|---|---|---|---|
| `n_agg_noanc` (winner) | 0.924 | - | 0.805 | 0.981 | 0.993 | 0.960 | - |
| `s_subspace` (maxT) | 0.460 | -0.444 [-0.490, -0.402] | 0.00 | 0.78 | 0.75 | 0.00 | reject |
| `s_subspace_ct95` | 0.502 | -0.406 [-0.453, -0.365] | 0.003 | 0.82 | 0.86 | 0.00 | reject |
| `n_2stage` | 0.924 | 0 [0, 0] | 0.805 | 0.981 | 0.993 | 0.960 | uninformative |
| `n_cauchy` | 0.931 | +0.006 [-0.003, +0.015] | 0.852 | 0.967 | 0.981 | 0.946 | reject (n.s.) |
| `n_agg_eb` | 0.924 | +0.0001 [0, 0.0004] | 0.805 | 0.981 | 0.993 | 0.960 | reject (no-op) |

## Findings

1. **`s_subspace` / `s_subspace_ct95` (WS0, joint-subspace sum-of-projectors).** Misses nearly
   all narrow (R20) and context programs (miss rate ~1.0) while recovering broad ones decently
   (R50/R100 0.75-0.86): the aggregate spectrum is dominated by the largest shared subspace
   directions, so a narrow program contributing one donor-coherent direction among ~15 factors
   per cell type never moves it. The liberal ct95 gate does not fix it (power problem, not error
   control) and costs param guardrail headroom (0.067 vs 0). Also ~47x slower (469 vs ~10
   s/dataset). Upstream recall is identical to the winner (0.9976), so the deficit is purely the
   statistic.
2. **`n_2stage` (WS1, two-stage pair-gate -> node test).** Exact algebraic identity of the
   winner: with the liberal gate at |T| >= 2 and the winner's tau = 2, every dropped edge
   contributes max(0, |T| - 2)^2 = 0 anyway. A meaningful two-stage test needs gate_t > tau.
3. **`n_cauchy` (WS1, Cauchy-combination node statistic).** The only arm that moved the needle,
   and in the right direction on narrow programs (R20 +0.047) -- exactly the WS2 gap -- but the
   overall paired CI includes 0, so it is not adopted. Guardrails <= 0.10 (semi genes-only false
   0.067 vs winner 0). Revisit only via a new prereg powered for R20 (e.g. n_perm = 999).
4. **`n_agg_eb` (WS2.1, EB score shrinkage by loading posterior SDs).** No-op under the
   rank-first prep: per-column rank-gauss equalizes score scale across donors, so
   SD-proportional shrinkage has nothing to bite on (median score change 0.45%; per-donor
   shrinkage tracks donor noise, r = 0.94). Uncertainty propagation is neutralized by the prep,
   not by the statistic; a raw-score prep would expose it but loses everywhere else.

Winner unchanged: `n_agg_noanc`. No held-out runs spent (nothing cleared dev).
