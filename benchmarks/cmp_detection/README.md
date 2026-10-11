# cmp-bench: benchmark for coordinated multicellular program detection

Benchmark-first development of `find_coordinated_programs()` (tusoskill workflow). Two pre-registered confirmations:

- `PREREG_cmp-bench-v1.1.md`: selected the caller `n_agg_noanc` (rank-first scores, ambient filter, node-aggregate max-T). Test seeds 1001-1050: F1 0.577, all 11 checks passed (`results/test_*.csv`).
- `PREREG_cmp-bench-v1.2.md`: confirmed `max_factors = 30` upstream of the frozen caller. Held-out seeds 1051-1100: all 7 criteria passed (`results/heldout_v1.2_*.csv`).

## Design

- 8 scenarios (null, cmp_distinct, cmp_mixed, context, ambient, lr, weak, genes_only) x 2 bases. The parametric base uses simulated programs (amplitude 0.45). The semi-synthetic base uses Stephenson COVID-19 PBMC pseudobulk with planted donor signal (12 cell types; amplitude 0.875).
- Planted CMPs span 20%, 50% or 100% of cell types, plus a context-restricted CMP.
- Primary endpoint: membership F1 of called CMPs (marginal mode).
- Guardrails per base: null FWER, ambient decoy rate, within-context false calls and genes-only false calls, each <= 0.10.

## v1.2 held-out result

| Run | Pooled F1 | Parametric | Semi |
|---|---|---|---|
| `max_factors = 30` + `n_agg_noanc` | 0.705 | 0.475 | 0.934 |
| `max_factors = 10` + `n_agg_noanc` | 0.571 | 0.472 | 0.671 |
| node correlation baseline | 0.349 | 0.018 | 0.679 |
| principal angles baseline | 0.355 | 0.310 | 0.400 |

- Paired difference vs cap 10: +0.127 (95% CI 0.113-0.142); semi +0.251, parametric +0.003.
- Winner guardrails (parametric / semi): null FWER 0.04 / 0.00; ambient decoy 0.04 / 0.02; within-context 0 / 0; genes-only 0.06 / 0.
- Mechanism: at cap 10, real cell types often never get a program for the CMP member. Upstream recall on the semi base rises from 0.73-0.83 to >= 0.995.

The package default stays `max_factors = 10`; see the README section "Coordinated multicellular programs" for the reason.

`results/tradeoffs.csv` lists variants that gained F1 only at or past a guardrail. None was adopted.

## Scripts (`scripts/`)

Paths are hard-coded to the original sandbox (`/mnt/shared-workspace/cmp/...`); edit `config.R` and the path constants before rerunning.

- `sim_cmp.R`: scenario generator (both bases).
- `fit_cache_v2.R`: per-cell-type flashier fits by arm (`main` = cap 10, `K20`, `K30`, `K40`, `featall`, `hvg4000`, `pnormal`, `noadj`).
- `cand_*.R`, `candidates.R`, `registry.R`: candidate callers.
- `run_cand_v2.R`: scores callers on a split (`dev`, `test`, `heldout2`, `nullx`).
- `eval_cmp.R`: scoring.
- `v2_score.R`, `prereg_check.R` (v1.1), `prereg_check_v12.R`: scorecards and pre-registered checks.
- `real_data.R`, `real_summary.R`, `cross_cohort_v2.R`, `split_half.R`, `real_compare_K.R`: real-data runs (Stephenson, COMBAT), site replication, cross-cohort matching and cap-vs-cap CMP comparison.
- `figs_v2.R`: held-out figures.
- `joint_fit.R`, `mofa_helper.py`: MOFA+ pilot (rejected).
- `lr_eval.R`, `annot_lr.R`, `annot_theme.R`: annotation benchmarks.
