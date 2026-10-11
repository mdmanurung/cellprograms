# Pre-registration: cmp-bench-v1.2 (fresh held-out confirmation of the v2 winner)

Written 2026-10-10, before any data for seeds 1051-1100 are simulated, fitted, or evaluated.

## Candidate under test (frozen)
- **v2 winner `K30 + n_agg_noanc`**: per-cell-type flashier fits identical to the v1 finalist except `max_factors = 30` (was 10).
  - Caller: the unchanged v1 finalist `n_agg_noanc` (rank-first, ambient filter amb_thr 0.3, tau 2, link_threshold 3, n_perm 199, alpha 0.05, min_shared 20, no anchors).
- **Selection basis:** dev screen only (seeds 1-15 K20; seeds 1-8 K30; full dev 1-30 and null calibration 31-130 to be completed before unblinding).
  - Choice of 30 over 40: 30 already reaches upstream recall 1.00 on the Stephenson-based base, and its fits are about 20% faster. The 40-program screen is recorded for information and does not change the candidate.
- **Simulated-base fits** at any cap >= 20 are identical to the cap-20 fits whenever no cell type reaches 20 programs. This was verified bitwise on one dataset (lr s11, 13 programs); the observed maximum is 14. Simulated-base held-out fits are run at cap 30.

## Comparators
- **v1 finalist:** `n_agg_noanc` on cap-10 fits (main arm), same datasets.
- **Repo baselines:** `b_node_cor` and `b_pa_cc` on cap-10 fits.
- **Sanity check:** `s_shuffled` on cap-30 fits.

## Data
- Seeds 1051-1100 (50 per base), all 8 scenarios, both bases. Amplitudes param 0.45, semi 0.875. Generator and seed offsets are unchanged from v1.1.
- Seeds 1001-1050 (v1.1 test) are not reused.

## Primary endpoint
- `membership_F1_pooled`: mean target F1 over positive scenarios in marginal mode.

## Criteria (all must hold to adopt cap 30 as the package default)
1. The winner beats the v1 finalist: paired bootstrap 95% CI of the per-dataset F1 difference > 0 (2000 resamples).
2. The winner beats both repo baselines: paired CI > 0.
3. Null FWER <= 0.10 per base (null scenario, 50 datasets per base).
4. Ambient decoy <= 0.10 per base.
5. Within-context false call rate <= 0.10 per base.
6. genes_only false call rate <= 0.10 per base.
7. `s_shuffled` F1 < 0.02.

## Prohibited
- No parameter is changed after unblinding.
- If any criterion fails, the v1 finalist (cap 10) remains the benchmarked default, and cap 30 is reported as dev-only evidence.
