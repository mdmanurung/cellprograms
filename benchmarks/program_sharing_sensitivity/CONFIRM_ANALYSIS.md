# CP-4 program-sharing sensitivity — confirmation arm

Date: 2026-10-18 | rows: 379 (0 errors) | reps per scenario x arm: 20

| Scenario | Arm | n | shared-promotion [95% CI] | artifact frac (promoted) | Gate |
|---|---|---|---|---|---|
| P1_shared_biology | full | 20 | 1.000 [0.832, 1.000] | 0.00 | recovery 1.000 [0.832, 1.000] (gate >= 0.80): PASS |
| P1_shared_biology | prespecified | 20 | 1.000 [0.832, 1.000] | 0.00 | recovery 1.000 [0.832, 1.000] (gate >= 0.80): PASS |
| P1_shared_biology | posthoc | 20 | 1.000 [0.832, 1.000] | 0.00 | recovery 1.000 [0.832, 1.000] (gate >= 0.80): PASS |
| N0_unrelated | full | 20 | 0.000 [0.000, 0.168] | NaN | promotion 0.000 [0.000, 0.168] (defect if CI lo > 0.05): PASS |
| N0_unrelated | prespecified | 20 | 0.000 [0.000, 0.168] | NaN | promotion 0.000 [0.000, 0.168] (defect if CI lo > 0.05): PASS |
| N0_unrelated | posthoc | 20 | 0.000 [0.000, 0.168] | NaN | promotion 0.000 [0.000, 0.168] (defect if CI lo > 0.05): PASS |
| N1_ambient_ig | full | 20 | 1.000 [0.832, 1.000] | 1.00 | artifact flag on promoted groups = 1.00 (gate = 1.00): PASS |
| N1_ambient_ig | prespecified | 20 | 0.000 [0.000, 0.168] | NaN | prespecified promotion 0.000 [0.000, 0.168] (gate = no promotion): PASS |
| N1_ambient_ig | posthoc | 20 | 0.600 [0.361, 0.809] | 1.00 | posthoc descriptive: promotion 0.600 [0.361, 0.809], mean posthoc_n_drop=15.5 |
| N2_vdj_lineage | full | 20 | 0.000 [0.000, 0.168] | NaN | full-arm promotion 0.000 [0.000, 0.168] (gate = 0): PASS |
| N2_vdj_lineage | prespecified | 20 | 0.000 [0.000, 0.168] | NaN | prespecified arm promotion 0.000 [0.000, 0.168] (documented cost) |
| N2_vdj_lineage | posthoc | 20 | 0.000 [0.000, 0.168] | NaN | posthoc arm promotion 0.000 [0.000, 0.168] (documented cost) |
| N3_sex_shared | full | 20 | 0.000 [0.000, 0.168] | NaN | no full-arm promotion; artifact-flag check vacuous (prereg expected promotion+flag; support dilution per Amendment 1): DOCUMENTED |
| N3_sex_shared | prespecified | 20 | 0.000 [0.000, 0.168] | NaN | prespecified promotion 0.000 [0.000, 0.168] (gate = no promotion): PASS |
| N3_sex_shared | posthoc | 20 | 0.000 [0.000, 0.168] | NaN | posthoc descriptive: promotion 0.000 [0.000, 0.168], mean posthoc_n_drop=11.4 |
| N4_genotype_shared | full | 20 | 0.000 [0.000, 0.168] | NaN | no full-arm promotion; artifact-flag check vacuous (prereg expected promotion+flag; support dilution per Amendment 1): DOCUMENTED |
| N4_genotype_shared | prespecified | 20 | 0.000 [0.000, 0.168] | NaN | prespecified promotion 0.000 [0.000, 0.168] (gate = no promotion): PASS |
| N4_genotype_shared | posthoc | 20 | 0.000 [0.000, 0.168] | NaN | posthoc descriptive: promotion 0.000 [0.000, 0.168], mean posthoc_n_drop=11.8 |

## Q1 two-axis sub-fixture

- sharing fires in both regimes: TRUE
- CMP fires only with correlated donor scores: TRUE
- Q1 gate: PASS
