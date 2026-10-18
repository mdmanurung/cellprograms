# CP-3A CMP calibration — confirmation arm

Date: 2026-10-18 | replicates/scenario: 500 | n_perm = 499 | alpha = 0.05
Valid rows: 5500 | Determinism vs interrupted run: 4000/4000 replicates identical (any_sig, n_sig)

| Scenario | n | FWER [95% CI] | null-node FPR [95% CI] | plant power [95% CI] | Verdict |
|---|---|---|---|---|---|
| A0_realistic | 500 | 0.054 [0.036, 0.078] | 0.002 [0.001, 0.003] | — | CALIBRATED (CI covers alpha) |
| A0_sentinel | 500 | 0.050 [0.033, 0.073] | 0.013 [0.010, 0.018] | — | CALIBRATED (CI covers alpha) |
| A1_realistic | 500 | 0.044 [0.028, 0.066] | 0.002 [0.001, 0.002] | — | CALIBRATED (CI covers alpha) |
| A1_sentinel | 500 | 0.044 [0.028, 0.066] | 0.012 [0.009, 0.016] | — | CALIBRATED (CI covers alpha) |
| B1_realistic | 500 | 1.000 [0.993, 1.000] | 0.001 [0.001, 0.002] | 1.000 [0.998, 1.000] | STRONG FWER CONTROL (null-node FPR CI <= alpha) |
| C1_realistic | 500 | 0.058 [0.039, 0.082] | 0.002 [0.001, 0.003] | — | FWER estimate 0.058 [0.039, 0.082] vs alpha 0.05 -> document |
| C2_realistic | 500 | 0.040 [0.025, 0.061] | 0.001 [0.001, 0.002] | — | FWER estimate 0.040 [0.025, 0.061] vs alpha 0.05 -> document |
| C3_realistic | 500 | 0.040 [0.025, 0.061] | 0.001 [0.001, 0.002] | — | FWER estimate 0.040 [0.025, 0.061] vs alpha 0.05 -> document |
| C4_realistic | 500 | 0.050 [0.033, 0.073] | 0.002 [0.001, 0.002] | — | FWER estimate 0.050 [0.033, 0.073] vs alpha 0.05 -> document |
| D1_realistic | 500 | 1.000 [0.993, 1.000] | 0.002 [0.001, 0.002] | 1.000 [0.998, 1.000] | NA |
| D2_realistic | 500 | 1.000 [0.993, 1.000] | 0.002 [0.001, 0.003] | 0.825 [0.805, 0.844] | NA |
