# Deviations from PREREG.md (after tag prereg-v1)

Each entry: what changed, why, whether any benchmark result had been seen. Only pilot data (discarded reps 9001-9020) had been seen for all entries below.

## D1. Amplitude rule (2026-10-08)
- **Grid extended** for `semi` S1 and S2 with amp 0.6, 0.7, 0.8, 0.9: base power jumped from about 0 at 0.5 to 1.0 (S1) / 0.75 (S2) at 1.0, so no grid point was near 0.5.
- **Tie rule:** equal distance from 0.5 goes to the larger amplitude (the registered rule did not say; the smaller one gave power 0).
- **S3 amplitude = S2 amplitude** (S7 = S1, as registered). Reason: the pilot showed `base` is mostly not callable on `semi` S3 (callable in 10-20% of reps; C-D share about 31 donors with K = 10 + 10, chance cosine above 1), so base power on S3 cannot be tuned to 0.5 by amplitude. The registered rule would have picked the largest amplitude for a reason unrelated to signal strength. Using the S2 amplitude is neutral to the arms.
- Consequence to report with the results: on `semi` S3, `base` power is capped by callability; the K rules are expected to help mainly through callability.

## D2. Chosen amplitudes (2026-10-08, from discarded pilot reps 9001-9020)
param: S1 0.5 (pilot power 0.60), S2 1.0 (0.60); semi: S1 0.7 (0.45), S2 0.8 (0.50). S3 uses the S2 value, S7 the S1 value (D1). Stored in `amplitudes.csv`, produced by `pilot_amp.R`.

## D3. Untested sharing tests are recorded, not scored as negative (2026-10-09)
- `arms.R` `.eval_fit` now writes `tested` and `fail_reason` per pair. Before, a `sharing_spectrum` error set `ss = NULL` and every pair scored `shared = FALSE`, so a failing arm looked perfectly calibrated (F = 0).
- `tested = FALSE` only when `sharing_spectrum` errored or the pair is absent (a cell type kept no program). Underpowered pairs are tested but not callable.
- `shared`, `callable` and `F` keep their registered definitions; no registered quantity changes. Touches a file frozen at `prereg-v1`. No benchmark result had been seen (smoke data only).

## D4. Coverage guardrail `frac_untested <= 0.05` (2026-10-09; user decision 2026-10-08)
- `summarize.R` adds a guardrail: the worst-scenario share of untested (arm, pair) rows must be <= 0.05, and `frac_untested` / `frac_untested_base` columns appear in `decisions.csv`. A missing or NA value fails the guardrail (results written before D3 have no `tested` column).
- Reason: with D3 an arm that fails to test is visible, but F is deliberately left as registered, so coverage needs its own gate. A broken arm is rejected on coverage, not rewarded with F = 0.
- Additive and conservative: it can only turn an adoption into a rejection. Touches a file frozen at `prereg-v1`. No benchmark result had been seen.

## D5. Named error when fewer than 2 cell types keep a program (2026-10-09)
- `R/principal-angles.R` `sharing_spectrum` now stops with "needs at least 2 cell types that retained a program; found N" instead of the opaque `n < m` from `combn`. Message only; no numeric effect. Test added in `tests/testthat/test-principal-angles.R`.
- The benchmark records this message in `fail_reason` (D3).

## D6. Main run drops `k_stab` and `w_S` (2026-10-09; user decision 2026-10-08)
- Both exceed the registered 3x time guard in the full-scale runtime pilot (1 dataset per base): `k_stab` 9.73x, `w_S` 4.51x; `k_20` and `vt_12` about 1.25x, `k_n` and `genes_after` under 1x. They could not be adopted whatever their effect. `k_stab` also prunes to K = 0 at the registered amplitudes.
- Main run: `ARMS=base,strata,none_strata,none,k_n,k_20,vt_12,genes_after`, via the submit-time `ARMS` override in `run.sbatch` / `run_rep.R`; no frozen file changes. About 324 CPU-h instead of about 1006.
- The dropping criterion (runtime) is independent of the endpoints. Only pilot data had been seen. `K_ARMS` and `R_ARMS` in `summarize.R` still list the two arms; they report NA and are not adopted.

## D7. Main run submitted in array chunks (2026-10-09)
- The cluster has `MaxArraySize = 125`, so `--array=1-450` is rejected. `run.sbatch` gains an `OFFSET` env (task line = OFFSET + array index) and each 450-line task file is submitted as 5 chunks (OFFSET 0, 100, 200, 300, 400). Same tasks, same seeds, same arms; scheduling only.
