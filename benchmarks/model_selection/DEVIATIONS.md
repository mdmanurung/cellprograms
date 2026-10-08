# Deviations from PREREG.md (after tag prereg-v1)

Each entry: what changed, why, whether any benchmark result had been seen. Only pilot data (discarded reps 9001-9020) had been seen for all entries below.

## D1. Amplitude rule (2026-10-08)
- **Grid extended** for `semi` S1 and S2 with amp 0.6, 0.7, 0.8, 0.9: base power jumped from about 0 at 0.5 to 1.0 (S1) / 0.75 (S2) at 1.0, so no grid point was near 0.5.
- **Tie rule:** equal distance from 0.5 goes to the larger amplitude (the registered rule did not say; the smaller one gave power 0).
- **S3 amplitude = S2 amplitude** (S7 = S1, as registered). Reason: the pilot showed `base` is mostly not callable on `semi` S3 (callable in 10-20% of reps; C-D share about 31 donors with K = 10 + 10, chance cosine above 1), so base power on S3 cannot be tuned to 0.5 by amplitude. The registered rule would have picked the largest amplitude for a reason unrelated to signal strength. Using the S2 amplitude is neutral to the arms.
- Consequence to report with the results: on `semi` S3, `base` power is capped by callability; the K rules are expected to help mainly through callability.

## D2. Chosen amplitudes (2026-10-08, from discarded pilot reps 9001-9020)
param: S1 0.5 (pilot power 0.60), S2 1.0 (0.60); semi: S1 0.7 (0.45), S2 0.8 (0.50). S3 uses the S2 value, S7 the S1 value (D1). Stored in `amplitudes.csv`, produced by `pilot_amp.R`.
