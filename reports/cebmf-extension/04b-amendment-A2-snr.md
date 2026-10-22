# 04b — Amendment A2 (2026-10-17, pre-run): synthetic generator SNR fix + no-CP10K sensitivity arm

Registered BEFORE any M3 battery result exists (only the S1 pilot ran). Gates unchanged.

## Pilot finding (rep 1, rate 0, beta 0)

S1 from-scratch fit (10 s, converged, k_eff=10 with no pruning — consistent with the
documented A5 limitation) recovered planted programs at matched |cor| 0.44-0.76,
top-50 Jaccard ~0 on broad programs. Diagnosis: NOT a fit or wrapper artifact —
(a) a 3x longer budget (maxit 300, patience 10) reproduced the identical solution;
(b) SVD itself, on the registered shifted-log recipe, caps recovery at 0.66-0.80 per
planted program; (c) SVD on plain log1p (no CP10K) reaches 0.87-0.92. The registered
A1 bar (>=0.9) is unreachable for ANY method on the generator as originally specified:
the base noise (sd 0.5 in log space) plus compositional attenuation under CP10K set a
ceiling below the gate. The rare program (C2 bar 0.7) sat exactly at the ceiling.

## Changes (generator only; gates, arms, metrics, protocol unchanged)

1. Base expression noise: sd 0.5 -> 0.25 (log space).
2. Per-factor usage noise: sd 0.3 -> 0.2 (E-design inheritance).
3. Expected effect: per-marker signal correlation ~0.98; SVD ceiling ~0.97-0.99 —
   the gates become reachable while the task stays non-trivial.
4. NEW sensitivity arm S1-nocp: S1 on a log1p-only recipe (no CP10K), rate 0, 10 reps —
   quantifies the compositional-normalization cost on factor recovery. Reported as a
   sensitivity, not a gate; the registered preprocessing contract is unchanged for all
   gated arms.

Original pilot numbers are preserved in the ledger entry for 2026-10-17. The nested
dropout design, donor structure, truth, seeds (100000+rep), and all gates
(A1-A5, B1-B3, C1-C4) are unchanged.
