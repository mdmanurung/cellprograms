# W5 SEACells-builder battery: summary and gate verdicts

30 combos (3 dropout rates x 10 seeds), synthetic truth (5000 cells x 1500
genes, 50 donors, 3 broad 150-gene blocks + 1 rare 80-gene block at 5% + 2
nulls, Poisson UMI, usage noise 0.2, nested dropout masks). Builder: SEACells
0.3.3 (53 metacells per combo, all converged, align_mismatch=0, ~0 metacells
under 10 cells, 0 zero-count metacells). Fitter: cebmf_torch at the registered
operating point. Harness certified by the k-means cross-check (1.7e-5).

## Gate verdicts (registered pre-run; all FAIL, honestly reported)

- **G1 rare-lift transfer: FAIL at all rates** (bar: rare_abs >= 0.7x per-rep
  top-10-SVD ceiling in >=8/10 seeds). Rate 0.0: 3/10 pass (mean ratio 0.658,
  range 0.573-0.760 — near-miss). Rates 0.1/0.3: 0/10 (mean ratios 0.092/0.053)
  — the rare program's metacell-level recovery collapses under dropout as its
  5% cell mass dilutes in ~71-cell pooled metacells.
- **G2 builder concordance: FAIL at all rates** (bar: mean matched |cor| >= 0.8
  over all K=10 Hungarian-matched pairs). Measured 0.427/0.411/0.330. The
  matched-pair detail shows WHY: the 5 real programs match at 0.67-0.86 (max
  0.86), the 5 noise columns at 0.01-0.09 — the mean over K=10 pairs is
  diluted by noise-noise matches (same K>p dilution as the CD8T route
  concordance). The registered metric is harsh for K=10 on 5 programs; the
  real-program concordance itself is moderate (0.67-0.86).
- **G4 soft-projection: FAIL at all rates** (bar: >=0.99 on strong factors).
  strong_mean 0.365, strong_min 0.011 (rate 0.0). The computation was audited:
  alignment verified (align_mismatch=0), soft rows sum to 1, per-factor
  correlations are real (not NaN; one noise column is NaN by zero variance —
  the known missing nan_to_num, conservative for the verdict). Mechanism
  hypothesis (Jensen compression): metacell factors are fit on SUM-aggregated
  counts while the exact projection is per-cell — log(1+sum) vs log(1+cell)
  are different nonlinear transforms, so soft-reconstructed and exactly
  projected cell scores disagree systematically. The 0.99 bar assumed
  aggregation-scale compatibility; the assumption fails. This is the
  quantitative form of the sum-then-log1p replay concern registered in the
  plan.
- **G3: NOT DIRECTLY MEASURED in this battery** — mcRigor dubious rates run in
  the W1/W2 smoke path (PASSED there); the battery records ZeroRate proxies
  only (per_combo_diagnostics.csv).

## Interpretation for M6

The SEACells builder itself is healthy (convergence, sizes, alignment all
clean), and the S4-seacells refinement arm is usable (s4_seacells.csv). But
the registered transfer gates fail: rare-program lift does not survive
metacell pooling at 5% mass under dropout, all-pairs concordance is diluted
below the bar, and soft-projection consistency is bounded away from 0.99 by
the aggregation-scale mismatch. These are the quantitative costs the metacell
plan anticipated; they motivate the plan's transparency requirements
(ZeroRate/K_eff reporting, fixed-K + pruning) and bound what the builder can
claim: metacells trade rare-state fidelity and projection exactness for
denoising/scale.

## Files

s3_seacells.csv, s4_seacells.csv, per_combo_diagnostics.csv,
g1_per_dropout.csv, g2_per_combo.csv, g4_per_combo.csv, runtime_per_rate.csv,
gate_verdicts.md. Raw JSONs on the shared store (w5_battery_out/, 30 files).
Note: the original aggregate_w5.py was lost to a worker-2 workspace snapshot
restore; this aggregation was reconstructed from the JSON schema
(aggregate_w5_v2.py) and reproduces the registered outputs.
