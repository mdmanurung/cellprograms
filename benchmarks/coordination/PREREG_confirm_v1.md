# prereg-coord-v1a — confirmatory amendment: block_E_cal

**Tag:** `prereg-coord-v1a` · **Date:** 2026-10-22 · **Status:** FILED before the
confirmatory run is launched. Amends `PREREG.md` (prereg-coord-v1) under its §8
execution-scale note and the dev-scale decision table
(`results/DECISION_prereg_coord_v1.md`). Everything not amended here inherits the
frozen prereg (simulator, adjustment position, stage-1, scoring, matching rule).

## A. Question

Does the permutation-calibrated block-mass threshold (`block_E_cal`) — the fix
indicated by the dev-scale decision — produce an adoptable coordination arm?
Dev scale: `block_E` clears both adoption margins (paired recall +0.130
[+0.004, +0.255]; precision +0.068) and fails only specificity: null FPR 1.000
(16/16 datasets, mean 15.9 spurious calls) and priv_corr trap rate 1.000 (24/24).

## B. Arms (focused)

| Arm | Role | Settings |
|---|---|---|
| `pa_cc` | baseline | unchanged (n_perm = 1999) |
| `block_E` | dev-scale reference | unchanged frozen settings (§5 of the prereg) |
| `block_E_cal` | the calibrated challenger | block_E + FWER layer: per called cluster, member activity v_ct = Y_used_ct %*% w_ct (w_ct = the shared factor's loadings restricted to ct's block); stat = mean pairwise \|cor\| across member types over shared donors; null = independent within-strata permutation of every member's donor rows (n_perm = 999, α = 0.05, min_shared = 10); max-T across called clusters → FWER p; clusters with fwer_p > α dropped, ids/activity columns renumbered; uncalibratable clusters (< 10 shared donors) retained with fwer_p = NA |
| `oracle` | detectability floor | unchanged (node-based arms only) |

`eb_cp` and the attribution arms (`subspace_sum`, `gene_match`, `node_cor`) are
frozen at their dev-scale verdicts (not adopted) and do not rerun.

Implementation note: `block_E` and `block_E_cal` share one `fit_joint_block` fit
per replicate (the joint fit dominates runtime; the calibration layer is matrix
products only). Both arms therefore consume identical block fits — the comparison
is the calibration layer alone.

## C. Scale and grid

- Primary: {R50, R20} × `act_only`, amp 0.50 (the preregistered final amplitude,
  §4), act_cor = 1 — **96 replicates each** (the prereg §8 confirmatory scale) → 192 datasets.
- Gates: `null` × 48; `priv_corr|R20` × 48.
- Total 288 datasets. Seeds from `make_grid.R`'s frozen generator
  (seed = 1000 + 37·k + r); the first 32 R50 replicates share seeds with the dev run.
- Estimated cost ≈ 1150 s/dataset (stage-1 212 s + pa_cc 13 s + block_E 908 s
  shared with block_E_cal + calibration ≈ 10 s + oracle 6 s) ≈ 92 core-h
  ≈ 6–8 h wall on 16 parallel shards.
- `out/amp_choice.csv` (consumed by `make_grid.R`) was lost with the recycled
  scratch volume and is reconstructed here at its frozen value: amp = 0.50 for
  R100/R50/R20/R20small (prereg §4 "Final amplitude = 0.50 for all breadths").
  This is a reconstruction of a frozen value, not a new choice.

## D. Endpoints

Unchanged from prereg §6, at the new scale:

- **Primary:** membership recall/precision/breadth on R50 + R20 each-alone,
  `act_only`, act_cor = 1, paired vs pa_cc (n = 192 paired replicate–program pairs).
- **Null FPR:** per-dataset spurious rate on the null family (48 replicates).
- **Trap:** call containing both members of a priv_corr decoy pair (48 replicates).
- **Secondary:** runtime per arm; calibration retention — how many block_E calls
  block_E_cal drops, by family (measures what the layer removes vs preserves).

Summarized with the frozen `summarize_main.R` (`OUT_SUBDIR=confirm`), unmodified.

## E. Decision rule (amended)

Adopt `block_E_cal` as the package default for coordination calls iff **all** of:

1. Paired recall vs pa_cc ≥ **+0.10** on the primary endpoint.
2. Paired precision not worse by > **0.05**.
3. Null FPR gate: P(≥ 1 spurious call) ≤ **0.10 absolute** AND ≤ pa_cc + 0.03.
4. Trap gate: trap rate ≤ **pa_cc + 0.03** (relative only — no absolute ceiling:
   the baseline itself traps at 0.458, and correlated-private trapping is an
   identifiability limit of activity-level evidence, not a calibration defect).

`block_E` (uncalibrated) is re-scored on the same four gates for reference; its
dev-scale null-FPR failure is expected to reproduce.

## F. Preregistered expectation (stated in advance)

The calibration layer tests cross-type activity alignment against a
donor-permutation null. Correlated-private decoy pairs (cor 0.7) **are aligned by
construction**, so block_E_cal is expected to retain them: trap rate ≈ block_E's
1.000, failing gate 4. The informative outcomes are (i) whether the null FPR
gate is met — the defect the fix targets — and (ii) whether the recall margin
survives calibration. If gates 1–3 pass and gate 4 fails, the preregistered
disposition is: *block_E_cal is null-safe but trap-vulnerable; adoption requires
gene-level or model-based evidence beyond activity alignment. This closes the
block_E program's re-adoption question; no further block_E arm iterations
without a new registration.*

**Addendum (pre-launch, 2026-10-22, before any registered-grid replicate):** a
harness-validation smoke on ONE null replicate (rep 999, seed 424242, outside the
registered grid) reproduced the block_E null defect (17/17 spurious calls) and
showed block_E_cal retaining 2 of those 17 calls (fwer_p ≤ 0.05). Mechanistic
expectation this raises, recorded in advance: the calibration null (row
permutation with the loadings fixed) does not capture the joint fit's
double-dipping bias — w_ct was fit on Y_ct, so the strongest factors' block
activities align by construction and survive the permutation test. The layer
therefore reduces (17 → 2 on the smoke) but may not eliminate null spurious
calls, and gate 3 may fail. **The registered rule stands unchanged**; this
addendum exists so the smoke observation cannot be mistaken for a post-hoc
finding. The smoke replicate's harness validations: all four arms produce
contract-conformant rows with no errors; the shared block fit works (block_E
442 s, block_E_cal 3 s); per-dataset cost ≈ 620 s (faster than the §C estimate).
