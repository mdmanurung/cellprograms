# cebmf_torch capability audit (source-grounded)

Repo: `github.com/william-denault/cebmf_torch` @ commit `dc5f30b` (HEAD of main at audit time).
Source root: `/workspace/cebmf_torch/src/cebmf_torch/`. All `file:line` refs are relative to repo root unless noted.
Audit performed on machine `worker-0`; smoke test run with the installed editable package.

---

## (a) Verified capabilities table

| # | Capability | Verdict | Evidence (file:line) |
|---|---|---|---|
| 1 | **Objective = EBMF variational (negative ELBO), not plain Frobenius** | Verified | `src/cebmf_torch/cebmf/cebmf.py:418-428` (`_cal_obj`: `loss = (-ll + KL).item()`); data term = expected normal-means loglik (`_compute_constant_loglik` :431-436, `_compute_elementwise_loglik` :438-442); KL terms `kl_l/kl_f` = `(-prior.loss) - normal_means_loglik(...)` (:370, :415); `normal_means_loglik` defined :681-753 |
| 2 | EM update loop | Verified | `iter_once` `cebmf.py:258-268`: for k in K → `_update_factors(k)` (:306-326) → `update_tau()` (:272-303) → `_backfit()` (:444-458) → `_cal_obj()` (:418). Coordinate-wise: update L[:,k] then F[:,k] against residual with factor k added back (`Rk = R + outer(L_k,F_k)`, :309-311, :318-320) |
| 3 | Factor update = weighted least squares + EBNM prior fit | Verified | `_update_L_factor` `cebmf.py:328-371`: `lhat = (W·Rk)@Fk / ((W)@Fk2)` (:337-346), `se_l = sqrt(1/(tau·denom))` (:341,:347); prior fit via `prior_L_fn.fit(...)` (:352-360); write-back `L[:,k]=post_mean`, `L2[:,k]=post_mean2` (:367-368). Mirror for F at :373-415 |
| 4 | Initialization: svd / random / zero / user-provided | Verified | `cebmf/_initialisation.py:155-159` (`INIT_STRATEGIES`), `svd_initialise` :62-92 (column-mean NaN imputation :34-51, zero-pad if K>rank :85-88), `random_initialise` :96-118 (scale `RANDOM_INIT_SCALE=0.01`, :7), `zero_initialise` :122-142, `user_provided_factors` :163-198. `initialise_factors` dispatch `cebmf.py:215-255`; `fit()` warns + SVD-inits if not initialised (:203-204) |
| 5 | Backfitting / factor pruning | Verified | `_backfit` `cebmf.py:444-458` (drops factors, never keeps 0 — keeps one if all prunable :449-451); `_should_prune_factor` :619-632: prune if **min π0 over coordinates on either side ≥ `prune_thresh`** (default `1-1e-3`, :17); `_prune_indices` :634-651 slices L, L2, F, F2, kl_l, kl_f, model states, π0 lists and sets `model.K = len(keep)`; resets `obj=[]` |
| 6 | Effective K detectable post-fit | Verified | `m.model.K` after fit (pruned count); per-factor π0 in `m.pi0_L` / `m.pi0_F` (initialized :516-518, set :371, :416). Smoke test: K stayed 3 with π0 ≈ 1e-4–1e-5 |
| 7 | Stopping rule: **fixed `maxit` only — NO objective-based convergence check** | Verified | `fit` `cebmf.py:193-212`: `for _ in range(maxit): self.iter_once()` (:210); `self.obj` history is recorded (:428) but never inspected for convergence anywhere in `src/` |
| 8 | Numerical safeguards | Verified | `NUMERICAL_EPS=1e-12` `cebmf.py:16`; denominator clamps :336,:349,:594; `R2.clamp_min(0)` :477; `R.nan_to_num_` :253,:679; stable loglik with eps + validity mask `normal_means_loglik` :681-753; LBFGS loss `nan_to_num(→1e30)` + line-search-failure fallback `ebnm/point_laplace.py:196-216`; `s.clamp_min(1e-12)` `ebnm/ash.py:231`; spike-only shortcut `point_laplace.py:262-272` |
| 9 | Prior registry: **11 names** | Verified | Registered names (runtime-confirmed): `norm, exp, laplace, gbinary, cash, cgb, cgb_sharp, emdn, lcash, po_lcash, spiked_emdn`. Registration loops `priors/registry.py:72-81`; `MixturePriorType` {NORM, EXP} `priors/mixture.py:17-19`; `PointPriorType` {LAPLACE, EXP, GBINARY} `priors/point.py:19-23`; `LearnedPriorType` {CASH, CGB, CGB_SHARP, EMDN, LCASH, PO_LCASH, SPIKED_EMDN} `priors/learned.py:22-42` |
| 10 | Prior → solver mapping | Verified | `norm`/`exp`(ASH) → `ash()` `ebnm/ash.py:164-241`; `laplace` → `ebnm_point_laplace` `ebnm/point_laplace.py:66-285` (spike+Laplace-slab, LBFGS on (w, a, μ), μ fixed 0 by default :69); `exp`(point) → `ebnm_point_exp` `ebnm/point_exp.py:118`; `gbinary` → `ebnm_gb` `ebnm/generalized_binary.py:34` (EM E-step :64, μ-optimization :95); `cash` → `cash_posterior_means` `cebnm/cash_solver.py`; `cgb`/`cgb_sharp` → `cebnm/cov_gb_prior.py:191`, `cebnm/cov_sharp_gb_prior.py`; `emdn` → `cebnm/emdn.py:162` (MDN); `lcash`/`po_lcash` → `cebnm/lcash.py:597,715` (softmax / proportional-odds nets); `spiked_emdn` → `cebnm/spiked_emdn.py:209` |
| 11 | Prior semantics | Verified | `norm`/`exp` = ASH grid mixtures (normal scale grid `utils/mixture.py:165`; exponential grid with spike at 0 `utils/mixture.py:207`); `laplace`/`exp`(point)/`gbinary` = spike-and-slab point-mass priors; `cash` = covariate-assisted spike-and-slab; `cgb`/`cgb_sharp` = covariate (sharp) generalized-binary; `emdn`/`spiked_emdn` = (spiked) mixture-density-network priors conditioned on covariates X; `lcash` = linear covariate ASH (multinomial-logistic π(x)); `po_lcash` = proportional-odds variant (`cebnm/lcash.py:73-147`) |
| 12 | **All priors return posterior 2nd moments** (post_mean2), not just means | Verified | `Prior` dataclass requires `post_mean`, `post_mean2` `priors/base.py:19-33`; every solver populates both (e.g. `point_laplace.py:273-275`, `utils/posterior.py:239-328`, `cebnm/lcash.py:464-466`, `cebnm/cash_solver.py:260-273`, `cebnm/emdn.py:291-307`) |
| 13 | Noise models CONSTANT / ROW_WISE / COLUMN_WISE / KNOWN | Verified | `NoiseType` `cebmf.py:43-47`; `update_tau` :272-303 (match → dim None/0/1; KNOWN = no-op :286-287); `_update_tau` :588-608: `tau = 1/mean(R2)` with mask-count normalization `m = mask.sum(dim).clamp_min(1)` :590-594. CONSTANT → scalar tau + broadcast `tau_map` view (:599-604); ROW/COLUMN → per-row/col tau expanded to (N,P) view (:605-608) |
| 14 | Known-SE mode (`S=`) | Verified | `S` forces `NoiseType.KNOWN` + warning `cebmf.py:156-165`; `_setup_known_variance` :523-586 broadcasts S to (N,P), `tau_map = 1/S²` :586; NaN/inf/≤0 S entries folded into mask with warning :565-583 |
| 15 | Missing data (NaN in Y) | Verified | `mask = (~isnan(Y)).float()`, `Y0 = nan_to_num(Y, 0)` `_initialise_tensors` `cebmf.py:499-500`; residual always `⊙ mask` (:253, :675-679); weighted updates use `W = tau_map·mask` (:339-345, :385-390); loglik restricted to observed entries (:438-442); SVD init imputes column means `_initialisation.py:34-51` |
| 16 | Covariates (external + self-covariate) | Verified | `CovariateParams` `cebmf.py:63-68`; `_build_covariate_matrix` :653-672 (hstack external + other factors; intercept fallback at K=1); passed as `X` to prior fit :352-356, :393-397 |
| 17 | Device selection CPU/CUDA/MPS | Verified | `utils/device.py:6-27` (`get_device`: CUDA > MPS > CPU); default `device or get_device()` `cebmf.py:151`; covariates moved :171-172; all core tensors created with `device=self.device` (:501-515); learned priors move nets `.to(device)` (`cebnm/cov_gb_prior.py:261`, `emdn.py:242`, `lcash.py:527`, `spiked_emdn.py:299`) |
| 18 | dtype: **hardcoded float32 for Y** | Verified | `self.Y = self.data.to(self.device).float()` `cebmf.py:185` — float64 input is silently downcast; mask float32 :499; truncnorm helpers keep caller dtype under `precision="auto"` (float32 default) `utils/maths.py:293-318`; L-BFGS mixture path upcasts to float64 `utils/mixture.py:305-338` |
| 19 | Mini-batching scope: **mixture-weight EM only** | Verified | `optimize_pi_logL(..., batch_size, shuffle, seed)` `utils/mixture.py:7-120` — mini-batch EM over rows of the **precomputed (J, K_grid) log-likelihood matrix**; batches accumulate expected counts `n_k` (:100-118). Wired via `AshConfig.batch_size=128` (default) `ebnm/ash.py:38` and `ash(batch_size=128)` :174. **Factor/residual/τ updates are always full-batch** (`cebmf.py:306-326` uses complete Rk, lhat, fhat) |
| 20 | Sparse / AnnData / chunked / streaming | **Absent** | `grep -rn sparse src/` → 0 hits; no anndata/scanpy/h5py/dask/zarr/chunk anywhere in `src/`. Dense `torch.Tensor` (N,P) only; `data` indexed as `self.Y.shape` `cebmf.py:186` |
| 21 | Tests | Verified | 30 files in `tests/`: rank-1 recovery (`test_basic_rank1_torch.py:16-27`), ASH posteriors (norm/exp), point-laplace/point-exp, generalized binary, CASH (loss stability, posterior), LC-ASH (component selection, marginal-loglik, warm start, state integration), GCBF, EMDN, known variance, NaN handling (`test_with_NA_torch.py`), GPU-residency regressions (`test_gpu_residency_improvements.py`), initializer, mixture optimizer, prior-builder isolation, pruning direction, truncnorm edge cases, doc-example generator |
| 22 | License MIT | Verified | `LICENCE.md:3` "Copyright 2025 William R.P. Denault" + MIT text (file spelled `LICENCE.md`, not `LICENSE`); `pyproject.toml:12` `license = { text = "MIT" }` |
| 23 | Dependencies | Verified | `pyproject.toml:20-22`: runtime = `torch>=2.0`, `tqdm` only; dev group :24-30 (matplotlib, numpy, pytest, pytest-cov, ruff); examples group :31-40 (jupyter, seaborn, torchvision, transformers, `uni`, …); `requires-python >= 3.11` :11 |

## (b) Limitations / gaps

1. **No posterior uncertainty in the result object (explicit gap).** `CEBMFResult` holds only `L, F, tau, history_obj` (`cebmf.py:21-39`). Per-entry posterior second moments **are** computed and maintained internally (`self.L2`, `self.F2`, `cebmf.py:502,504`, updated :368,:414), so `Var(L) ≈ L2 − L²` is recoverable from the model object, but it is not returned, and per-entry posterior **sd** from the EBNM solvers (e.g. `EBNMLaplaceResult.post_sd`, `point_laplace.py:60`) is dropped at the `Prior` boundary (`priors/base.py:19-33` has no `post_sd` field). **No LFSR / local-false-sign-rate anywhere in `src/`** (grep: 0 hits). Any downstream uncertainty/lfsr feature must be added.
2. **No objective-based stopping.** `fit` runs exactly `maxit` EM iterations (`cebmf.py:210`); the ELBO history (`self.obj`) is recorded but never tested. Also `obj` is reset to `[]` on pruning (:651), so history is discontinuous after backfitting.
3. **No sparse input, no AnnData/chunking/streaming, no mini-batching of the factorization.** The (N,P) matrix must fit in device memory as a dense tensor; every iteration materializes multiple (N,P) temporaries (see (e)).
4. **`exp` name collision.** Both `MixturePriorType.EXP` (ASH exponential mixture) and `PointPriorType.EXP` (point exponential) register the name `"exp"`; because point priors register second (`registry.py:76-79`), `cEBMF(prior='exp')` resolves to the **point-exponential** builder (runtime-verified: `get_builder('exp') → PointBuilder`). The ASH exponential mixture is reachable only via direct `ash(prior='exp')`.
5. **Hardcoded float32.** `.float()` at `cebmf.py:185` discards float64 input precision; there is no dtype option. Truncnorm internals stay float32 under `"auto"` (`utils/maths.py:313-318`).
6. **Pruning semantics are one-sided and conservative.** A factor is pruned when the *minimum* π0 across coordinates on **either** side ≥ `prune_thresh` (`cebmf.py:619-632`); with covariate priors π0 is a per-coordinate vector (`learned.py:154-163`), so a single strongly-null coordinate can trigger pruning of a factor that is non-null elsewhere. `prune_thresh` is a π0 (slab-weight) threshold, **not** an effect-size or variance-explained threshold.
7. **GPU path untested here** (worker-0 is CPU-only). Code review finds full device coverage (all tensors created on `self.device`; no `.cpu()`/`.item()` host syncs in the hot loop — the code deliberately keeps scalars as 0-d tensors, e.g. `point_laplace.py:53-55`, `utils/mixture.py:110-115`), but `torch.linalg.svd` on MPS (`_initialisation.py:79`) and CUDA LBFGS behavior were not executed.
8. **Determinism holds only with fixed seeds and default settings.** SVD init is deterministic; ASH mini-batch EM is deterministic by default (`shuffle=False`, `seed=None` → sequential batches, `utils/mixture.py:92-97`); point-laplace LBFGS is deterministic. `shuffle=True` without `seed` would break reproducibility (`ash.py:38-39`).

## (c) Untested claims (README claims not verified in this audit)

- **"GPU support for all operations" / "GPU-accelerated"** (`README.md:31,41`): plausible from source (see (a)#17) but not executed — no CUDA device available on the audit machine.
- **"Mini-batch EM and Adam optimizers for mixture weights"** (`README.md:44`) and **"Mini-batch EM for `pi` is implemented via Adam on logits (recommended) or online EM"** (`README.md:118`): **inaccurate vs. source.** The mixture-weight optimizer is EM (`optimizer="em"` default) or **L-BFGS** on softmax logits (`ebnm/ash.py:38-39`, `utils/mixture.py:264-351`); Adam appears only as an optional 8-step warm-start inside point-Laplace (`point_laplace.py:73-77`, default off). There is no "Adam on logits" path for ASH mixture weights.
- **"Posterior inference: compute posterior means and variances for all supported models"** (`README.md:45`): true at the EBNM-solver level (`post_mean2` everywhere), but variances are **not surfaced** through `cEBMF.fit()` (see gap 1) — the claim overstates the user-facing API.
- Docker GPU image, codecov badge, readthedocs: not verified (external).
- MPS (Apple Silicon) path: selected by `get_device` (`device.py:23-24`) but never tested in CI tests reviewed.

## (d) Smoke-test results (worker-0, CPU, torch float32)

Setup: 200×50 rank-3 matrix (`Ltrue@Ftrue.T + 0.1·N(0,1)`), `cEBMF(Y, K=3, prior_L='laplace', prior_F='laplace', noise_type=NoiseType.COLUMN_WISE)`, `fit(maxit=20)`.

- **Runs**: yes, 5.08 s wall (includes SVD-init warning; `fit()` auto-initializes via SVD, `cebmf.py:203-204`).
- **Result fields**: exactly `L (200,3) float32 cpu`, `F (50,3)`, `tau` — shape `(200,50)` for COLUMN_WISE (expanded `tau_map` view, `cebmf.py:605-608`), `history_obj` (len 20).
- **Objective**: `history_obj` decreases −12203.4 → −12125.2 and plateaus (last 3 values within 0.02) — behaving like the negative ELBO.
- **Pruning**: effective K stayed 3; π0 per factor ≈ 7.6e-5 / 2.7e-4 / 7.5e-6 (far below `prune_thresh = 0.999`), as expected for true rank 3.
- **Determinism**: two identical runs (same `torch.manual_seed(42)`) → **bit-identical L and F** (`torch.equal == True`, max abs diff 0.0).

## (e) Per-operation complexity (N rows, P cols, K factors, G = ASH grid size, E = LBFGS iters)

Per `iter_once` (`cebmf.py:258-268`), for each factor k:

| Operation | Where | Cost |
|---|---|---|
| Residual recompute `R = (Y0 − L Fᵀ) ⊙ mask` | `cebmf.py:675-679`, called 3×/factor (:309,:317,:323) | O(NPK) each → O(NP·K²) per iteration |
| `Rk = R + outer(L_k, F_k)` | :309-311, :318-320 | O(NP) |
| L/F least-squares stats (`mask@Fk2`, `Rk@Fk`) | :337-347, :383-392 | O(NP) |
| EBNM prior fit (ASH) | `ebnm/ash.py:70-106` | O(J·G) loglik + O(iters·J·G) EM, J = N or P |
| EBNM prior fit (point-Laplace) | `ebnm/point_laplace.py:66-216` | O(J·E) LBFGS over 3 scalars |
| τ update | :272-303, :588-608 | O(NP) |
| ELBO `_expected_residuals_squared` (2 matmuls) + loglik | :461-477, :418-428 | O(NPK) |

**Total ≈ O(NP·K²) flops per EM iteration** (dominated by the 3K full residual matmuls) + O(K·(N+P)·G) for prior fits; memory traffic ≈ 10–15 × (N,P) tensor passes per iteration.

**Resident memory, dense N×P float32 (4 bytes/entry), persistent tensors:** `Y` (:185), `Y0` (:500), `mask` (:499), `R` (:253/:679), `tau_map` (free expanded **view** for CONSTANT/ROW/COLUMN, :599-608; full (N,P) materialization only in KNOWN mode, :586), plus `S` (N,P) in KNOWN mode → **≈ 4–5 × 4NP bytes persistent** (≈ 16–20 bytes/entry). Transient per-iteration peaks add ~5–8 more (N,P) temporaries (`Yfit`, `resid_mean_sq`, `first_moment_sq`, `second_moment`, `R2`, `recon`, `Rk`, `W`) → **peak ≈ 45–55 bytes/entry**. Factor tensors are negligible: 4·K·(N+P)·4 B. Example: N=5M, P=1k (5e9 entries) → ≈ 20 GB persistent, ≈ 250 GB peak on this implementation — dense-only, single-device; **this is the main blocker for 5M-cell scale** (would require sparse input, streaming/chunked residual updates, or out-of-core factors, none of which exist in the current source).
