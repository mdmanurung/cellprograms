# Pre-registration: which model modifications earn their complexity

Registered by git tag `prereg-v1`. This file, the harness (`sim_param.R`, `sim_semi.R`, `arms.R`, `run_rep.R`, `run.sbatch`, `make_tasks.sh`) and `summarize.R` are frozen at the tag. Any later change goes into `DEVIATIONS.md` with a reason. Analyses not listed here are labelled exploratory.

## Question
Does any of the candidate modifications below improve program recovery or sharing-test performance over the current defaults, by more than the margin stated for it, without making calibration, recovery or runtime worse? A candidate that does not clear its rule is **not** ported into `R/`; the current default stays.

## Arms
All arms: `covariates = "institute"`, `covariate_mode = "residualize"`, `point_laplace`, `var_type = 1`, 2000 variable genes, `max_factors = 10`, seed 1, scores-space `sharing_spectrum` with `n_perm = 999` (floor 0.001, below the BH level of one true pair among 6), unless the arm states otherwise.

| Arm | Change vs `base` |
|---|---|
| `base` | none |
| `k_n` | K_ct = min(10, floor(n_ct / 10)), fit per cell type |
| `k_20` | `max_factors = 20` |
| `k_stab` | `base` fit, then keep programs with `recovery_freq >= 0.7` over 10 subsamples of 80% of donors (`assess_program_stability`) |
| `genes_after` | residualize on institute first, then take the top 2000 genes by residual variance |
| `vt_12` | `var_type = c(1, 2)` |
| `w_S` | `flash_control = list(S = sqrt(c / n_cells_i))`, `var_type = 2`; per cell type c is set so median S^2 = 0.5 x median per-gene residual variance of the `base` fit |
| `strata` | `base` fit; the sharing null permutes within institute (`strata=`) |
| `none`, `none_strata` | no covariates, without / with `strata=`. Diagnostic only (confounder check); never adopted. |

## Data
Four cell types A, B, C, D, six pairs. Institute is the only covariate and the only batch.

**`param` (parametric).** 4000 genes per cell type (so that gene selection of 2000 is non-trivial). A, B: 100 donors; C: d001-d050; D: d021-d070 (C and D share 30 donors). Institute has 3 levels (0.6 / 0.25 / 0.15), batch loadings are shared by all cell types (as in `08_sim_adjustment.R`). Each cell type has 2 private programs (50 genes, loading sd 1). n_cells is lognormal (median 500, sd_log 0.8) and per-donor noise sd is sqrt(median(n)/n_i), exactly the `w_S` model, so `w_S` is advantaged on this base and agreement in sign there is not evidence for it.

**`semi` (semi-synthetic on COMBAT).** A = CD4 (121 donors), B = cMono (120), C = GDT (50), D = DP (57); C and D share about 31 donors. Per replicate, donors are permuted within Institute independently per cell type (rows keep their real n_cells). This removes real cross-cell-type biology and keeps Institute structure, noise and heteroscedasticity. Programs are planted on 50 genes drawn from the full matrix before any arm selects genes, with loading sd = `amp` x median gene SD (log-CPM).

## Scenarios and replicates
| ID | Truth | Reps |
|---|---|---|
| S0 | null (private programs / real structure, batch) | 200 |
| S1 | shared z, pair A-B (large-large) | 50 |
| S2 | shared z, pair A-C (large-small) | 50 |
| S3 | shared z, pair C-D (small-small, about 30 shared donors) | 100 |
| S7 | same genes in A and B, independent z (no sharing) | 50 |

Every arm runs on the identical dataset in each replicate (paired). Seeds: `base_offset(0 | 5e6) + 1e5 * scenario_index + rep`; pilot reps 9001-9020; confirmation reps 10001+.

**Amplitude.** Chosen per base and scenario (S1, S2, S3; S7 uses the S1 value) in a pilot (reps 9001-9020, results discarded) over `amp` in {0.1, 0.15, 0.25, 0.5, 1, 2}: the value where `base` power is closest to 0.5 and inside [0.3, 0.7]; if none is inside, the closest. The chosen values are stored in `amplitudes.csv` and listed in `DEVIATIONS.md`.

## Metrics
- **R (recovery).** For each cell type in the true pair, max |cor(true z, score column)| minus that arm's chance level (mean over all replicates of max |cor| between 50 random z and the arm's score columns, per cell type). R is the mean over the two cell types.
- **P (power).** The true pair has `shared = TRUE` (BH q < 0.05 over callable pairs).
- **F (false positive).** Any of the 6 pairs has `shared = TRUE`, in S0 or S7. False calls on non-true pairs in S1-S3 (`F_other`) are reported, not gated.
- **Time.** Median wall seconds per dataset (for `k_stab`, including the base fit).
- Secondary, reported only: callability, `excess_frac`, `z_pair`. No AUROC.

## Decision rules
Effects are paired over replicates; 95% CIs by paired bootstrap (4000 resamples, `summarize.R`). "Pooled" means replicates of both listed scenarios in one vector.

| Arm(s) | Primary endpoint | Effect passes if |
|---|---|---|
| `k_n`, `k_20`, `k_stab` | dP on S3 | estimate >= +0.10 and CI lower bound > 0 |
| `genes_after`, `vt_12`, `w_S` | dR pooled over S1 and S2 | estimate >= +0.02 and CI lower bound > 0 |
| `strata` | dF on S0 | `base` F on S0 > 0.05, the CI upper bound of dF < 0, and mean P cost over S1-S3 <= 0.05 |

**Guardrails** (all must hold, on both bases):
1. F(arm) <= F(`base`) + 0.03 in each of S0 and S7.
2. dR pooled over S1 and S2 >= -0.02 (not applied to `strata`; it does not change the fit).
3. Median time <= 3 x `base`.

**Adopt** an arm iff its effect passes on `semi`, its estimate has the same sign on `param` (for `strata`: dF <= 0), and the guardrails hold on both bases. Ties among adopted arms: K rules prefer `k_n` > `k_20` > `k_stab`; weights prefer `vt_12` > `w_S`. `none` and `none_strata` are never adopted.

**Reported, not gated.** Absolute calibration of `base` (F on S0 and S7, Clopper-Pearson 95% CI against 0.05; the earlier hard-adjusted simulation rejected 0.15). The confounder check: `none` must show F inflated on S0, otherwise the null contains no usable batch confounder and the strata result is uninformative.

## Detectable effects
With 100 paired replicates, a paired difference in a binary endpoint has SE <= sqrt(discordance/100). At discordance 0.4 (worst case) the SE is 0.063 and the 80%-power minimum detectable effect about 0.18; the +0.10 margin is detectable only if arms disagree on at most about 13% of datasets. A rejection therefore reads "no gain above the detectable effect", not "no effect". For R the MDE depends on the replicate-to-replicate SD of dR and is reported by `summarize.R` as the bootstrap CI width.

## After the main run
If two or more arms are adopted, the combination is tested against `base` on confirmation seeds (10001+) under the same rules before anything is ported. Ported arguments get one test each; defaults flip only for confirmed winners.
