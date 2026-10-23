# Pre-registration: mr2-cnmf-v1 — cNMF as the stage-1 factorization backend

**Tag:** `mr2-cnmf-v1` · **Date:** 2026-10-22 · **Status:** FILED before any
benchmark dataset is simulated, fitted, or evaluated for this comparison.
Amends the MR-2 draft rule in `reports/STAGE_B_ROADMAP.md` (which is a draft
sketch, not a filing). Repo law: this file is committed before the run starts.

## 1. Question

Does cNMF (consensus NMF; Kotliar et al. 2019) replace the flashier EBMF backend
as the package's stage-1 factorization on the same pseudobulk inputs? The
registered Stage A operating point is the cmp-bench-v1.2 winner: per-cell-type
flashier fits at `max_factors = 30` (point-laplace loadings, var_type 1, 2000
HVGs, institute residualization) upstream of the frozen caller `n_agg_noanc`.

## 2. Arms

| Arm | Stage-1 | Caller |
|---|---|---|
| `stageA_K30` | flashier, the registered Stage A recipe (v1.2 winner settings, frozen in `PREREG_cmp-bench-v1.2.md`) | `n_agg_noanc` (unchanged) |
| `cnmf` | cNMF per cell type on the SAME processed matrices (same HVG panel, same institute residualization, same centering), K grid {5, 10, 20, 30}, K selected per cell type by cNMF's own consensus-stability criterion (largest K whose mean stability is within 0.05 of the grid maximum; ties → smaller K) | `n_agg_noanc` (unchanged) |

Both arms consume identical simulated datasets (same seeds, same planted truth);
the comparison is the stage-1 backend alone. cNMF scores are non-negative
(usage); the cNMF arm skips `canonicalize_programs` sign/scale alignment (no
signed factors exist) — recorded here as a structural difference, not a
deviation. The cNMF adapter emits the same fit-object fields the caller
consumes (`scores`, `loadings`, `cell_types`, `observation_ids`,
`sample_metadata`).

## 3. Data

The cmp-bench-v1.2 held-out harness: simulator `sim_cmp.R` (frozen), 8
scenarios, planted-CMP engine unchanged. **Deviations from the as-run v1.2
harness, all forced by the 2026-10-22 loss of the `/mnt/shared-workspace/cmp`
scratch volume, all disclosed:**

1. **Semi base rebuilt from source.** The original `semi_base.rds` (Stephenson
   COVID-19 PBMC pseudobulk) is lost and rebuilt from ArrayExpress
   **E-MTAB-10026** (the Stephenson et al. 2021 accession; the
   `mofacellular_comparison/README.md` reference to "E-MTAB-8731" is an
   accession error — that record is an unrelated mouse miRNA study). Rebuild:
   donor × gene count pseudobulk per cell type from
   `covid_portal_210320_with_raw.h5ad` + `annotation_200112.csv`; 12 major cell
   types (mapping frozen in the rebuild script `01_rebuild_semi_base.R`, which
   also reports per-type donor/cell counts); donors with < 50 cells of a type
   dropped for that type; per-type `Y` (counts) and `cpm_mean` profiles;
   `donor_meta` with `donor_id`, `site`, `sex`,
   `Status_on_day_collection_summary`. The rebuilt substrate differs from the
   lost original in unspecified preprocessing details of the in-session build;
   both arms consume the SAME rebuilt substrate, so the paired comparison is
   internally valid. The committed v1.2 numbers are a reference, not a pairing
   target.
2. **Gene-set resources reconstructed.** `res/hallmark.csv` (MSigDB Hallmark,
   official GMT download) and `res/lr_consensus.csv` (CellChatDB consensus
   ligand–receptor pairs, secreted + receptor–ligand interactions, reduced to
   `ligand`/`receptor` gene-symbol columns) are rebuilt from public sources and
   committed under `benchmarks/cmp_detection/res/`. The original in-session
   files are lost; the `lr` scenario's planted gene pairs therefore differ from
   the as-run v1.2. Both arms consume the same files.
3. **Harness paths adapted** to the repo layout (env-var overrides; no scoring
   or planting logic touched).

## 4. Design (reduced from v1.2 scale, registered)

- **Base: semi only.** The v1.2 result showed the parametric base barely
  differentiates arms (K30-vs-K10 paired diff +0.003 param vs +0.251 semi);
  the backend question concerns real-data structure, which the semi base
  carries.
- **Seeds 1051–1075** (25 of the 50 held-out seeds), 8 scenarios → 200 datasets
  per arm, 400 stage-1 fits total. Amplitude: semi 0.875 (v1.2 frozen).
- Rationale: the lost fit cache forces Stage A reruns (~25 min/dataset at K30);
  200 paired datasets give the paired-difference SE ≈ 0.01–0.015, adequate for
  the 2-point margin below. The remaining held-out seeds 1076–1100 stay
  untouched for a possible scale-up.
- Paired comparison: per-dataset F1 (marginal mode) on common datasets, paired
  bootstrap 2000 resamples (seed 1), as in `prereg_check_v12.R`.

## 5. Endpoint and decision rule

- **Endpoint:** pooled membership F1 over positive scenarios (marginal mode),
  as in v1.2 (`membership_F1_pooled`).
- **Adoption rule (from the MR-2 draft, made precise):** cNMF replaces the
  flashier backend as the registered stage-1 default iff
  1. paired F1 difference (cnmf − stageA_K30) ≥ **+0.02 absolute**, AND
  2. the paired bootstrap 95% CI excludes 0 (non-overlapping with the null),
  AND
  3. guardrails on the cNMF arm do not degrade beyond +0.05 absolute vs the
     Stage A arm on any of: null FWER, ambient decoy rate, within-context
     false-call rate, genes-only false-call rate (each ≤ 0.10 absolute
     regardless, inheriting the v1.2 guardrail levels).
- Otherwise the flashier backend remains the default and the result is ledgered
  as the MR-2 disposition either way.

## 6. Cost and execution

- Stage A: 200 fits × ~25 min ≈ 83 CPU-h; cNMF: 200 × ~2 min ≈ 7 CPU-h; caller
  eval + scoring ≈ 1 CPU-h. Total ≈ 91 CPU-h ≈ 6 h wall on 16 cores.
- Execution: sharded, per-dataset cache checkpoints (resumable), dedicated
  machine; the run starts only after this file and the rebuilt substrate are in
  place.

## Addendum (2026-10-22, before the run; no rule changes)

1. **Rebuild script language.** The registered "01_rebuild_semi_base.R" is
   implemented as `scripts/01_rebuild_semi_base.py` (h5ad access is python-side)
   + `scripts/02_assemble_semi_base.R` (assembles semi_base.rds). The frozen
   12-major mapping lives in the python script verbatim; B_malignant excluded;
   ASDC assigned to DC. Naming deviation only; the registered recipe is
   unchanged.
2. **cpm_mean operationalized** as the type-level mean CPM profile (gene
   vector, size factor 1e6), stored per type alongside Y.
3. **Stability rule operationalized** (cnmf_fit.py): stability_k = 1 -
   local_density_k / max(local_density over the K grid) from cNMF's
   k-selection output; chosen K = largest k with stability_k >=
   max(stability) - 0.05 (ties -> smaller k), per the registered rule.
4. **cNMF run settings**: seed 14, t=1, num_iter=100, consensus
   density_threshold 0.2 / local_density_threshold 0.05 (cNMF defaults).

## Addendum A (2026-10-22) — cNMF input resolution

The prereg §2 specified cNMF on the "SAME processed matrices" as K30, but those
matrices are continuous (planted signal is centered; min < 0) while cNMF requires
non-negative integer counts. **Resolved:** the cNMF arm consumes the planted
donor × gene matrix rounded + clipped to ≥ 0 (its native count input); cNMF does
its own HVG selection and TPM normalization internally. K30 consumes the
centered/institute-residualized matrix as registered. Each backend is thus tested
on its natural input — the methodologically fair comparison. Recorded as a
registered deviation from the §2 "same matrices" wording, not a silent change.
Also: cnmf 1.7.1 reports stability as consensus *silhouette* (higher = stabler),
replacing the older local_density; the §2 K-selection rule (largest K within
0.05 of max stability, ties → smaller K) is applied to silhouette unchanged.

## Addendum B (2026-10-22) — Y scale correction in the rebuilt substrate

The first rebuilt semi_base stored RAW COUNTS as Y; the sim_cmp.R contract is
Y = log2(CPM + 1) (the ambient scenario inverts it via 2^Y - 1, and the
planting amplitude msd is a log-scale sd). The raw-count base made planted
signals vanishingly small in relative terms and crashed the ambient transform
(2^counts = Inf). Fixed in the rebuild script (Y.mtx.gz now stores
log2(CPM+1)); all fits from the raw-count base are DISCARDED. The cNMF input
(Addendum A) is amended accordingly: cNMF consumes round(max(2^Y - 1, 0)) —
CPM-scale pseudo-counts recovered from the planted log-space matrix — and
runs its own normalization internally. Also: HSC (15 donors) and RBC (17) are
dropped from the rebuilt base (MIN_DONORS = 32): the C4 residual-df guard
(added 2026-10-16, after the v1.2 run) aborts any type with n < K + 2 = 32 at
the registered max_factors = 30.
