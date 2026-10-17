# Amendment A4 — Normalization-sensitivity arms (PFlog, CLR) — PRE-RUN

Status: REGISTERED 2026-10-17, before any PFlog/CLR run. Parent: `04-experiment-plan.md`
(E-M3.1), as amended by A1 (DIVAS primary), A2 (SNR), A3 (Poisson generator +
ceiling-relative gates). No gate is changed by this amendment.

## Motivation

User request (2026-10-17): "Try pflog clr as alternative for log10 normalization
for singlecell or metacell ebmf." The registered recipe is the shifted logarithm
`y = log1p(x / (alpha * s_i))`, alpha=1, library-size size factors, followed by
gene center+scale (flashier single-cell vignette recipe; Ahlmann-Eltze & Huber
2023). PFlog and CLR are published alternatives from the depth-normalization
literature (Booeshaghi et al. 2022, bioRxiv 2022.05.06.490859; "PFlog" /
"shifted CLR" family), not present in cellprograms, flashier, or cebmf_torch.

## Registered transform definitions

All arms share the registered gene center+scale stage (fit on TRAIN donors,
replayed on test chunks), so arms differ ONLY in count-scale/centering geometry:

| Arm       | Transform (per cell row)                                        |
|-----------|------------------------------------------------------------------|
| shiftedlog (control) | `zscore_gene( log1p(x / s_c) )`, s_c = lib size / mean |
| pflog     | `zscore_gene( rowcenter( log1p(4 * ahat * x / s_c) ) )`          |
| clr       | `zscore_gene( rowcenter( log1p(x) ) )`                           |

- `ahat` (PFlog): pooled negative-binomial overdispersion by moment matching on
  depth-normalized TRAIN counts (`var = mu + ahat*mu^2`), median over expressed
  genes (mean > 1% of max mean), floored at 0.25 so `4*ahat >= 1` (Poisson-limit
  fallback = conventional pseudocount-1 scale). Recorded per run.
- rowcenter = subtract the cell's own row mean (CLR geometry; this is the depth
  mechanism for clr, which uses no size factors).
- Implementation: `m3_norm.py` (`PFLog`, `CLRLog`), drop-in compatible with
  `cebmf_ext.preproc.ShiftedLog` (`fit_transform` / `transform`), consumed
  unchanged by `project_scores`.

## Registered arms (additive; existing battery unchanged)

- Synthetic battery, S1 (single-cell) and S3 (metacell) strategies only,
  norm in {pflog, clr}, dropout rate in {0, 0.1, 0.3}, reps 1-10:
  2 x 2 x 3 x 10 = 120 additional runs. Control = the already-registered
  shiftedlog arms.
- Real-data sensitivity (secondary, pilot-grade): UCSF CD8 singlecell and
  metacell granularity discovery fits (workstream F recipe: top-2000 HVGs,
  K=10, laplace/laplace, 5 seeds) under pflog and clr, compared to the
  completed shiftedlog fits by K_eff, factor concordance with the control fit,
  and runtime. No acceptance gate; reported as sensitivity.

## Gates (UNCHANGED, A3 operationalization)

- A1 (broad recovery) and C2 (rare recovery) remain RELATIVE to the per-rep SVD
  ceiling computed on the SAME transformed train matrix (m3_norm arms recompute
  the ceiling on their own Y_train, so the comparison is ceiling-fair).
- C1, C3, C4 untouched. Normalization arms are SENSITIVITY analyses: they
  qualify or contextualize the S1/S3 strategy conclusions; a normalization that
  changes the S1-vs-S3 ranking is reported as a material caveat to C1/C2
  interpretation, not as a new acceptance bar.

## Known caveat (registered)

The A3 synthetic generator is Poisson by construction (`counts ~ Poisson(lambda)`),
so the NB moment-matched `ahat` is expected to hit the 0.25 floor on synthetic
data (PFlog's variance-stabilization premise — NB overdispersion — is absent
there). The PFlog synthetic arm therefore tests the CLR-centering geometry at
the pseudocount-1 scale, not the Anscombe calibration; the UCSF real-data arms
carry the calibration. This is inherent to the registered generator and is
reported alongside results.
