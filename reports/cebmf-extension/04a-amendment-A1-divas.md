# 04a — Amendment A1 (2026-10-17, pre-run): DIVAS COVID-19 multi-omics as the primary real-data example

User-suggested amendment to `04-experiment-plan.md`, registered BEFORE any M3 run
executes. No registered criterion (A1-A5, B1-B3, C1-C4, leakage rules B-L1..L3) is
altered; the DIVAS arms ADD evidence and give gate C3 a second, stronger reference.

## Rationale

The DIVAS COVID-19 case study (Su et al. 2020; data + precomputed DIVAS results public:
GitHub ByronSyun/DIVAS_COVID19_CaseStudy, Zenodo 17430294) is already wired into this
repo (`benchmarks/divas_comparison/`: download script, cellprograms fits at DIVAS
parity, published comparison REPORT). It offers: (a) a published reference
decomposition (DIVAS components, nsim=500) for concordance gating; (b) in-session
cellprograms/flashier reference fits (packages verified installed); (c) real sample-level
covariates (WHO ordinal score wos, severity, sex, age) — a far stronger approach-B test
bed than the negative UCSF age pilot; (d) 6 blocks (4 gene blocks, 8,634 shared genes +
proteomics 481 + metabolomics 763; 114 samples) exercising multi-block robustness.

## Structural note (binding)

DIVAS gene blocks are SAMPLE-level pseudobulks (114 samples); there is no cell-level
matrix to project onto. The projection-based strategies (S2 pseudobulk-init+project,
S3 metacell-init+project) therefore do NOT apply to DIVAS. They remain evaluated on:
(a) the synthetic battery (unchanged; gates C1/C2 require planted truth, which DIVAS
cannot provide), and (b) the UCSF granularities as SECONDARY real cell-level arms
(cebmf fits reusable from workstream F).

## DIVAS arms (new, registered)

Per block (CD4+ T, CD8+ T, CD14+ monocytes, NK, proteomics, metabolomics), 5 seeds:
- **D1 from-scratch**: cebmf fit at DIVAS parity — K=30, center=TRUE, scale=FALSE,
  laplace (point-Laplace) priors, var_type=1, backfitting + pruning, seed family
  20261016 — matching the registered cellprograms DIVAS-parity recipe
  (`max_factors=30, features="all", loading_prior="point_laplace"`).
- **D2 anchored refinement**: user-init passthrough from (i) published DIVAS block
  loadings and (ii) in-session cellprograms program loadings; refine (full EM, maxit 50,
  ELBO-plateau stop). Reports ELBO vs from-scratch and concordance of the refined fit to
  its anchor (realizes init+refine on real data; MR-4 overlap noted, pilot-grade).
- **D3 sensitivity**: registered operating point (K=10) on the CD8+ T block.
- **Metrics** (divas_comparison conventions): matched |Spearman| of scores across the
  114 samples vs DIVAS components and vs cellprograms programs; top-100 loading Jaccard
  for matched gene-block programs; wos (Spearman) and severity (Kruskal) associations,
  BH-adjusted; seed stability; K_eff; ZeroRate (C4); runtime + peak RSS.
- **Gate C3 (DIVAS analogue)**: cebmf-vs-DIVAS concordance >= cellprograms-vs-DIVAS
  concordance on the same blocks, same metrics (published REPORT.md numbers are the
  registered reference point; in-session refits confirm them).

## Data status and mitigation (2026-10-17)

The 8 input CSVs (6 blocks + sample metadata + summary) are downloaded and verified on
the shared workspace (`G_m3/divas/input/`, 44 MB). The precomputed DIVAS results RDS
(419 MB, Zenodo 17430294) is TEMPORARILY unavailable (HTTP 403, Zenodo rate limit
"unusual traffic"; GitHub mirror carries only a metadata stub). Mitigation: (i) the
synthetic battery and the in-session cellprograms reference fits start immediately (no
RDS dependency); (ii) D2 arm (i) (DIVAS-loadings init) and the cebmf-vs-DIVAS concordance
metric wait for a Zenodo retry; (iii) if Zenodo remains blocked, the C3 DIVAS analogue
falls back to the published REPORT.md concordance numbers as the registered reference,
with the RDS retry documented. No criterion is evaluated against missing data.

## E-M3.2 amendment (approach B)

PRIMARY covariate experiments move to DIVAS: wos (continuous) and severity (ordinal) as
moderation covariates on the 4 gene blocks, lcash prior on L with X_l = covariate,
evaluated on held-out PATIENTS. Binding split rule: 85 patients, T1/T2 paired — splits
are BY PATIENT, never by sample (B-L2; the divas_comparison pairing sensitivity is the
precedent). Arms: B-fit vs A-fit vs A + post-hoc regression (B3), 5 patient splits x 2
folds. The synthetic beta-grid power experiment (E-M3.2 as originally registered) is
RETAINED unchanged — it answers the smallest-detectable-beta question the real data
cannot. UCSF age arm demoted to secondary pilot (completed, negative at this operating
point).

## Compute impact

DIVAS blocks are small (114 x 8634 max): D1 5 seeds x 6 blocks ~ minutes-to-tens of
minutes; D2 ~ same; D3 trivial. E-M3.2 DIVAS arms ~ 1 h. Total plan estimate unchanged
(~9-13 h). Execution: worker-1 (R references + cebmf arms via staged adapter) and
worker-0/3 (synthetic battery).
