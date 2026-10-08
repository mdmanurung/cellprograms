# cellprograms Real-Data Benchmarks — Phase 1 Plan

**Date**: 2026-10-06
**Status**: Proposed (awaiting approval)
**Supersedes**: v2 implementation plan (completed; preserved at `PLAN_v2_completed.md`)

## Objective

Validate the cellprograms pipeline (per-cell-type EBMF + the L1/L2/L3-block borrowing ladder, all simulation-validated) on five real cohorts, against sanity/cell-type-aware baselines and distributional comparators, answering the v1 benchmarking plan's three questions:

1. **Sample-biology preservation** — do sample-level representations recover clinical/biological covariates?
2. **Program stability & interpretability** — are discovered programs reproducible and biologically coherent?
3. **Local-first correctness** — does the borrowing ladder recover shared/private structure without false sharing on real data?

## Locked decisions (user, 2026-10-06)

- **Comparators**: full suite — cheap baselines (random, composition, CLR-composition, global pseudobulk, grouped pseudobulk, per-CT PCA) **plus GloScope and PILOT**. MOFA-FLEX/mc-ASTRA remain excluded per the approved v2 plan (literature context only).
- **Atlas scale**: capped subsets — HLCA core (~584k cells), OneK1K capped at ~300 donors.
- **IBD dataset**: ~~Smillie et al. 2019 UC colon atlas (via CELLxGENE)~~ **DROPPED 2026-10-06 (user decision)** — the CELLxGENE deposit contains only 12 healthy donors' epithelial cells (no UC samples); Martin anti-TNF (11 donors, all normal) and Kong Crohn's (5 case donors) are also inadequate. Full Smillie atlas requires Broad SCP259 auth. Four cohorts this phase: COMBAT, Stephenson, HLCA-core, OneK1K-300.
- **Sequencing**: COMBAT end-to-end first as the template, then roll out to the other four datasets.

## Datasets and access routes (verified 2026-10-06)

All five cohorts are openly reachable; `cellxgene_census` 1.17.0 + scanpy 1.11.4 are installed in this environment.

| Dataset | Benchmark role | Access route | Scale |
|---|---|---|---|
| COMBAT (Cell 2022) | Primary case study (clinical: disease/source/severity/outcome/mortality; technical: institute/pool) | CELLxGENE Census, collection `8f126edf-5405-4731-8374-b5ce11f53e82` | 836,148 PBMC cells |
| Stephenson (Nat Med 2021) | Severity + trajectory, PBMC | CELLxGENE Census, collection `ddfad306-714d-4cc0-9985-d9072820c530` | ~781k cells / 130 donors |
| OneK1K (Science 2022) | Large-N scalability, many cell types | CELLxGENE Census, collection `dde06e0f-ab3b-46be-96a2-a8082383c4a1` | capped ~300 donors (~390k cells) |
| HLCA (Nat Med 2023) | Heterogeneous-tissue stress test (tissue/disease/smoking/assay) | CELLxGENE Census; core atlas; collection/dataset ID resolved programmatically at execution | ~584k cells (core) |
| Smillie UC colon (Cell 2019) | IBD biological comparison (epithelial + immune programs) | CELLxGENE Census; dataset ID resolved programmatically at execution | ~366k cells |

Correction to earlier notes: GSE158055 is the Ren et al. 2021 cohort, **not** COMBAT; the Census route is used regardless. OneK1K raw data are also in GEO/SRA (GSE196830) but Census is simpler.

## Execution flow (per dataset, COMBAT first)

1. **Acquire + pseudobulk** (Python, background job on a ≥32 GB worker): pin `CENSUS_VERSION`; query Census for the dataset; verify disease/tissue/donor labels before analysis; build donor × cell-type pseudobulk (skill `single-cell-census-query` scripts: `enumerate_labels.py`, `build_pseudobulk.py`); per-CT HVG selection capped at ~3,000 genes/CT (matches validated EV-BIDIFAC feasibility envelope); log-normalize; save pseudobulk matrices + sample metadata to `/mnt/shared-workspace/`.
2. **Stage-1 fits** (R): per-CT EBMF with flashier, identical settings to the validated simulations.
3. **Borrowing ladder** (R): L1 shared priors (global + groups modes), L2 initialization borrowing, L3 `block` decomposition (default backend), PVA + `match_programs`; evb-S scores-level cross-check as a lightweight confirmation.
4. **Cheap baselines** (R, existing env): random embedding, composition, CLR-composition, global pseudobulk, grouped pseudobulk, per-CT PCA.
5. **Distributional comparators**: GloScope (Bioconductor; R ≥4.4 compatible) and PILOT (`pip install pilotpy`, v2.0.15). Both consume single-cell embeddings, not pseudobulk → run with per-sample cell caps (≤2,000 cells/sample) to produce sample×sample distance matrices → MDS embedding → evaluated with the same Panel A/B metrics. Reported under fairness regime B (author-recommended input), since matched-input (regime A) does not apply.
6. **Panels A–E** per the v1 benchmarking plan: A sample biology (KNN macro-F1 for categorical covariates; Spearman for continuous), B technical retention (reported separately; **no single aggregate score**), C representation quality (program-level R², held-out reconstruction, subsample stability using the package's existing machinery), D interpretability (loading concentration, pathway coherence, cell-type localization, technical contamination), E practicality (runtime/memory logging).
7. **Cross-cohort reproducibility**: COMBAT ↔ Stephenson (both COVID-19 PBMC) program matching via loading similarity (`match_programs`).

## Compute estimates

- Pseudobulk streaming: ~65 min per ~1M cells on a 32 GB worker → COMBAT ~55 min, Stephenson ~50 min, HLCA-core ~40 min, Smillie ~25 min, OneK1K-300 ~25 min. Run as tracked background jobs (max 2 concurrent; keep foreground light to avoid BLAS contention).
- Stage-1 EBMF + ladder: minutes per dataset (validated on sims).
- GloScope/PILOT: heaviest comparator step (pairwise sample distances); estimate firmed up after the COMBAT pilot; cell caps keep it tractable.
- Expected total: roughly one day of wall-clock with backgrounding.

## Acceptance criteria (definition of done, this phase)

1. COMBAT end-to-end: QC'd pseudobulk, stage-1 + ladder fits, all baselines + GloScope + PILOT, panels A–E, per-dataset report.
2. Template rolled out to Stephenson, Smillie-IBD, HLCA-core, OneK1K-300, each with a per-dataset report.
3. Cross-dataset summary including COMBAT↔Stephenson program matching.
4. All deliverables under `/mnt/results/benchmarks/<dataset>/` (figures PNG); computational record in the execution-trace notebook.

## Risks and assumptions

- HLCA and Smillie Census IDs are resolved programmatically at execution (query collections by name); fallback is direct h5ad download from CELLxGENE Discover.
- GloScope/PILOT use single-cell input (regime B); comparisons are interpreted accordingly.
- OneK1K is healthy donors only → used for scalability/stability (Panels C/E), not disease-metric Panel A.
- Technical-covariate availability varies by dataset; Panel B is adapted per dataset.
- Pseudobulk is the only heavy step; if a dataset exceeds the 32 GB worker, provision a larger worker via ManageMachine rather than HPC.
