# cellprograms — Package Development Plan v2

## Summary

`cellprograms` is a lightweight R package for discovering cell-type-specific transcriptional programs from pseudobulk single-cell RNA-seq data using empirical-Bayes matrix factorization (`flashier`), with a principled ladder of mechanisms for borrowing information across cell types. Version 2 of the plan makes three changes relative to the original scope document:

1. **MOFA-FLEX is removed entirely** — no backend, no benchmark comparator. The joint-method role is filled by an internal Level-3 race.
2. **Level 3 borrowing is an empirical race between two candidates**, each run at two levels (expression and program scores): a block-structured joint EBMF built from `flashier` primitives versus an EV-BIDIFAC-style linked decomposition. The benchmark decides which becomes the package's structured-integration backend.
3. **Level 4 borrowing is a fully developed principal-angle framework** for aligning and comparing program subspaces across cell types, cohorts, and bootstrap refits — diagnostics, matching, null calibration, and a lightweight stage-2 integration method.

The core philosophy is unchanged: **factorize locally first; integrate later only when scientifically useful.** Sharing is an output, not an assumption. Every downstream quantity remains traceable: systemic state → cell-type program → gene loadings → original pseudobulk matrix.

---

## 1. Scientific objective and design principles

After pseudobulking, the data are matrices \(Y_c \in \mathbb{R}^{N_c \times G_c}\) with cell type \(c\), sample-timepoint rows, gene columns. Stage 1 discovers programs independently per cell type:

\[
Y_c \approx Z_c W_c^\top + E_c
\]

and the combined sample × cellular-program matrix \(Z_{\mathrm{all}} = [Z_B \mid Z_{CD4} \mid \cdots]\) remains the central downstream object.

Design principles carried from v1, restated with the borrowing ladder:

- **P1 — Orchestration, not invention.** Compose published machinery (`flashier`/`ebnm`, ICA, BIDIFAC-family, principal-angle linear algebra). No custom variational inference or novel factor model — with one gated exception (EV-BIDIFAC reimplementation, §4.3).
- **P2 — Cell-type autonomy at discovery.** Never concatenate expression across cell types before stage-1 factorization.
- **P3 — Sharing is an output.** Cross-cell-type coordination is measured (L4) or modeled (L1–L3) after, or carefully around, local discovery.
- **P4 — First-stage programs are never replaced.** Integrated states are optional summaries layered on top of `B_F1`, `Mono_F3`, etc.
- **P5 — Full traceability** from any integrated state back to genes and pseudobulk.
- **P6 — Unsupervised discovery.** Outcomes never influence factorization, feature selection, priors, or alignment.
- **P7 — No hidden imputation.** Missing sample × cell-type coverage is carried explicitly.

### The borrowing ladder

| Level | Mechanism | What is shared | False-sharing risk |
|---|---|---|---|
| L1 | Shared EBNM prior hyperparameters | Sparsity/scale of loadings | Minimal |
| L2 | Initialization borrowing | Where to look (donor loadings) | Low (shrink-to-zero) |
| L3 | Structured joint decomposition | Factor identity (global/partial/private) | Real — must be benchmarked |
| L4 | Principal-angle alignment | Nothing in the fits; post-hoc geometry | None in fits; confounds in interpretation |

L1–L2 trade almost nothing for efficiency. L3 trades false-splitting for false-sharing — the trade the benchmark must quantify. L4 changes no fits at all.

---

## 2. Input data contract and preprocessing

Unchanged from v1, condensed:

- Canonical input: named list of pseudobulk matrices (rows = `observation_id`, typically participant × timepoint; columns = genes), plus sample metadata (`observation_id`, `participant_id`, `cohort`, `timepoint`, `treatment`, `batch`, `cell_type`, …). Outcomes live in metadata only.
- Expect approximately Gaussian continuous expression (log-CPM or equivalent). Default `center = TRUE`, `scale = FALSE`, per cell type.
- Feature selection is per cell type: `"all"`, `"variable"`, or a named gene list. Never outcome-guided.
- Missing cell types: each cell-type fit uses available observations; scores align by `observation_id`; missingness is never silently zero-filled. Backend-specific rules (ICA requires completeness or explicit user-requested imputation; EV-BIDIFAC uses its native missing-data machinery).

---

## 3. Stage 1 — per-cell-type EBMF (core)

Unchanged from v1:

```r
flashier::flash(
    Y,
    ebnm_fn = list(ebnm::ebnm_normal, ebnm::ebnm_point_normal),
    var_type = 2,
    greedy_Kmax = 30,
    backfit = TRUE,
    nullcheck = TRUE
)
```

- Sample-side factors: Normal prior. Gene-side loadings: sparse signed point-Normal (sensitivity: point-Laplace, unimodal adaptive shrinkage).
- Rank: `max_factors = 30` permits up to 30; report candidates attempted / retained after greedy / after backfit / after nullcheck.
- Canonicalization: unit-RMS sample scores with scale absorbed into loadings; sign fixed by largest absolute loading (deterministic tie-break); verify reconstruction invariance.
- Stability: bootstrap over biological samples; match refit factors to the original (now via L4 subspace-aware matching, §4.4); report recovery frequency, score correlation, loading cosine similarity, sign consistency, top-gene recovery.

---

## 4. Borrowing architecture

### 4.1 Level 1 — shared EBNM priors

**Idea.** Keep separate factorizations; pool the *prior hyperparameters* of the gene-loading EBNM across cell types (sparsity scale, slab variance). Each cell type keeps its own programs, rank, and fit.

**API.**

```r
fit <- fit_celltype_programs(
    x,
    share_prior = c("none", "global", "groups"),
    prior_groups = NULL,        # e.g. lineage map: list(T = c("CD4","CD8"), ...)
    prior_max_iter = 5,
    prior_tol = 1e-3
)
```

**Algorithm.** Outer empirical-Bayes loop: (1) fit all cell types independently; (2) extract each fit's fitted loading prior `g_c`; (3) compute a pooled prior (global, or within `prior_groups`) by pooling the loading estimates across cell types and re-estimating the EBNM; (4) refit each cell type with the pooled prior as initialization; (5) iterate until prior parameters move less than `prior_tol` or `prior_max_iter` is hit.

**Diagnostics.** Per-cell-type vs pooled prior parameter table; convergence trace; warning when a cell type's individually estimated prior deviates strongly from the pooled one (masking genuine biology) — in that case recommend `share_prior = "groups"` or `"none"`.

**Failure modes.** Pooling across genuinely different sparsity regimes (e.g., a near-noise cell type with a program-rich one) can over-sparsify the rich type. Mitigation: the deviation warning, and the benchmark's N-sweep to show where pooling pays off.

### 4.2 Level 2 — initialization borrowing

**Idea.** Use programs from data-rich cell types (or a reference fit) to initialize fits in data-poor cell types. `flashier` supports initializing from supplied loadings; borrowed loadings are free to shrink to zero during backfit, so the borrowing is a suggestion, not a constraint.

**API.**

```r
fit <- fit_celltype_programs(
    x,
    init_from = list(NK = donor_spec(fit_ref, celltype = "CD8", programs = c("F001","F004")))
)
```

**Rules.** Donor and target loadings are matched on the union gene representation (§4.4, gene-universe decision); genes absent from the donor contribute zero. Borrowed factors enter as additional greedy initializations, never as fixed structure, unless the user explicitly requests `fix_borrowed = TRUE` (reference-projection mode).

**Diagnostics.** Borrowed-loading survival rate (fraction of borrowed programs retained after backfit/nullcheck); correlation between borrowed initialization and final loading; warning when survival is ~0 (donor inappropriate) or ~1 with near-zero target-specific refitting (target data uninformative — flag rather than celebrate).

### 4.3 Level 3 — structured joint decomposition: the race

Two candidates, each run at **two levels**:

- **Expression level (E):** input is the column-linked stack of cell-type expression matrices (shared sample rows, cell-type column blocks).
- **Score level (S):** input is \(Z_{\mathrm{all}}\), the concatenated stage-1 program-score matrix (column blocks = cell types). Blocks are scaled to unit Frobenius norm before decomposition so no cell type dominates by scale.

#### Candidate A — block-structured joint EBMF (`"block"`)

Greedy two-pass strategy built from `flashier` primitives — no new inference:

1. **Shared pass.** Fit EBMF on the stacked matrix (expression) or on \(Z_{\mathrm{all}}\) (scores). Factors whose loadings concentrate in one block are reassigned as private candidates; the rest are shared structure.
2. **Private pass.** For each cell type, backfit EBMF on its residual after removing the shared reconstruction, recovering private programs the shared pass diluted.
3. **Order ablation.** Also implement the reverse order (private-then-shared: fit per-cell-type first, then fit shared structure on the concatenated *residual* activities). Greedy sequential decompositions are order-sensitive; the benchmark must report both.

Output is annotated factor-by-factor as private-to-c / shared-across-subset / global, using block-wise loading concentration (e.g., max block mass vs entropy over blocks).

#### Candidate B — EV-BIDIFAC-style decomposition (`"ev_bidifac"`)

EV-BIDIFAC (Empirical Bayes Linked Matrix Decomposition; arXiv:2408.00237) is an empirical variational Bayes linked factorization that decomposes linked matrices into modules shared across arbitrary row/column subsets — globally shared, partially shared, and specific — is tuning-parameter-free, and natively handles entry-wise and block-wise missingness via EM-style imputation. This maps directly onto the package's private / partially shared / broadly coordinated taxonomy.

**Backend gate (Phase B0).** The Lock lab distributes BIDIFAC/BIDIFAC+ as loose R scripts on GitHub (`lockEF/bidifac`), not a CRAN package, and EV-BIDIFAC code availability/licensing must be verified. Decision gate:

1. If EV-BIDIFAC code is available under a compatible license → thin wrapper backend.
2. Else → implement the published variational algorithm in-package from the manuscript (it is a fully specified, tuning-free mean-field VB scheme; this is a *port*, not a new model — the single sanctioned exception to P1), and validate against the manuscript's published simulation settings and the authors' simulation gallery before any use.

**Unidimensional vs bidimensional.** The base design links matrices through shared samples (columns) — the unidimensional case. The bidimensional capability (e.g., cell types × cohorts simultaneously) is noted as a v2+ extension, not part of the race.

#### Race design and decision criteria

Four configurations: `block-E`, `block-S`, `evb-E`, `evb-S`. Evaluated on the simulation suite (§9) with pre-registered criteria:

| Criterion | Metric | Scenarios |
|---|---|---|
| Shared-state recovery | score correlation to true shared activity | 2, 4, 6 |
| False sharing | rate of merging independent programs | 1, 5, 7 |
| Private recovery | private-program score/loading recovery | 1, 4 |
| Partial sharing | subset-membership precision/recall | 3, 4 |
| Reconstruction | held-out RMSE | all |
| Stability | bootstrap recovery of module structure | all |
| Practicality | runtime, memory, failure rate | all |

**Decision rule.** The winner (per criterion profile, not a single aggregate) becomes the default structured backend of `integrate_programs()`. The loser remains shipped as a non-default method through v1.0 so the benchmark is reproducible. If the race is a wash, `ev_bidifac`-on-scores wins by tie-break: it preserves local-first stage 1 and is computationally trivial.

### 4.4 Level 4 — principal-angle framework

**Math.** For cell types \(a, b\) with loading matrices \(W_a, W_b\) (genes × programs): orthonormalize to subspace bases \(U_a, U_b\) (QR); form the cross-Gram \(M = U_a^\top U_b\); SVD \(M = P \Sigma Q^\top\). Singular values \(\sigma_i = \cos\theta_i\) are the principal angles; columns of \(U_a P\), \(U_b Q\) are the principal vectors — the maximally aligned direction pairs. The spectrum is invariant to rotations within either subspace, which is exactly the invariance factor matching needs. (Same machinery as SpaGE's cross-dataset alignment; classical Björck–Golub principal angles.)

**Gene-universe decision (delegated, decided).** Loadings are represented on the **union** gene set with zero-padding. Key fact: padded coordinates contribute zero to every inner product, so the cross-Gram — and therefore all angles — is *identical* to intersect-only computation. Union representation is therefore strictly preferable: identical geometry, plus full-length principal vectors for gene traceability. What actually matters is the **effective dimension** of the comparison = size of the gene intersection, which (i) is reported per pair as a diagnostic with a warning below a threshold (default < 200 genes), and (ii) sets the ambient dimension of the null calibration. Disjoint feature sets correctly yield all angles = 90° (no measurable sharing).

**Null calibration.** Sharing is only claimed against an explicit null:

- *Asymptotic null:* for random \(k_a, k_b\)-dimensional subspaces in ambient dimension \(G_{\cap}\), squared cosines follow the Jacobi/MANOVA ensemble; largest-angle approximation \(\cos\theta_1 \approx (\sqrt{k_a}+\sqrt{k_b})/\sqrt{G_{\cap}}\).
- *Permutation null:* shuffle gene labels of one cell type's loadings (destroys gene-wise correspondence, preserves sparsity/scale), recompute spectrum; report empirical percentiles per angle. Default 200 permutations.

A pair's sharing summary: number of angles exceeding the 95th percentile of the permutation null, plus the full spectrum.

**Functions.**

```r
principal_angles(fit, ct_a, ct_b, weight = c("none", "variance"), n_perm = 200)
sharing_spectrum(fit)                    # all pairs → matrix of spectra + summaries
plot_sharing_spectrum(fit)               # heatmap of shared-dimension counts / angle spectra
align_programs(fit, reference = "CD8")   # principal vectors + aligned coordinates
match_programs(fit_a, fit_b, method = c("subspace", "correlation"))
```

`weight = "variance"` scales each subspace direction by explained variance before the Gram (sensitivity option; default `"none"`).

**Four uses.**

1. **Sharing diagnostic (primary).** Pairwise sharing spectra across cell types — a graded, assumption-free map of cross-cell-type structure. Becomes a standard output alongside `plot_program_correlations()`.
2. **Matching upgrade.** `match_programs()` gains `method = "subspace"`: principal-vector rotation initializes Hungarian matching on projected loadings. Used for bootstrap stability matching and cross-cohort reproducibility. Subspace agreement is reported alongside factor-level matches, so false splitting shows up as subspace agreement with failed one-to-one matches rather than as instability.
3. **Lightweight stage-2 integration.** `integrate_programs(fit, method = "pva")`: project each cell type's programs onto shared principal vectors → multicellular axes with explicit gene traceability via \(g_c = \sum_j a_j w_j\). A near-free competitor to ICA in the stage-2 slot.
4. **Cross-cohort reproducibility metric.** Angle spectra between discovery and validation refits complement loading-similarity matching.

**Caveats encoded as diagnostics, not hidden.**

- *Correlated-but-distinct confound* (Scenario 5): small angles do not prove shared biology. The sharing spectrum is always reported next to the permutation null, and documentation states plainly that PVA measures geometric overlap, not identity.
- *Rank sensitivity:* trailing noise factors inflate angles. Report spectra as a function of included rank \(k\) (sensitivity curve) when a cell type retains many factors.
- *Sparsity loss:* principal vectors are dense combinations of sparse loadings. Gene-level interpretation of aligned axes goes back through the sparse \(w_j\) combination weights, never through thresholding the dense vector directly.

---

## 5. Combined representation and integration API

\(Z_{\mathrm{all}}\) remains the central export. Integration methods:

```r
integrate_programs(fit, method = c("none", "ica", "pva", "block", "ev_bidifac"),
                   level = c("scores", "expression"))   # level only for block/ev_bidifac
```

- `"ica"`: multi-start ICA on \(Z_{\mathrm{all}}\) with component matching and stability (unchanged from v1).
- `"pva"`: principal-vector alignment (§4.4).
- `"block"` / `"ev_bidifac"`: L3 candidates; `method` default set to the race winner after benchmarking.
- `trace_state(fit, integ, state)` maps any integrated state back to per-cell-type gene summaries \(g_c = \sum_{j \in c} a_j w_j\).

---

## 6. Diagnostics and QC

Carried from v1 (per cell type: N, G, missingness, factors retained, variance reconstructed, residual variance, sparsity, max factor correlation, leverage, stability, technical-covariate associations; warnings for small N, sparse pseudobulk, zero-variance genes, duplicate IDs, single-observation-dominated factors, library-size association, unstable programs, rank saturation) **plus**:

- L1: prior deviation table/warnings.
- L2: borrowed-loading survival diagnostics.
- L3: module classification table (private/partial/global ranks), variance decomposition by module class, block-scale report, convergence.
- L4: intersection sizes, angle spectra + null percentiles, rank-sensitivity curves.

---

## 7. Object structure and provenance

```r
cell_program_fit <- list(
    input = ..., preprocessing = ...,
    fits = list(B = flash_object, ...),
    programs = ..., scores = ..., loadings = ...,
    canonicalization = ..., stability = ...,
    borrowing = list(
        level1 = list(share_prior, pooled_priors, convergence),
        level2 = list(init_map, survival)
    ),
    alignment = list(            # L4
        spectra = ..., nulls = ..., matches = ...
    ),
    integration = NULL,          # ica / pva / block / ev_bidifac fit + metadata
    provenance = list(package_versions, backend_versions, seed, timestamp, parameters)
)
```

Standard export gains `sharing_spectrum.tsv`, `alignment/`, and `borrowing.tsv` alongside the v1 layout.

---

## 8. Public API summary

Core (v1, unchanged): `as_cell_program_data`, `validate_cell_program_data`, `fit_celltype_programs`, `canonicalize_programs`, `assess_program_stability`, `program_scores`, `program_loadings`, `top_genes`, `plot_program*`.

New/changed:

- `fit_celltype_programs(..., share_prior, prior_groups, init_from)` — L1/L2
- `principal_angles`, `sharing_spectrum`, `plot_sharing_spectrum`, `align_programs` — L4
- `match_programs(..., method = c("subspace", "correlation"))` — upgraded
- `integrate_programs(..., method = c("none","ica","pva","block","ev_bidifac"), level)` — extended
- `trace_state`, `integration_summary`, `program_summary` — extended for module taxonomy

---

## 9. Simulation and benchmark integration

The v1 simulation suite (Scenarios 1–7) now maps onto the borrowing levels:

| Scenario | Stress-tests |
|---|---|
| 1 — purely private | L3/L4 false-sharing rate; L4 null calibration |
| 2 — one global program | L3 shared recovery; L4 should show near-zero angles everywhere |
| 3 — lineage-restricted sharing | L3 partial-sharing membership; L4 spectrum shape |
| 4 — mixed private/shared/global (**principal**) | everything |
| 5 — correlated but distinct | L3 merging risk; L4 interpretation confound |
| 6 — same state, different genes | L3-S and PVA reconnection of stage-1 factors |
| 7 — same genes, independent activities (**negative control**) | L3/L4 must not infer coordination |

New components:

- **L4 null calibration study.** Under Scenario 1, empirical angle spectra must sit inside the permutation-null envelope; under Scenario 7, sharing claims must stay at null level despite overlapping gene identities.
- **L3 race** (§4.3): four configurations × scenarios × N-sweep (N = 20…1000) × missingness sweep. Decision criteria pre-registered as above.
- **Ablation ladder** (extends v1 §16): grouped pseudobulk → per-CT PCA → per-CT EBMF → +L1 → +L2 → +stability filtering → +L3 winner → +L4/PVA integration. Each rung isolates one source of gain.
- **Real-data benchmarks** (COMBAT, Stephenson COVID-19, OneK1K, HLCA, IBD) proceed per the v1 benchmarking plan with one change: MOFA-FLEX/mc-ASTRA is removed as a comparator; the joint-method evidence comes from the L3 race, with published mc-ASTRA results cited as literature context rather than re-run.

---

## 10. Implementation phases and roadmap

| Phase | Content | Version |
|---|---|---|
| B0 | Feasibility gates: EV-BIDIFAC code/license audit; principal-angle prototype + union≡intersection unit proof | — |
| B1 | L4 core: `principal_angles`, `sharing_spectrum`, permutation nulls, subspace `match_programs`; wire into stability + cross-cohort matching | v0.2–v0.3 |
| B2 | PVA as stage-2 method (`integrate_programs(method="pva")`) + ICA (as v1) | v0.3 |
| B3 | L1 shared priors + L2 initialization borrowing | v0.4 |
| B4 | L3 candidates: `block` (both orders) and `ev_bidifac` (wrapper or port), both levels | v0.5 |
| B5 | Simulation suite + race + null calibration → backend decision | v0.5–v0.6 |
| B6 | Real-data benchmark integration, vignettes, docs | v0.6–v1.0 |

Rationale for ordering: L4 is pure R linear algebra with no dependencies and immediately strengthens existing stability/matching machinery; L1/L2 are small `flashier`-level changes; L3 is the largest and is gated on B0 and on the simulation harness built in B1–B3.

---

## 11. Testing and acceptance criteria

**Unit tests.**

- Principal angles are invariant to orthogonal rotation of either fit's factors and to factor sign/scale canonicalization.
- Union zero-padded vs intersect-only computation give identical angles (the §4.4 equivalence, tested explicitly).
- Disjoint feature sets → all angles 90°.
- L1: pooled-prior refit changes reconstruction R² by < tolerance on simulated data; prior iteration converges.
- L2: borrowed initialization with irrelevant donor shrinks to zero (survival diagnostic fires).
- Canonicalization preserves reconstruction; all ID alignment exact; no silent imputation (v1 carried).

**Simulation acceptance.**

- Scenario 1: L4 sharing claims at null level; L3 false-sharing rate reported and below pre-set tolerance for the winner.
- Scenario 2: L3 recovers the global factor (score correlation > 0.9 at adequate N); L4 spectra near zero angle across all pairs.
- Scenario 6: stage-1 + PVA/EV-BIDIFAC-on-scores reconnects the shared state.
- Scenario 7: no method infers coordination from gene overlap alone.
- Race decision criteria (§4.3) computed and reported per criterion, no single aggregate.

---

## 12. Compute and execution targets

- **Package development, unit tests, simulation suite:** sandbox default machine (16 CPU / 64 GB) is ample. Simulation matrices are small (N ≤ 1000, G ≤ 5000, C ≤ 8 cell types); the full race grid is embarrassingly parallel across replicates and will be chunked with checkpoints to shared storage. Estimated hours, not days.
- **EV-BIDIFAC at expression level:** cost driven by SVDs of the stacked matrix; at benchmark scale (N ~ 500, ~3000 genes/cell type, 5–10 cell types) this is sandbox-feasible; HLCA-scale runs may need gene filtering or a larger worker.
- **Real-data pseudobulk construction** from large atlases (HLCA ~2.4M cells) is the only heavy step and belongs to the later benchmark phase: larger worker or HPC when reached; not needed for B0–B5.

---

## 13. Explicit assumptions and open questions

**Assumptions (decided).**

1. Gene universe for L4: union representation with zero-padding; angles identical to intersection-only; intersection size reported and used for null calibration. (Delegated by user; decided here with the equivalence argument.)
2. MOFA-FLEX removed everywhere, including benchmarks (user decision); joint-method evidence rests on the L3 race plus literature-reported mc-ASTRA results.
3. Block-structured candidate = greedy shared-then-residual (and reverse-order ablation) using `flashier` primitives; a simultaneous constrained EBMF is excluded as custom inference (P1).
4. The L3 loser remains shipped as a non-default method through v1.0 for benchmark reproducibility.
5. PVA is a first-class stage-2 integration method alongside ICA.
6. Tie-break if the race is inconclusive: `ev_bidifac` on scores, for local-first fidelity and computational cost.

**Open questions (gates, not blockers).**

1. EV-BIDIFAC code availability and license → B0 gate; fallback is the validated in-package port.
2. Bidimensional extension (cell types × cohorts) — scoped out of the race; candidate for v2+.
3. Whether L1 prior pooling should propagate prior uncertainty into score uncertainty — deferred; v1 treats scores as estimated summaries (per v1 §24).
