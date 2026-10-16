# 01 — Scientific specification-to-implementation matrix (Phase 1)

Reconstruction of the methodology cellprograms implements, mapped to primary references, with deviation classification. Runtime verification evidence is collected by workstreams A-D (see TASKS.md); this document states what the implementation *is* and what each piece *claims*, and classifies the relationship to the underlying method.

## 1. References (verified this audit)

Core machinery:

- **flashier**: Willwerscheid J, Carbonetto P, Wang W, Stephens M. *flashier: Empirical Bayes Matrix Factorization.* R package 1.0.60, commit `fd44811`, https://github.com/willwerscheid/flashier. [Verified via `citation()` on the pinned install.]
- **ebnm**: Willwerscheid J, Carbonetto P, Stephens M (2025). *ebnm: An R Package for Solving the Empirical Bayes Normal Means Problem Using a Variety of Prior Families.* Journal of Statistical Software 114(3):1-33. doi:10.18637/jss.v114.i03. [Verified via `citation()`.]
- **mashr**: Urbut SM, Wang G, Carbonetto P, Stephens M (2019). *Flexible statistical methods for estimating and testing effects in genomic studies with multiple conditions.* Nature Genetics 51(1):187-195. doi:10.1038/s41588-018-0268-8. [Verified via `citation()`.]
- **cNMF**: Kotliar D et al. (2019). *Identifying gene expression programs of cell-type identity and cellular activity with single-cell RNA-Seq.* eLife 8:e43803. doi:10.7554/eLife.43803.
- **MOFA+**: Argelaguet R et al. (2020). *MOFA+: a statistical framework for comprehensive integration of multi-modal single-cell data.* Genome Biology 21:11. doi:10.1186/s13059-020-02015-1.

Comparators (cited by the package's benchmarks/vignettes):

- **GeneNMF**: Andreatta M, Carmona S. *GeneNMF: Non-Negative Matrix Factorization for Single-Cell Omics.* CRAN v0.9.2 (2025-09-11), GPL-3; NMF via RcppML (DeBruine et al. 2021, doi:10.1101/2021.09.01.458620); metaprogram concept traces to Barkley et al. (2022) and Gavish et al. (2023). [Verified from CRAN README/refman.]
- **MOFAcellulaR / multicellular factor analysis**: Ramirez Flores RO, Lanzer JD, Dimitrov D, Velten B, Saez-Rodriguez J (2023). *Multicellular factor analysis of single-cell data for a tissue-centric understanding of disease.* eLife 12:e93161. doi:10.7554/eLife.93161. [Verified.]
- **DIVAS**: Sun Y, Marron JS, Lê Cao K-A, Mao J (2026). *DIVAS: an R package for identifying shared and individual variations of multiomics data.* Bioinformatics, btag615. doi:10.1093/bioinformatics/btag615; method paper Prothero J et al. (2022). *Data Integration Via Analysis of Subspaces.* arXiv:2212.00703. doi:10.48550/arXiv.2212.00703. [Verified.]
- **mcRigor**: Liu P, Li JJ (2025). *mcRigor: a statistical method to enhance the rigor of metacell partitioning in single-cell data analysis.* Nature Communications. doi:10.1038/s41467-025-63626-5. [Verified.]
- **Spectra** (supervised comparator family): Kunes RZ, Walle T, Land M, Nawy T, Pe'er D (2023). *Supervised discovery of interpretable gene programs from single-cell data.* Nature Biotechnology. doi:10.1038/s41587-023-01940-3. [Verified.]
- **decoupler/CollecTRI convention**: implemented in `prepare_collectri` (OmniPath interactions endpoint, `datasets=collectri&genesymbols=yes`; `mor` = -1 only for consensus-inhibition-only; NFKB/AP1 complex merging). Primary decoupler/CollecTRI references to be pinned in Phase 9.

## 2. Reconstructed methodology, component by component

Notation: for cell type c, Y_c is the n_c x p_c pseudobulk matrix (donors x genes); the EBMF model is Y_c = L_c F_c^T + E_c with L_c donor scores, F_c gene loadings.

### 2.1 EBMF core — `R/ebmf.R` -> flashier

| Element | Specification | Implementation | Classification |
|---|---|---|---|
| Objective | EBMF ELBO with per-loading empirical-Bayes priors (flashier) | Delegated entirely to `flashier::flash_init -> flash_greedy(greedy_Kmax = max_factors) -> flash_backfit -> nullcheck` | **Faithful use** (no reimplementation of the factorization) |
| Priors | Loading prior families from ebnm | `loading_prior` in {point_laplace (default), point_normal, unimodal} -> `ebnm::ebnm_point_laplace/_point_normal/_unimodal`; scores use `ebnm::ebnm_normal` | **Faithful use** |
| Residual variance | var_type pass-through | 0 = constant; 1 = one variance per row/donor (default); 2 = per-gene (documented "unstable on log-normalized pseudobulk"); c(1,2) = Kronecker rank-one | **Faithful use** with documented operational guidance (Phase 3C will test the var_type=2 claim) |
| Fallback ladder | Robustness wrapper | point_laplace -> point_normal (ebnm_fn swap) -> no_nullcheck variants; non-finite ELBO treated as failure (`.fit_ok`); every attempt re-seeded identically; fallback warned and recorded in `fits[[ct]]$fallback` | **Deliberate extension** (sound; requires that fallbacks are surfaced to users — verified in code, runtime check in child A) |
| Residual-df guard | K capped by residual degrees of freedom (n - q - 1) | Present (referenced in tests) | **Deliberate extension** |
| Feature selection | top-N variable genes / preset artifact blacklist opt-in | `features = "variable"`, `n_variable_genes`; blacklist never consumes variable slots; explicit `exclude_genes` kept | **Deliberate extension** |
| Identifiability | Rotation/sign/scale ambiguity inherent to EBMF | Handled downstream by canonicalization; factor order arbitrary | Documented |

### 2.2 Covariate adjustment — two paths

| Path | Specification | Implementation | Classification |
|---|---|---|---|
| Hard residualization | Per-gene OLS on covariate design X, then EBMF on residuals | Rank-guarded model matrix (redundant columns dropped); per cell type | **Standard extension** (two-stage regression + EBMF) |
| Fixed mode | Joint fit with covariate factors fixed | `.flash_fixed`: fixed factors initialized at OLS fitted values (`flash_factors_init`, normal prior mode 0, scale 10), free factors added greedily, backfit over free factors | **Faithful use** of flashier fixed-factor machinery; **suspected defect**: no rank guard on `solve(crossprod(Xc))` (CAND-02, child A) |
| Interpretation | Residualized scores are orthogonal to the covariate span by construction | Documented in `?fit_celltype_programs` | Correct caveat; association tests on residualized scores cannot detect covariate-associated biology (by design, not defect) |

### 2.3 Canonicalization — `R/extract.R`

- Specification: simultaneous sign flip and RMS scaling, Z* = Z/s_k, W* = s_k W_k, preserving the reconstruction exactly; PSD multipliers stored and applied downstream rather than transforming PSDs in place.
- Hand verification: tcrossprod(S/s, sW) = tcrossprod(S, W) for any nonzero s; sign flip cancels. Zero-norm guard needed (RMS = 0).
- Classification: **equivalent implementation** of a standard convention. Runtime invariance + upstream flashier.utils parity claim ("verified to 0 difference") pending (child A, CAND-11).

### 2.4 Stability — `R/stability.R`

- Specification (package-defined estimand, no external reference): B donor subsamples; refit with replayed feature sets (exclusions carried inside `pre$feature_sets`), covariates, and seed + b; match programs to the reference fit by loading cosine; report recovery frequency, loading/score similarities, sign consistency, top-gene overlap.
- Classification: **deliberate extension** (package-defined stability notion). Open issues: failed-refit denominator (CAND-01) and non-replay of `preprocessing_spec` (CAND-12).

### 2.5 Principal angles / sharing spectrum — `R/principal-angles.R`

- Specification: SVD principal angles between program subspaces (scores or loadings space); permutation null by permuting donor labels within strata (`.perm_within`); pair-level statistic = sum of squared cosines; graded effect sizes (z_pair, excess_frac).
- Asymptotic reference: chance-level largest cosine = (sqrt(k_a) + sqrt(k_b))/sqrt(p_eff) — the Marchenko-Pastur spectral edge for two independent random subspaces of dimensions k_a, k_b in p_eff dimensions. **Theoretically justified** (standard random-matrix result); `n_shared_05` requires both permutation p < 0.05 and cosine above this reference (effect-size guard); `underpowered` guard when the reference >= 1 (pair never called); BH across callable pairs with an n_perm floor warning.
- `match_programs`: greedy loading-cosine matching; docstring says "Hungarian-style" — **documentation defect** (CAND-09); greedy-vs-optimal behavior quantified by child D.
- Classification: **theoretically grounded extension**; permutation validity of `.perm_within` to be exercised by child D/B evidence.

### 2.6 Program sharing — `R/program-sharing.R`

- Specification: for matched program pairs across cell types, classify shared/partial/private gene-loading structure using (i) a gene-label permutation null on loading correlation and (ii) mashr-style empirical-Bayes classification with loading posterior means + SDs, using a CorShrink-style shrunk null-gene correlation matrix (`.estimate_V`) as the mashr covariance; PSD handled via `F_psd * scale_multiplier`; tiny-SE floor and broken-gene dropping documented in code.
- Relationship to mashr (Urbut 2019): mashr estimates condition-specific effects with correlated residual covariance across conditions. Here it is **adapted** from condition effects to cross-cell-type program loadings — a nonstandard use whose key assumption (posterior SDs from the factorization are valid measurement SEs, conditional on the fitted model and the matching/orientation steps) is the central statistical question. Classification: **deliberate extension**; validity partially open (CAND-07 permutation-null exchangeability; SD-validity question in Phase 8).
- The shared-loadings vs coordinated-activity distinction is maintained structurally: `assess_program_sharing` (loadings) vs `find_coordinated_programs` (donor activities). Documentation consistency check in Phase 8.

### 2.7 CMP detection — `R/coordination.R` (novel procedure; internal mathematical specification)

- Estimand: cross-cell-type coordination of donor-level program activities, beyond within-cell-type structure.
- Pipeline: rank-normalize scores; residualize on strata (rank-guarded); pairwise donor-intersection correlations; T = atanh(r)*sqrt(n-3) with |r| clipped at 0.999999; node statistic A_v = sum_c max(0, max|T| - tau)^2 (tau = 2); max-T permutation calibration (donor labels permuted within strata, independently per cell type, all program rows of a cell type permuted jointly — preserves within-type program correlation); p = (1 + #{max >= a})/(B+1); significant nodes joined by a greedy graph with edge admission |T| >= link_threshold (3) and min_shared/min_donors floors; consensus scores from member activities.
- Null hypothesis tested: **complete null** (no cross-type coordination anywhere) for the max-T calibration — family-wise error control over node rejections holds by construction under exchangeability. Strong control under partial alternatives is NOT implied by the construction and is not claimed in code; membership precision of the greedy graph is a separate, uncalibrated quantity (CAND-04).
- Classification: **novel extension** (no external reference implementation). Its validation is entirely empirical — child B's calibration study is the evidence base; the benchmark artifacts (cmp-bench-v1.2) are development evidence, not calibration proof.

### 2.8 Interpretation layer — `R/cmp-interpret.R`, `R/cmp-annotate.R`

- TF footprints: decoupler-style ULM (competitive) + self-contained sign-flip gate (`p_self`) to prevent "depleted regulon" artifacts; gene-permutation max-|t| for p_fwer; in-code documentation states the competitive null is anticonservative under correlated gene blocks and that footprints are resemblance, not causality. Classification: **documented heuristic with explicit limitations** (rare and commendable).
- Ligand-receptor: degree-preserving resource permutation (ligand labels among ligands, receptor labels among receptors), max-over-CMP null; FWER or pooled-FDR control; optional MOFAtalk-style donor co-variation check; "compatible with, not evidence of causation" documented. Classification: **documented heuristic**.
- Themes (ULM on |loadings|), gene-overlap classes (hypergeometric, BH), context association (OLS on consensus scores, constant covariates dropped): standard, documented.
- Nuisance screens (`cmp_nuisance_genes`): documented screening heuristics; VDJ regex misses joining genes (CAND-08).

### 2.9 Frozen preprocessing + projection — `R/preprocess-freeze.R`, `R/cell-projection.R`

- Specification: content-addressed (sha256) preprocessing spec computed before the digest field is added; projection replays the frozen spec fail-closed (refuses legacy/fixed fits, undeclared feature_scale, design mismatch); sparse inputs unsupported v1 (documented).
- Classification: **deliberate extension** with fail-closed semantics (the conservative choice). Numerical + biological validation pending (Phase 10).

### 2.10 Simulation — `R/simulation.R`

- Linear-Gaussian pseudobulk simulator; scenarios private / mixed / same_genes_independent / global with planted shared-activity vs shared-genes controls (the controls needed to validate CMP and sharing truth claims).
- Classification: **documented limitation** (CAND-10): no count-level (NB/dropout/ambient) process, so count-level pseudobulk behavior cannot be validated by this simulator; Phase 12's independent simulation matrix must state this scope explicitly.

## 3. Deviation classification summary

| Component | Classification |
|---|---|
| EBMF fitting, priors, var_type | Faithful use (delegation) |
| Fallback ladder, residual-df guard, feature selection | Deliberate extension |
| Hard residualization | Standard extension |
| Fixed-mode adjustment | Faithful use + suspected defect (CAND-02) |
| Canonicalization | Equivalent implementation |
| Stability | Deliberate extension (package-defined estimand) |
| Principal angles / sharing spectrum | Theoretically grounded extension |
| Program sharing | Deliberate extension (mashr adaptation; assumptions open) |
| CMP detection | Novel extension (empirical validation required) |
| Interpretation layer | Documented heuristics |
| Freeze + projection | Deliberate extension (fail-closed) |
| Simulation | Documented limitation |

## 4. Open specification questions routed to later phases

1. Does var_type = 1 (per-donor variance) absorb biological heterogeneity on log1p pseudobulk, and does the var_type = 2 instability claim reproduce? (Phase 3C; child A-5)
2. Are loading posterior SDs valid measurement SEs for the mashr adaptation, given factor selection, matching, sign orientation, and data reuse? (Phase 8; child D)
3. Does the CMP max-T procedure maintain strong control under partial alternatives, and what error rate does the greedy graph + link_threshold = 3 admission actually control? (Phase 7; child B)
4. Is the gene-label permutation null for sharing anticonservative under gene-gene covariance? (Phase 8; child D-1)
5. What exactly does projection estimate for new donors/cells, and is its uncertainty calibrated? (Phase 10)
6. What evidence supports the K = 10 vs K = 30 operating points, and does rank selection respond to signal rather than to reconstruction alone? (Phase 3D; model_selection benchmark)
