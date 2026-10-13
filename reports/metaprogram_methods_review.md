# Metaprogram detection across separate per-cell-type EBMFs: method evaluation

Date 2026-10-09 (revised after two adversarial review rounds, see section 9). Goal (user): find programs that recur, are shared or are coordinated across cell types, simultaneously across all cell types, from separate per-cell-type EBMF fits. The pairwise principal-angle output is not the target. Each metaprogram carries a label `coordinated` (activity covaries across donors), `gene_reuse` (loadings overlap) or `both`, reported separately.

**Evidence labels:** `[F]` fetched this session (web pages were read through a summarizer, so details may be thin). `[R]` repo file read. `[M]` my recollection, not verified. `[I]` my inference. Unlabelled statements are design choices or definitions.

## 1. The structural problem every method must solve
- Each cell type has its own program set and its own donor set. Programs are not aligned across cell types (program 3 in B is unrelated to program 3 in T).
- **Donor coverage is partial.** In Biomni's 10-cell-type benchmark only 9 of 124 COMBAT donors have all 10 cell types `[R biomni/REVIEW.md]`. COMBAT has 16 cell types `[R README.md]`, so complete coverage across all of them is rarer still `[I]`. A complete-case method collapses to a few cell types, or to pairs.
- A method therefore needs (a) a way to use incomplete donor coverage, (b) a null that respects donor-level confounders (Institute, sex), (c) a program-level output with breadth (how many cell types), not a cell-type-pair verdict.
- Zero-filling absent donors leaks the missingness pattern into the signal `[R biomni/REVIEW.md finding 2]`.
- **Output unit matters for scoring.** The draft benchmark scores recall, precision and breadth per (cell type, program) node. Methods that do not output nodes (`pa_cc`, MOFA factors, DIVAS components, GMDF programs) need an adapter or a second, cell-type-level metric tier (Biomni review, finding 8).

## 2. Summary table

| Method | Input | Output unit | Simultaneous | Partial sharing | Missing donors | Strata null affordable? | Label | Verdict |
|---|---|---|---|---|---|---|---|---|
| `node_cor` (repo) | EBMF scores | node pairs, merged into groups | pairwise, then merged | via clustering | observed donors per pair | yes, cheap | coordinated | keep |
| `subspace_sum` (repo) | EBMF scores | donor direction + member cell types; node by best `abs(cor)` | yes | membership gate (uncorrected 95% quantile) | weighted, heuristic | yes, cheap | coordinated | keep, fix step-down |
| `gene_match` (repo) | EBMF loadings | node pairs, merged | pairwise, then merged | via clustering | n/a (genes) | gene-label null, no strata | gene_reuse | keep |
| `pa_cc` (repo) | EBMF scores | cell-type groups, no node | pairwise | no | pairwise | yes | coordinated (cell-type level) | baseline |
| Stacked two-stage EBMF (`ebmf2`) | EBMF scores (and loadings on the union panel for gene reuse) | stage-2 factors with sparse node loadings, member nodes directly | yes | yes, via sparse node loadings | `NA` rows, native | needs a stage-2 refit per permutation, small | coordinated (gene_reuse on loadings) | add as in-repo arm |
| DIALOGUE-style multi-CCA | EBMF scores | one canonical variate per cell type, a mixture of programs | yes | no: all blocks present | complete-case only | moderate (refit per permutation) | coordinated | add |
| DIVAS | scores or pseudobulk (unverified) | joint / partial / individual components per block set | yes | **yes, core design** | unverified | no strata option documented | coordinated | add after reading methods |
| MOFA / MOFAcellulaR | **pseudobulk** `[R R/mofacellular.R]` | factors with per-view weights | yes | via per-view weights `[M]` | native `[M]` | infeasible (refit per permutation) | coordinated | add, power and null data only |
| GMDF | pseudobulk | gene programs, shared + context-specific | yes | only pre-annotated contexts | rows need not match `[I]` | n/a (gene space) | gene_reuse | optional |
| nnTensor jNMF / siNMF (non-negative) | pseudobulk counts | shared + specific gene programs (genes as rows) or activity (donors as rows) | yes | only via penalised `H_k` rows | genes-as-rows: tolerant; donors-as-rows: no | n/a in gene space | gene_reuse (coordinated if donors as rows) | add as gene-reuse comparator, preferred over GMDF |
| RaJIVE / aJIVE | EBMF scores | joint, individual, residual per block | yes | **no (joint = all blocks)** | complete samples | cheap (theoretical rank bounds) | coordinated, all-block only | add as comparator and negative control |
| Tensor CP / Tucker | pseudobulk | components with a cell-type mode | yes | via cell-type mode | needs masking | infeasible | coordinated | optional |
| PARAFAC2 | pseudobulk, shared genes | per-slice donor scores | yes | unclear | rows may differ | n/a | gene_reuse-like `[M]` | optional |
| `multiblock` shared-sample (GCA, DISCO, SCA...) | EBMF scores | block and common scores | yes | DISCO local components `[M]` | complete rows `[I]` | moderate | coordinated | optional |
| `multiblock` shared-variable (JIVE, STATIS, HOGSVD) | pseudobulk | joint and individual gene structure | yes | JIVE: no `[M]` | rows need not match `[I]` | n/a | gene_reuse | optional |
| CRAN `multiCCA` | n/a | n/a | n/a | n/a | n/a | n/a | n/a | reject |
| `mbpls` | n/a | n/a | n/a | n/a | n/a | n/a | n/a | reject |

## 3. Per-framework evaluation

### 3.1 cNMF consensus (`cnmf.py`, `consensus`) `[F]`
- **What it does:** pools spectra from many replicate NMF runs per K, L2-normalizes them, drops outliers by local density (mean distance to the nearest `0.30 * replicates / K` neighbors, drop if >= 0.5), runs k-means with K clusters, takes the per-gene median of each cluster, refits usages against the fixed medians, and picks K by silhouette stability against prediction error.
- **Borrow:**
  - A density filter on programs pooled over replicate fits: keep a program only if it has near neighbors, so unstable programs never seed a metaprogram.
  - The per-gene median of a metaprogram's members as its signature (after sign alignment with `canonicalize_programs`).
  - Refit on fixed signatures, to score low-count donors or a new cohort. It cannot score a donor absent from a cell type, which has no data there.
- **Fit:** a shared pre-filter on nodes before any cross-cell-type step. It generalizes `recovery_freq >= 0.7` to a continuous score. Note the replicates differ: cNMF varies NMF initialization, ours would vary donor subsamples, so the density measures sampling stability.
- **Pros:** cheap to add, improves every arm.
- **Cons:** fixed-K k-means does not fit metaprograms (unknown count, variable breadth, at most one node per cell type). Many replicates per K is costly: our `k_stab` arm was 9.7x over the time guard with 10 subsamples `[R HANDOFF.md]`.
- **Verdict:** borrow the filter and the median; do not copy the clustering.

### 3.2 Tensor factorization (CP, Tucker, PARAFAC2; scITD and Tensor-cell2cell style) `[M]`
- **Borrow:** the cell-type mode reads out which cell types hold a component, so breadth is a direct output.
- **Fit (corrected 2026-10-09; earlier text said "poor" too strongly):** a tensor decomposition is feasible as a comparator on pseudobulk once a common gene panel exists (section 10) and absent donors are masked (TensorLy has a `mask` option for CP `[M]`; the TensorLy user-guide page I fetched shows only CP, Tucker and tensor-train basics and does not cover masks, PARAFAC2 or sparsity `[F]`). The limits are what the model assumes: a CP component is donor activity x gene signature x cell-type weights, so every member cell type shares **both** the same activity and the same genes, scaled by a cell-type weight. It represents only the `both` case with activity correlation 1. A program that is coordinated across cell types but uses different genes in each cannot be one CP component; it splits into several components with unrelated gene patterns and near-duplicate donor vectors. Tucker's core tensor relaxes this but gives dense, rotation-ambiguous cell-type factors and no direct breadth. It also replaces the per-cell-type EBMF with a joint model, so coordination is assumed instead of discovered, and private programs have to be extra components. A donor x program x cell type tensor stays meaningless because programs are not aligned across cell types.
- **PARAFAC2 differs:** it lets rows differ per slice, so it tolerates different donor sets per cell type. For the same reason it does not link donors across cell types, so it behaves like GMDF: gene reuse, not donor coordination.
- **Pros:** one joint model with breadth as an output.
- **Cons:** alignment and missingness, no strata null, no node-level output.
- **Verdict:** optional comparator on raw pseudobulk only; not in the package.

### 3.3 DIVAS, bioRxiv 10.64898/2026.01.12.698985 (Sun, Marron, Le Cao, Mao, Jan 2026) `[F: abstract and README only]`
- **What it does:** R package for jointly shared, **partially shared** and individual variation across data blocks. Angle-based subspace analysis, hierarchical search over combinations of blocks, rotational bootstrap inference, COVID-19 case study. AGPL-3.0, 78 commits, 12 stars, repo `ByronSyun/DIVAS_Develop`. Main call `DIVASmain(list_of_blocks)`.
- **Borrow:** the joint / partial / individual split, which is the R100 / R50 / R20 structure, and its inference.
- **Fit:** blocks = cell types. Our principal angles is the pairwise version of the same angle idea `[I]`.
- **Unverified:** whether blocks must share all samples, how missing samples are handled, rank selection, runtime. For 10 cell types the search can cover up to 1023 block subsets `[I]`.
- **Concern (my inference):** DIVAS is built for high-dimensional blocks with noise. On blocks of at most 10 program scores its rank estimation has little to do. It may be better suited to pseudobulk, which bypasses the EBMFs.
- **Pros:** closest published design to partial sharing.
- **Cons:** AGPL-3.0 cannot be bundled in an MIT package, so use it as an external comparator only. No strata-based null is documented. Young package.
- **Verdict:** add as a comparator arm after reading the full methods.

### 3.4 CRAN `multiCCA` (Gorecki, v0.1.0, 2026-03-23) `[F]`
- Multiple **kernel** CCA and multiple **functional** CCA for repeated measures. Sole author, v0.1.0, no documented missing-data or permutation support.
- Very likely not the multi-CCA DIALOGUE uses (its package name is unconfirmed, section 8).
- **Borrow:** nothing. **Fit:** none. **Pros:** none for this task. **Cons:** wrong problem, immature.
- **Verdict:** reject.

### 3.5 DIALOGUE multi-CCA (`DIALOGUE.main.R`, `DIALOGUE1.PMD`) `[F]`
- **What it does:** per cell type, reduce to sample-level medians, ANOVA-filter features (BH 0.05), center and scale. Rows are matched across cell types by sample name, and **only samples present in every cell type are kept** (stops below 5). It calls `MultiCCA(...)` with `ncomponents=k`, a penalty from `MultiCCA.permute`, `niter=100` (package name not shown; PMA is likely `[I]`). Per-cell-type score is `X_c %*% w_c`.
- **Significance:** shuffles each feature column within each cell type, reruns, compares correlations with a rank-sum test, and keeps pairs with `emp.p < 0.1` (a lenient threshold). Signatures come from gene correlations with the residualized scores, then mixed models and NNLS.
- **Borrow:** one canonical variate per cell type per program (our "one node per cell type"); sparse penalties; correlation-based signatures.
- **Fit:** run on EBMF scores (at most 10 columns per cell type) instead of genes.
- **Pros:** simultaneous across cell types; cheap on scores; sparsity selects which programs take part.
- **Cons:**
  - **Output unit:** a canonical variate on scores is a mixture of a cell type's programs, not one node. Assign it to the dominant weight, with sparsity helping.
  - **Complete-case only:** with partial coverage it needs cell-type subsets, and enumerating subsets turns it into DIVAS's combinatorial search or into pairwise runs.
  - **Null:** the feature-shuffle null misses donor-level confounders (my inference for why the repo's PMD run called 45 of 45 pairs shared; unconfirmed). The penalty search uses the same null, so it needs strata-aware permutation too.
  - **Not the PMD question:** a run on scores tests a different model than DIALOGUE on genes. It does not settle whether the earlier PMD result was caused by the null.
- **Verdict:** add as a comparator arm, with a donor-level within-strata null for both testing and penalty choice.

### 3.6 `mbpls` (Python, Multiblock PLS) `[F]`
- Supervised: `MBPLS(n_components).fit(X, y)` predicts a response from blocks. The docs page does not describe unsupervised use. It says missing data is handled without detail, and gives no version or maintenance statement.
- **Borrow:** nothing for discovery. **Fit:** none, there is no response. **Pros:** possible later use predicting severity from metaprogram scores, a different task. **Cons:** supervised, Python bridge.
- **Verdict:** reject for discovery.

### 3.7 R `multiblock` unsupervised methods `[F: names and call shapes only]`
- The vignette lists SCA, GCA, GPA, MFA, PCA-GCA, DISCO, HPCA, MCOA as shared-sample methods and JIVE, STATIS, HOGSVD as shared-variable methods. Input is a named list of blocks, output has `scores`, `loadings`, `blockScores`, `blockLoadings`. It does not explain algorithms or missing-data support.
- **Shared-sample group (donors = rows):** from memory `[M]` GCA solves a generalized CCA and DISCO separates common, local (some blocks) and distinct components. It appears to need complete rows `[I]`. DISCO is the only one aimed at partial sharing.
- **Shared-variable group (genes = columns):** these take datasets sharing the gene panel, so donors need not match. They are gene-reuse comparators like GMDF `[I]`. JIVE gives joint plus individual with no partial sharing `[M]`.
- **Borrow:** DISCO's local components. **Pros:** one package, common interface. **Cons:** missing-data behaviour unknown; algorithms undocumented in the vignette.
- **Verdict:** do not adopt now; DISCO is an optional extra on complete-donor subsets.

### 3.8 MOFA2 / MOFAcellulaR `[M]`
- Cell types as views, donors as samples. MOFA handles missing samples per view natively and has per-view sparsity, so factor weights show which views a factor loads on `[M]`.
- **Fit:** the repo exports **raw pseudobulk** for it `[R R/mofacellular.R:1-8]`, so it bypasses the EBMFs, like GMDF.
- **Borrow:** the joint model with native missing-view handling.
- **Pros:** the one candidate here with documented-by-memory native missing-donor handling; export already exists; the simulator already returns pseudobulk `[R sim_coord.R:83]`, so no new simulator work.
- **Cons:** a strata null needs a full refit per permutation, infeasible over 200+ replicates and several scenarios, so it can be scored on power and null data only without a calibrated FPR. Factors need a mapping to nodes before node-level scoring.
- **Verdict:** add as a comparator, evaluated without a calibrated null.

### 3.9 GMDF, Generalized Matrix Decomposition Framework (`livnatje/GMDF`) `[F: README only]`
- **What it does:** input is a list of expression matrices `E`, one per dataset, an annotation matrix `a` (datasets x contexts), `k` shared programs, `k1` context-specific programs per context, and `N1` runs to combine. Each dataset is modelled as `Hw_i x W` (shared programs `W`, dataset-specific usages) plus context-specific programs weighted by `a`. Depends on `plyr` and `rliger`. With `N1 > 1` it combines runs and clusters programs (`W.multi.run`, `W.clusters`, `Wf`, `sig`), a consensus step like cNMF.
- **What it detects:** shared gene programs (`gene_reuse`), not donor coordination. Each dataset has its own usage matrix, so donors need not match `[I]`, which also makes it indifferent to missing donors.
- **Mapping:** dataset i = cell type i. Partial sharing only through contexts you annotate in `a`, so a random 5-of-10 or 2-of-10 program is detectable only if an annotation covers that subset.
- **Fit:** factorizes the data, not our EBMFs. Use as an independent gene-reuse comparator on the same pseudobulk, then compare `W` with EBMF node loadings. Inside the package it would replace `fit_celltype_programs`, which the repo's earlier decision rules out (ROADMAP section E).
- **Unknown:** nonnegativity (an `rliger`-style NMF would need a shift for log-normalized pseudobulk), rank choice, significance testing. No citation, license, benchmarks, 7 commits, 2 stars.
- **Pros:** missing-donor indifference; borrows strength across cell types for gene programs.
- **Cons:** gene reuse only; annotation-dependent partial sharing; research-script maturity; an unlicensed repo cannot be bundled.
- **Verdict:** optional gene-reuse arm, below `gene_match`.

### 3.10 Stacked two-stage EBMF ("EBMF of EBMFs") `[I]`, toy check run 2026-10-09
- **Idea:** stage 1 is the existing per-cell-type EBMF. Stage 2 runs `flashier` again on the donor x node matrix formed by stacking every cell type's program scores (columns unit-RMS scaled, a donor's rows absent from a cell type left as `NA`). Each stage-2 factor is a metaprogram. Its loading over nodes is sparse under a point-Laplace prior, so the nonzero nodes give the member cell types and breadth directly. A factor with one nonzero node is a private program. For the gene-reuse label, run the same stage on loadings over the common gene panel (route B of section 10).
- **Why it fits:** no change to the per-cell-type fits; no pairwise tests; missing donors enter as `NA` (flashier accepts missing entries `[M]`, and the toy below ran with 22% missing); sparsity gives variable breadth and one-node-per-cell-type is checkable from the loadings. Node score uncertainty could be passed through flashier's `S` standard errors `[M]`, as in the `w_S` arm.
- **Toy check (feasibility only, not a benchmark):** 6 cell types, 4 programs each, 100 donors, 22% missing, a metaprogram A in 3 cell types and a metaprogram B in 2, signal sd 0.5. Over 12 seeds flashier (point-Laplace, greedy, backfit) recovered A (all 3 member nodes) in 10 of 12 seeds, **never recovered B**, and returned about one factor per run with no spurious factors under my loose scoring. So the stage works mechanically, but power for a narrow 2-cell-type program looks weak with default settings. Scripts were in the session scratchpad; no repo files were changed.
- **Follow-up toy (complete donors, 20 seeds per setting, scores scored differently from the run above, so not directly comparable):** stage-2 EBMF recovered the donor scores of both A and B in 90-100% of seeds and the exact member-node set in 80-100%. Plain PCA recovered the scores equally well but the exact node set for B in only 50% (dense loadings). `ica::icafast` with 6 components recovered scores in only 30% of seeds when the latent z was Gaussian and 50-70% when z was skewed, with the exact node set in 0-50%. So the earlier miss of B appears to come mainly from the 22% missing donors, not from the narrow breadth alone. One toy setting, arbitrary thresholds and component count; not a benchmark.
- **Risks:** stage 2 captures any donor-level covariance among nodes, including Institute or sex, so it needs the same within-strata permutation null or residualization as the other arms; stage-1 score error is ignored unless `S` is used; K and the number of stage-2 factors are choices; permutation nulls need a stage-2 refit per permutation, small here (about 100 x 100) but unmeasured.
- **Relation to rejected ideas:** ROADMAP section E rejects shared priors (L1) and initialization borrowing (L2) inside the per-cell-type fit. This stage leaves the fits untouched, so it does not revisit them.
- **Joint single-model alternative (not recommended):** one EBMF on donors x (genes x cell types) with shared donor scores. Flashier's priors are element-wise, so it has no native group-level sparsity to make a factor load on only some cell-type blocks `[M]`, and the earlier concatenation result at K >= 30 overfit `[R README.md]`.
- **Verdict:** add as an in-repo comparator arm (`ebmf2`) in Phase 1 or 2, after the common-gene-panel work. It is cheap to build on existing code but is unvalidated.

### 3.11 nnTensor: jNMF, siNMF and NTF (CRAN vignettes nnTensor-2 and nnTensor-3) `[F: partial vignette text]`
- **Models:** siNMF `X_k ~ W H_k` and jNMF `X_k ~ (W + V_k) H_k`, all factors non-negative, with `X_k` of size N x M_k. The **rows N are shared** across matrices and the columns M_k may differ (`jNMF(X_list, algorithm = "KL", J = 3)`; outputs `W`, list `V`, list `H`, error traces). NTF is non-negative CP; the vignette says NTF does not work on Tucker-structured toy data, consistent with the CP assumption in section 3.2. Missing-value masks, penalties and rank selection are not shown in the text I read.
- **Two ways to orient our data:**
  - **Rows = genes (union panel), columns = donors per cell type.** `W` = shared gene programs, `V_k` = cell-type-specific variant of the same program, `H_k` = donor usages that need no matching across cell types. Gene reuse (`gene_reuse`), tolerant of missing donors. Same family as GMDF but documented and on CRAN, so **prefer it over GMDF** for the gene-reuse arm `[I]`.
  - **Rows = donors, columns = genes per cell type.** `W` = shared donor activity, `V_k` = cell-type-specific activity, `H_k` = cell-type-specific gene loadings. This can represent coordinated activity with different genes per cell type, which CP cannot. Costs: rows must be shared, so absent donors are a problem (no mask shown); partial sharing needs `H_k` rows driven to zero by a penalty not shown; the common/specific split is ambiguous.
- **Non-negativity is the main barrier.** EBMF scores are signed, so jNMF cannot be applied to them; it can only run on non-negative pseudobulk counts, bypassing the EBMFs, and covariate residualization is impossible. Making EBMF non-negative (an exponential-type prior in `flashier` `[M]`) would change the fits and need its own benchmark arm.
- **Verdict:** add jNMF (genes-as-rows orientation, KL on pseudobulk counts, union panel) as the gene-reuse comparator; donors-as-rows only as an optional coordination comparator on complete-donor subsets.

### 3.12 RaJIVE (`ericaponzi/RaJIVE`, robust angle-based JIVE) `[F: README only]`
- R package on CRAN and GitHub (25 commits, 0 stars, licence type not named on the page). `Rajive(data, initial_signal_ranks)` splits a list of blocks into joint, individual and residual parts, with `get_joint_rank()` and `get_individual_rank()`. Blocks share the sample dimension (inferred from the example). Missing values are not mentioned.
- **Relation to jNMF:** the same idea of common plus block-specific parts, but signed, SVD-based, with the individual part separate from the joint part, instead of non-negative with `V_k` tied to the same `H_k`.
- **Fit to the EBMF framework:** it can take each cell type's EBMF `K_c` as `initial_signal_ranks` and run on score blocks. As I recall `[M]`, angle-based JIVE stacks the per-block signal subspaces and takes an SVD, which is the same computation as `subspace_sum`, with theoretical thresholds for joint rank in place of a permutation null `[I]`.
- **Partial sharing:** "joint" means shared by all blocks. The DIVAS abstract states that earlier integration methods find shared-by-all or unique variation and often miss partially shared variation `[F]`. Expect RaJIVE to find R100 and miss or mislabel R50 and R20.
- **Verdict:** add as a cheap comparator and negative control on complete-donor subsets; its joint-rank estimate is also a useful calibration check for `subspace_sum`.

## 4. Relation to what is already coded
- `subspace_sum` (`benchmarks/coordination/methods.R:146-195` `[R]`) takes the top eigenvectors of the donor-weighted sum of projections `sum_c U_c U_c'` over the score subspaces. With **complete data** that is the MAXVAR generalized CCA objective `[I]`. With missing donors the `D^{-1/2}` weighting is a heuristic with no calibration evidence, so "multi-CCA on scores is already implemented" holds only in the complete-data case.
- DIALOGUE's multi-CCA optimizes pairwise correlations with sparsity, which is a different objective from MAXVAR `[I]`. What it adds beyond `subspace_sum`: sparsity (which programs per cell type) and permutation-chosen penalties.
- Missing for every in-repo method: cluster-level error control, a coherence statistic, and a repaired step-down (`cummax` plus an early `break` is a sequential stop rule, not a formal step-down, per the Biomni review).

## 5. Recommended plan with a stopping rule
**Phase 1, in-repo and cheap (no new dependency):**
1. Repair `subspace_sum`'s step-down; add a cluster-level coherence statistic and cluster-level test.
2. Add the cNMF-style density filter on nodes, shared by all arms.
3. **Common gene panel for the gene-reuse label (user request 2026-10-09: union of the per-cell-type HVGs).** `gene_match` and any loading-based comparison need loadings over the same genes. Today `features = "variable"` picks the top `n_variable_genes` by variance independently per cell type (`R/ebmf.R:304-330` `[R]`), and `shared_genes` is the intersection across cell types (`R/data-model.R:81` `[R]`). Two ways to get a union panel:
   - **(B) Refit loadings on fixed scores (recommended).** Keep each cell type's fit on its own HVGs, then compute loadings for the union panel by regressing each union gene's pseudobulk on the fixed program scores, restricted to genes measured in that cell type. No change to the fits or scores, which stay the primary analysis. This is a new post-fit helper, so it needs an `R/` edit and a test. It follows the same reasoning as cNMF's refit step (section 3.1) `[I]`.
   - **(A) Fit on the union.** Already possible without code changes: `features = list(<ct> = union_panel_available_in_<ct>, ...)`. The list branch errors if a requested gene is absent from a cell type, so the union must be intersected with each cell type's measured genes first. It changes every fit, so by the repo rule it needs its own benchmark arm (for example `genes_union`) before it becomes a default.
   - Dependency: M2's `genes_after` arm (HVG selection after covariate adjustment) decides whether the HVGs are chosen before or after adjustment, so define the union only after that result.
   - **HVG rule (user decision 2026-10-09): Seurat v3 (`vst`).** The package's own rule is plain variance, and the earlier COMBAT tuning used scran HVGs on TMM-normalized pseudobulk `[R README.md]`. Routes A and B are to be compared empirically, see section 10.
4. Reconcile the prereg method list with the code; apply the four blockers from the Biomni review (replicates and paired inference, node-to-truth matching rule, pinned EBMF config, within-stratum nuisance).

**Phase 2, one external arm with a calibrated null:** DIALOGUE-style multi-CCA on scores with a within-strata donor permutation for the test and the penalty search. Cell-type-subset runs are required at partial coverage. Compute is probably tractable (blocks of at most 10 columns) `[I]`.

**Phase 3, optional, bounded budget:** DIVAS (after reading its methods), MOFA (power and null datasets only), GMDF (gene reuse only, on the `genes_only` simulator mode against `gene_match`).

**Kill rule (pre-register):** drop an external arm if it cannot (a) produce node-level or cell-type-level output that can be scored, (b) run on at least 100 replicates per scenario within a stated CPU-hour budget, or (c) show a null-calibrated false-positive rate. Effort estimates are guesses until measured; none are given here.

## 6. Open design decisions (yours)
1. **`both` rule:** the `coordinated` and `gene_reuse` families each produce their own groups, and the groups can differ in membership. Is a metaprogram `both` when the groups intersect, when one contains the other, or by a membership overlap threshold?
2. **Scoring tiers:** one node-level tier for methods that output nodes, one cell-type-level tier for the rest, compared only within a tier.
3. **Null for external arms:** accept uncalibrated arms evaluated on power and null datasets only, or exclude arms that cannot be strata-calibrated?
4. **Subset handling:** how a complete-case method chooses cell-type subsets on real data.
5. **Union panel route:** settled as an empirical question (section 10): HVG rule is Seurat v3; routes A (fit on the union) and B (refit loadings on fixed scores) are tested head to head, with B as the default if neither wins by the registered margin.

## 7. Common requirements from the Biomni review of the draft prereg
Reconcile the method list, set replicates and paired inference, pre-register the node-to-truth matching rule, pin the EBMF fit config, add a within-stratum nuisance, align error control or report operating curves, add a null-calibration gate.

## 8. Open checks before building
- Read DIVAS full methods or vignette: sample overlap, missing-sample handling, rank selection, runtime.
- Confirm which package provides `MultiCCA` in DIALOGUE, and its license.
- Confirm the repo's PMD run used the feature-shuffle null, by reading `integrate-pmd.R` in the Biomni zip.
- Check licenses for GMDF and DIALOGUE before any code reuse.

## 9. Revision notes
- **Round 1 corrections:** removed "refit scores donors absent from a cell type" (cNMF); corrected the MOFA row to pseudobulk input and removed "only option with native missing-donor handling"; replaced "answers the PMD question properly"; scoped the MAXVAR claim to complete data; replaced "MIT-like" with an unchecked license; softened the claim about which package DIALOGUE uses; added multiblock shared-variable methods and the PARAFAC2 distinction; added compute and null-feasibility; merged the duplicated per-framework section into section 3.
- **Round 2 additions:** output-unit mismatch and the canonical-variate-as-mixture problem; COMBAT coverage scope (10 of 16 cell types); undefined `both` merge rule; lenient `emp.p < 0.1`; DIVAS fit concern on low-dimensional score blocks; a stopping rule and phase ordering that puts in-repo fixes before external arms.

## 10. Upstream HVG test plan: Seurat v3 `vst`, routes A and B (draft, not registered)
**Question:** when the gene-reuse analysis needs loadings over a common gene panel, is it better to (A) fit each cell type on the union of the per-cell-type HVGs, or (B) keep fits on each cell type's own HVGs and refit loadings over the union on fixed scores? Both start from the same upstream selection.

**Upstream selection (user decision):** Seurat v3 `vst`, applied per cell type. As I recall it `[M]`, `vst` fits a loess of log10 variance on log10 mean of raw counts, standardizes counts by the expected sd (clipped), and ranks genes by the variance of the standardized values. Raw counts are required. `Seurat` is **not installed** in the `cellprograms-r` environment `[R]`, so it must be added, or the vst rule re-implemented and checked against Seurat on a test matrix.
- **Level:** pseudobulk raw counts (samples as columns) is consistent with the earlier COMBAT tuning, which selected HVGs on pseudobulk. Cell-level selection is an optional variant: it is dominated by a few high-count donors. COMBAT has `pb_counts.parquet` `[R benchmarks/biomni_replication/data/combat]`; whether it is raw counts, and whether cell-level counts are available on this cluster, are open checks.
- **Size:** top 2000 per cell type, as in the earlier tuning. The union of 10 or more cell types can be several thousand genes `[I]`. An optional variant takes the top N genes by how many cell types call them variable, like Seurat's `SelectIntegrationFeatures` `[M]`.
- **Gene availability:** intersect the union with each cell type's measured genes before any fit (the `features` list branch errors on absent genes, `R/ebmf.R:310-315` `[R]`).

**Arms** (same pseudobulk, same seeds across arms):
- `own`: each cell type fit on its own HVGs; gene reuse scored on the intersection of own HVGs. This is today's practice and the baseline.
- `A_union`: each cell type fit on union ∩ measured genes. Needs no code change, only a `features` list.
- `B_refit`: fit as in `own`; then loadings over union ∩ measured genes by regression of each gene's pseudobulk on the fixed program scores. Needs a new post-fit helper plus a test that scores are bit-identical to `own`.
- Optional: `A_top`, the frequency-ranked top-N panel.

**Data:**
1. Semi-synthetic on real COMBAT pseudobulk, injecting gene programs on the count scale (the `semi` generator in `benchmarks/model_selection/sim_semi.R` is the template; it needs raw counts and a count-scale injection).
2. A count-level simulator (negative-binomial) with planted gene-reuse programs at 100 / 50 / 0 % gene overlap, extending `sim_coord.R`. The current simulator is Gaussian on a log scale `[R]`, where `vst` is not meaningful.
3. Real COMBAT, exploratory only: loading stability under donor subsampling for A vs B.

**Endpoints (pre-register margins after a pilot, as in M2):**
- Primary, gene reuse: Found rate of planted gene-reuse programs by `gene_match` on the common panel, against `own`; false-positive rate on null and nuisance-only data.
- Loading accuracy: correlation of estimated loadings with planted weights over panel genes. For `B_refit`, report separately on genes inside the cell type's own HVG set and on genes outside it.
- Scores (matters for A only, since B leaves them unchanged): planted-program recovery R as in M2 (best `abs(cor)` minus chance).
- Cost: fit time.

**Decision rule (draft):**
- Adopt `B_refit` over `own` if gene-reuse Found rate improves by the registered margin with false-positive rate no worse than `own` + 0.03. B cannot hurt scores by construction.
- Adopt `A_union` over `B_refit` only if it beats B on gene reuse by the registered margin **and** R is no worse than `own` − 0.02 **and** time is at most 3x `own`.
- If neither passes, keep `own` and report gene reuse on the intersection only.
- At least 200 paired replicates per scenario cell, paired bootstrap CIs, margins at least 2x the paired standard error (from the Biomni review).

**Sequence and timing:**
1. After M2 finishes and `R/` is safe to edit: add Seurat (or a checked vst re-implementation) to the environment; write `B_refit` helper + test; build the count-level simulator.
2. Pilot to set margins; register (tag); run; summarize.
3. Effort is a guess until measured: a day or two for the helper and simulator, plus run time. I have not estimated run time.

**Open checks:** raw vs normalized content of `pb_counts.parquet` and `mat/ct_*.csv`; availability of cell-level counts; Seurat installation or a validated re-implementation; whether HVGs are chosen before or after covariate adjustment, which depends on M2's `genes_after` arm.
