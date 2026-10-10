# cellprograms: selective coupled gene-program discovery and conditional hierarchical EBMF

Date: 2026-10-09  
Last updated: 2026-10-09  
Project: `/exports/para-lipg-hpc/mdmanurung/cellprograms`  
Repository: [mdmanurung/cellprograms](https://github.com/mdmanurung/cellprograms)  
Reviewed commit: `584c50236c78730f4b33fe2d82bfd211b71a5629`  
Corresponding Git tree: `f975b9ec1e4b074db762359f261db9628e6397bd`  
Status: reviewed research plan and progress tracker; implementation not started; optimizer blocked on GS-00  
Preregistration: not registered; numerical performance margins and selection rules are not frozen

## 1. Decision, objective, and authorization

Develop an experimental method for discovering fully shared, partially shared, and private transcriptional gene programs across scRNA-seq cell-type pseudobulk views. The primary conserved quantity is a signed gene-loading direction, up to reversal of the entire factor. Donor activities remain view-specific; donor matching and activity coordination are not required for discovery.

Start with a selective exact-core coupled factorization. Introduce flexible loading deviations and empirical prior learning only after simpler models demonstrate incremental benefit. A zero-shared-program result is valid.

**Primary hypothesis:** joint estimation improves recovery of weak recurrent loading directions over a strong independent-EBMF consensus method with fixed-template residual recruitment, without unacceptable false membership, false sharing, or loss of private programs.

The current user request authorizes saving this revised plan and updating `HANDOFF.md`. It does not authorize package implementation, dependency installation, scheduler submissions, benchmark execution, commits, pushes, or release. Later implementation requires a separate action request and completion of the scientific gates below. Read-only review and the limited validation recorded in section 12 have already occurred.

This is a separate gene-sharing workstream. Preserve the existing coordination workstream, independent fitting API, defaults, and registered model-selection benchmark. The earlier rejection of L1/L2 modifications to the independent backend remains in force; the new method is an explicitly experimental alternative, not an automatic reversal of those baseline decisions.

`HANDOFF.md` records an active model-selection main run using this checkout. That operational status was not refreshed during this documentation update. Do not modify files consumed by that run while it remains active. Documentation can be saved now; implementation should use an isolated checkout and immutable source snapshots for future jobs.

## 2. Review resolutions and remaining blockers

| Review finding | Revision and status |
|---|---|
| Shared contributions can be represented by private factors; scale/sign conventions do not establish identity | Addressed in the design: sharing is selected through a specified structural criterion. Exact penalties, rank budgets, ambiguity rules, and independent-only behavior remain **unresolved GS-00 blockers**. |
| Weighted error is not a complete learned-noise objective; posterior means do not define an EBMF ELBO | Full likelihood, precision behavior, penalty definitions, and update acceptance are required before optimization. Derivation remains **unresolved GS-00**, with a separate EB gate in GS-08. |
| Exact recurrence, approximate conservation, correlated loadings, and gene overlap lack an operational boundary | A loading-coordinate, practical-equivalence, signal, and eligibility contract is required. Numerical thresholds remain **unresolved GS-00/GS-04** and will be locked before confirmation. |
| Consensus-only initialization can miss programs detected in one source view; recruitment can be mistaken for evidence | Singleton source candidates, a fixed-template recruitment comparator, and independent target-view evidence are now required. Candidate and reporting rules must be specified in GS-00. |
| Comparator assumptions, adaptive calibration, failures, evaluation denominators, and package integration were underspecified | Explicit adapters, mixed-truth error metrics, full-pipeline calibration, truth-based evaluation, environment evidence, and gated package integration are now included. Their executable artifacts remain pending. |

The plan is ready to guide specification and feasibility work. It is not a frozen optimizer specification or a claim that any coupled method is validated.

## 3. Scientific and data contracts

### 3.1 Views and preprocessing

For view \(c\), let \(X_c\in\mathbb R^{n_c\times G}\) contain observations in rows and genes in columns. Observation sets and sizes may differ across views. Version 1 uses complete numeric matrices on an explicitly ordered common measured gene panel.

Keep `as_cell_program_data()` permissive about unequal panels for existing independent methods. Add joint-method panel validation at the coupled fitting boundary, rather than changing the general constructor.

Before fitting, record:

- Gene identifiers, duplicate-ID checks, panel construction, ordered genes, and genes excluded with reasons. Select features from training data only. A training-derived union of informative genes can be intersected with genes measured in all eligible views; do not silently reduce the panel to genes independently selected in every view.
- Transformation, nuisance adjustment, centering, and gene scaling. Store fitted transformation parameters for held-out projection. Preserve the original baseline's preprocessing behavior.
- View eligibility rules based on observation availability, precision, and supported input assumptions; thresholds must be fixed independently of whether local EBMF discovers a program.
- Handling of measured but constant/uninformative genes and completely degenerate views. Do not dynamically remove genes in different views and then pretend loadings share one panel. Never encode unavailable genes as expression zero.
- Sample/donor IDs, repeated measurements, cell counts when available, residual degrees of freedom, noise assumptions, data provenance, and source/environment versions.

Signed-loading conservation is defined in one common coordinate system. If a method uses view-specific gene scaling, map loadings back into that system before evaluation. A view-specific diagonal scaling changes loading direction and cannot generally be absorbed by one scalar amplitude. Weighting the likelihood by precision and changing the coordinate system are separate operations.

Unequal gene panels and entry-level missing features are deferred extensions. Version 1 should align compatible panels or reject unsupported inputs with a named error. Missing donor observations across views are allowed and are not zero-filled.

### 3.2 Program identity and reporting

| Quantity | Role |
|---|---|
| Orientation-invariant signed-vector cosine, \(|g^\top w|/(\|g\|\|w\|)\) | Primary loading recovery measure in the frozen common coordinates |
| Gene-support overlap and gene-weight rankings | Secondary; thresholds for support must be defined |
| Loading-subspace similarity | Secondary and primary fallback for explicitly ambiguous/split cases |
| Donor-score covariance or correlation | Post hoc; excluded from gene-sharing selection |

Before simulation results are used for method selection, define exact-core truth, practical approximate conservation, distinguishable correlated programs, minimum signal requirements, and an indeterminate region. A cosine of 0.8 cannot be labeled shared in one scenario and unrelated in another solely because of a hidden simulator label.

A reported shared program requires supported recurrence in at least two views. A single supported member is private; a singleton may still seed recruitment. Fully shared means supported in every eligible view under a prespecified denominator. Partially shared means supported in a proper subset of that denominator.

Report per-view states `supported`, `unsupported`, `indeterminate`, and `unobserved`, with reasons. An unsupported membership is not proof of biological absence. Computational failure has its own diagnostic state and cannot silently remove a view from the eligible denominator. Record both the planned view universe and the analyzable coverage.

Fitted model membership, reported evidence, and biological interpretation are distinct. A loading constrained to equal the common template does not supply independent evidence of recurrence. Common technical signatures can be statistically recurrent without being biological programs; qualify them using metadata and sensitivity analysis.

## 4. Model and objective specification

### 4.1 Exact-core selective coupled model

\[
X_c=\sum_{k=1}^{K_s}h_{ck}b_{ck}g_k^\top+U_cV_c^\top+E_c,
\qquad \|g_k\|_2=1,
\qquad h_{ck}\in\{0,1\}.
\]

Fit the view-specific activity vector \(b_{ck}\) directly in the initial model. For nonzero activities, amplitude and RMS-standardized scores may be reported afterward as

\[
a_{ck}=\sqrt{n_c^{-1}\sum_i b_{cki}^2},\qquad z_{ck}=b_{ck}/a_{ck}.
\]

For zero activities, use explicit inactive states; do not divide by zero. Apply a deterministic global sign convention to each template with the compensating change to its activity vectors. There is no donor correspondence constraint or cross-view activity covariance penalty.

Every shared term can be represented by a private factor with \(u_{ck}=h_{ck}b_{ck}\) and \(v_{ck}=g_k\). Private factors can also duplicate or cancel a shared contribution. Consequently, likelihood reconstruction alone does not determine shared identity. The structural prior/penalty, complexity rule, and ambiguity treatment must resolve or report these cases.

GS-00 must specify:

- Comparable costs for a common template, member activity vectors, private copies, memberships, and factor existence. All costs must have defined scale behavior and be included in the reported objective.
- Per-view total rank limits, including shared contributions and private factors, respecting residual rank and degrees of freedom. State any overcomplete representation explicitly rather than treating it as identifiable.
- Rules for candidate birth, deletion, merger, privatization, duplication, and cancellation; single-member candidates cannot remain reported shared programs.
- The behavior when all memberships are zero. Distinguish mathematical nesting of a private-only model, equality to the new backend's private-only implementation, and numerical equality to existing independent EBMF. Only claim the equalities actually implemented and checked.
- An ambiguity rule for observationally equivalent factorizations and a deterministic treatment of ties. Do not enforce universal shared/private orthogonality merely to obtain uniqueness.

### 4.2 Complete optimization criterion

If precisions are learned, the Gaussian likelihood must include its normalization terms. An example schematic negative log-likelihood is

\[
\mathcal L(\Theta,\tau)=
\frac12\sum_{c,i,j\ \mathrm{observed}}
\left\{\tau_{cij}\big[X_{cij}-\widehat X_{cij}(\Theta)\big]^2-\log\tau_{cij}\right\}
+\mathcal P_s+\mathcal P_p+\mathcal P_h.
\]

This is a schematic criterion, **not the frozen objective**. Define the precision structure, positivity/variance safeguards, nuisance model, view weighting, and every penalty before coding updates. If precision is fixed, record how it is obtained from training data and why the corresponding omitted terms are constant.

Covariate residualization affects residual degrees of freedom and can induce dependence; the simulator and calibration must follow the implemented preprocessing. Do not claim nominal inference from an unverified independent-noise approximation.

Fixed-penalty v1 is called selective coupled factorization. Using EBMF for initialization does not make its new objective empirical Bayes. Existing `flashier` ELBOs and the coupled criterion are not interchangeable ranking quantities.

### 4.3 Candidate search and optimizer acceptance

Generate candidates from coherent loading consensus and from individually detected source factors. A candidate can enter recruitment without already being a reported shared program. Specify the limitation when no view detects a program; residual-based births are permitted only if derived, budgeted, and calibrated as part of the complete search.

Test a source-derived template in target expression while allowing a competitive private explanation. A fully fitted private model can absorb the true target signal, so specify whether recruitment replaces/refits a private factor rather than relying only on a fixed residual. Refit affected parameters after membership changes.

For each update and structural move, define the criterion, acceptance tolerance, iteration limit, restart rule, and rollback behavior. Require a finite objective trace and safeguards against unacceptable decreases in a maximization criterion or increases in a minimization criterion. Heuristic proposals must be checked against the declared criterion before acceptance.

Statuses distinguish convergence, iteration exhaustion, numerical failure, and fallback. Never treat a crash or an untested selection as a successful zero-shared fit.

A row-stacked representation with common gene factors and view-restricted activity support provides a small joint feasibility comparator. Existing `flashier` initialization and fixed-entry facilities may help; verify their actual behavior at pinned versions. They do not automatically implement adaptive memberships, comparable private costs, independent view priors, or the complete coupled objective.

### 4.4 Conditional flexible and empirical Bayes extensions

Only after exact-core confirmation, consider

\[
w_{ck}=a_{ck}g_k+\delta_{ck},\qquad g_k^\top\delta_{ck}=0.
\]

Orthogonality separates amplitude from directional deviation but does not establish global factor identity. Gate the entire contribution by membership, and define inactive deviations. Begin with fixed shrinkage, keeping private exits and zero sharing possible.

Empirical prior learning is a later gate. Specify score, core, deviation, membership, and factor-existence priors together with their hyperparameters and full variational objective. If unit-norm constraints are used, derive how they interact with priors and EBNM updates; ordinary independent-prior updates cannot simply be followed by normalization without checking the objective. Actual posterior inclusion probabilities require calibration and are not interchangeable with selection stability or fixed-penalty scores.

## 5. Simulation and truth validation

### 5.1 Fidelity tiers

1. **Gaussian parametric pseudobulk:** exact loading directions, memberships, activity distributions, private factors, noise, and SNR. This is the first implementation and calibration tier.
2. **Count-based simulation:** reference-informed negative-binomial or another justified model, varying donor cell counts, gene means/dispersion, coverage, and cell-type prevalence; pseudobulk and preprocess through the actual pipeline. Choose one suitable simulator before adding multiple dependencies. Candidate packages include scDesign3 and Splatter.
3. **Semi-synthetic real backgrounds:** inject programs into COMBAT/Stephenson-derived data with documented injection coordinates and amplitudes. Donor permutation does not destroy real conserved gene-loading structure. Score planted recovery and transfer without declaring all additional discoveries false.

For count-based simulations, distinguish latent log-mean effects from loading directions after aggregation and normalization. Validate their relationship rather than assuming the planted coefficients remain exact Gaussian loading truth.

### 5.2 Scenario families

| ID | Loading structure | Primary activities | Evaluation |
|---|---|---|---|
| G0 | Private factors only, including overlapping supports and correlated backgrounds | Independent | Global structural-null false sharing; valid zero-shared results |
| G1 | Exact signed core in all eligible views | Independent | Global sharing and loading recovery |
| G2 | Exact signed core in a selected subset | Independent | Partial membership and false recruitment into nonmembers |
| G3 | Core plus controlled directional deviations | Independent | Exact-core misspecification initially; flexible-model recovery later |
| G4 | Same support with unrelated signs/weights | Independent | Separate support reuse from weighted conservation |
| G5 | Different gene directions | Coordinated in this deliberate sensitivity scenario | No gene-sharing call based only on score coordination |
| G6 | Correlated but distinguishable programs under the frozen identity rule | Independent | Overmerging and ambiguity handling |
| G7 | Common technical loading over private biology | Independent or coordinated | Statistical recurrence plus nuisance qualification |
| G8 | Split representations, overlapping programs, and selected rotational ambiguities | Independent | Split/merge and subspace recovery; unique identity only where identifiable |
| G9 | True shared program weak in a member view; strong source and genuine nonmembers included | Independent | Recruitment benefit and false target inclusion |

For G9 include a case with only one independently detected source factor. Do not restrict the test to programs already detected in two strong views. Include mixed datasets with true shared programs and private near-neighbors for membership calibration.

Initial design proposals: 4/8 views, 50/100 observations per view, 2,000 genes, 0-3 shared programs, 1-3 private programs per view, and 30-100 active genes. Stress tests include 2-20 views, 20-300 unbalanced observations, 1,000-10,000 genes, larger ranks/supports, correlated backgrounds, technical confounding, and varying precision. These are proposals, not a Cartesian grid or tuned defaults.

Specify missing donor profiles separately from entirely unavailable views and unavailable features. Version 1 supports unequal donor coverage; unsupported gene panels are boundary tests, not silently supported performance scenarios.

### 5.3 Simulator assertions

Before looking at method results, test exact membership/eligibility masks, normalized loading directions, support overlap, SNR, noise structure, panel order, missingness, and fixed-seed reproducibility. Assert no planted common technical direction in an unadjusted strict G0 null.

Generate primary activities independently across views, except G5 and explicit coordination sensitivity arms. Finite-sample accidental correlations are expected. Assess their distribution across replicates with sample-size-appropriate tolerances; do not reject seeds until every observed correlation is artificially small.

Existing `sim_param.R` S7 conserves signed weighted loadings, while `simulate_pseudobulk('same_genes_independent')` redraws weights on common support. Preserve their original labels and estimands. The existing parametric S0 includes a common batch loading and is not an unadjusted gene-sharing structural null.

## 6. Competitors, ablations, and common reporting

| Label | Method | Purpose and reporting requirement |
|---|---|---|
| B0 | Existing independent `flashier` EBMF | Preserve reference fitting. Local loading/private recovery is native; shared-call metrics require an explicit adapter or are not applicable. |
| B1 | Independent EBMF plus signed, coherent consensus | Strong custom post-hoc comparator; full-loading similarity, sign alignment, calibrated rejection, and unknown shared count. |
| B1-R | B1/source candidates plus fixed-template target recruitment | Essential comparator: isolate recruitment from jointly updating templates. Candidate generation/search budget should match B4 where possible. |
| Bstack | Small row-stacked EBMF feasibility model with explicit membership reporting | Signed Gaussian joint reference. Document structural assumptions and limitations; do not assume equivalence to B4. |
| B2 | LIGER iNMF | Nonnegative common plus dataset-specific weights; define extraction and membership from evidence, not merely existence of its common W. |
| B3 | MOFAcellulaR/MOFA2 in its usual donor-oriented representation | Different native estimand: shared donor activities. Define loading-recurrence adapter and use only compatible donor tracks. |
| B4 | Exact-core selective coupled factorization | Proposed joint estimator; compare with B1-R as well as B1. |
| B4-0 | Same new backend with shared memberships disabled | Check the new objective's private-only nesting and effect of joint structural selection. Does not automatically equal B0. |
| B5 | Flexible coupled factorization with fixed shrinkage | Conditional deviation ablation. |
| B6 | Flexible hierarchical empirical Bayes | Conditional adaptive-prior contribution after B5 succeeds. |
| Reference oracle | Memberships known, weights unknown | Measure subset-selection cost. Treat as a reference unless optimization/model conditions establish an upper bound. |

B1 is GeneNMF-inspired, not an unchanged run of GeneNMF. Specify deviations from its nonnegative workflow, fixed metaprogram count, specificity weighting, clustering, and normalization. Avoid giving one view extra consensus weight merely because it contributes more restarts or factors. Freeze whole-cluster coherence and consistent sign alignment; high absolute pairwise similarities alone do not define a coherent common template.

For LIGER distinguish the shared W from complete view-specific W+V_c, and define activity-supported membership. Signed Gaussian centered inputs are not directly compatible with its nonnegative likelihood. Use defensible method-specific/count-based inputs and map evaluation coordinates back; unsupported tracks are marked not applicable, not scored as failures or zero discoveries.

Standard MOFAcellulaR links observations through donor activities. Do not manufacture pairing in an unpaired benchmark. A transposed group-factor/MOFA representation with genes as the matched dimension is a possible bounded feasibility study, explicitly labeled an adaptation until audited and validated. Prefer a small relevant joint reference over proliferating external methods.

Across applicable methods use common biological sample filters, target measured genes, evaluation partitions, seeds, tuning budgets, rank accounting, and failure/coverage records. Include both a controlled common-input track where feasible and a clearly labeled method-appropriate preprocessing track. Runtime includes candidate generation, tuning, restarts, and selection, not just one final fit.

## 7. Evaluation, calibration, and prediction

### 7.1 Primary endpoints and matching

Primary efficacy is unconditional recovery of true weak shared programs/members, comparing B4 with B1-R on paired datasets. Separate loading recovery from membership recovery. Define any detectable-program stratum using an independent prespecified criterion or oracle; never define detectability by the tested method's successful discoveries.

Required reporting includes global-null false-sharing frequency, false member recruitment in mixed truth, shared-program recall, signed-loading recovery, private-program preservation, duplicate/split/merge rates, gene-support recovery, convergence/failure/coverage, runtime, and peak RAM. Numerical improvement and noninferiority margins are pending GS-04.

Use optimal global assignment with explicit dummy unmatched costs and admissible-match thresholds. Fix the truth-matching rule independently of method results. Define how a missed true program, an unmatched predicted program, and indeterminate membership affect each denominator. Do not compute membership accuracy only among recovered programs. Keep duplicate predictions as errors; do not allow all split fragments to receive full recall credit.

Evaluate ambiguous factors using the prespecified subspace metric. Unique factor identity should be scored only when the simulator makes it identifiable. Gene-support thresholds are required because a posterior-mean loading need not be exactly zero.

Truth matching belongs only to evaluation and may use evaluation truth. It must not influence fitting, tuning, candidate choice, or selection thresholds. This differs from training-derived alignment used for prediction or external projection.

### 7.2 Selection calibration

The proposed operating target is approximately 5% global structural-null probability of at least one reported shared program. This is neither a current guarantee nor a complete partial-membership error criterion.

Before confirmation, freeze the hypothesis family, nominal target, tolerated departure, confidence-bound acceptance rule, independent replicate counts, and evaluation coverage guardrail. Calibrate the complete adaptive procedure, including candidate search, recruitment, model selection, and reporting. Fresh confirmation datasets are required after selecting thresholds on a calibration pilot.

For planning, 500 independent null datasets at an event rate of 5% yield an approximate standard error of one percentage point. Determine replicate counts from desired interval width and paired efficacy power rather than adopting this example automatically. Do not pool stress scenarios in a way that hides a poorly calibrated subgroup.

Mixed-truth datasets separately measure incorrect recruitment into nonmembers, false new programs, private-to-shared misclassification, and overmerging. Posterior probabilities, if eventually implemented, need their own calibration; Bayesian probability thresholds and frequentist family-wise control are not identical guarantees.

Record attempted, succeeded, failed, and unevaluated runs with reasons. Report conditional estimates among evaluated runs plus coverage. Untested runs are neither true negatives nor silently discarded successes. A missing coverage result fails an adoption gate. Rerun policy, timeouts, and paired-analysis handling must be preregistered.

### 7.3 Leakage-resistant projection and transfer

Split by biological donor where observations repeat, keeping all measurements from one donor in one outer fold across views. Independent sampling across views does not justify putting the same donor in training and testing in different views.

Fit feature selection, nuisance regression, centering/scaling, precision estimates, candidate generation, priors, rank, and reporting thresholds within training data. Freeze those quantities for the outer test set.

Estimate held-out activities using an observed gene subset; evaluate reconstruction on a disjoint masked subset whose loading weights were learned from training donors. Private activities and nuisance estimates cannot use masked expression. Genes never observed in any training donor require an additional loading model and are outside ordinary masked-gene reconstruction.

For source-only transfer, learn the template without the target view, estimate target activities from observed target genes, and score held-out target expression. If target training data update the template, label that as target adaptation and evaluate it separately. Specify projection shrinkage and preserve the training scale convention.

Supplement reconstruction with independent recurrence evidence. An exact-core constraint forces participating loading directions to match; that fitted equality cannot validate the membership that imposed it.

## 8. Compatibility repairs and package boundaries

| Area | Planned change and acceptance |
|---|---|
| `R/mofacellular.R` | Align equal gene panels by identifiers or reject unsupported panels. Test reversed order, equal-size different IDs, unequal counts, partial overlap, missing observations, and raw-count versus transformed usage. Never pass transformed data through raw-count TMM. |
| `R/stability.R` | Preserve legacy unconditional recovery for compatibility; add attempted/succeeded counts, failure reasons, and successful-fit conditional recovery. All-failed conditional recovery is NA. Specify reference-seed and zero-factor behavior. |
| `R/utils.R` / `R/principal-angles.R` | Use optimal matching in the new evaluator with threshold/unmatched handling. Production replacement of legacy greedy stability matching requires its own decision and invariance review. Correct claims that greedy matching is Hungarian assignment. |
| `R/data-model.R` / `R/ebmf.R` | Keep permissive independent input panels and existing fit defaults. Coupled validation/preprocessing has a separate contract and provenance; do not reuse differently selected Y_used matrices as if already harmonized. |
| Coordination helpers and test infrastructure | New gene-sharing clustering has explicit coherence. Do not rewrite coordination estimands. Add `tests/testthat.R` and dependency-enabled package validation before release; separate dependency-free adapter tests from fit-dependent tests. |

Use the existing data object and extraction conventions where suitable. The proposed fitter can receive `cell_program_data` plus an optional independent fit for initialization; an input fit also carries pseudobulk, but the source matrix and reprocessing contract must be explicit. Do not restrict recruitment to genes retained by local feature selection when the declared joint panel includes others.

Use an experimental name such as `fit_coupled_programs(..., model = 'exact_core')` until the EB objective exists. Preserve current public functions. Zero sharing and private exits are required behavior, not optional switches that can be disabled to force a result. Delay exports/classes until the result contract and evidence gate pass.

Outputs retain ordered templates, per-view loading/activity matrices, private factors, computational memberships, reported evidence states/reasons, eligibility/masks, fitted transformations, objective traces, restart/selection diagnostics, conditional and unconditional stability, and data/software provenance. Use distinct global-template IDs and local factor IDs so current `<cell_type>_<k>` identifiers remain interpretable.

Keep Monte Carlo calibration/recovery tests in benchmark scripts. Package tests should check objective calculations, update acceptance, dimensions, reconstruction, masks, projection leakage, and deterministic edge cases; passing a few synthetic tests is not a calibration claim.

## 9. Progress summary and implementation steps

No implementation step is complete. GS-00 is the next bounded specification task. GS-05 and later depend on the completed design and evidence gates, not merely elapsed time or successful compilation.

| ID | Outcome | Status |
|---|---|---|
| GS-00 | Freeze mathematical, reporting, and selection contract | not started |
| GS-01 | Establish isolated baseline and interoperability evidence | not started |
| GS-02 | Implement and verify Tier 1 simulator and evaluator | not started |
| GS-03 | Implement strong consensus/recruitment and joint feasibility comparators | not started |
| GS-04 | Pilot feasibility; freeze new preregistration and confirmation rules | not started |
| GS-05 | Implement exact-core optimizer | blocked: GS-00 through GS-04 |
| GS-06 | Run independent exact-core confirmation and decide continuation | blocked: GS-05 |
| GS-07 | Add fixed-shrinkage loading deviations if warranted | blocked: GS-06 benefit gate |
| GS-08 | Derive adaptive hierarchical EBMF if warranted | blocked: GS-07 benefit gate |
| GS-09 | Validate count/semi-synthetic realism and external cohorts | blocked: viable simpler method |
| GS-10 | Integrate and document the validated package API | blocked: applicable scientific and software gates |

### [ ] GS-00 — Freeze the model and selection contract

- **Status:** not started.
- **Outcome:** a source-reviewed specification sufficient to implement a single objective and unambiguous evaluation/reporting rules.
- **Actions:** resolve section 2 blockers; write the full likelihood/penalties, coordinate and eligibility rules, rank limits, membership evidence, birth/recruitment search, private-only behavior, and optimizer acceptance rules. Define the operating identity distinction for G3/G6 and ambiguity handling for G8. Separate candidate generation from discoveries.
- **Dependencies:** this plan and existing baseline contracts; specification work requires no production code changes.
- **Affected areas:** this document or a linked project-local model contract; planned gene-sharing benchmark specification.
- **Validation:** work through zero sharing, one source view, zero-amplitude memberships, equivalent private copies, cancellation/duplication, and view/gene order permutations. Derive scale behavior and precision updates. No unresolved choice may be left for the optimizer author to invent.
- **Evidence:** pending reviewed contract, explicit remaining decisions, and acceptance checklist.

### [ ] GS-01 — Establish baseline and interoperability evidence

- **Status:** not started.
- **Outcome:** a reproducible source-loaded baseline plus safe coupled/competitor inputs.
- **Actions:** confirm whether existing arrays remain active before editing their source; create an isolated implementation checkout; pin dependency versions/source commits and R environment; run dependency-enabled tests; freeze a representative real-data or strongest-fixture baseline artifact. Implement the bounded adapter/stability repairs in section 8 after the runtime-source boundary is safe. Add a package-check test entry point.
- **Dependencies:** GS-00 data/reporting contract and a later implementation request. Do not install dependencies as part of this documentation-only task.
- **Affected areas:** `R/mofacellular.R`, `R/stability.R`, relevant tests, `tests/testthat.R`, environment/source manifest, baseline verification artifacts.
- **Validation:** adapter gene/value round trips, partial-panel rejection/alignment, all-failed stability diagnostics, dependency-enabled existing suite, and before/after baseline numerical comparison. Label fixture evidence if real-data comparison is unavailable.
- **Evidence:** pending test logs, source/environment hashes, baseline invariance receipt, and repair review.

### [ ] GS-02 — Build the Tier 1 simulator and evaluator

- **Status:** not started.
- **Outcome:** reusable Gaussian truth generation and method-independent scoring.
- **Actions:** implement G0-G9 as scenario families with explicit masks and independent primary activities; start with compact settings rather than a full grid. Add mixed-truth near-neighbor cases and truth self-checks. Implement optimal assignment, unmatched/duplicate penalties, membership denominators, private recovery, and ambiguity/subspace scoring.
- **Dependencies:** GS-00 definitions and GS-01 reproducible environment/input boundary.
- **Affected areas:** `benchmarks/gene_sharing/simulations/parametric.R`, `evaluate.R`, simulator/evaluator self-checks. Reuse existing utilities only where their semantics match.
- **Validation:** truth assertions before any model comparison; hand-calculated assignment cases including a greedy counterexample; missed/extra/split/merged prediction fixtures; permutation and sign invariance.
- **Evidence:** pending simulator manifest, invariant logs, evaluator fixtures, and reproducible seeds.

### [ ] GS-03 — Establish strong simpler comparators

- **Status:** not started.
- **Outcome:** B0/B1/B1-R and a bounded signed joint feasibility reference with common reporting.
- **Actions:** implement full-loading consensus, balanced view contributions, sign-consistent coherence, and no-sharing selection. Implement singleton-source fixed-template recruitment with independent target evidence. Check a small stacked-EBMF formulation before writing a custom optimizer. Specify external comparator adapters and unsupported tracks without yet requiring a full external benchmark matrix.
- **Dependencies:** GS-00 contract, GS-01 baseline, and GS-02 truth/evaluation checks.
- **Affected areas:** `benchmarks/gene_sharing/methods/`; only stable reused package helpers where separately authorized.
- **Validation:** strong positives, G0 no sharing, G4/G5 counterexamples, G9 weak recruitment, mixed-truth nonmember rejection, and coherence/duplicate cases. Validate candidate-search budgets and loading coordinate mapping.
- **Evidence:** pending comparator implementation/source manifest and pilot results explicitly labeled exploratory.

### [ ] GS-04 — Pilot and freeze the new benchmark

- **Status:** not started.
- **Outcome:** a separate tagged preregistration with executable acceptance rules and measured compute requirements.
- **Actions:** use discarded pilot/calibration seeds to choose a compact first design, tune selection, estimate paired variance/power, and measure cost. Freeze null tolerance/CI rule, membership/private noninferiority margins, efficacy margin, failure/coverage threshold, tuning/restart budget, confirmation seeds, and rerun policy. Reserve independent confirmation seeds before inspecting B4 confirmation results.
- **Dependencies:** GS-02/GS-03 feasible outputs; GS-00 fully resolved. B4 feasibility work after the complete objective may inform a documented pilot but cannot become confirmation evidence.
- **Affected areas:** `benchmarks/gene_sharing/PREREG.md`, `DEVIATIONS.md`, seed/task manifests, pilot cost/decision receipts.
- **Validation:** final rules are computable from evaluator outputs, including failed/missing runs; no TBD margin in registered endpoints. Review source snapshots and training-only preprocessing. Keep original `prereg-v1` and its deviations unchanged.
- **Evidence:** pending preregistration tag/commit, discarded-seed ledger, power/interval rationale, and measured resource estimate. Tagging or publishing requires an appropriate later action request.

### [ ] GS-05 — Implement exact-core selective coupled fitting

- **Status:** blocked until GS-00 through GS-04 are accepted.
- **Outcome:** objective-correct fitting with variable membership, private exits, and zero sharing.
- **Actions:** implement activity/core/private updates and accepted structural moves; recruit into missed views, refitting competing private contributions as specified; track objective and convergence; retain the independent-only candidate and multistart diagnostics. Use the smallest implementation that satisfies the contract.
- **Dependencies:** frozen objective and selection contract, validated simulator/evaluator, simpler comparators, and pilot/confirmation separation.
- **Affected areas:** initially benchmark-local coupled implementation and focused mathematical tests; public `R/` integration waits for GS-10.
- **Validation:** independently computed tiny objective, accepted-move behavior, reconstruction and scale/sign invariance, zero-sharing/private-only nesting as actually promised, candidate deletion, weak-view recruitment, duplicate/cancellation handling, and failure propagation. Check invariance to view order up to labels and declared numerical tolerances.
- **Evidence:** pending objective/update receipts and deterministic self-checks; none constitutes statistical calibration.

### [ ] GS-06 — Confirm incremental exact-core benefit

- **Status:** blocked until GS-05.
- **Outcome:** an adopt/stop decision on independent confirmation datasets.
- **Actions:** compare B4 against B1-R, B1, B0 recovery, B4-0, and the applicable joint reference. Prioritize G0/G1/G2/G4/G5/G9, with G6 overmerging checks and G3 misspecification sensitivity. Run the frozen compact design, paired analysis, subgroup reporting, and full cost accounting. Extend to applicable B2/B3 tracks when inputs/adapters are verified.
- **Dependencies:** GS-04 registered rules and GS-05 checks; scheduler execution requires a later authorized task.
- **Affected areas:** new benchmark runner, immutable manifests/results, `summarize.R`, report and decision receipts.
- **Validation:** complete attempted/success/failed coverage; null and mixed-membership gates; efficacy against B1-R; private preservation; runtime/RAM; no confirmation-seed tuning.
- **Evidence:** pending uncertainty intervals, coverage/cost table, and explicit continuation decision. If the exact-core model fails the practical-benefit gate, stop custom hierarchical expansion; report the simpler successful method or the negative result.

### [ ] GS-07 — Test fixed-shrinkage flexible loadings

- **Status:** blocked on a successful GS-06 gate.
- **Outcome:** incremental evidence for controlled loading deviations.
- **Actions:** derive the fixed-shrinkage objective and parameter conventions, including inactive deviations; implement B5 only after the derivation; evaluate G3/G6/G9 and null/mixed truth with fresh confirmation data. Treat heterogeneous panels as a separate model/input extension rather than enabling them incidentally.
- **Dependencies:** GS-06 benefit, exact-core reproducibility, and a separate frozen flexibility comparison.
- **Affected areas:** benchmark-local deviation model, mathematical tests, extension preregistration and report.
- **Validation:** objective-correct updates, amplitude/deviation orthogonality as defined, private exits, and improvement without unacceptable null or overmerging degradation.
- **Evidence:** pending fixed-shrinkage benefit and robustness receipt. Failure stops GS-08.

### [ ] GS-08 — Derive and evaluate adaptive hierarchical EBMF

- **Status:** blocked on a successful GS-07 gate.
- **Outcome:** a complete empirical Bayes model with incremental evidence beyond fixed shrinkage.
- **Actions:** define prior families/hyperparameters, variational moments, scale constraints, and factor/membership existence treatment; derive the coupled ELBO and updates; implement only the justified adaptive components. Evaluate probability calibration, prior sensitivity, and practical cost against B5.
- **Dependencies:** GS-07 benefit and an accepted full probabilistic derivation; exploratory hierarchical ideas alone are insufficient.
- **Affected areas:** experimental EB model/fit/selection code, objective tests, independent extension benchmark.
- **Validation:** ELBO and posterior-moment calculations against independent tiny cases; update acceptance; prior/scale behavior; probability calibration and incremental recovery/cost gates.
- **Evidence:** pending derivation, probability calibration, and comparison receipt. Until then do not label v1 fixed-penalty membership scores posterior probabilities or call B4 a validated EB extension.

### [ ] GS-09 — Establish realism and external validity

- **Status:** blocked until a simpler method is viable; does not require GS-08 to succeed.
- **Outcome:** count/semi-synthetic transfer evidence and independently replicated real-data programs.
- **Actions:** implement one justified count simulator; validate transformed truth; use real-background injections with limited truth claims. Analyze COMBAT/Stephenson exploratorily, then verify cohort overlap and data availability before treating HLCA/OneK1K or other cohorts as external validation. Freeze templates and mapping rules for projection; screen technical associations and small-view precision.
- **Dependencies:** viable GS-03 or GS-06 method, appropriate data access, preprocessing/projection contract, and later execution request.
- **Affected areas:** count/semi-synthetic generators, real-data scripts, manifests, projection tests, external validation report.
- **Validation:** leakage-resistant held-out and source-only transfer; technical-sensitivity results; cross-cohort loading/subspace recurrence; cohort/donor independence; measured resource requirements. Biological causality is not inferred from expression recurrence.
- **Evidence:** pending cohort provenance, replication/transfer intervals, nuisance analyses, negative results, and limits.

### [ ] GS-10 — Integrate the supported package API

- **Status:** blocked on applicable scientific and software gates.
- **Outcome:** a documented experimental or validated API with explicit supported assumptions and preserved baseline behavior.
- **Actions:** expose only the method variant whose gates passed; reuse data/extraction conventions; add exports, result class, error handling, projection helpers, provenance, examples, and documentation. Keep unsupported flexible/EB variants unexported. Add reproducible environment and package-check evidence.
- **Dependencies:** successful method-specific evidence, GS-01 compatibility repairs, stable result contract, and later implementation/release request.
- **Affected areas:** minimal `R/consensus-programs.R` / `R/coupled-programs.R` or existing modules as justified, `NAMESPACE`, `DESCRIPTION`, documentation, package tests and benchmark report. Split additional files only when implementation warrants it.
- **Validation:** source-loaded dependency-enabled tests, package-check test execution, baseline artifact comparison, zero-factor and masked/unobserved outputs, source-only projection, examples, and result serialization/provenance. Report skipped optional integrations explicitly.
- **Evidence:** pending package check, compatibility receipt, user documentation, and supported-method status. Production-ready claims require the complete applicable gate set; experimental release labels remain explicit.

## 10. Expected files and verification commands

Create files only as their steps require them. Do not scaffold the original six hierarchical modules before the model has passed its gates.

```text
docs/plans/2026-10-09-selective-coupled-ebmf-gene-sharing.md
benchmarks/gene_sharing/
  PREREG.md
  DEVIATIONS.md
  simulations/parametric.R
  methods/                    # only implemented/verified competitors
  evaluate.R
  run_rep.R
  summarize.R
  report.R
  test_simulator.R
  test_evaluator.R
tests/testthat.R              # package-check entry point
tests/testthat/
  test-mofacellular.R
  test-extract-stability.R
  test-coupled-objective.R    # only after GS-05
  test-coupled-projection.R   # only after projection implementation
```

Count/semi-synthetic generators and scheduler/workflow files are added later when authorized and needed. Benchmark results need immutable task/source/seed manifests, a reproducible completion check, and receipts that distinguish pilot, confirmation, software tests, and scientific evidence. Use the project's established runner where appropriate; if a new dependency workflow is needed, follow the user's Snakemake preference rather than building an ad hoc orchestrator.

Use R4_51 for R work. The review located:

```text
/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/Rscript
```

Current source-loaded package-test command:

```bash
/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/Rscript --vanilla -e 'testthat::test_local(".", reporter="summary")'
```

Current dependency-free coordination checks, useful as compatibility smoke evidence:

```bash
/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/Rscript --vanilla benchmarks/coordination/test_sim_coord.R
/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/Rscript --vanilla benchmarks/coordination/test_methods.R
```

Planned new simulator/evaluator checks are not runnable until their files exist. Package build/check commands and required dependency versions must be recorded in GS-01. Do not count a zero exit status from a test run configured with `stop_on_failure = FALSE` as a passing suite.

## 11. Final acceptance and exclusions

Acceptance requires mathematical validity, zero sharing as a supported successful outcome, calibrated structural-null and mixed-membership behavior, meaningful improvement over B1-R, private preservation, adequate failure/evaluation coverage, robustness, and justified compute cost. External and package-readiness claims require their additional gates.

Adversarial checks include unsupported initialized shared candidates; a view with no identifiable program; singleton/zero-amplitude memberships; duplicate and cancelling cores; factor sign reversal; view and gene-order permutations; missing views and unsupported panels; weak source/target recruitment; coordinated activities with unrelated loadings; a shared nuisance signature; ambiguous subspaces; and optimizer failure. Their expected outputs must be specified before execution.

Exclude supervised outcome-based discovery, mandatory donor correspondence, coordination/covariance optimization, arbitrary cross-omics likelihoods, exhaustive subset enumeration at large view counts, unmeasured-gene zero filling, and retroactive modification of the existing registered estimand. Do not add flexibility solely because it can be implemented.

If B4 fails to improve over calibrated fixed-template recruitment, stop before full hierarchical optimization. A negative result is a completed scientific outcome. If B4 succeeds, describe it as selective coupled factorization; only successful GS-08 evidence supports the hierarchical empirical Bayes claim.

## 12. Dated evidence and decision log

- **2026-10-09 — Review:** current checkout commit and tree verified. Package functions, relevant tests, model-selection/coordination simulators, preregistration/deviation documents, and primary method sources inspected. This is source evidence, not numerical validation of a coupled method.
- **2026-10-09 — Dependency-free R4_51 checks:** reproduced reversed gene-order corruption in the MOFA adapter; a unit-vector example where greedy assignment totals 1.113341 versus optimal 1.732051; chain merging without an A-C edge; and all failed stability refits recorded as zero recovery without success counts. Verified an exact shared contribution has an identical private-only reconstruction.
- **2026-10-09 — Existing smoke tests:** `benchmarks/coordination/test_sim_coord.R` and `test_methods.R` completed successfully. The latter warned that 199 permutations are insufficient for its 45-pair sparse-signal BH resolution. These are smoke/regression checks, not FWER calibration.
- **2026-10-09 — Package test attempt:** R4_51 4.5.1 had `pkgload` and `testthat`, but lacked `flashier` and `ebnm`. `testthat::test_local(..., stop_on_failure = FALSE)` produced 15 skips and one dependency error in the degenerate-view fitting test. No full passing suite or baseline-fit invariance was established. The testthat build-version warning and two small-view/zero-variance warnings were also observed.
- **2026-10-09 — Plan saved:** review incorporated into this project-local progress document and linked in `HANDOFF.md`. No package/benchmark code changed, dependencies installed, jobs submitted, or scientific stages executed in this documentation task.

Append future entries with date, step ID, exact command/source versions, artifact path/hash, result, coverage, and remaining limitations. Record deviations without rewriting pilot history or completed evidence. Historical handoff queue counts and earlier package-test reports are not current validation for this plan.

## 13. Primary sources and version requirements

Local compatibility and estimand sources: [EBMF backend](../../R/ebmf.R), [data model](../../R/data-model.R), [MOFA adapter](../../R/mofacellular.R), [stability](../../R/stability.R), [matching and geometry](../../R/principal-angles.R), [utilities](../../R/utils.R), [package simulator](../../R/simulation.R), [model-selection preregistration](../../benchmarks/model_selection/PREREG.md), [deviations](../../benchmarks/model_selection/DEVIATIONS.md), [parametric generator](../../benchmarks/model_selection/sim_param.R), [coordination workstream](../../benchmarks/coordination/PREREG.md), and [roadmap](../../ROADMAP.md).

Method sources:

- [Wang and Stephens, Empirical Bayes Matrix Factorization, JMLR 2021](https://jmlr.org/papers/v22/20-589.html); [flashier objective](https://github.com/willwerscheid/flashier/blob/master/R/objective.R) and [fixed factor entries](https://willwerscheid.github.io/flashier/reference/flash_factors_fix.html).
- [GeneNMF source, master branch](https://github.com/carmonalab/GeneNMF/blob/master/R/main.R). The original plan's `main` link did not resolve during review.
- [LIGER iNMF objective and input contract](https://welch-lab.github.io/liger/reference/runINMF.html).
- [MOFAcellulaR multicellular vignette](https://github.com/saezlab/MOFAcellulaR/blob/main/vignettes/get-started.Rmd) and [MOFA2 FAQ](https://biofam.github.io/MOFA2/faq.html).
- [Bayesian group factor analysis with structured sparsity](https://www.jmlr.org/papers/v17/14-472.html); count-simulator candidates [scDesign3](https://github.com/SONGDONGYUAN1994/scDesign3) and [Splatter](https://github.com/Oshlack/splatter).

Live upstream references describe reviewed behavior and are not immutable dependency pins. GS-01/GS-04 must record exact commits/package versions actually used. Count simulators and transposed group-factor adaptations remain candidates, not source-audited implemented comparators.
