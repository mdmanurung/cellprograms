# C1: coordinated program detection across cell types

Status: prospective specification; registered only when `prereg-coord-v1` and the source manifest exist. The original unregistered text is retained in PREREG_initial_draft.md. The independently reviewed implementation plan is docs/C1_PLAN.md; review resolutions are docs/C1_PLAN_REVIEW.md. No main-evaluation outcome has been inspected.

## Target and inferential unit
Donors, not cells, are independent biological units. A coordinated program is a set of locally estimated factor nodes, at most one node per cell type, with shared or oppositely directed donor activity beyond the measured categorical nuisance model. Activity coordination is primary; gene reuse is secondary. This benchmark does not establish causality, general validity under unmeasured confounding, or external cohort performance.

## Candidate selection, pilot and references
The selected candidate is `node_cor`: the package's `find_coordinated_programs`, with joint donor max-T testing followed by complete, sign-consistent linkage. Every member pair must be eligible and called; incompatible cell-type/sign merges are rejected. Its called edges and all tested edges are retained. Edge p-values are not cluster p-values. Existing `pa_cc`, `subspace_sum`, and `gene_match` are descriptive comparators; no main-data method selection is allowed. The legacy node_cor helper in methods.R is retained only for existing regression tests; the runner explicitly uses the package implementation.

`pa_cc` tests cell-type subspaces and produces cell-type components with no identifiable program nodes. It is scored only on the membership tier, with node metrics unavailable. Node-level methods also receive membership-tier scores for a fair breadth comparison. Gene matching is a gene-usage comparator, not an activity detector; legacy inferential labels are not interpreted as validated cluster error control.

Pilot seeds are 90,000,000 + 10,000 * pilot-cell-index + replicate. Initial cells have three replicates (2–4), plus the local mid_small replicate 1. The extended pilot compares amplitudes 1 and 2 in both disjoint-gene and same-gene small-R20 cells, ten replicates each. These observations are discarded for validation. Choose the smaller level meeting pilot Found >= .80 for R20 and >= .90 for R100/R50 in both extended activity cells, pooled node precision >= .90 and spurious incidence <= .10; if none qualifies, extend the discarded pilot prospectively. The selected main amplitude, weak/strong amplitudes and runtime are recorded in PILOT.md and main_cells.csv before tagging. Numerical main gates below do not depend on main outcomes.

Ideal truth-score references use the same donor sets, categorical adjustment, significance level and candidate clustering, including private scores. They are algorithm-specific ideal references, not universal optimal detectability bounds. Their power is not tuned by expression amplitude. genes_only has no coordinated-activity truth; any activity call is false. Oracle runs bypass factor estimation and therefore separate local fitting losses from detection losses.

Prospective pilot extension, independently reviewed before any main run: scalar amplitudes 1 and 2 failed the same-gene/small-R20 selection thresholds. Test profiles (R100,R50,R20)=(1,1,4), then (2,2,4), each in both gene modes with ten fresh datasets per cell. Choose the first profile meeting the original per-mode Found, pooled-precision and spurious thresholds; do not pool modes. No detector, backend, matching or qualification gate changes. If neither qualifies, retain the failures and review another prospective change. Any eventual support is restricted to the selected signal regime, including its R20/R100 amplitude ratio. Original scalar 1 and 2 same-gene/small-R20 settings remain separate descriptive stress cells, in addition to weak signal. Pilot task parser failures after finalized exploratory outputs remain recorded as failures; structurally validated outputs may guide pilot selection but never count as main evidence.

## Data and fits
Ten cell types, 100 donors; randomly assigned four 100-, three 75-, and three 50-donor cell types. Small cell types share a randomly selected 30-donor core plus 20 independent additional donors. R20 members are explicitly selected from the small or large tiers. Missingness is independent of activity in this simulator; disease-related missingness is not covered.

2,000 genes; 40 genes per planted program; R100 in ten, R50 in five random, R20 in two specified-tier cell types; 1–2 private programs per cell type. Same, half-overlap and disjoint-gene variants are explicit. All activity and loading truths, including private and measured nuisance, are retained. At activity correlation .7, the shared latent coefficient is sqrt(.7), so pairwise population correlation is .7.

Measured nuisance is the standardized sum of Institute-level offsets and a second categorical donor variable varying within Institute. It has its own 40-gene loading support and amplitude 1.5. Both variables are residualized; permutations are blocked by their cross-classification. The raw unadjusted/unblocked null is a positive nuisance-detectability control. Primary residual noise is homoscedastic; heterogeneous donor scales are a descriptive stress condition with its own negative controls, not a qualified operating regime by default.

Fits are pinned: all genes; max_factors=10; point_laplace; var_type=1; center=TRUE, scale=FALSE, backfit/nullcheck=TRUE; covariate_mode=residualize; covariates Institute and nuisance_group in adjusted cells. Per-dataset fit seed=simulation seed+2,000. Every method uses the same fitted object. Fallbacks, retained ranks, runtime, fit settings, R/package versions and frozen source manifest are retained. Runtime is R4_51 (R 4.5.1); installed flashier/ebnm versions are saved in receipts.

## Inference and scope
Alpha=.05; 999 permutations; minimum 20 shared donors. Statistic is |atanh(clipped score correlation)| * sqrt(n_shared-3). Per draw, independently permute every cell type's donor rows within its named blocks, moving all factor columns together and retaining original donor identifiers; use that one dataset draw for all node pairs. Max-T adjusted edge probability is (1 + #null maxima >= observed statistic)/(1,000). The simulation and permutation seeds are distinct. The public API preserves the caller's RNG state.

Validity requires score exchangeability within blocks conditional on the observed donor sets. Categorical hard residualization commutes with those block permutations. Arbitrary continuous adjustment, soft/fixed adjustment and immobile blocks are rejected by the public interface. Heterogeneous residual scales and missingness related to activity need independent calibration. Global-null calibration is assessed empirically; strong mixed-null familywise control is not established by this benchmark. Mixed-alternative spurious calls are separately gated. Local-factor estimation uncertainty is not propagated.

Global constant/zero factors and pair-intersection-constant columns are ineligible. Too-small donor intersections are explicitly reported. Empty memberships are valid; unsupported inputs or failed computation raise errors, never become clean negative results.

## Truth matching and metrics
Within each cell type, match fitted scores to adjusted recurring/private truth scores on observed donors by absolute correlation >= .5. Choose a one-to-one assignment maximizing total eligible correlation; unmatched nodes are allowed. Truth columns and fitted IDs are sorted; a deterministic dynamic program chooses the first optimum on ties. Activity identity takes precedence over gene identity. Only unmatched adjusted nodes can match the removed nuisance by loading cosine >= .5, one nuisance node per cell type. Raw fits match raw nuisance activity in the same score assignment.

Assign called clusters one-to-one to activity-true R100/R50/R20 to maximize correctly matched member count; deterministic cluster and truth order resolves ties. pa_cc uses cell-type membership overlap, without invented node identities. Node metrics use correctly mapped member nodes. Recall=correct members/true breadth; precision=correct members/called cluster size; breadth error=called size minus true breadth. Found requires recall >= .5, precision >= .5 and at least two correctly matched members. Missing true nodes remain in recall denominators.

A cluster not assigned a successful Found match is conservatively spurious; unmatched duplicate fragments also count. Fragmentation is reported separately. A cluster containing at least two nuisance-matched nodes is a decoy. Report dataset incidence, not independent edge counts. Pooled node precision per activity cell is all correctly assigned called nodes / all called nodes; an entirely no-call cell has precision zero. Bootstrap complete datasets, not individual nodes.

## Main cells, seeds and uncertainty
Primary: null and genes_only, 200 replicates each; act+genes and act_only crossed with R20-small/large, 200 each; raw nuisance control, 200. Secondary: half-overlap, weak/strong signal, activity correlation .7, heterogeneous activity and heterogeneous null/genes_only controls, 100 each. All breadths are planted together in each activity dataset. Exact cells and amplitudes are in main_cells.csv; task manifest is tasks_main.txt.

Main seed=100,000,000 + 10,000 * main-cell-index + replicate. These seeds are disjoint from the pilot and M2. Permutation seed=main seed+3,000. Outputs are immutable; records are finalized by rename and existing records are never overwritten. Startup verifies frozen executable hashes before sourcing, each main record retains that manifest, and summarization rejects mismatches or missing/duplicate datasets.

4,000 dataset bootstrap draws, report seed 8,102,026. Fixed equal-weight activity-cell macro means use stratified resampling within each cell, paired across methods. Bootstrap intervals are two-sided percentile 95%. Report two-sided Clopper–Pearson intervals and the explicitly gated one-sided 95% bounds below. Both R50/R20 gates must pass; no best-cell choice or pooling of negative controls is allowed.

## Qualification: every gate must pass
- Each primary null and genes_only cell separately: dataset spurious-or-decoy incidence point <= .05 and one-sided exact 95% upper <= .10.
- Raw nuisance control: nuisance detection point >= .80 and one-sided exact 95% lower >= .70.
- Both R50 and R20: equal-weight macro Found gain over pa_cc membership tier >= .10 and paired lower > 0 across the four primary activity cells.
- R100 macro Found non-loss: paired lower >= -.05.
- Each primary activity cell: pooled node precision >= .80 and dataset-bootstrap lower >= .70.
- Both R50 and R20: ideal candidate reference macro Found >= .80; fitted macro Found >= .60 and lower >= .50. Hard replicates are never excluded on their oracle result.
- Each primary activity cell: spurious incidence <= .10 and paired upper of candidate-minus-pa_cc spurious incidence <= .03.

Only the prospectively selected detector can qualify. Descriptive comparator/stress results do not expand its support. If any gate fails, retain the negative result, independently review a prospective revision and use a fresh held-out seed block. Never weaken these gates after inspecting main data.

## Sources and completion
Permutation plus-one: Phipson & Smyth, 2010, https://pubmed.ncbi.nlm.nih.gov/21044043/. Max-T/strong-control assumptions: Westfall & Troendle, 2008, https://pmc.ncbi.nlm.nih.gov/articles/PMC3117234/. These support the stated limitations, not a proof for estimated pseudobulk factors.

C1 is done only after complete valid outputs, all qualification gates, usable package API/extraction/docs/examples, required checks, independent final source/evidence review, compatibility evidence and an accurate completion audit. A source tag or a passing synthetic unit test alone is not scientific qualification.
