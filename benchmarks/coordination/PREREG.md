# PREREG `prereg-coord-v1` — DRAFT (not registered): recurring programs across cell types

Goal: after EBMF per cell type on pseudobulk, find programs that **recur** in all cell types, in a fraction of them, or in only a few, and tell them apart from private programs and from donor-level nuisance. Not registered until the pilot margin is filled in and the tag `prereg-coord-v1` exists.

Status: **DRAFT** for user review (registration is a scientific commitment). The four
blockers recorded in HANDOFF 2026-10-11 are resolved here: method list (below),
replicates + paired inference (below), node-to-truth rule (`scoring.R`, ported from
`reports/biomni_transfer/part1_callers_scoring.md` section 3), and the pinned
post-Stage-A EBMF config (below). The one remaining TBD is the adoption margin,
fixed after the pilot and before the main run.

## Simulation (`sim_coord.R`; extends `R/simulation.R`)
- 10 cell types (so 20% = 2 and 50% = 5 are exact), 100 donors, per-cell-type donor sets of 50-100 with random missingness (as in COMBAT, 50-121).
- Planted per dataset:
  - R100: one program active in all 10 cell types.
  - R50: one program in 5 cell types, chosen at random each replicate.
  - R20: one program in 2 cell types, chosen at random.
  - P: 1-2 private programs per cell type.
  - N: a donor-level nuisance (Institute) acting on every cell type, shared genes. Decoy: measured and adjustable.
- Activity z is shared across the member cell types of a program (corr 1; robustness arm 0.7). Cell types outside the program carry no signal for it.
- **Factor 1, breadth:** 100%, 50%, 20% (above), all in the same dataset and also each alone.
- **Factor 2, gene reuse:** member cell types use the same genes / 50% overlapping genes / disjoint genes. Activity coordination is the primary target; same-gene reuse is the secondary one.
- **Factor 3, signal:** amplitude set per breadth class by a pilot so a pairwise oracle (true z correlation test) has power about 0.5; also one weaker and one stronger level.
- **Factor 4, small cell types:** members drawn from the small cell types (about 50 donors) vs the large ones, reported separately.

## Pinned post-Stage-A EBMF config (blocker 4)

Fits on simulated pseudobulk use exactly the Stage-A CMP-detection settings, never
changed between compared arms ("never change two factors" rule):

- `fit_celltype_programs(x, loading_prior = "point_laplace", max_factors = 30,
  center = TRUE, scale = FALSE, features = "all", var_type = 1,
  backfit = TRUE, nullcheck = TRUE, covariates = ~institute,
  covariate_mode = "residualize", seed = <replicate seed>)` then
  `canonicalize_programs()`.
- `features = "all"`: simulated genes are exchangeable; HVG selection adds no
  realism here and would confound the gene-reuse factor.
- The C4 residual-df guard never binds (100 donors, K = 30, q = 3).
- flashier/ebnm versions recorded in `fit$provenance` per run.

## Methods compared (`methods.R`; the registered list is the implemented list — blocker 1)

- `pa_cc`: current pairwise principal-angle test (`sharing_spectrum`, scores space, `strata`) then connected components on called pairs. Cell-type-level (`program_id = NA`); the lenient baseline.
- `node_cor`: node-by-node score correlation (studentized, |T| = |atanh(r)|·sqrt(n−3)) with a within-strata permutation max-T over all cross-cell-type edges, then greedy single linkage.
- `node_agg`: node-level evidence aggregation (A_v = Σ_{c′} (max_u |T_vu| − τ)₊², τ = 2, on rank-first scores), node max-T, edges kept between significant nodes with |T| ≥ 3. Ported from the cmp_detection finalist `n_agg_noanc` (part1 §1.2). This is the **cluster-level coherence arm**: it pools weak cross-cell-type evidence that per-edge max-T cannot.
- `node_cor_den` / `node_agg_den`: the two node-level methods with the **node density filter** (`min_density`): a cluster is kept only if its called edges cover at least that fraction of the possible cross-cell-type member pairs (drops single-linkage chains). `min_density` is fixed by the pilot (smallest value that removes ≥ 90% of null-scenario clusters without removing planted ones) and pinned before the main run.
- `subspace_sum`: sum-of-subspaces generalized eigenproblem with fixed-sequence testing (secondary arm).
- `gene_match`: loading-cosine matching with gene-label permutation (secondary arm; targets gene reuse, not activity coordination).

## Replicates and paired inference (blocker 2)

- **n_rep = 50 replicate datasets per (mode × amplitude) cell** of the factor
  design, seeds `1000·cell + rep`; every method runs on the identical replicate
  (paired by construction). The null scenario adds 50 extra seeds.
- **Primary endpoint:** `membership_F1_pooled` — mean over planted-program rows of
  the per-program best-cluster F1 (`scoring.R: eval_calls_coord`), positive modes,
  marginal strata.
- **Paired comparison:** for each candidate method vs the reference `pa_cc`,
  per-replicate differences of mean F1 over target signals, percentile bootstrap
  CI (2000 resamples of replicates, `set.seed(1)`), reported with the point
  estimate. Same datasets enter both sides (common-replicate intersection).
- **Guardrails (dataset-level binary, then averaged):** null FWER (any call on null
  seeds; Clopper-Pearson 95% upper limit reported), decoy calls on N, spurious
  calls in positive modes. Bound: each ≤ 0.10.
- **Oracle floor:** the same metrics for an oracle that sees the true z of every
  node; methods are judged against the floor, not against 1.

## Node-to-truth rule (blocker 3; `scoring.R`)

Ported from part1 §3.5 (`node_assign` / `eval_calls`), adapted to `sim_coord.R`'s
truth: each node is assigned to the planted program whose true z it correlates
with (|cor| ≥ 0.4, member cell types only) or to SITE (multiple R² of
score ~ strata); a cluster's label is the majority signal over its member nodes
(≥ 50%), else spurious; SITE-majority clusters are decoy calls. Per planted
program: recall = member cell types with a node assigned to it in the best
cluster / |members|; precision = fraction of cluster nodes assigned to it;
F1 = harmonic mean; **found** = recall ≥ 0.5 and precision ≥ 0.5 with ≥ 2 member
nodes. Private programs are not signals: nodes matching only a private program
are unassigned and their clusters spurious.

## What is recovered (per planted program)
A method outputs clusters of per-cell-type programs (nodes) with a cell-type membership set.
- **Recall:** fraction of true member cell types whose true-program node is in one cluster.
- **Precision:** fraction of cluster members that are true members.
- **Breadth error:** called breadth minus true breadth.
- **Found:** recall >= 0.5 and precision >= 0.5 with at least 2 members.
- **Spurious call:** a cluster that matches no planted program, and **decoy call:** a cluster aligned with N after adjustment. Both reported per dataset as FPR.
- **Detectability floor:** the same metrics for an oracle that sees the true z of every node. Methods are judged against the floor, not against 1.

## Decision rules
To be written after the pilot, before tagging. Placeholder: a method is preferred over `pa_cc` only if recall on R50 and R20 improves by at least a margin (TBD after the pilot; candidate +0.05 membership F1 with the paired CI excluding 0) with decoy and spurious FPR no worse than `pa_cc` + 0.03.

## Artifact decoy arms (`sim_coord.R`, `artifacts = c("vdj", "deletion")`)

Two artifact-shaped signals are planted as **decoy classes** (never targets;
clusters aligned with them are artifact calls, reported like decoys, never in
the F1):

- `VDJ` (Ig/TCR-like): a donor-composition latent expressed on a 40-gene
  receptor pool (renamed `IGLV*`/`TRBV*` so the package's `artifact_genes()`
  regex flags every gene) in a random 3-5 cell-type subset, same genes across
  members — the real shape of Ig/TCR programs.
- `DEL` (donor-deletion-like): 1-2 donors shifted on a shared 100-gene set in
  every cell type; the "activity" is the donor indicator — the exact shape of
  a low-quality-donor artifact.

The count-level variant `gen_coord_counts()` realizes the same latents as
Poisson counts (library size `ncell`, relative expression exp(Y)) and emits
raw counts plus the log2-CPM fit input. Stated limit: Poisson aggregation
only, no NB overdispersion arm.

## Known limits to state up front
- A 20% program in two small cell types sharing about 30 donors may be undetectable for any method; that is why the oracle floor is reported.
- Nuisance here is measured. An unmeasured donor-level nuisance is out of scope for this registration.
