# Lessons from the Biomni CMP benchmark (`sess_a9d08850c561`)

**Status 2026-10-11: partial.** Part A ("ours") is from this repo. Part B ("theirs") now includes their Part 1 reply (callers, scoring, with code; verbatim copy in `docs/biomni_transfer/part1_callers_scoring.md`, message `msg_2e6a1f26d5229c06`). Parts 2-4 (simulator, prereg and screen results, real data) were requested but the agent run paused at its 750-credit limit before answering; the real-data summary it produced instead is in B1c. Their code was run through our `methods.R` helpers (`.perm_rows`, `.cluster_edges`), so it is directly portable.

Sources: Biomni final reply (message `msg_8143fe9c581ee72a`) and its file names under `tusoskill/cmp_v1`, `cmp_v2`. Biomni terms: CMP = cross-cell-type metaprogram, "cap" = max programs (K) per cell type.

## Part A. Ours (what the repo has)

### A1. Callers (`benchmarks/coordination/methods.R`)
| Arm | Input | Statistic | Null | Output | Lines |
|---|---|---|---|---|---|
| `node_cor` | scores | Fisher-z score correlation per node pair, shared donors >= 20 | rows permuted within strata (`.perm_rows`), max-T over all pairs | edges, then greedy single-linkage clusters | `:81` |
| `gene_match` | loadings | cosine over shared genes x sqrt(#genes), genes >= 100 | gene-label permutation, max-T | edges, then clusters | `:105` |
| `pa_cc` | scores | principal-angle test (`sharing_spectrum`), connected components of called cell-type pairs | strata permutation, BH inside `sharing_spectrum` | cell-type level only (`program_id` NA) | `:132` |
| `subspace_sum` | scores | top eigenvalue of `sum_c U_c U_c'` with donor weights `D` | deflated sequential permutation (`10c9b55`), fixed-sequence stop | cell-type sets per direction | `:156` |
- Shared pieces: node list `.nodes` (`:7`), within-strata row permutation `.perm_rows` (`:18`), `.cluster_edges` (`:26`: strongest edge first, skip a merge that puts two nodes of one cell type in a cluster), `.maxT_edges` (`:59`).
- Toy check (`test_methods.R`, 10 seeds, 999 perms; R100/R50/R20 recovery): `node_cor` 1/1/1, `gene_match` 1/1/1, `pa_cc` 1/0/0, `subspace_sum` 1/1/0.3. Null gives 0 clusters for all.

### A2. Simulator (`benchmarks/coordination/sim_coord.R`)
- 10 cell types, 100 donors, donors per cell type 100/75/50 with a 30-donor core shared by the small ones (`:21-35`); Institute nuisance on every cell type (`:19-20`, `:77`).
- Planted R100, R50, R20 (`COORD_BREADTH`, `:11`), 1-2 private programs per cell type, 5 modes crossing activity sharing x gene sharing (`:4-8`). Gaussian data, so no count structure, no ambient RNA, no Ig/TCR genes.

### A3. Scoring and prereg (`benchmarks/coordination/PREREG.md`, draft, not registered)
- Metrics defined: recall, precision, breadth error, "found" (recall and precision >= 0.5, >= 2 members), spurious and decoy FPR, oracle floor.
- Open blockers (HANDOFF): method list vs `methods.R`; replicates and paired inference; node-to-truth matching rule (not written); pinned EBMF config; decision-rule margins (TBD).

### A4. Process that worked (`benchmarks/model_selection/`)
- Prereg tag `prereg-v1`, numbered deviations (`DEVIATIONS.md` D1-D8), discarded pilot reps for tuning amplitudes, guardrail columns in `decisions.csv`, array-chunked SLURM runs.

### A5. Known gaps in ours
1. No cluster-level coherence statistic or test (HANDOFF (b)).
2. No node density filter (cNMF-style; HANDOFF (c)).
3. Simulator has no real-data nuisances (ambient, Ig/TCR, donor repertoire).
4. No real-data CMP table; no site-replication check.
5. Guard design flaw from M2: absolute `frac_untested`, should be relative to base.

## Part B. Theirs: stated in the reply, plus what to retrieve

### B1. Known from the reply
- **Screen:** simulations based on Stephenson data, seeds 1-8; F1 0.9325 at cap 30 and cap 40; true-member recall 1.00; four false-call checks all 0. Cap 40 rejected for cost. Cap 30 is the pre-registered candidate.
- **Structure:** prereg v1.2, a development seed range (1-30 so far), a held-out range (1051-1076, 1100), 100 extra null datasets, separate "nullx" null screen, ledger backup.
- **Real data (Stephenson, COMBAT), cap 10 vs 40:** CMPs 13 -> 31 and 11 -> 30; Ig/TCR-dominated CMPs 4 -> 6 and 1 -> 9; ambient CMPs flat (4 -> 3, 4 -> 4); severity-linked (p < 0.001) 3 -> 5 and 3 -> 3. Four cell types hit 40 at cap 40.
- **Site replication in Stephenson:** concordant sign 0.92 -> 0.90 (Cambridge), 0.98 -> 0.97 (Newcastle); also nominally significant 0.62 -> 0.54, 0.77 -> 0.77.
- **Confound they flagged:** COMBAT cap 10 vs 40 also changed sex adjustment.
- **Ops:** out-of-memory when scoring ran next to 16 parallel fits; reduced to 8 workers; Stephenson-based held-out fits 6-8 h.

### B1b. From their Part 1 (callers and scoring)
Stated by them, not re-run by us.
**Caller finalist `n_agg_noanc`** (all v2 tables call it `main`; upstream arms K20/K30/K40 only change the fit):
1. **Rank first, then residualize on strata.** Rank-Gauss each score column, then remove strata means. The reverse order let single-donor spikes shift a stratum mean, so residual ranks encoded strata. Ours (`node_cor`, `methods.R:81`) uses raw scores.
2. **Ambient filter before testing.** Drop nodes whose loading correlates > 0.3 with (log2 CPM profile of another cell type minus own). Needs per-cell-type mean profiles.
3. **Node statistic, not edge statistic.** `A_v = sum over partner cell types c' of (max_u |T_vu| - tau)_+^2`, tau = 2, `T = atanh(r) sqrt(n-3)` on shared donors (>= 20). It pools weak evidence across partner cell types, so a program coupled to 3-4 cell types is found when no single edge passes edge max-T. Null: 199 within-strata row permutations, single-step max of `A_v` over nodes gives node FWER.
4. **Membership:** keep edges with both end nodes significant and |T| >= 3, then the same greedy single linkage as ours (`.cluster_edges`).
5. **What did not work:** adding edge "anchors" (edge max-T and node test, each at alpha, not jointly controlled) gave null FWER 0.14; a minP over both families was tried to repair it. Adaptive top-k (`c1_adapt`) uses minP across k from the same permutations.
6. **Context variable (C2):** site-residualized context vs every program, max-|T| permuted within strata, Bonferroni split of alpha with the node test; extra "context CMPs" grouped by sign.

**Scoring (their `eval_cmp.R`)**
- Node-to-truth rule: a node is assigned to the best planted signal among those that include its cell type (plus a SITE pseudo-signal, sqrt(R^2) of score ~ site) if |cor(score, true z)| >= 0.4, else unassigned. CMP label = majority signal if > 50% of its nodes, else "spurious". Signs ignored.
- F1 per target signal = max over called CMPs of F1(recall = distinct member cell types with an assigned node / members; precision = assigned nodes / nodes). No one-to-one matching. Cell-type-level callers (`pa_cc`) score more leniently.
- Primary endpoint: mean F1 over signal rows in positive scenarios (marginal mode). Paired difference vs reference per dataset, 2000-resample percentile bootstrap over datasets. Breadth strata R100/R50/R20 from signal names.
- `oracle` = upstream recall (member cell types having any node assigned to the signal), a caller-independent ceiling.
- Four dataset-level false-call checks, guardrail <= 0.10 each per base: null FWER (130 null datasets, Clopper-Pearson upper bound), ambient decoy, within-context false call, genes_only false call. Extras: `site_any`, `spurious_pos`.
- Their own caveats: the 0.4 cut-off drives F1 and recall; loadings are not scored; no one-to-one matching.

### B1c. From their latest final reply (real data at cap 30; message `msg_71e605ac08c700f2`, run paused at its 750-credit limit)
Stated by them, not re-run by us.
| Cohort | Cap 10 / 30 / 40 | CMPs | Narrow (2-4 cts) | V(D)J-dominated | Genotype-marked | Severity p<1e-3 | Cross-cohort matches |
|---|---|---|---|---|---|---|---|
| Stephenson | | 13 / 27 / 31 | 2 / 5 / 7 | 4 / 2 / 6 | 1 / 8 / 9 | 3 / 2 / 5 | 4 / 7 / 5 |
| COMBAT | | 11 / 25 / 30 | 3 / 10 / 15 | 1 / 6 / 9 | 2 / 6 / 11 | 3 / 5 / 3 | 4 / 7 / 5 |
- **New flag: genotype-marked CMPs.** A CMP is marked if >= 1 gene from a fixed list (donor deletion, copy number, HLA type, sex chromosome: GSTM1, HLA-DQ/DR, CCL3L3, ...) is among its 12 most recurrent top genes. A donor lacking GSTM1 lacks it in every cell type, so the caller is right that these co-vary; they are not coordinated programs. 6 of 13 new CMPs per cohort were marked; the GSTM1 CMP even matches across cohorts. Only 1 of 26 new CMPs associates with severity (p < 0.05).
- Cross-cohort match: Stephenson CMP whose best COMBAT match has loading |r| >= 0.5 over >= 3 shared cell types.
- Existing cap-10 multi-cell-type CMPs all survive at cap 30; one severity association weakened (p 2e-4 -> 0.065). Y-chromosome genes (RPS4Y1, DDX3Y) appear in one new CMP despite the sex covariate; they say to check it against donor sex labels.
- Cap 30 chosen over 40; global `max_factors` default stays 10, logged as a deviation from their prereg v1.2.
- **Real-data calls need a triage step** (ambient, V(D)J, genotype genes, severity/context association, site replication). Their simulator has no V(D)J or genotype structure, so its false-call checks do not cover these.
- **Side effect to know about:** their agent built one commit (`feat/coordinated-programs`, 57 files, on `b0c6986`, 42 tests pass) with a bundle, patch and `PR_DESCRIPTION.md` under `pr/`; the GitHub push got 403. Not applied here; files are only reachable through the Biomni web app. Review before applying anything to `main`.

### B2. Still to retrieve (Parts 2-4)
| Topic | Question | Ours to compare |
|---|---|---|
| Simulator | How Stephenson-based "semi" and "param" data are built; planting; amplitude tuning; nullx; any Ig/TCR signal | A2 |
| Prereg v1.1/v1.2 | Decision rules, guardrails, adoption rule, seed splits, deviations | A4, `PREREG.md` margins |
| Screen scorecard | Per-caller rows and failure modes | A1 toy check |
| Real-data flags | Ig/TCR, severity, site-replication definitions; covariates; K cap choice | A5.4 |
| Ops and redesign | What breaks; what they would change | all |

## Part C. Adoption candidates so far (provisional; finalize after Part B)
Do now (no Part B needed):
1. Register the untested/coverage guard relative to `base` in `PREREG.md` (A5.5, D8).
2. Add a real-data nuisance mode to `sim_coord.R`: an Ig/TCR-like donor-specific block of genes confined to 2-4 cell types. Without it, no arm can be scored on the failure Biomni saw at high K.
3. Add a "K cap exceeded" diagnostic to real runs (count cell types that hit `K`); report it with every real-data result.
4. Real-data stage: copy their CMP table (size, Ig/TCR share, ambient share, severity, site replication) and their dev / held-out seed discipline.
5. Never change two factors between compared runs; write that as a prereg line.

From Part 1 (stated, not yet tested by us):
6. **Node aggregation statistic `A_v`** is the best candidate for our missing cluster-level coherence test (HANDOFF (b)): it addresses R50/R20 power where edge max-T and `subspace_sum` are weak. Port as a new arm `node_agg` and compare to `node_cor` on the toy check and the registered simulator.
7. **Rank-first-then-residualize** score preparation and the **ambient filter** as shared pre-steps for all arms (needs profiles).
8. **Closes a PREREG blocker:** their node-to-truth rule (|cor| >= 0.4, SITE pseudo-signal, majority label > 50%, F1 = max over CMPs, paired bootstrap, four guardrails <= 0.10) can be adopted as our matching rule, with their caveat (no one-to-one) fixed by a Hungarian match.
9b. **Add a genotype-marker flag** (fixed gene list: GSTM1, HLA class II, CCL3L3, sex-chromosome genes) to the real-data CMP table, next to Ig/TCR and ambient; and plant a donor-deletion-like signal (a gene block absent in a donor subset across all cell types) in the simulator, because a method that "wins" on coordination will call it.
9. **Do not mix test families** (their anchors raised null FWER to 0.14): keep one family per call or control them jointly.

Later (after retrieval):
- Their null screen ("nullx") design for our null-calibration gate.
- Their caller details, if any beat `node_cor` on the scorecard.
