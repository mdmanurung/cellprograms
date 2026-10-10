# Lessons from the Biomni CMP benchmark (`sess_a9d08850c561`)

**Status 2026-10-10: partial.** Part A ("ours") is written from this repo. Part B ("theirs") is filled only with what the session's final reply stated; their code, simulator, scoring and prereg are not yet retrieved (the Biomni `send_message` returned HTTP 402; a retry is scheduled for 2026-10-11 00:30, or the files can be downloaded by hand). Every "theirs: not stated" below is a question to ask, not a finding.

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

### B2. Still to retrieve (one row each; fill when the retry returns)
| Topic | Question for them | Ours to compare |
|---|---|---|
| Caller | Exact algorithm, null, strata, multiplicity, cluster/merge, membership rule, code | A1 |
| Density filter / noise nodes | Is there a node filter before calling? | A5.2 |
| Simulator | How Stephenson-based "semi" data are built; planted breadth/genes/amplitude; "nullx" | A2 |
| Scoring | F1 definition; node-to-truth matching; "member share"; the four false-call checks | A3 blockers |
| Prereg v1.2 | Decision rules, guardrails, adoption rule, deviations | A4, PREREG.md margins |
| Real-data flags | How ambient, Ig/TCR, severity, site replication are computed; covariates | A5.4 |
| Failure modes | What failed in the screen; what they would change | all |

## Part C. Adoption candidates so far (provisional; finalize after Part B)
Do now (no Part B needed):
1. Register the untested/coverage guard relative to `base` in `PREREG.md` (A5.5, D8).
2. Add a real-data nuisance mode to `sim_coord.R`: an Ig/TCR-like donor-specific block of genes confined to 2-4 cell types. Without it, no arm can be scored on the failure Biomni saw at high K.
3. Add a "K cap exceeded" diagnostic to real runs (count cell types that hit `K`); report it with every real-data result.
4. Real-data stage: copy their CMP table (size, Ig/TCR share, ambient share, severity, site replication) and their dev / held-out seed discipline.
5. Never change two factors between compared runs; write that as a prereg line.

Later (after retrieval):
- Their F1 and node-to-truth matching, if clearly specified, to close the PREREG blocker.
- Their null screen ("nullx") design for our null-calibration gate.
- Their caller details, if any beat `node_cor` on the scorecard.
