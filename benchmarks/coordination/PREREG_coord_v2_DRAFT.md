# prereg-coord-v2 — DRAFT (unregistered)

**Status:** DRAFT, 2026-11-17. Not registered, no tag, nothing run. `prereg-coord-v1`
(`PREREG.md`, tag `prereg-coord-v1`, FILLED 2026-10-21) is frozen and untouched.
**Purpose:** one place for the changes the benchmark needs before a second registered run, and for
every choice that still needs your value. Each row says where its value comes from.

Row sources: `user` (you decided it), `inherited: <path>` (v1 set it; needs a one-line confirm
before registration, because the earlier stage's code does not settle it), `open` (no value yet).

## 1. What changes against v1

| # | Change | Code | Status |
|---|---|---|---|
| 1 | New arm `node_agg` (statistic `A_v = sum_c' (max_u \|T_vu\| - tau)_+^2` on rank-first scores, node max-T, greedy linkage). Restored from `df0edaa`; it was in the pre-freeze suite and is not one of the v1 arms. Same statistic as the package default caller (`n_agg_noanc`, `DECISION_default_caller.md`). | `methods.R::method_node_agg`, tested in `test_methods.R` | code done |
| 2 | New trap scenario `donor_deletion`: a gene set absent in a fraction of donors, in every cell type. | `sim_coord.R`, tested in `test_sim_coord.R` | code done, values below |
| 3 | New trap scenario `ig_tcr`: a clonal gene set, large in a few donors, independently per affected cell type. | same | code done, values below |
| 4 | Comparison family named in the decision table; p-values unadjusted. | `family.R`, `summarize_main.R`, tested in `test_summary.R` | code done |
| 5 | Relative untested-pair guard (arm <= base + 0.05), opt-in `UNTESTED_GUARD=relative`. | `benchmarks/model_selection/summarize.R` | code done; default stays absolute |
| 6 | Confounder table (F on S0 for `base`, `none`, `none_strata`, `strata`, Clopper-Pearson bounds, no verdict). | same | code done |
| 7 | Stage-2 EBMF probe (`ebmf2`), continuous metrics only. | `probe_ebmf2_stage2.R` | script done, smoke-run on one seed, no result |

## 2. Rows decided by you

| Row | Value | Source |
|---|---|---|
| Multiplicity across arms | none; report unadjusted p and name the family and its size | user, 2026-11-17 |
| Untested-pair guard, relative form | arm untested fraction <= base + 0.05 | user, 2026-11-17 |
| `donor_deletion` fraction of deleted donors | 0.5 | user, 2026-11-17 |
| `donor_deletion` number of deleted genes | 5 | user, 2026-11-17 |
| `ig_tcr` fraction of clonal donors per affected cell type | 0.1 | user, 2026-11-17 |
| `ig_tcr` clonal amplitude | 3 x `amp` (gene weights with sd 3 * amp) | user, 2026-11-17 |
| Three labels (`coordinated`, `gene_reuse`, `both`) | reported separately | user, 2026-10-09 (HANDOFF) |
| Method set excludes MOFA-on-scores and GFA | MOFAcellulaR only as optional comparator | user, 2026-10-09 (HANDOFF) |

## 3. Open rows (need your value; no defaults)

| Row | Why it is open |
|---|---|
| `donor_deletion` depth (`del_depth`, the downward shift in the deleted donors) | `gen_coord()` requires it; no value chosen |
| `ig_tcr` number of affected cell types (`n_clone_cts`) and genes (`n_clone_genes`) | required arguments; no value chosen |
| Base amplitude `amp` for the trap families | v1 used 0.50 (`PREREG.md`, amplitude pilot); whether the trap families reuse it is open |
| Replicates per trap family, and which arms run on them | not set |
| `both` label merge rule | none assumed: labels stay separate until you set a rule |
| How much `none` F on S0 counts as "inflated" in the confounder check | `PREREG.md` (model-selection) does not say. A Clopper-Pearson lower bound above 0.05 would call semi `none` (F = 0.09 [0.054, 0.139]) inflated, which `RESULTS.md` records as not inflated |
| Node-level scorer | v1 scores membership-set F1 only. `reports/historical/c1/code/evaluate.R` is a node-to-truth scorer that does not run on this simulator (it reads `truth$scores`, `truth$loadings`, a `nuisance_group` covariate); porting is open |
| Scoring tier for arms without node output (`pa_cc`, `gene_match`) | not set |
| Strata null for external comparators | not set |
| Meaning of `act_cor` | `sim_coord.R` sets each member's activity to `act_cor * z + sqrt(1 - act_cor^2) * noise`, so the pairwise correlation is `act_cor^2` (0.49 at 0.7). The v1 primary endpoint uses `act_cor = 1`, where both readings agree. Keep as is, or redefine it as pairwise |
| `ebmf2` threshold for "recovered" | the probe reports `share` and `n_ct` only; a call needs a threshold, and no figure plan is approved |

## 4. Rows inherited from v1 (confirm each before registration)

All inherited from `benchmarks/coordination/PREREG.md`; none is settled by the code that uses it.

| Row | v1 value | Section |
|---|---|---|
| Hierarchy of arms | `eb_cp` > `block_E` > attribution arms vs `pa_cc` | 7 |
| Adoption margin | paired recall gain >= +0.10, precision not worse by > 0.05 | 7 |
| FPR gate | <= 0.10 absolute and <= `pa_cc` + 0.03 | 7 |
| Found rule | recall >= 0.5 and precision >= 0.5 | 6 |
| Decoy rule | activity cor with nuisance >= 0.5 | 6 |
| Primary family | act_only, R50 + R20, `act_cor` = 1 | 6, 7 |
| Replicates | 32 per primary family (development scale); 96 for the confirmatory scale | 8 |
| Gene pool and program size | G = 2000, 50 genes per program | 4 |
| Permutations | `pa_cc` 1999; other arms 999 | 5 |
| `node_agg` parameters | tau = 2, link threshold 3, alpha 0.05, min shared donors 20 (as in `n_agg_noanc`) | method code |
| Seeds | shared dataset seeds across arms; the seed list | 8 |

## 5. Not planned here

No SLURM job is submitted for any of this. A run needs: this draft turned into a registered
`prereg-coord-v2` with every open row filled, and your approval of the compute (the v1 main run
was 272 datasets, about 5.4 h on 12 workers).
