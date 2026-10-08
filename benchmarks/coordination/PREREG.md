# DRAFT (not registered): recurring programs across cell types

Goal: after EBMF per cell type on pseudobulk, find programs that **recur** in all cell types, in a fraction of them, or in only a few, and tell them apart from private programs and from donor-level nuisance. Not registered until margins are filled in and the tag `prereg-coord-v1` exists.

## Simulation (`sim_coord.R`, to write; extends `R/simulation.R`)
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

## What is recovered (per planted program)
A method outputs clusters of per-cell-type programs (nodes) with a cell-type membership set.
- **Recall:** fraction of true member cell types whose true-program node is in one cluster.
- **Precision:** fraction of cluster members that are true members.
- **Breadth error:** called breadth minus true breadth.
- **Found:** recall >= 0.5 and precision >= 0.5 with at least 2 members.
- **Spurious call:** a cluster that matches no planted program, and **decoy call:** a cluster aligned with N after adjustment. Both reported per dataset as FPR.
- **Detectability floor:** the same metrics for an oracle that sees the true z of every node. Methods are judged against the floor, not against 1.

## Methods compared
- `pa_cc`: current pairwise principal-angle test (`sharing_spectrum`, scores space, `strata`) then connected components on called pairs.
- `pa_prog`: program-level calls from the principal vectors, then clustering.
- `score_cor`: node-by-node score correlation with a permutation null (strata), then clustering.
- (no PMD arm, per the earlier review; reconsider only if all of the above fail)

## Decision rules
To be written after the pilot, before tagging. Placeholders: a method is preferred over `pa_cc` only if recall on R50 and R20 improves by at least a margin (TBD) with decoy and spurious FPR no worse than `pa_cc` + 0.03.

## Known limits to state up front
- A 20% program in two small cell types sharing about 30 donors may be undetectable for any method; that is why the oracle floor is reported.
- Nuisance here is measured. An unmeasured donor-level nuisance is out of scope for this registration.
