# Plan: SEACells + mcRigor metacell builder for cellprograms (audit-side, feeds M6)

## Summary

Borrow the *metacell construction layer* of hdWGCNA without adopting its network model:
**SEACells** (GAT-based metacells with soft assignments; Persad et al. 2023 Nat Biotechnol,
doi:10.1038/s41587-023-01716-9) becomes the primary builder, governed by **mcRigor**
(Liu & Li 2025 Nat Commun, doi:10.1038/s41467-025-63626-5) — mcDiv dubious-metacell
detection, DubRate/ZeroRate Score optimization of the granularity γ, and two-step
re-partitioning for rare states. From hdWGCNA (Morabito et al. 2023, GPLv3) we borrow
*principles only* (donor×celltype stratification, k=20–75 size targets, max_shared
overlap control, min_cells gates, metacell QC wrappers) and use the package two ways:
(1) as an **input-preparation wrapper** (MetacellsByGroups → NormalizeMetacells →
metacell matrix → cellprograms/cebmf — the same borrow-as-wrapper pattern as the
MOFAcellular functions), and (2) as a **5th benchmark arm** (full hdWGCNA pipeline,
eigengenes as programs) so the "why EBMF rather than WGCNA" claim becomes empirical.
License rule (binding): hdWGCNA is GPLv3, cellprograms is MIT — **dependency/wrapper
only, never source copying**. Program model unchanged: cebmf adapter (metacell fit →
chunked projection → optional refine + divergence guard). Audit-side only; all results
fold into the M6 adopt/revise/abandon evidence; no package code changes in this plan.

## Workstreams

- **W1 — SEACells builder module** (Python, worker-0): `build_metacells_seacells()`
  taking the UCSF CD8 AnnData (53,474 clean cells, top-2000 HVGs, donor column),
  donor-stratified per MetacellsByGroups principle. Outputs a versioned contract:
  metacell matrix (metacells × genes), **soft-assignment matrix** (cells × metacells),
  n_cells weights, and a QC table. Doublet-removal status of the input verified before
  running (mcRigor precondition).
- **W2 — mcRigor γ optimization + QC**: sweep γ ∈ {25, 50, 100} (cells per metacell),
  compute DubRate/ZeroRate/Score (w=0.5) per γ and per donor; pick argmax Score;
  run mcDiv dubious-metacell detection (double permutation, size-specific thresholds)
  on the chosen partition; two-step re-partitioning of dubious metacells (mcRigor's
  rare-state extension). If the mcRigor implementation is not publicly runnable,
  implement mcDiv + double permutation exactly per the Methods (fully specified).
- **W3 — hdWGCNA input-prep wrapper** (R, worker-1): pin hdWGCNA 0.4.12 (GitHub dev),
  wrap MetacellsByGroups(group.by = donor × celltype-if-available, k=25, max_shared=10,
  min_cells) → NormalizeMetacells → export genes×metacells CSV + cell-count weights →
  feed the cebmf adapter. Wrapper only — no WGCNA network steps in this path.
- **W4 — hdWGCNA benchmark arm** (R, worker-1): full pipeline (TestSoftPowers →
  ConstructNetwork → ModuleEigengenes → kME) on the same three granularities;
  eigengenes exported as programs and compared with the existing Hungarian/Jaccard
  machinery in compare_F.py (13-arm rerun + refreshed heatmap).
- **W5 — Battery extension** (worker-0): new S3-seacells arm (SEACells-metacell-init +
  chunked project) and S4-seacells (init + refine), 10 reps × dropout {0, 0.1, 0.3},
  identical generator and gates as the registered battery — tests whether the
  rare-lift-above-ceiling survives a different builder (kNN-soft vs k-means).
- **W6 — Downstream tooling borrow**: kME-style loading-based hub ranking as the
  interim top_genes metric (fills the lfsr gap without upstream changes); WGCNA
  module-preservation (Zsummary) as a second stability metric beside Hungarian
  concordance; UCell signature scores as the registered B3 post-hoc comparator.
- **W7 — Projection-route test**: program-on-cells via SEACells soft assignments
  (cells ≈ weighted metacell combinations) vs the exact chunked L-update (M2 path).

## Registered gates (committed before any run; no post-hoc relaxation)

- **G1 rare-lift transfer**: S3-seacells rare recovery ≥ 0.7× per-rep SVD ceiling in
  ≥ 8/10 seeds (C2 analogue). Failure = builder-dependent lift → documented limitation.
- **G2 builder concordance**: cebmf programs from SEACells-builder vs k-means-builder
  matched |cor| ≥ 0.8 on shared programs. Failure = builder choice is load-bearing.
- **G3 mcRigor QC**: chosen partition's dubious rate < 10% of metacells; ZeroRate of
  the SEACell matrix ≤ the k-means metacell matrix ZeroRate at the same γ.
- **G4 soft-projection agreement**: soft-assignment projection vs chunked L-update
  matched |cor| ≥ 0.99 on strong factors (M2 bars), reported for all factors.
- **G5 hdWGCNA arm**: descriptive row (no pass bar); rare-program recovery vs ceiling
  and cross-granularity stability reported either way.

## Compute & execution

- SEACells on CPU (no CUDA in session): pilot on the largest single donor first to
  measure throughput, then full 53k-cell cohort (est. 1–3 h, GAT training CPU-bound);
  mcRigor sweep adds ~4 partitions × QC compute (~1 h). Worker-0, background with
  checkpoints to the shared store.
- hdWGCNA install (worker-1, R 4.4.2, Seurat present): WGCNA + harmony + UCell deps;
  pipeline runs are minutes at metacell scale (2,138 metacells) — the 30k single-cell
  granularity metacell build is the long pole (~1 h est. from hdWGCNA k=25 guidance).
- Total estimate ~6–10 h across both machines; no HPC; all fits seeded (20261016
  family), revisions pinned (hdWGCNA 0.4.12, SEACells GitHub HEAD recorded).

## Deliverables & sequence

1. Plan doc committed pre-run: `docs/cebmf-extension/06-metacell-plan.md` + ledger
   entry (this document's content), then implementation starts.
2. W1+W2 → builder contract + QC tables (G3).
3. W5 → battery S3-seacells/S4-seacells arms (G1, G2).
4. W3+W4 → wrapper + 13-arm comparison (G5).
5. W6+W7 → tooling + projection test (G4).
6. compare_F.py rerun, heatmap, vignette update (metacell-plan section), ledger +
   commits per milestone; M6 evidence table gains: builder-independence of the rare
   lift, mcRigor-optimized γ, hdWGCNA empirical comparison, soft-projection route.

## Explicit assumptions

- UCSF cohort is doublet-cleaned (53,474 "clean cells"); verified at execution.
- mcRigor code availability checked at execution; fallback is faithful re-implementation
  of mcDiv + double permutation from the Methods.
- SEACells CPU training is feasible at 53k cells; if not, fallback per mcRigor's own
  method comparison: MetaCell (its top-scored method) with the same mcRigor governance.
- No cell-type labels in the UCSF granularities yet (pilot note) — stratification is
  donor-only until h5ad re-extraction adds subtypes.
