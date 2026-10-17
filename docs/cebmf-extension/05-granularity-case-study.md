# Case study: pseudobulk vs metacell vs single-cell programs with cellprograms-style EBMF

**Status**: pilot (workstream F + G-M3 synthetic battery), UCSF CD8+ T cell cohort
(53,474 clean cells, 99 donors, top-2,000 HVGs, seed family 20261017). FINAL at
10/12 arms: the two GeneNMF fine-granularity arms return zero metaprograms at
registered defaults (documented below with mechanism). Generated 2026-10-17;
commits `1b40d77`..`89a4064` + ledger.

## What was compared

Three input granularities built from the same cells (builder: `build_granularity.py`,
gamma=25, MAX_CELLS=30000):

| Granularity | Matrix | Construction |
|---|---|---|
| pseudobulk | 99 x 2000 | mean counts per donor |
| metacell | 2,138 x 2000 | k-means (k=5) within donor on shifted-log PCs |
| single-cell | 29,998 x 2000 | cells as observations (capped) |

Four methods x three granularities = 12 arms: **ebmf** (flashier, R, the cellprograms
backend), **cebmf_torch** (PyTorch, laplace/laplace, column-wise), **cNMF** (100
restarts, density-silhouette K selection), **GeneNMF** (K grid 2-15, metaprogram
consensus). All arms share the shifted-log recipe (alpha=1, library-size size factors)
and the same HVG panel.

## Arm summary (K, runtime, primary seed)

| Method | pseudobulk | metacell | single-cell |
|---|---|---|---|
| ebmf | K=10, 45 s | K=10, 115 s | K=10, 640 s |
| cebmf_torch | K_eff=6, 28 s | K_eff=10, 75 s | K_eff=4, 982 s |
| cNMF | k=6 (sil 0.96), 5 s | k=2 (sil 1.00), 2,637 s | k=2 (sil 1.00), 49,532 s |
| GeneNMF | 10 MPs, 7 s | 0 MPs (consensus collapse) | 0 MPs (consensus collapse) |

Seed stability (ebmf, 3 seeds, greedy-matched loadings): 1.00 in all granularities.
cebmf determinism: bit-identical refit at fixed seed (`F_cebmf_determinism.csv`).

## Finding 1 — metacell is the most stable granularity for cebmf

Within-method, across-granularity concordance (mean Hungarian-matched |cor| of gene
loadings; top-50 Jaccard in parens):

| Pair | cebmf_torch | ebmf (flashier) | cNMF |
|---|---|---|---|
| pseudobulk <-> metacell | 0.663 (0.33) | 0.340 (0.09) | 0.728 (0.18) |
| pseudobulk <-> single-cell | 0.558 (0.29) | 0.192 (0.04) | 0.710 (0.21) |
| metacell <-> single-cell | **0.822 (0.40)** | 0.507 (0.19) | 0.954 (0.49) |

cebmf_torch's metacell and single-cell fits agree at 0.82 mean matched |cor| — the
highest cross-granularity stability observed in this study. flashier's arms agree far
less across granularities (0.19-0.51), i.e. its recovered programs are more
granularity-dependent. cNMF's K selection is unstable across granularities (k=6 at
pseudobulk vs k=2 at metacell AND single-cell, the latter two with degenerate perfect
silhouettes of 1.00 — the density-silhouette criterion collapses to the coarsest split
whenever observations are abundant), so its cross-granularity numbers reflect the
single 2-program axis; on that axis cNMF is self-consistent (mc<->sc 0.954).

## Finding 2 — the two EBMF implementations converge best at single-cell granularity

Same-granularity, cross-method concordance:

| Pair | pseudobulk | metacell | single-cell |
|---|---|---|---|
| ebmf <-> cebmf_torch | 0.460 | 0.610 | **0.754** |
| cNMF <-> cebmf_torch | 0.294 | 0.651 | 0.600 |
| GeneNMF <-> cebmf_torch | 0.474 | 0 MPs | 0 MPs |
| ebmf <-> cNMF | 0.282 | 0.570 | 0.498 |

At single-cell granularity the R and PyTorch EBMF stacks recover substantially the
same program space (0.75); at pseudobulk granularity they diverge most (0.46) —
consistent with cebmf's pseudobulk arm pruning to K_eff=6 while flashier retains
K=10 (the extra flashier factors are the weaker, less reproducible ones).

## Finding 3 — programs are biologically coherent in every arm, and one age axis is universal

Best-program marker enrichment (mean rank percentile of canonical CD8 T marker sets —
cytotoxicity, IFN response, cell cycle, naive/memory, terminal effector, heat shock,
antigen presentation, NFkB): 0.94-1.00 in every loaded arm; metacell arms reach
0.997 (cebmf and ebmf). Metacell aggregation loses nothing in marker coherence.

Donor-score vs age (Spearman, 99 donors): every arm contains a dominant age-associated
program with |rho| 0.34-0.49 (cebmf pseudobulk -0.459, ebmf single-cell +0.479, cNMF
metacell -0.493/+0.493 as a sign-paired axis). The same age axis survives across all
four methods and all three granularities — the strongest cross-arm biological
invariant in this study.

## Finding 4 — metacell aggregation lifts rare-program recovery above the single-cell SVD ceiling (synthetic truth)

G-M3 synthetic battery (5,000 cells x 1,500 genes, 50 donors, 3 broad + 1 rare (5%)
+ 2 null programs, Poisson UMI generator, nested dropout 0/10/30%, 10 reps; gates
evaluated relative to the per-rep SVD ceiling per Amendment A3):

- **S3 (metacell-init + project)** recovers the rare program at mean 0.840 matched
  |cor| over 10 reps (range 0.744-0.964), exceeding the cell-level SVD-10 ceiling
  (mean 0.780) in 8/10 reps (max +0.18). S1 (from-scratch single-cell) sits exactly
  AT the ceiling (0.777 vs 0.780) — the lift is specific to the metacell input.
- **Mechanism (tested by controls, not assumed)**: aggregation does NOT lift the
  linear ceiling itself — SVD on the metacell matrix has rare ceiling 0.776 (no
  lift), and projecting cells onto metacell-estimated SVD directions recovers rare
  at only 0.487 (worse than cell-level SVD). The lift arises because k-means
  metacells CONCENTRATE rare states (rare-high cells cluster into dedicated
  metacells, raising their usage fraction from ~5% of cells to majorities of a few
  metacells) while sqrt(n) averaging suppresses sampling noise; that makes the rare
  direction identifiable for cebmf's OBLIQUE sparse factorization, which is not
  bound by the orthogonality constraint that caps SVD. S1 vs S3 differ only in the
  source of F (same projection path), so the lift is in factor estimation.
- S3 held-out prediction 0.798 vs S1 0.855: metacell init trades a little held-out
  fidelity for rare-program recovery and a ~40% faster fit.
- S2 (donor pseudobulk-init) collapses (matched |cor| 0.03-0.38): donor averaging
  discards the cell-level usage variance the model needs.
- S4mc (metacell-init + refine) matches S1 held-out (0.855) — initialization is
  safe but, at 5k cells, costs 1.6x more wall-clock than from-scratch; the C1
  cost claim is a large-N claim (tested at M4/M5).
- cebmf battery (270 cells: S1-S4 x shifted-log + S1/S3 x pflog/clr x 3 dropout
  rates x 10 reps): COMPLETE and aggregated. S1/S3 pass the A1 broad gate in
  nearly all conditions and C2 rare gate almost everywhere; S2 (donor
  pseudobulk-init) is dead (0/10 gates at every rate); 67/270 cells (24.8%)
  entered the silent-divergence regime and were all rescued by the adapter's
  loud-raise + K-halving fallback (~2x wall-clock cost); CLR-family normalization
  plus heavy zeros destroys rare recovery at metacell granularity (C2 0/10 at
  30% dropout, rare |cor| 0.048).
- cNMF anchor battery (S5, 30 cells): RUNNING; anchor-arm conclusions deferred
  until it lands.

## Finding 4b — the same ceiling-relative logic on DIVAS: anchoring beats the linear ceiling

DIVAS blocks are sample-level pseudobulks (114 samples; no cell-level matrix), so
metacell strategies do not apply — but the ceiling-relative logic does, against the
in-session cellprograms reference fits (DIVAS-parity recipe, seed=1). D1 =
from-scratch cebmf (K=30, 5 seeds, bit-stable); D2 = anchored refinement
(init_F = reference loadings, init_L = least-squares, maxit=50):

| Block | SVD-30 ceiling | D1 from-scratch | D2 anchored |
|---|---|---|---|
| CD4 T | 0.409 | 0.278 | 0.636 |
| CD8 T | 0.364 | 0.224 | 0.681 |
| CD14 monocytes | 0.425 | 0.277 | 0.661 |
| NK | 0.333 | 0.209 | 0.507 |
| proteomics (481) | 0.489 | 0.422 | 1.000 |
| metabolomics (763) | 0.630 | 0.474 | 1.000 |

Three observations. (1) D1 sits BELOW the linear ceiling on every block: from-scratch
K=30 on 114 samples is heavily underdetermined and converges to an optimum that
reproduces the cellprograms program space no better than orthogonal SVD. (2) D2
EXCEEDS the ceiling on every block (0.51-1.0 vs 0.33-0.63) — the same oblique-factor
phenomenon as the synthetic rare program: prior-shrunk oblique factors match the
reference better than any orthogonal 30-dim subspace. (3) On proteomics/metabolomics
the anchor is a fixed point (retention 1.000) — small blocks refine onto the anchor
exactly; gene blocks retain 27-29 of 30 programs. Anchoring multiplies concordance
with the cellprograms reference ~3x on gene blocks, ~2x on small blocks. (DIVAS-vs-
DIVAS concordance, gate C3, still awaits the blocked Zenodo RDS.)

## Finding 5 — normalization sensitivity is material (Amendment A4)

PFlog/CLR (shifted-CLR family, Booeshaghi et al. 2022) vs the registered shifted-log
recipe, UCSF real data (alpha-hat calibration active; concordance = greedy-matched
|cor| vs the shifted-log fit, K=10 everywhere):

| Arm | PFlog | CLR |
|---|---|---|
| single-cell | 0.580 | 0.562 |
| metacell | 0.745 | 0.705 |

Normalization choice moves recovered programs as much as the method choice does at
single-cell granularity. On the Poisson synthetic, PFlog == CLR (alpha-hat hits the
registered floor); the calibration only differentiates on overdispersed real data.

**DIVAS extension (registered pre-run)**: pflog/clr vs the center-only parity
control on the DIVAS blocks (from-scratch, D1 protocol, 5 seeds, 50 fits, zero
divergence, seed-stable). Proteomics excluded (1,703 negative values make log1p
undefined; no ad-hoc shift invented). Result: pflog/clr sit BELOW center-only on
every block (pflog 0.154-0.172, clr 0.158-0.216 vs D1 0.209-0.278 gene blocks /
0.474 metabolomics); the largest drop is metabolomics (3x) — the one count-like
block where the transforms apply as intended. Gene-block values are already
log-scale, so the transforms act as monotone compressions there (calibration
premise void, as on the Poisson synthetic). Interpretation: normalization choice
materially moves DIVAS programs away from the cellprograms reference; center-only
parity remains the registered DIVAS recipe; cross-normalization program
comparisons are not interchangeable.

## Status of metacell extensions in cellprograms (audit view)

- **Exists in-session, not in the package**: the granularity builder (donor
  pseudobulk, k-means-within-donor metacells, single-cell) is an audit script
  (`build_granularity.py`); cellprograms/flashier ship no metacell builder.
- **Metacell is where cebmf_torch wins**: highest cross-granularity stability (0.82),
  rare-program recovery above the single-cell SVD ceiling, and 13x faster fits than
  single-cell (75 s vs 982 s at this cohort size).
- **Known hazards**: cNMF's density-silhouette K selection degenerates whenever
  observations are abundant (k=2, silhouette 1.00 at metacell AND single-cell);
  GeneNMF metaprogram consensus collapses at metacell/single-cell — mechanism
  established by reconstruction: confidence = fraction of member patterns whose
  top-gene set contains a gene, and the max is <= 0.2 in every cross-K cluster
  (< min.confidence 0.5; still 0 MPs at a min.confidence=0.2 sensitivity), i.e.
  per-K NMF patterns on fine-grained input are too inconsistent for any cross-K
  consensus; ebmf-style variance collapse under entrywise dropout (workstream E)
  argues for metacell construction BEFORE fitting when zeros are mean-dependent.
- **mcRigor link (MR-3)**: the ZeroRate budget that drives the ebmf dropout collapse
  is a per-granularity design constraint for any metacell builder cellprograms adopts.
- **Pending**: cNMF anchor battery S5 (running); M4 1M-cell chunked run; M5 500k
  shard scaling; M6 adopt/revise/abandon gate; upstream frozen-prior request
  (weak-factor projection fidelity, M2).

## Reproduction

Scripts: `python/benchmarks/m3/` (synthetic battery, A4 norms, DIVAS D1/D2 +
normalization extension), audit workstream F (`build_granularity.py`, `ebmf_F.R`,
`cebmf_F.py`, `cnmf_F.py`, `genenmf_F.R`, `genenmf_conf02.R`, `compare_F.py`).
Outputs: `F_comparison_pairs.csv`, `F_concordance_heatmap.png`,
`F_marker_enrichment.csv`, `F_age_assoc.csv`, `F_ebmf_norm_sensitivity.csv`;
benchmark artifact tables committed in-repo under `reports/audit/G_m3/`;
per-cell battery JSONs on the shared store. Nothing in this vignette is a
registered acceptance criterion; gates C1-C4 are as pre-registered
(04-experiment-plan.md, amendments A1-A4).
