# Package index

## Data & simulation

Build and check `cell_program_data` objects, simulate pseudobulk ground
truth, and bridge to MOFAcellular.

- [`simulate_pseudobulk()`](https://mdmanurung.github.io/cellprograms/reference/simulate_pseudobulk.md)
  : Simulate cell-type pseudobulk matrices with planted programs
- [`as_cell_program_data()`](https://mdmanurung.github.io/cellprograms/reference/as_cell_program_data.md)
  : Create a cell program data object
- [`validate_cell_program_data()`](https://mdmanurung.github.io/cellprograms/reference/validate_cell_program_data.md)
  : Validate a cell program data object
- [`artifact_gene_sets`](https://mdmanurung.github.io/cellprograms/reference/artifact_gene_sets.md)
  : Artifact-gene sets
- [`artifact_genes()`](https://mdmanurung.github.io/cellprograms/reference/artifact_genes.md)
  : Expand artifact-gene sets against a gene universe
- [`as_mofacellular_input()`](https://mdmanurung.github.io/cellprograms/reference/as_mofacellular_input.md)
  : Export pseudobulk profiles to MOFAcellulaR input format

## Fit & canonicalize

Fit per-cell-type programs (flashier R backend or cebmf_torch Python
backend), canonicalize sign/scale, freeze preprocessing, and flag
artifact-carrying programs.

- [`fit_celltype_programs()`](https://mdmanurung.github.io/cellprograms/reference/fit_celltype_programs.md)
  : Fit cell-type-specific programs with flashier
- [`canonicalize_programs()`](https://mdmanurung.github.io/cellprograms/reference/canonicalize_programs.md)
  : Canonicalize program scale and sign
- [`freeze_preprocessing()`](https://mdmanurung.github.io/cellprograms/reference/freeze_preprocessing.md)
  : Freeze a preprocessing specification
- [`flag_artifact_genes()`](https://mdmanurung.github.io/cellprograms/reference/flag_artifact_genes.md)
  : Flag programs or CMPs carrying artifact genes

## Extract

Summaries, scores, loadings and top genes from a canonicalized fit.

- [`program_summary()`](https://mdmanurung.github.io/cellprograms/reference/program_summary.md)
  : Summary of programs per cell type (n, K, R2 of the reconstruction,
  ELBO)
- [`program_scores()`](https://mdmanurung.github.io/cellprograms/reference/program_scores.md)
  : Program score matrix across all cell types
- [`program_loadings()`](https://mdmanurung.github.io/cellprograms/reference/program_loadings.md)
  : Program loading matrices
- [`top_genes()`](https://mdmanurung.github.io/cellprograms/reference/top_genes.md)
  : Top genes of a program
- [`project_programs_to_cells()`](https://mdmanurung.github.io/cellprograms/reference/project_programs_to_cells.md)
  : Project fitted programs onto individual cells

## Compare, align & share

Program stability, cross-fit matching, within-fit alignment, and
donor-level sharing structure.

- [`assess_program_stability()`](https://mdmanurung.github.io/cellprograms/reference/assess_program_stability.md)
  : Assess program stability by donor subsampling
- [`match_programs()`](https://mdmanurung.github.io/cellprograms/reference/match_programs.md)
  : Match programs across two fits (bootstrap refits or cohorts)
- [`align_programs()`](https://mdmanurung.github.io/cellprograms/reference/align_programs.md)
  : Align all cell types' program subspaces to a reference cell type
- [`program_sharing()`](https://mdmanurung.github.io/cellprograms/reference/program_sharing.md)
  : Post-hoc shared-vs-specific assessment of fitted programs (mashr)
- [`sharing_spectrum()`](https://mdmanurung.github.io/cellprograms/reference/sharing_spectrum.md)
  : Pairwise sharing spectra across all cell types
- [`principal_angles()`](https://mdmanurung.github.io/cellprograms/reference/principal_angles.md)
  : Principal angles between the program subspaces of two cell types

## Coordination & annotation

Cross-cell-type coordinated programs (CMPs): detection, triage, and
annotation against TF regulons, ligand-receptor pairs and gene sets.

- [`find_coordinated_programs()`](https://mdmanurung.github.io/cellprograms/reference/find_coordinated_programs.md)
  : Find coordinated multicellular programs
- [`cmp_triage()`](https://mdmanurung.github.io/cellprograms/reference/cmp_triage.md)
  : One-row-per-CMP triage table
- [`annotate_cmp_tf()`](https://mdmanurung.github.io/cellprograms/reference/annotate_cmp_tf.md)
  : Transcription-factor footprints of coordinated programs
- [`annotate_cmp_lr()`](https://mdmanurung.github.io/cellprograms/reference/annotate_cmp_lr.md)
  : Ligand-receptor compatibility of coordinated programs
- [`annotate_cmp_themes()`](https://mdmanurung.github.io/cellprograms/reference/annotate_cmp_themes.md)
  : Themes of coordinated programs
- [`assess_cmp_stability()`](https://mdmanurung.github.io/cellprograms/reference/assess_cmp_stability.md)
  : Assess CMP stability under synchronized donor resampling
- [`prepare_collectri()`](https://mdmanurung.github.io/cellprograms/reference/prepare_collectri.md)
  : Format an OmniPath CollecTRI table as a signed TF-target network
- [`ambient_program_scores()`](https://mdmanurung.github.io/cellprograms/reference/ambient_program_scores.md)
  : Ambient-likeness of programs
- [`cmp_context_association()`](https://mdmanurung.github.io/cellprograms/reference/cmp_context_association.md)
  : Association of coordinated programs with a donor variable

## Joint & tensor factorization

Post-fit joint factorization across cell types: empirical-Bayes CP
tensor factorization (genes x donors x cell types) and block-structured
joint EBMF at expression or scores level, with benchmark-contract
coordination calls. See `vignette("coordination_joint")` and the
preregistered benchmark under `benchmarks/coordination/`.

- [`fit_eb_cp()`](https://mdmanurung.github.io/cellprograms/reference/fit_eb_cp.md)
  : EB-CP tensor factorization across cell types
- [`eb_cp_from_fit()`](https://mdmanurung.github.io/cellprograms/reference/eb_cp_from_fit.md)
  : EB-CP coordination programs from a fitted \`cell_program_fit\`
- [`eb_cp_programs()`](https://mdmanurung.github.io/cellprograms/reference/eb_cp_programs.md)
  : Coordination programs from an EB-CP fit
- [`fit_joint_block()`](https://mdmanurung.github.io/cellprograms/reference/fit_joint_block.md)
  : Block-structured joint EBMF (Level 3 candidate A)
- [`block_clusters()`](https://mdmanurung.github.io/cellprograms/reference/block_clusters.md)
  : Coordination calls from a block_joint fit (benchmark output
  contract)

## Metacells (SuperCell + mcRigor)

Optional upstream aggregation of single cells into metacells before
program discovery. The default backend is SuperCell 2.0 (pure R, no
Python) with mcRigor QC and granularity optimization. The SEACells
backend remains available and additionally requires basilisk (the
SEACells Python environment is built on first use).

- [`MetacellsByGroupsSuperCell()`](https://mdmanurung.github.io/cellprograms/reference/MetacellsByGroupsSuperCell.md)
  : Construct metacells by groups with SuperCell and mcRigor QC
- [`OptimizeGammaSuperCell()`](https://mdmanurung.github.io/cellprograms/reference/OptimizeGammaSuperCell.md)
  : Optimize SuperCell granularity with mcRigor
- [`GetMetacellSuperCell()`](https://mdmanurung.github.io/cellprograms/reference/GetMetacellSuperCell.md)
  : Access metacell results stored by MetacellsByGroupsSuperCell
- [`MetacellsByGroupsSEACells()`](https://mdmanurung.github.io/cellprograms/reference/MetacellsByGroupsSEACells.md)
  : Construct metacells by groups with SEACells and mcRigor QC
- [`OptimizeGammaSEACells()`](https://mdmanurung.github.io/cellprograms/reference/OptimizeGammaSEACells.md)
  : Optimize SEACells granularity with mcRigor
- [`GetMetacellSEACells()`](https://mdmanurung.github.io/cellprograms/reference/GetMetacellSEACells.md)
  : Access metacell results stored by MetacellsByGroupsSEACells

## Program plots

Fit-level visualisations.

- [`plot_program_scores()`](https://mdmanurung.github.io/cellprograms/reference/plot_program_scores.md)
  : Plot program scores across donors
- [`plot_top_genes()`](https://mdmanurung.github.io/cellprograms/reference/plot_top_genes.md)
  : Plot top-loaded genes per program
- [`plot_variance_explained()`](https://mdmanurung.github.io/cellprograms/reference/plot_variance_explained.md)
  : Plot variance explained per program
- [`plot_artifact_flags()`](https://mdmanurung.github.io/cellprograms/reference/plot_artifact_flags.md)
  : Plot artifact-gene flags
- [`plot_stability()`](https://mdmanurung.github.io/cellprograms/reference/plot_stability.md)
  : Plot program stability
- [`plot_alignment()`](https://mdmanurung.github.io/cellprograms/reference/plot_alignment.md)
  : Plot program alignment (principal angles)

## Coordination & sharing plots

CMP and sharing visualisations.

- [`plot_cmp_network()`](https://mdmanurung.github.io/cellprograms/reference/plot_cmp_network.md)
  : Plot the coordinated-program network
- [`plot_cmp_nodes()`](https://mdmanurung.github.io/cellprograms/reference/plot_cmp_nodes.md)
  : Plot CMP node statistics
- [`plot_cmp_membership()`](https://mdmanurung.github.io/cellprograms/reference/plot_cmp_membership.md)
  : Plot CMP membership across settings
- [`plot_context_association()`](https://mdmanurung.github.io/cellprograms/reference/plot_context_association.md)
  : Plot CMP context associations
- [`plot_sharing()`](https://mdmanurung.github.io/cellprograms/reference/plot_sharing.md)
  : Plot program sharing results
- [`plot_sharing_spectrum()`](https://mdmanurung.github.io/cellprograms/reference/plot_sharing_spectrum.md)
  : Plot the sharing spectrum

## Internal

Internal helpers, S3 methods and package overview, not usually called
directly.

- [`cellprograms`](https://mdmanurung.github.io/cellprograms/reference/cellprograms-package.md)
  [`cellprograms-package`](https://mdmanurung.github.io/cellprograms/reference/cellprograms-package.md)
  : cellprograms: Cell-Type-Specific Transcriptional Programs via EBMF
- [`cmp_gene_overlap()`](https://mdmanurung.github.io/cellprograms/reference/cmp_gene_overlap.md)
  : Gene-overlap class of coordinated programs
- [`cmp_nuisance_genes()`](https://mdmanurung.github.io/cellprograms/reference/cmp_nuisance_genes.md)
  : Default nuisance-gene patterns for CMP triage
- [`cmp_scores()`](https://mdmanurung.github.io/cellprograms/reference/cmp_scores.md)
  : Consensus donor scores of coordinated programs
- [`.degenerate_programs()`](https://mdmanurung.github.io/cellprograms/reference/dot-degenerate_programs.md)
  : Numerically degenerate programs, two kinds: 1. flashier can retain
  factors whose loadings underflow to numerical noise (median
  \|loading\| many orders of magnitude below the cell type's strongest
  program; e.g. sparse factors whose support is a few genes and whose
  remaining entries sit at ~1e-16). Their posterior SDs underflow as
  well, producing garbage z-scores, and correlation-based matching is
  scale-invariant, so they spuriously match every strong program. 2.
  Fixed, never-updated factors: real-scale loadings whose posterior SD
  column collapsed to machine precision (observed for the smallest cell
  type in both real cohorts). Point estimates without usable
  uncertainty: mashr cannot assess them. A program is degenerate when
  median \|loading\| \< \`rel_tol\` x the cell type's largest, or when
  its median posterior SD (scaled to the loadings) is \<= \`sd_floor\`.
- [`.estimate_V()`](https://mdmanurung.github.io/cellprograms/reference/dot-estimate_V.md)
  : CorShrink-style null-correlation estimate for a matched program set
- [`.match_programs_loadings()`](https://mdmanurung.github.io/cellprograms/reference/dot-match_programs_loadings.md)
  : Candidate cross-cell-type program matches from loading correlation
- [`print(`*`<cell_program_data>`*`)`](https://mdmanurung.github.io/cellprograms/reference/print.cell_program_data.md)
  : Print a cell program data object
- [`print(`*`<cmp_result>`*`)`](https://mdmanurung.github.io/cellprograms/reference/print.cmp_result.md)
  : Print a CMP detection result
- [`print(`*`<cmp_tf>`*`)`](https://mdmanurung.github.io/cellprograms/reference/print.cmp_tf.md)
  : Print a CMP transcription-factor triage result
- [`print(`*`<principal_angles>`*`)`](https://mdmanurung.github.io/cellprograms/reference/print.principal_angles.md)
  : Print a principal-angles summary
- [`print(`*`<program_sharing>`*`)`](https://mdmanurung.github.io/cellprograms/reference/print.program_sharing.md)
  : Print a program-sharing summary
- [`program_metadata()`](https://mdmanurung.github.io/cellprograms/reference/program_metadata.md)
  : Program metadata table (one row per program)
- [`summary(`*`<cell_program_data>`*`)`](https://mdmanurung.github.io/cellprograms/reference/summary.cell_program_data.md)
  : Summarize a cell program data object
- [`union_panel_loadings()`](https://mdmanurung.github.io/cellprograms/reference/union_panel_loadings.md)
  : Refit program loadings on a common (union) gene panel
