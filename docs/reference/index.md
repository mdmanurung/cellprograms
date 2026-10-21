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

## Metacells (SEACells + mcRigor)

Optional upstream aggregation of single cells into metacells before
program discovery. Requires the suggested packages Seurat, mcRigor and
basilisk (the SEACells Python environment is built on first use).

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
