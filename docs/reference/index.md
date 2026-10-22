<div class="page-header">

# Package index

</div>

<div class="section level2">

## Data & simulation<a href="#data-simulation" class="anchor" aria-label="anchor"></a>

<div class="section-desc">

Build and check `cell_program_data` objects, simulate pseudobulk ground truth, and bridge to MOFAcellular.

</div>

</div>

<div class="section level2">

[`simulate_pseudobulk()`](simulate_pseudobulk.html)<br>
Simulate cell-type pseudobulk matrices with planted programs

[`as_cell_program_data()`](as_cell_program_data.html)<br>
Create a cell program data object

[`validate_cell_program_data()`](validate_cell_program_data.html)<br>
Validate a cell program data object

[`artifact_gene_sets`](artifact_gene_sets.html)<br>
Artifact-gene sets

[`artifact_genes()`](artifact_genes.html)<br>
Expand artifact-gene sets against a gene universe

[`as_mofacellular_input()`](as_mofacellular_input.html)<br>
Export pseudobulk profiles to MOFAcellulaR input format

</div>

<div class="section level2">

## Fit & canonicalize<a href="#fit-canonicalize" class="anchor" aria-label="anchor"></a>

<div class="section-desc">

Fit per-cell-type programs (flashier R backend or cebmf_torch Python backend), canonicalize sign/scale, freeze preprocessing, and flag artifact-carrying programs.

</div>

</div>

<div class="section level2">

[`fit_celltype_programs()`](fit_celltype_programs.html)<br>
Fit cell-type-specific programs with flashier

[`canonicalize_programs()`](canonicalize_programs.html)<br>
Canonicalize program scale and sign

[`freeze_preprocessing()`](freeze_preprocessing.html)<br>
Freeze a preprocessing specification

[`flag_artifact_genes()`](flag_artifact_genes.html)<br>
Flag programs or CMPs carrying artifact genes

</div>

<div class="section level2">

## Extract<a href="#extract" class="anchor" aria-label="anchor"></a>

<div class="section-desc">

Summaries, scores, loadings and top genes from a canonicalized fit.

</div>

</div>

<div class="section level2">

[`program_summary()`](program_summary.html)<br>
Summary of programs per cell type (n, K, R2 of the reconstruction, ELBO)

[`program_scores()`](program_scores.html)<br>
Program score matrix across all cell types

[`program_loadings()`](program_loadings.html)<br>
Program loading matrices

[`top_genes()`](top_genes.html)<br>
Top genes of a program

[`project_programs_to_cells()`](project_programs_to_cells.html)<br>
Project fitted programs onto individual cells

</div>

<div class="section level2">

## Compare, align & share<a href="#compare-align-share" class="anchor" aria-label="anchor"></a>

<div class="section-desc">

Program stability, cross-fit matching, within-fit alignment, and donor-level sharing structure.

</div>

</div>

<div class="section level2">

[`assess_program_stability()`](assess_program_stability.html)<br>
Assess program stability by donor subsampling

[`match_programs()`](match_programs.html)<br>
Match programs across two fits (bootstrap refits or cohorts)

[`align_programs()`](align_programs.html)<br>
Align all cell types' program subspaces to a reference cell type

[`span_alignment()`](span_alignment.html)<br>
Principal-angle alignment between two loading subspaces

[`program_sharing()`](program_sharing.html)<br>
Post-hoc shared-vs-specific assessment of fitted programs (mashr)

[`sharing_spectrum()`](sharing_spectrum.html)<br>
Pairwise sharing spectra across all cell types

[`principal_angles()`](principal_angles.html)<br>
Principal angles between the program subspaces of two cell types

</div>

<div class="section level2">

## Coordination & annotation<a href="#coordination-annotation" class="anchor" aria-label="anchor"></a>

<div class="section-desc">

Cross-cell-type coordinated programs (CMPs): detection, triage, and annotation against TF regulons, ligand-receptor pairs and gene sets.

</div>

</div>

<div class="section level2">

[`find_coordinated_programs()`](find_coordinated_programs.html)<br>
Find coordinated multicellular programs

[`cmp_triage()`](cmp_triage.html)<br>
One-row-per-CMP triage table

[`annotate_cmp_tf()`](annotate_cmp_tf.html)<br>
Transcription-factor footprints of coordinated programs

[`annotate_cmp_lr()`](annotate_cmp_lr.html)<br>
Ligand-receptor compatibility of coordinated programs

[`annotate_cmp_themes()`](annotate_cmp_themes.html)<br>
Themes of coordinated programs

[`assess_cmp_stability()`](assess_cmp_stability.html)<br>
Assess CMP stability under synchronized donor resampling

[`prepare_collectri()`](prepare_collectri.html)<br>
Format an OmniPath CollecTRI table as a signed TF-target network

[`ambient_program_scores()`](ambient_program_scores.html)<br>
Ambient-likeness of programs

[`cmp_context_association()`](cmp_context_association.html)<br>
Association of coordinated programs with a donor variable

</div>

<div class="section level2">

## Joint & tensor factorization<a href="#joint-tensor-factorization" class="anchor" aria-label="anchor"></a>

<div class="section-desc">

Post-fit joint factorization across cell types: empirical-Bayes CP tensor factorization (genes x donors x cell types) and block-structured joint EBMF at expression or scores level, with benchmark-contract coordination calls. See `vignette("coordination_joint")` and the preregistered benchmark under `benchmarks/coordination/`.

</div>

</div>

<div class="section level2">

[`fit_eb_cp()`](fit_eb_cp.html)<br>
EB-CP tensor factorization across cell types

[`eb_cp_from_fit()`](eb_cp_from_fit.html)<br>
EB-CP coordination programs from a fitted \`cell_program_fit\`

[`eb_cp_programs()`](eb_cp_programs.html)<br>
Coordination programs from an EB-CP fit

[`fit_joint_block()`](fit_joint_block.html)<br>
Block-structured joint EBMF (Level 3 candidate A)

[`block_clusters()`](block_clusters.html)<br>
Coordination calls from a block_joint fit (benchmark output contract)

</div>

<div class="section level2">

## Metacells (SuperCell + mcRigor)<a href="#metacells-supercell-mcrigor-" class="anchor" aria-label="anchor"></a>

<div class="section-desc">

Optional upstream aggregation of single cells into metacells before program discovery. The default backend is SuperCell 2.0 (pure R, no Python) with mcRigor QC and granularity optimization. The SEACells backend remains available and additionally requires basilisk (the SEACells Python environment is built on first use).

</div>

</div>

<div class="section level2">

[`MetacellsByGroupsSuperCell()`](MetacellsByGroupsSuperCell.html)<br>
Construct metacells by groups with SuperCell and mcRigor QC

[`OptimizeGammaSuperCell()`](OptimizeGammaSuperCell.html)<br>
Optimize SuperCell granularity with mcRigor

[`GetMetacellSuperCell()`](GetMetacellSuperCell.html)<br>
Access metacell results stored by MetacellsByGroupsSuperCell

[`MetacellsByGroupsSEACells()`](MetacellsByGroupsSEACells.html)<br>
Construct metacells by groups with SEACells and mcRigor QC

[`OptimizeGammaSEACells()`](OptimizeGammaSEACells.html)<br>
Optimize SEACells granularity with mcRigor

[`GetMetacellSEACells()`](GetMetacellSEACells.html)<br>
Access metacell results stored by MetacellsByGroupsSEACells

</div>

<div class="section level2">

## Program plots<a href="#program-plots" class="anchor" aria-label="anchor"></a>

<div class="section-desc">

Fit-level visualisations.

</div>

</div>

<div class="section level2">

[`plot_program_scores()`](plot_program_scores.html)<br>
Plot program scores across donors

[`plot_top_genes()`](plot_top_genes.html)<br>
Plot top-loaded genes per program

[`plot_variance_explained()`](plot_variance_explained.html)<br>
Plot variance explained per program

[`plot_artifact_flags()`](plot_artifact_flags.html)<br>
Plot artifact-gene flags

[`plot_stability()`](plot_stability.html)<br>
Plot program stability

[`plot_alignment()`](plot_alignment.html)<br>
Plot program alignment (principal angles)

</div>

<div class="section level2">

## Coordination & sharing plots<a href="#coordination-sharing-plots" class="anchor" aria-label="anchor"></a>

<div class="section-desc">

CMP and sharing visualisations.

</div>

</div>

<div class="section level2">

[`plot_cmp_network()`](plot_cmp_network.html)<br>
Plot the coordinated-program network

[`plot_cmp_nodes()`](plot_cmp_nodes.html)<br>
Plot CMP node statistics

[`plot_cmp_membership()`](plot_cmp_membership.html)<br>
Plot CMP membership across settings

[`plot_context_association()`](plot_context_association.html)<br>
Plot CMP context associations

[`plot_sharing()`](plot_sharing.html)<br>
Plot program sharing results

[`plot_sharing_spectrum()`](plot_sharing_spectrum.html)<br>
Plot the sharing spectrum

</div>

<div class="section level2">

## Internal<a href="#internal" class="anchor" aria-label="anchor"></a>

<div class="section-desc">

Internal helpers, S3 methods and package overview, not usually called directly.

</div>

</div>

<div class="section level2">

[`cellprograms`](cellprograms-package.html) [`cellprograms-package`](cellprograms-package.html)<br>
cellprograms: Cell-Type-Specific Transcriptional Programs via EBMF

[`cmp_gene_overlap()`](cmp_gene_overlap.html)<br>
Gene-overlap class of coordinated programs

[`cmp_nuisance_genes()`](cmp_nuisance_genes.html)<br>
Default nuisance-gene patterns for CMP triage

[`cmp_scores()`](cmp_scores.html)<br>
Consensus donor scores of coordinated programs

[`.degenerate_programs()`](dot-degenerate_programs.html)<br>
Numerically degenerate programs, two kinds: 1. flashier can retain factors whose loadings underflow to numerical noise (median \|loading\| many orders of magnitude below the cell type's strongest program; e.g. sparse factors whose support is a few genes and whose remaining entries sit at ~1e-16). Their posterior SDs underflow as well, producing garbage z-scores, and correlation-based matching is scale-invariant, so they spuriously match every strong program. 2. Fixed, never-updated factors: real-scale loadings whose posterior SD column collapsed to machine precision (observed for the smallest cell type in both real cohorts). Point estimates without usable uncertainty: mashr cannot assess them. A program is degenerate when median \|loading\| \< \`rel_tol\` x the cell type's largest, or when its median posterior SD (scaled to the loadings) is \<= \`sd_floor\`.

[`.estimate_V()`](dot-estimate_V.html)<br>
CorShrink-style null-correlation estimate for a matched program set

[`.match_programs_loadings()`](dot-match_programs_loadings.html)<br>
Candidate cross-cell-type program matches from loading correlation

[`print(`*`<cell_program_data>`*`)`](print.cell_program_data.html)<br>
Print a cell program data object

[`print(`*`<cmp_result>`*`)`](print.cmp_result.html)<br>
Print a CMP detection result

[`print(`*`<cmp_tf>`*`)`](print.cmp_tf.html)<br>
Print a CMP transcription-factor triage result

[`print(`*`<principal_angles>`*`)`](print.principal_angles.html)<br>
Print a principal-angles summary

[`print(`*`<program_sharing>`*`)`](print.program_sharing.html)<br>
Print a program-sharing summary

[`program_metadata()`](program_metadata.html)<br>
Program metadata table (one row per program)

[`summary(`*`<cell_program_data>`*`)`](summary.cell_program_data.html)<br>
Summarize a cell program data object

[`union_panel_loadings()`](union_panel_loadings.html)<br>
Refit program loadings on a common (union) gene panel

</div>
