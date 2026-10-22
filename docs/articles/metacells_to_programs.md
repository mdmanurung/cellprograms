# From metacells to programs: the SuperCell + mcRigor workflow

This article runs the **full `cellprograms` pipeline downstream of
metacell construction**. The starting point is a single-cell Seurat
object; metacells are built with the default SuperCell 2.0 + mcRigor
backend (pure R — no Python; see “Metacells: SuperCell (default) and
SEACells” for construction details and the benchmark behind the
default), and the metacell counts then feed program discovery, stability
assessment, coordinated multicellular program (CMP) detection, sharing
assessment, and cell-level projection.

The fixture plants known truth so every downstream call can be checked:
two cell types across 24 donors, each type with a **private program**
(constant across donors) and a **coordinated program** whose donor
activities are correlated across the two cell types.

## Simulate a single-cell dataset with planted programs

24 donors x two cell types x 60 cells, 800 genes. Genes 1:40 are a
T-cell-private marker block; genes 81:120 a B-cell-private block; genes
41:80 form the **coordinated block**, expressed in both cell types with
a donor activity gradient (lambda rising from 0.5 to 8.5 across donors)
that is *shared* across types — the planted CMP. Everything else is
background Poisson noise.

\
[`library`](https://rdrr.io/r/base/library.html)`(``cellprograms``)`\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``7``)`\
\
`n_donor`` ``<-`` ``24``; ``n_cell`` ``<-`` ``60``; ``n_gene`` ``<-`` ``800`\
`types`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``"Tcells"``, ``"Bcells"``)`\
`donors`` ``<-`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"d"``, `[`seq_len`](https://rdrr.io/r/base/seq.html)`(``n_donor``)``)`\
`a_grad`` ``<-`` `[`seq_len`](https://rdrr.io/r/base/seq.html)`(``n_donor``)``                ``# planted donor activity gradient`\
\
`groups`` ``<-`` `[`unlist`](https://rdrr.io/r/base/unlist.html)`(`[`lapply`](https://rdrr.io/r/base/lapply.html)`(``donors``, ``function``(``d``)`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``d``, ``"_"``, `[`rep`](https://rdrr.io/r/base/rep.html)`(``types``, each ``=`` ``n_cell``)``)``)``,`\
`                 use.names ``=`` ``FALSE``)`\
`lam`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(``0.5``, nrow ``=`` ``n_gene``, ncol ``=`` `[`length`](https://rdrr.io/r/base/length.html)`(``groups``)``)`\
`for`` ``(``d`` ``in`` `[`seq_len`](https://rdrr.io/r/base/seq.html)`(``n_donor``)``)`` ``{`\
`  ``for`` ``(``t`` ``in`` `[`seq_along`](https://rdrr.io/r/base/seq.html)`(``types``)``)`` ``{`\
`    ``cols`` ``<-`` `[`which`](https://rdrr.io/r/base/which.html)`(``groups`` ``==`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``donors``[``d``]``, ``"_"``, ``types``[``t``]``)``)`\
`    ``private`` ``<-`` ``if`` ``(``types``[``t``]`` ``==`` ``"Tcells"``)`` ``1``:``40`` ``else`` ``81``:``120`\
`    ``lam``[``private``, ``cols``]`` ``<-`` ``8``                          ``# private block, constant`\
`    ``lam``[``41``:``80``, ``cols``]`` ``<-`` ``0.5`` ``+`` ``8`` ``*`` ``(``d`` ``-`` ``1``)`` ``/`` ``(``n_donor`` ``-`` ``1``)``  ``# coordinated gradient`\
`  ``}`\
`}`\
`counts`` ``<-`` ``Matrix``::`[`Matrix`](https://rdrr.io/pkg/Matrix/man/Matrix.html)`(`[`rpois`](https://rdrr.io/r/stats/Poisson.html)`(``n_gene`` ``*`` `[`length`](https://rdrr.io/r/base/length.html)`(``groups``)``, `[`as.vector`](https://rdrr.io/r/base/vector.html)`(``lam``)``)``,`\
`                         nrow ``=`` ``n_gene``, ncol ``=`` `[`length`](https://rdrr.io/r/base/length.html)`(``groups``)``, sparse ``=`` ``TRUE``)`\
[`rownames`](https://rdrr.io/r/base/colnames.html)`(``counts``)`` ``<-`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"gene"``, `[`seq_len`](https://rdrr.io/r/base/seq.html)`(``n_gene``)``)`\
[`colnames`](https://rdrr.io/r/base/colnames.html)`(``counts``)`` ``<-`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"cell"``, `[`seq_len`](https://rdrr.io/r/base/seq.html)`(`[`length`](https://rdrr.io/r/base/length.html)`(``groups``)``)``)`\
\
`meta`` ``<-`` `[`data.frame`](https://rdrr.io/r/base/data.frame.html)`(`\
`  donor ``=`` `[`sub`](https://rdrr.io/r/base/grep.html)`(``"_.*"``, ``""``, ``groups``)``,`\
`  celltype ``=`` `[`sub`](https://rdrr.io/r/base/grep.html)`(``".*_"``, ``""``, ``groups``)``,`\
`  row.names ``=`` `[`colnames`](https://rdrr.io/r/base/colnames.html)`(``counts``)`\
`)`

## Build metacells (SuperCell 2.0 + mcRigor)

Metacells are built per donor x cell type so they never cross either
variable. `cells_per_metacell = 15` turns each 60-cell group into ~4
metacells; mcRigor QC runs on every group.

\
`obj`` ``<-`` ``Seurat``::`[`CreateSeuratObject`](https://satijalab.github.io/seurat-object/reference/CreateSeuratObject.html)`(``counts ``=`` ``counts``, meta.data ``=`` ``meta``)`\
`obj`` ``<-`` ``Seurat``::`[`NormalizeData`](https://satijalab.org/seurat/reference/NormalizeData.html)`(``obj``, verbose ``=`` ``FALSE``)`\
`obj`` ``<-`` ``Seurat``::`[`FindVariableFeatures`](https://satijalab.org/seurat/reference/FindVariableFeatures.html)`(``obj``, nfeatures ``=`` ``200``, verbose ``=`` ``FALSE``)`\
`obj`` ``<-`` ``Seurat``::`[`ScaleData`](https://satijalab.org/seurat/reference/ScaleData.html)`(``obj``, verbose ``=`` ``FALSE``)`\
`obj`` ``<-`` ``Seurat``::`[`RunPCA`](https://satijalab.org/seurat/reference/RunPCA.html)`(``obj``, npcs ``=`` ``10``, verbose ``=`` ``FALSE``)`\
\
`obj`` ``<-`` `[`MetacellsByGroupsSuperCell`](https://mdmanurung.github.io/cellprograms/reference/MetacellsByGroupsSuperCell.md)`(``obj``,`\
`                                  group.by ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"donor"``, ``"celltype"``)``,`\
`                                  cells_per_metacell ``=`` ``15``,`\
`                                  min_cells ``=`` ``50``,`\
`                                  n_top_genes ``=`` ``200``,`\
`                                  n_pc ``=`` ``10``,`\
`                                  mcrigor_args ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``Nrep ``=`` ``5``)``,`\
`                                  seed ``=`` ``7``, verbose ``=`` ``FALSE``)`\
`res`` ``<-`` `[`GetMetacellSuperCell`](https://mdmanurung.github.io/cellprograms/reference/GetMetacellSuperCell.md)`(``obj``)`\
`mc`` ``<-`` ``res``$``metacell_object`\
[`ncol`](https://rdrr.io/r/base/nrow.html)`(``mc``)``                          ``# total metacells`\
[`head`](https://rdrr.io/r/utils/head.html)`(``res``$``qc``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"metacell"``, ``"mcDiv"``, ``"dubious"``)``]``)`

[`OptimizeGammaSuperCell()`](https://mdmanurung.github.io/cellprograms/reference/OptimizeGammaSuperCell.md)
can replace the fixed `cells_per_metacell` when you want the granularity
chosen by the mcRigor score (see the metacells article); the downstream
workflow is identical either way.

## Metacell counts to donor pseudobulk

With `mode = "sum"` (the default), summing metacell counts within donor
x cell type reproduces the donor pseudobulk exactly — the metacell layer
is lossless at the donor level and simply denoises the cells the
pseudobulk is built from. We log-normalize the aggregated counts with
the standard `log1p(count / library size x 10^4)` recipe and declare the
fit’s measurement units as `"logcounts"` (this is what later lets us
project the programs back onto single cells).

\
`counts_mc`` ``<-`` ``SeuratObject``::`[`LayerData`](https://satijalab.github.io/seurat-object/reference/Layers.html)`(``mc``[[``"RNA"``]``]``, layer ``=`` ``"counts"``)``  ``# genes x metacells`\
`mc_donor`` ``<-`` ``mc``$``donor``                       ``# group columns are copied to the`\
`mc_type`` ``<-`` ``mc``$``celltype``                     ``# metacell object's metadata`\
\
`pb`` ``<-`` `[`lapply`](https://rdrr.io/r/base/lapply.html)`(``types``, ``function``(``t``)`` ``{`\
`  ``mcs`` ``<-`` `[`colnames`](https://rdrr.io/r/base/colnames.html)`(``mc``)``[``mc_type`` ``==`` ``t``]`\
`  ``agg`` ``<-`` `[`rowsum`](https://rdrr.io/r/base/rowsum.html)`(`[`t`](https://rdrr.io/r/base/t.html)`(`[`as.matrix`](https://rdrr.io/r/base/matrix.html)`(``counts_mc``[``, ``mcs``, drop ``=`` ``FALSE``]``)``)``,`\
`                group ``=`` ``mc_donor``[``mcs``]``)``              ``# donors x genes (counts)`\
`  ``agg`` ``<-`` ``agg``[``donors``, , drop ``=`` ``FALSE``]`\
`  ``lib`` ``<-`` `[`rowSums`](https://rdrr.io/r/base/colSums.html)`(``agg``)`\
`  ``out`` ``<-`` `[`log1p`](https://rdrr.io/r/base/Log.html)`(`[`sweep`](https://rdrr.io/r/base/sweep.html)`(``agg``, ``1``, ``lib``, ``"/"``)`` ``*`` ``1e4``)``       ``# donors x genes (logcounts)`\
`  `[`rownames`](https://rdrr.io/r/base/colnames.html)`(``out``)`` ``<-`` ``donors``   ``# same donor id across types: how CMP pairing`\
`  ``out``                       ``# knows which donors are shared`\
`}``)`\
[`names`](https://rdrr.io/r/base/names.html)`(``pb``)`` ``<-`` ``types`\
[`lapply`](https://rdrr.io/r/base/lapply.html)`(``pb``, ``dim``)`\
\
`sample_metadata`` ``<-`` `[`data.frame`](https://rdrr.io/r/base/data.frame.html)`(`\
`  observation_id ``=`` ``donors``,`\
`  planted_activity ``=`` ``a_grad`\
`)`\
`x`` ``<-`` `[`as_cell_program_data`](https://mdmanurung.github.io/cellprograms/reference/as_cell_program_data.md)`(``pb``, sample_metadata ``=`` ``sample_metadata``)`

## Fit and canonicalize programs

One empirical-Bayes matrix factorization per cell type; the coordinated
block’s donor-varying activity and the private blocks’ constant activity
should both be recovered as programs.

\
`fit`` ``<-`` `[`fit_celltype_programs`](https://mdmanurung.github.io/cellprograms/reference/fit_celltype_programs.md)`(``x``, loading_prior ``=`` ``"point_laplace"``,`\
`                             max_factors ``=`` ``3``, features ``=`` ``"variable"``,`\
`                             n_variable_genes ``=`` ``400``,`\
`                             feature_scale ``=`` ``"logcounts"``)`\
`fit`` ``<-`` `[`canonicalize_programs`](https://mdmanurung.github.io/cellprograms/reference/canonicalize_programs.md)`(``fit``)`\
[`program_summary`](https://mdmanurung.github.io/cellprograms/reference/program_summary.md)`(``fit``)``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"cell_type"``, ``"n_programs"``, ``"R2"``, ``"elbo"``)``]`

\
[`top_genes`](https://mdmanurung.github.io/cellprograms/reference/top_genes.md)`(``fit``, ``"Tcells_1"``, n ``=`` ``6``)`\
[`top_genes`](https://mdmanurung.github.io/cellprograms/reference/top_genes.md)`(``fit``, ``"Bcells_1"``, n ``=`` ``6``)`

The donor x program score matrix carries the planted gradient: the
coordinated programs’ scores track the donor index, the private programs
are flat.

\
`Z`` ``<-`` `[`program_scores`](https://mdmanurung.github.io/cellprograms/reference/program_scores.md)`(``fit``)`\
[`plot_program_scores`](https://mdmanurung.github.io/cellprograms/reference/plot_program_scores.md)`(``fit``, annotation ``=`` ``"planted_activity"``)`\
[`plot_variance_explained`](https://mdmanurung.github.io/cellprograms/reference/plot_variance_explained.md)`(``fit``)`

## Stability

[`assess_program_stability()`](https://mdmanurung.github.io/cellprograms/reference/assess_program_stability.md)
refits on donor subsamples and reports how often each program is
recovered — the same guarantee as in the pseudobulk workflow, now on
metacell-derived inputs.

\
`fit`` ``<-`` `[`assess_program_stability`](https://mdmanurung.github.io/cellprograms/reference/assess_program_stability.md)`(``fit``, n_boot ``=`` ``4``, resample_frac ``=`` ``0.7``,`\
`                                seed ``=`` ``7``, progress ``=`` ``FALSE``)`\
`fit``$``stability``$``per_program``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"program_id"``, ``"recovery_freq"``)``]`

## Coordinated multicellular programs

[`find_coordinated_programs()`](https://mdmanurung.github.io/cellprograms/reference/find_coordinated_programs.md)
tests whether programs from *different* cell types co-vary across
donors. The planted truth: the two programs of the coordinated block
(one per type) share donor activity; the private programs do not.

\
`cmp`` ``<-`` `[`find_coordinated_programs`](https://mdmanurung.github.io/cellprograms/reference/find_coordinated_programs.md)`(``fit``, n_perm ``=`` ``99``, seed ``=`` ``7``)`\
`cmp`\
[`plot_cmp_network`](https://mdmanurung.github.io/cellprograms/reference/plot_cmp_network.md)`(``cmp``)`

The coordinated-block pair is called; the private programs (constant
donor activity, independent noise) are not.
[`annotate_cmp_tf()`](https://mdmanurung.github.io/cellprograms/reference/annotate_cmp_tf.md)
/
[`cmp_triage()`](https://mdmanurung.github.io/cellprograms/reference/cmp_triage.md)
add TF-footprint annotation and triage on real data (external CollecTRI
download, not run here).

## Shared vs specific programs

[`program_sharing()`](https://mdmanurung.github.io/cellprograms/reference/program_sharing.md)
asks whether matched programs share their gene *effects* — here the
coordinated block is literally the same genes in both types, so a shared
call on the coordinated pair is the planted truth.

\
`if`` ``(`[`requireNamespace`](https://rdrr.io/r/base/ns-load.html)`(``"mashr"``, quietly ``=`` ``TRUE``)``)`` ``{`\
`  ``shr`` ``<-`` `[`program_sharing`](https://mdmanurung.github.io/cellprograms/reference/program_sharing.md)`(``fit``, n_perm ``=`` ``99``, seed ``=`` ``7``)`\
`  ``shr``$``matches`\
`  `[`plot_sharing`](https://mdmanurung.github.io/cellprograms/reference/plot_sharing.md)`(``shr``)`\
`}`

## Project programs back onto single cells

Because the fit declared `feature_scale = "logcounts"`, the fitted
programs can be projected onto the original single cells
([`project_programs_to_cells()`](https://mdmanurung.github.io/cellprograms/reference/project_programs_to_cells.md)
replays the fit’s preprocessing and infers per-cell activities against
the frozen loadings). The per-cell activity of a coordinated program
should recover the planted donor gradient — this is the bridge from
donor-level discovery to cell-level readout.

\
`logcounts_sc`` ``<-`` `[`as.matrix`](https://rdrr.io/r/base/matrix.html)`(``SeuratObject``::`[`LayerData`](https://satijalab.github.io/seurat-object/reference/Layers.html)`(``obj``[[``"RNA"``]``]``, layer ``=`` ``"data"``)``)`\
`logcounts_sc`` ``<-`` `[`t`](https://rdrr.io/r/base/t.html)`(``logcounts_sc``)``                   ``# cells x genes (dense)`\
`proj`` ``<-`` `[`project_programs_to_cells`](https://mdmanurung.github.io/cellprograms/reference/project_programs_to_cells.md)`(``fit``, expression ``=`` ``logcounts_sc``,`\
`                                  cell_types ``=`` `[`setNames`](https://rdrr.io/r/stats/setNames.html)`(``meta``$``celltype``,`\
`                                                        `[`rownames`](https://rdrr.io/r/base/colnames.html)`(``meta``)``)``)`\
`act`` ``<-`` ``proj``$``activity`\
`coord_prog`` ``<-`` `[`grep`](https://rdrr.io/r/base/grep.html)`(``"^Tcells"``, `[`colnames`](https://rdrr.io/r/base/colnames.html)`(``act``)``, value ``=`` ``TRUE``)``[``1``]`\
`# Activity is type-specific: cells outside the program's source type get NA,`\
`# so donor means are computed within T cells only.`\
`is_t`` ``<-`` ``meta``$``celltype`` ``==`` ``"Tcells"`\
`donor_mean`` ``<-`` `[`tapply`](https://rdrr.io/r/base/tapply.html)`(``act``[``is_t``, ``coord_prog``]``, ``meta``$``donor``[``is_t``]``, ``mean``)`\
[`round`](https://rdrr.io/r/base/Round.html)`(`[`cor`](https://rdrr.io/r/stats/cor.html)`(``donor_mean``[``donors``]``, ``a_grad``, method ``=`` ``"spearman"``)``, ``3``)`

The high rank correlation between the projected per-cell activity
(averaged within donor) and the planted gradient confirms the metacell
-\> donor fit -\> cell projection round trip preserves the biology.
Absolute activity levels live on the fit’s scale; the donor ordering is
the readout.

## Session info
