# Metacells with SEACells and mcRigor gamma optimization

When the starting point is single-cell data, `cellprograms` can
aggregate cells into **metacells** — transcriptomically homogeneous
small groups — before pseudobulk factorization. Metacells are built with
[SEACells](https://github.com/dpeerlab/SEACells) (run in an isolated,
pinned basilisk Python environment built on first use) and QC’d with
mcRigor, stratified so that **metacells never cross the grouping
variables** (donor x cell type by default).

This vignette runs the full workflow on a simulated dataset with planted
cell-type identity, then feeds the metacell counts into
[`fit_celltype_programs()`](https://mdmanurung.github.io/cellprograms/reference/fit_celltype_programs.md).

## Simulate a small single-cell dataset

Two donors x two cell types, 80 cells each, 200 genes. Each cell type
gets a private marker block; everything else is background Poisson
noise.

\
[`library`](https://rdrr.io/r/base/library.html)`(``cellprograms``)`\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``11``)`\
\
`n_per`` ``<-`` ``80`\
`types`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``"Tcells"``, ``"Bcells"``)`\
`donors`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``"donor1"``, ``"donor2"``)`\
`groups`` ``<-`` `[`rep`](https://rdrr.io/r/base/rep.html)`(`[`paste0`](https://rdrr.io/r/base/paste.html)`(`[`rep`](https://rdrr.io/r/base/rep.html)`(``donors``, each ``=`` `[`length`](https://rdrr.io/r/base/length.html)`(``types``)``)``, ``"_"``, `[`rep`](https://rdrr.io/r/base/rep.html)`(``types``, `[`length`](https://rdrr.io/r/base/length.html)`(``donors``)``)``)``,`\
`              each ``=`` ``n_per``)`\
\
`lam`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(``0.5``, nrow ``=`` ``200``, ncol ``=`` `[`length`](https://rdrr.io/r/base/length.html)`(``groups``)``)`\
`markers`` ``<-`` `[`list`](https://rdrr.io/r/base/list.html)`(``Tcells ``=`` ``1``:``12``, Bcells ``=`` ``13``:``24``)`\
`for`` ``(``g`` ``in`` `[`unique`](https://rdrr.io/r/base/unique.html)`(``groups``)``)`` ``{`\
`  ``type`` ``<-`` `[`sub`](https://rdrr.io/r/base/grep.html)`(``".*_"``, ``""``, ``g``)`\
`  ``cols`` ``<-`` `[`which`](https://rdrr.io/r/base/which.html)`(``groups`` ``==`` ``g``)`\
`  ``lam``[``markers``[[``type``]``]``, ``cols``]`` ``<-`` ``8`\
`}`\
`counts`` ``<-`` ``Matrix``::`[`Matrix`](https://rdrr.io/pkg/Matrix/man/Matrix.html)`(`[`rpois`](https://rdrr.io/r/stats/Poisson.html)`(``200`` ``*`` `[`length`](https://rdrr.io/r/base/length.html)`(``groups``)``, `[`as.vector`](https://rdrr.io/r/base/vector.html)`(``lam``)``)``,`\
`                         nrow ``=`` ``200``, ncol ``=`` `[`length`](https://rdrr.io/r/base/length.html)`(``groups``)``, sparse ``=`` ``TRUE``)`\
[`rownames`](https://rdrr.io/r/base/colnames.html)`(``counts``)`` ``<-`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"gene"``, `[`seq_len`](https://rdrr.io/r/base/seq.html)`(``200``)``)`\
[`colnames`](https://rdrr.io/r/base/colnames.html)`(``counts``)`` ``<-`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"cell"``, `[`seq_len`](https://rdrr.io/r/base/seq.html)`(`[`length`](https://rdrr.io/r/base/length.html)`(``groups``)``)``)`\
\
`meta`` ``<-`` `[`data.frame`](https://rdrr.io/r/base/data.frame.html)`(`\
`  donor ``=`` `[`sub`](https://rdrr.io/r/base/grep.html)`(``"_.*"``, ``""``, ``groups``)``,`\
`  celltype ``=`` `[`sub`](https://rdrr.io/r/base/grep.html)`(``".*_"``, ``""``, ``groups``)``,`\
`  row.names ``=`` `[`colnames`](https://rdrr.io/r/base/colnames.html)`(``counts``)`\
`)`\
\
`obj`` ``<-`` ``Seurat``::`[`CreateSeuratObject`](https://satijalab.github.io/seurat-object/reference/CreateSeuratObject.html)`(``counts ``=`` ``counts``, meta.data ``=`` ``meta``)`\
`obj`` ``<-`` ``Seurat``::`[`NormalizeData`](https://satijalab.org/seurat/reference/NormalizeData.html)`(``obj``, verbose ``=`` ``FALSE``)`\
`obj`` ``<-`` ``Seurat``::`[`FindVariableFeatures`](https://satijalab.org/seurat/reference/FindVariableFeatures.html)`(``obj``, nfeatures ``=`` ``60``, verbose ``=`` ``FALSE``)`\
`obj`` ``<-`` ``Seurat``::`[`ScaleData`](https://satijalab.org/seurat/reference/ScaleData.html)`(``obj``, verbose ``=`` ``FALSE``)`\
`obj`` ``<-`` ``Seurat``::`[`RunPCA`](https://satijalab.org/seurat/reference/RunPCA.html)`(``obj``, npcs ``=`` ``10``, verbose ``=`` ``FALSE``)`

## Build metacells per group

[`MetacellsByGroupsSEACells()`](https://mdmanurung.github.io/cellprograms/reference/MetacellsByGroupsSEACells.md)
fits SEACells separately within each `group.by` stratum and stores
everything under the Seurat object’s `misc`:

\
`obj`` ``<-`` `[`MetacellsByGroupsSEACells`](https://mdmanurung.github.io/cellprograms/reference/MetacellsByGroupsSEACells.md)`(``obj``,`\
`                                 group.by ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"donor"``, ``"celltype"``)``,`\
`                                 cells_per_metacell ``=`` ``20``,`\
`                                 min_cells ``=`` ``50``,`\
`                                 n_top_genes ``=`` ``60``,`\
`                                 n_pc ``=`` ``10``,`\
`                                 seed ``=`` ``11``, verbose ``=`` ``FALSE``)`\
`res`` ``<-`` `[`GetMetacellSEACells`](https://mdmanurung.github.io/cellprograms/reference/GetMetacellSEACells.md)`(``obj``)`\
[`names`](https://rdrr.io/r/base/names.html)`(``res``)`

The result carries the metacell Seurat object, the soft (cells x
archetypes) assignment matrix, per-metacell mcRigor QC, fit statistics,
the parameters used, and any skipped groups:

\
[`head`](https://rdrr.io/r/utils/head.html)`(``res``$``qc``)``                 ``# mcDiv, dubious flag, ZeroRate per metacell`\
`res``$``skipped_groups``           ``# groups dropped (e.g. too few cells) - usually empty`

Groups with too few cells (below `min_cells`) are skipped with a warning
rather than failing the call, and degenerate SEACells fits (all hard
assignments collapsing to one archetype) drop that group’s soft
assignments gracefully — check `skipped_groups` when working with small
strata.

## Optimize the number of metacells (gamma)

SEACells’ `gamma` is the number of archetypes — the granularity.
[`OptimizeGammaSEACells()`](https://mdmanurung.github.io/cellprograms/reference/OptimizeGammaSEACells.md)
refits SEACells across candidate gammas per group and picks the winner
by the mcRigor optimization score:

\
`opt`` ``<-`` `[`OptimizeGammaSEACells`](https://mdmanurung.github.io/cellprograms/reference/OptimizeGammaSEACells.md)`(``obj``,`\
`                             group.by ``=`` ``"donor"``,`\
`                             gammas ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``10``, ``20``)``,`\
`                             n_top_genes ``=`` ``60``,`\
`                             n_pc ``=`` ``10``,`\
`                             min_cells ``=`` ``50``,`\
`                             mcrigor_args ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``Nrep ``=`` ``5``)``,`\
`                             seed ``=`` ``11``, verbose ``=`` ``FALSE``)`\
`opt``$``optimal_gamma``            ``# winning gamma per donor group`\
`opt``$``scores``                   ``# mcRigor score per gamma x group`

Gammas whose membership labels degenerate (fewer than two distinct
labels after Seurat v5 aggregation) are dropped with a warning and
excluded from the pick rather than crashing the run.

## Metacell counts feed the program fit directly

\
`mc`` ``<-`` ``res``$``metacell_object`\
`counts_mc`` ``<-`` ``SeuratObject``::`[`LayerData`](https://satijalab.github.io/seurat-object/reference/Layers.html)`(``mc``[[``"RNA"``]``]``, layer ``=`` ``"counts"``)``  ``# genes x metacells`\
`by_type`` ``<-`` `[`split`](https://rdrr.io/r/base/split.html)`(`[`colnames`](https://rdrr.io/r/base/colnames.html)`(``mc``)``, ``mc``$``celltype``)`\
`pb`` ``<-`` `[`lapply`](https://rdrr.io/r/base/lapply.html)`(``by_type``, ``function``(``cells``)`` `[`t`](https://rdrr.io/r/base/t.html)`(`[`as.matrix`](https://rdrr.io/r/base/matrix.html)`(``counts_mc``[``, ``cells``, drop ``=`` ``FALSE``]``)``)``)`\
\
`x`` ``<-`` `[`as_cell_program_data`](https://mdmanurung.github.io/cellprograms/reference/as_cell_program_data.md)`(``pb``)`\
`fit`` ``<-`` `[`fit_celltype_programs`](https://mdmanurung.github.io/cellprograms/reference/fit_celltype_programs.md)`(``x``, loading_prior ``=`` ``"point_laplace"``,`\
`                             max_factors ``=`` ``4``, features ``=`` ``"variable"``,`\
`                             n_variable_genes ``=`` ``100``)`\
`fit`` ``<-`` `[`canonicalize_programs`](https://mdmanurung.github.io/cellprograms/reference/canonicalize_programs.md)`(``fit``)`\
[`program_summary`](https://mdmanurung.github.io/cellprograms/reference/program_summary.md)`(``fit``)``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"cell_type"``, ``"n_programs"``, ``"R2"``, ``"elbo"``)``]`

The planted structure survives the aggregation: each cell type’s
programs are dominated by its private marker block.

\
`marker_expr`` ``<-`` `[`as.matrix`](https://rdrr.io/r/base/matrix.html)`(``counts_mc``[`[`unlist`](https://rdrr.io/r/base/unlist.html)`(``markers``)``, , drop ``=`` ``FALSE``]``)`\
`marker_expr`` ``<-`` `[`t`](https://rdrr.io/r/base/t.html)`(`[`scale`](https://rdrr.io/r/base/scale.html)`(`[`t`](https://rdrr.io/r/base/t.html)`(``marker_expr``)``)``)`\
`df`` ``<-`` `[`data.frame`](https://rdrr.io/r/base/data.frame.html)`(`\
`  metacell ``=`` `[`factor`](https://rdrr.io/r/base/factor.html)`(`[`colnames`](https://rdrr.io/r/base/colnames.html)`(``marker_expr``)``, levels ``=`` `[`colnames`](https://rdrr.io/r/base/colnames.html)`(``marker_expr``)``)``,`\
`  gene ``=`` `[`rep`](https://rdrr.io/r/base/rep.html)`(`[`rownames`](https://rdrr.io/r/base/colnames.html)`(``marker_expr``)``, times ``=`` `[`ncol`](https://rdrr.io/r/base/nrow.html)`(``marker_expr``)``)``,`\
`  z ``=`` `[`as.vector`](https://rdrr.io/r/base/vector.html)`(``marker_expr``)`\
`)`\
`ggplot2``::`[`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)`(``df``, ``ggplot2``::`[`aes`](https://ggplot2.tidyverse.org/reference/aes.html)`(``metacell``, ``gene``, fill ``=`` ``z``)``)`` ``+`\
`  ``ggplot2``::`[`geom_tile`](https://ggplot2.tidyverse.org/reference/geom_tile.html)`(``)`` ``+`\
`  ``ggplot2``::`[`scale_fill_gradient2`](https://ggplot2.tidyverse.org/reference/scale_gradient.html)`(``)`` ``+`\
`  ``ggplot2``::`[`labs`](https://ggplot2.tidyverse.org/reference/labs.html)`(``x ``=`` ``"metacell"``, y ``=`` ``"marker gene"``)`` ``+`\
`  ``ggplot2``::`[`theme`](https://ggplot2.tidyverse.org/reference/theme.html)`(``axis.text.x ``=`` ``ggplot2``::`[`element_blank`](https://ggplot2.tidyverse.org/reference/element.html)`(``)``,`\
`                 axis.text.y ``=`` ``ggplot2``::`[`element_text`](https://ggplot2.tidyverse.org/reference/element.html)`(``size ``=`` ``6``)``)`

## Session info

    #> R version 4.4.2 (2024-10-31)
    #> Platform: x86_64-conda-linux-gnu
    #> Running under: Ubuntu 24.04.5 LTS
    #> 
    #> Matrix products: default
    #> BLAS/LAPACK: /opt/conda/lib/libopenblasp-r0.3.34.so;  LAPACK version 3.12.0
    #> 
    #> locale:
    #>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    #>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    #>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    #> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    #> 
    #> time zone: Etc/UTC
    #> tzcode source: system (glibc)
    #> 
    #> attached base packages:
    #> [1] stats     graphics  grDevices utils     datasets  methods   base     
    #> 
    #> other attached packages:
    #> [1] cellprograms_0.2.0.9001
    #> 
    #> loaded via a namespace (and not attached):
    #>   [1] deldir_2.0-4           pbapply_1.7-4          gridExtra_2.3.1       
    #>   [4] rlang_1.3.0            magrittr_2.0.5         RcppAnnoy_0.0.23      
    #>   [7] otel_0.2.0             spatstat.geom_3.8-2    matrixStats_1.5.0     
    #>  [10] ggridges_0.5.7         compiler_4.4.2         png_0.1-9             
    #>  [13] systemfonts_1.3.2      vctrs_0.7.3            reshape2_1.4.5        
    #>  [16] stringr_1.6.0          pkgconfig_2.0.3        fastmap_1.2.0         
    #>  [19] promises_1.5.0         rmarkdown_2.32         ragg_1.5.2            
    #>  [22] purrr_1.2.2            xfun_0.61              cachem_1.1.0          
    #>  [25] jsonlite_2.0.0         goftest_1.2-3          later_1.4.8           
    #>  [28] spatstat.utils_3.2-5   irlba_2.4.1            parallel_4.4.2        
    #>  [31] cluster_2.1.8.3        R6_2.6.1               ica_1.0-3             
    #>  [34] spatstat.data_3.1-9    bslib_0.12.0           stringi_1.8.7         
    #>  [37] RColorBrewer_1.1-3     reticulate_1.47.0      spatstat.univar_3.2-0 
    #>  [40] parallelly_1.48.0      lmtest_0.9-40          jquerylib_0.1.4       
    #>  [43] scattermore_1.2        Rcpp_1.1.2             knitr_1.52            
    #>  [46] tensor_1.5.1           future.apply_1.21.0    zoo_1.9-0             
    #>  [49] sctransform_0.4.3      httpuv_1.6.17          Matrix_1.7-6          
    #>  [52] splines_4.4.2          igraph_2.1.4           tidyselect_1.2.1      
    #>  [55] abind_1.4-8            yaml_2.3.12            spatstat.random_3.5-1 
    #>  [58] spatstat.explore_3.8-2 codetools_0.2-20       miniUI_0.1.2          
    #>  [61] listenv_1.1.0          lattice_0.23-1         tibble_3.3.1          
    #>  [64] plyr_1.8.9             shiny_1.14.0           S7_0.2.2              
    #>  [67] ROCR_1.0-12            evaluate_1.0.5         Rtsne_0.17            
    #>  [70] future_1.80.0          fastDummies_1.7.6      desc_1.4.3            
    #>  [73] survival_3.8-12        polyclip_1.10-7        fitdistrplus_1.2-6    
    #>  [76] pillar_1.11.1          Seurat_5.5.1           KernSmooth_2.23-27    
    #>  [79] plotly_4.12.1          generics_0.1.4         RcppHNSW_0.7.0        
    #>  [82] sp_2.2-3               ggplot2_4.0.3          scales_1.4.0          
    #>  [85] globals_0.20.0         xtable_1.8-8           glue_1.8.1            
    #>  [88] tools_4.4.2            data.table_1.18.6.1    RSpectra_0.16-2       
    #>  [91] RANN_2.6.2             fs_2.1.0               dotCall64_1.2         
    #>  [94] cowplot_1.2.0          grid_4.4.2             tidyr_1.3.2           
    #>  [97] nlme_3.1-170           patchwork_1.3.2        cli_3.6.6             
    #> [100] spatstat.sparse_3.2-0  textshaping_0.4.0      spam_2.11-4           
    #> [103] viridisLite_0.4.3      dplyr_1.2.1            uwot_0.2.5            
    #> [106] gtable_0.3.6           sass_0.4.10            digest_0.6.39         
    #> [109] progressr_1.1.0        ggrepel_0.9.8          htmlwidgets_1.6.4     
    #> [112] SeuratObject_5.4.0     farver_2.1.2           htmltools_0.5.9       
    #> [115] pkgdown_2.2.1          lifecycle_1.0.5        httr_1.4.9            
    #> [118] mime_0.13              MASS_7.3-66
