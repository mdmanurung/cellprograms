<div class="page-header">

# From metacells to programs: the SuperCell + mcRigor workflow

<div class="d-none name">

`metacells_to_programs.Rmd`

</div>

</div>

This article runs the **full `cellprograms` pipeline downstream of metacell construction**. The starting point is a single-cell Seurat object; metacells are built with the default SuperCell 2.0 + mcRigor backend (pure R — no Python; see “Metacells: SuperCell (default) and SEACells” for construction details and the benchmark behind the default), and the metacell counts then feed program discovery, stability assessment, coordinated multicellular program (CMP) detection, sharing assessment, and cell-level projection.

The fixture plants known truth so every downstream call can be checked: two cell types across 24 donors, each type with a **private program** (constant across donors) and a **coordinated program** whose donor activities are correlated across the two cell types.

<div class="section level2">

## Simulate a single-cell dataset with planted programs<a href="#simulate-a-single-cell-dataset-with-planted-programs" class="anchor" aria-label="anchor"></a>

24 donors x two cell types x 60 cells, 800 genes. Genes 1:40 are a T-cell-private marker block; genes 81:120 a B-cell-private block; genes 41:80 form the **coordinated block**, expressed in both cell types with a donor activity gradient (lambda rising from 0.5 to 8.5 across donors) that is *shared* across types — the planted CMP. Everything else is background Poisson noise.

<div id="cb1" class="sourceCode">

\
<span class="kw"><a href="https://rdrr.io/r/base/library.html" class="external-link"><code class="sourceCode r R">library</code></a></span><span class="op">`(`</span><span class="va">`cellprograms`</span><span class="op">`)`</span>\
<span class="fu"><a href="https://rdrr.io/r/base/Random.html" class="external-link"><code class="sourceCode r R">set.seed</code></a></span><span class="op">`(`</span><span class="fl">`7`</span><span class="op">`)`</span>\
\
<span class="va">`n_donor`</span>` `<span class="op">`<-`</span>` `<span class="fl">`24`</span>`; `<span class="va">`n_cell`</span>` `<span class="op">`<-`</span>` `<span class="fl">`60`</span>`; `<span class="va">`n_gene`</span>` `<span class="op">`<-`</span>` `<span class="fl">`800`</span>\
<span class="va">`types`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link"><code class="sourceCode r R">c</code></a></span><span class="op">`(`</span><span class="st">`"Tcells"`</span>`, `<span class="st">`"Bcells"`</span><span class="op">`)`</span>\
<span class="va">`donors`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/paste.html" class="external-link"><code class="sourceCode r R">paste0</code></a></span><span class="op">`(`</span><span class="st">`"d"`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/seq.html" class="external-link"><code class="sourceCode r R">seq_len</code></a></span><span class="op">`(`</span><span class="va">`n_donor`</span><span class="op">`)`</span><span class="op">`)`</span>\
<span class="va">`a_grad`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/seq.html" class="external-link"><code class="sourceCode r R">seq_len</code></a></span><span class="op">`(`</span><span class="va">`n_donor`</span><span class="op">`)`</span>`                `<span class="co">`# planted donor activity gradient`</span>\
\
<span class="va">`groups`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/unlist.html" class="external-link"><code class="sourceCode r R">unlist</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/lapply.html" class="external-link"><code class="sourceCode r R">lapply</code></a></span><span class="op">`(`</span><span class="va">`donors`</span>`, `<span class="kw">`function`</span><span class="op">`(`</span><span class="va">`d`</span><span class="op">`)`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/paste.html" class="external-link"><code class="sourceCode r R">paste0</code></a></span><span class="op">`(`</span><span class="va">`d`</span>`, `<span class="st">`"_"`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/rep.html" class="external-link"><code class="sourceCode r R">rep</code></a></span><span class="op">`(`</span><span class="va">`types`</span>`, each `<span class="op">`=`</span>` `<span class="va">`n_cell`</span><span class="op">`)`</span><span class="op">`)`</span><span class="op">`)`</span>`,`\
`                 use.names `<span class="op">`=`</span>` `<span class="cn">`FALSE`</span><span class="op">`)`</span>\
<span class="va">`lam`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/matrix.html" class="external-link"><code class="sourceCode r R">matrix</code></a></span><span class="op">`(`</span><span class="fl">`0.5`</span>`, nrow `<span class="op">`=`</span>` `<span class="va">`n_gene`</span>`, ncol `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/length.html" class="external-link"><code class="sourceCode r R">length</code></a></span><span class="op">`(`</span><span class="va">`groups`</span><span class="op">`)`</span><span class="op">`)`</span>\
<span class="kw">`for`</span>` `<span class="op">`(`</span><span class="va">`d`</span>` `<span class="kw">`in`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/seq.html" class="external-link"><code class="sourceCode r R">seq_len</code></a></span><span class="op">`(`</span><span class="va">`n_donor`</span><span class="op">`)`</span><span class="op">`)`</span>` `<span class="op">`{`</span>\
`  `<span class="kw">`for`</span>` `<span class="op">`(`</span><span class="va">`t`</span>` `<span class="kw">`in`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/seq.html" class="external-link"><code class="sourceCode r R">seq_along</code></a></span><span class="op">`(`</span><span class="va">`types`</span><span class="op">`)`</span><span class="op">`)`</span>` `<span class="op">`{`</span>\
`    `<span class="va">`cols`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/which.html" class="external-link"><code class="sourceCode r R">which</code></a></span><span class="op">`(`</span><span class="va">`groups`</span>` `<span class="op">`==`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/paste.html" class="external-link"><code class="sourceCode r R">paste0</code></a></span><span class="op">`(`</span><span class="va">`donors`</span><span class="op">`[`</span><span class="va">`d`</span><span class="op">`]`</span>`, `<span class="st">`"_"`</span>`, `<span class="va">`types`</span><span class="op">`[`</span><span class="va">`t`</span><span class="op">`]`</span><span class="op">`)`</span><span class="op">`)`</span>\
`    `<span class="va">`private`</span>` `<span class="op">`<-`</span>` `<span class="kw">`if`</span>` `<span class="op">`(`</span><span class="va">`types`</span><span class="op">`[`</span><span class="va">`t`</span><span class="op">`]`</span>` `<span class="op">`==`</span>` `<span class="st">`"Tcells"`</span><span class="op">`)`</span>` `<span class="fl">`1`</span><span class="op">`:`</span><span class="fl">`40`</span>` `<span class="kw">`else`</span>` `<span class="fl">`81`</span><span class="op">`:`</span><span class="fl">`120`</span>\
`    `<span class="va">`lam`</span><span class="op">`[`</span><span class="va">`private`</span>`, `<span class="va">`cols`</span><span class="op">`]`</span>` `<span class="op">`<-`</span>` `<span class="fl">`8`</span>`                          `<span class="co">`# private block, constant`</span>\
`    `<span class="va">`lam`</span><span class="op">`[`</span><span class="fl">`41`</span><span class="op">`:`</span><span class="fl">`80`</span>`, `<span class="va">`cols`</span><span class="op">`]`</span>` `<span class="op">`<-`</span>` `<span class="fl">`0.5`</span>` `<span class="op">`+`</span>` `<span class="fl">`8`</span>` `<span class="op">`*`</span>` `<span class="op">`(`</span><span class="va">`d`</span>` `<span class="op">`-`</span>` `<span class="fl">`1`</span><span class="op">`)`</span>` `<span class="op">`/`</span>` `<span class="op">`(`</span><span class="va">`n_donor`</span>` `<span class="op">`-`</span>` `<span class="fl">`1`</span><span class="op">`)`</span>`  `<span class="co">`# coordinated gradient`</span>\
`  `<span class="op">`}`</span>\
<span class="op">`}`</span>\
<span class="va">`counts`</span>` `<span class="op">`<-`</span>` `<span class="fu">`Matrix`</span><span class="fu">`::`</span><span class="fu"><a href="https://rdrr.io/pkg/Matrix/man/Matrix.html" class="external-link"><code class="sourceCode r R">Matrix</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/stats/Poisson.html" class="external-link"><code class="sourceCode r R">rpois</code></a></span><span class="op">`(`</span><span class="va">`n_gene`</span>` `<span class="op">`*`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/length.html" class="external-link"><code class="sourceCode r R">length</code></a></span><span class="op">`(`</span><span class="va">`groups`</span><span class="op">`)`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/vector.html" class="external-link"><code class="sourceCode r R">as.vector</code></a></span><span class="op">`(`</span><span class="va">`lam`</span><span class="op">`)`</span><span class="op">`)`</span>`,`\
`                         nrow `<span class="op">`=`</span>` `<span class="va">`n_gene`</span>`, ncol `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/length.html" class="external-link"><code class="sourceCode r R">length</code></a></span><span class="op">`(`</span><span class="va">`groups`</span><span class="op">`)`</span>`, sparse `<span class="op">`=`</span>` `<span class="cn">`TRUE`</span><span class="op">`)`</span>\
<span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">rownames</code></a></span><span class="op">`(`</span><span class="va">`counts`</span><span class="op">`)`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/paste.html" class="external-link"><code class="sourceCode r R">paste0</code></a></span><span class="op">`(`</span><span class="st">`"gene"`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/seq.html" class="external-link"><code class="sourceCode r R">seq_len</code></a></span><span class="op">`(`</span><span class="va">`n_gene`</span><span class="op">`)`</span><span class="op">`)`</span>\
<span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">colnames</code></a></span><span class="op">`(`</span><span class="va">`counts`</span><span class="op">`)`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/paste.html" class="external-link"><code class="sourceCode r R">paste0</code></a></span><span class="op">`(`</span><span class="st">`"cell"`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/seq.html" class="external-link"><code class="sourceCode r R">seq_len</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/length.html" class="external-link"><code class="sourceCode r R">length</code></a></span><span class="op">`(`</span><span class="va">`groups`</span><span class="op">`)`</span><span class="op">`)`</span><span class="op">`)`</span>\
\
<span class="va">`meta`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/data.frame.html" class="external-link"><code class="sourceCode r R">data.frame</code></a></span><span class="op">`(`</span>\
`  donor `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/grep.html" class="external-link"><code class="sourceCode r R">sub</code></a></span><span class="op">`(`</span><span class="st">`"_.*"`</span>`, `<span class="st">`""`</span>`, `<span class="va">`groups`</span><span class="op">`)`</span>`,`\
`  celltype `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/grep.html" class="external-link"><code class="sourceCode r R">sub</code></a></span><span class="op">`(`</span><span class="st">`".*_"`</span>`, `<span class="st">`""`</span>`, `<span class="va">`groups`</span><span class="op">`)`</span>`,`\
`  row.names `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">colnames</code></a></span><span class="op">`(`</span><span class="va">`counts`</span><span class="op">`)`</span>\
<span class="op">`)`</span>

</div>

</div>

<div class="section level2">

## Build metacells (SuperCell 2.0 + mcRigor)<a href="#build-metacells-supercell-2-0-mcrigor" class="anchor" aria-label="anchor"></a>

Metacells are built per donor x cell type so they never cross either variable. `cells_per_metacell = 15` turns each 60-cell group into ~4 metacells; mcRigor QC runs on every group.

<div id="cb2" class="sourceCode">

\
<span class="va">`obj`</span>` `<span class="op">`<-`</span>` `<span class="fu">`Seurat`</span><span class="fu">`::`</span><span class="fu"><a href="https://satijalab.github.io/seurat-object/reference/CreateSeuratObject.html" class="external-link"><code class="sourceCode r R">CreateSeuratObject</code></a></span><span class="op">`(`</span>`counts `<span class="op">`=`</span>` `<span class="va">`counts`</span>`, meta.data `<span class="op">`=`</span>` `<span class="va">`meta`</span><span class="op">`)`</span>\
<span class="va">`obj`</span>` `<span class="op">`<-`</span>` `<span class="fu">`Seurat`</span><span class="fu">`::`</span><span class="fu"><a href="https://satijalab.org/seurat/reference/NormalizeData.html" class="external-link"><code class="sourceCode r R">NormalizeData</code></a></span><span class="op">`(`</span><span class="va">`obj`</span>`, verbose `<span class="op">`=`</span>` `<span class="cn">`FALSE`</span><span class="op">`)`</span>\
<span class="va">`obj`</span>` `<span class="op">`<-`</span>` `<span class="fu">`Seurat`</span><span class="fu">`::`</span><span class="fu"><a href="https://satijalab.org/seurat/reference/FindVariableFeatures.html" class="external-link"><code class="sourceCode r R">FindVariableFeatures</code></a></span><span class="op">`(`</span><span class="va">`obj`</span>`, nfeatures `<span class="op">`=`</span>` `<span class="fl">`200`</span>`, verbose `<span class="op">`=`</span>` `<span class="cn">`FALSE`</span><span class="op">`)`</span>\
<span class="va">`obj`</span>` `<span class="op">`<-`</span>` `<span class="fu">`Seurat`</span><span class="fu">`::`</span><span class="fu"><a href="https://satijalab.org/seurat/reference/ScaleData.html" class="external-link"><code class="sourceCode r R">ScaleData</code></a></span><span class="op">`(`</span><span class="va">`obj`</span>`, verbose `<span class="op">`=`</span>` `<span class="cn">`FALSE`</span><span class="op">`)`</span>\
<span class="va">`obj`</span>` `<span class="op">`<-`</span>` `<span class="fu">`Seurat`</span><span class="fu">`::`</span><span class="fu"><a href="https://satijalab.org/seurat/reference/RunPCA.html" class="external-link"><code class="sourceCode r R">RunPCA</code></a></span><span class="op">`(`</span><span class="va">`obj`</span>`, npcs `<span class="op">`=`</span>` `<span class="fl">`10`</span>`, verbose `<span class="op">`=`</span>` `<span class="cn">`FALSE`</span><span class="op">`)`</span>\
\
<span class="va">`obj`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`MetacellsByGroupsSuperCell`](../reference/MetacellsByGroupsSuperCell.html)</span><span class="op">`(`</span><span class="va">`obj`</span>`,`\
`                                  group.by `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link"><code class="sourceCode r R">c</code></a></span><span class="op">`(`</span><span class="st">`"donor"`</span>`, `<span class="st">`"celltype"`</span><span class="op">`)`</span>`,`\
`                                  cells_per_metacell `<span class="op">`=`</span>` `<span class="fl">`15`</span>`,`\
`                                  min_cells `<span class="op">`=`</span>` `<span class="fl">`50`</span>`,`\
`                                  n_top_genes `<span class="op">`=`</span>` `<span class="fl">`200`</span>`,`\
`                                  n_pc `<span class="op">`=`</span>` `<span class="fl">`10`</span>`,`\
`                                  mcrigor_args `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/list.html" class="external-link"><code class="sourceCode r R">list</code></a></span><span class="op">`(`</span>`Nrep `<span class="op">`=`</span>` `<span class="fl">`5`</span><span class="op">`)`</span>`,`\
`                                  seed `<span class="op">`=`</span>` `<span class="fl">`7`</span>`, verbose `<span class="op">`=`</span>` `<span class="cn">`FALSE`</span><span class="op">`)`</span>\
<span class="va">`res`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`GetMetacellSuperCell`](../reference/GetMetacellSuperCell.html)</span><span class="op">`(`</span><span class="va">`obj`</span><span class="op">`)`</span>\
<span class="va">`mc`</span>` `<span class="op">`<-`</span>` `<span class="va">`res`</span><span class="op">`$`</span><span class="va">`metacell_object`</span>\
<span class="fu"><a href="https://rdrr.io/r/base/nrow.html" class="external-link"><code class="sourceCode r R">ncol</code></a></span><span class="op">`(`</span><span class="va">`mc`</span><span class="op">`)`</span>`                          `<span class="co">`# total metacells`</span>\
<span class="fu"><a href="https://rdrr.io/r/utils/head.html" class="external-link"><code class="sourceCode r R">head</code></a></span><span class="op">`(`</span><span class="va">`res`</span><span class="op">`$`</span><span class="va">`qc`</span><span class="op">`[`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link"><code class="sourceCode r R">c</code></a></span><span class="op">`(`</span><span class="st">`"metacell"`</span>`, `<span class="st">`"mcDiv"`</span>`, `<span class="st">`"dubious"`</span><span class="op">`)`</span><span class="op">`]`</span><span class="op">`)`</span>

</div>

[`OptimizeGammaSuperCell()`](../reference/OptimizeGammaSuperCell.html) can replace the fixed `cells_per_metacell` when you want the granularity chosen by the mcRigor score (see the metacells article); the downstream workflow is identical either way.

</div>

<div class="section level2">

## Metacell counts to donor pseudobulk<a href="#metacell-counts-to-donor-pseudobulk" class="anchor" aria-label="anchor"></a>

With `mode = "sum"` (the default), summing metacell counts within donor x cell type reproduces the donor pseudobulk exactly — the metacell layer is lossless at the donor level and simply denoises the cells the pseudobulk is built from. We log-normalize the aggregated counts with the standard `log1p(count / library size x 10^4)` recipe and declare the fit’s measurement units as `"logcounts"` (this is what later lets us project the programs back onto single cells).

<div id="cb3" class="sourceCode">

\
<span class="va">`counts_mc`</span>` `<span class="op">`<-`</span>` `<span class="fu">`SeuratObject`</span><span class="fu">`::`</span><span class="fu"><a href="https://satijalab.github.io/seurat-object/reference/Layers.html" class="external-link"><code class="sourceCode r R">LayerData</code></a></span><span class="op">`(`</span><span class="va">`mc`</span><span class="op">`[[`</span><span class="st">`"RNA"`</span><span class="op">`]`</span><span class="op">`]`</span>`, layer `<span class="op">`=`</span>` `<span class="st">`"counts"`</span><span class="op">`)`</span>`  `<span class="co">`# genes x metacells`</span>\
<span class="va">`mc_donor`</span>` `<span class="op">`<-`</span>` `<span class="va">`mc`</span><span class="op">`$`</span><span class="va">`donor`</span>`                       `<span class="co">`# group columns are copied to the`</span>\
<span class="va">`mc_type`</span>` `<span class="op">`<-`</span>` `<span class="va">`mc`</span><span class="op">`$`</span><span class="va">`celltype`</span>`                     `<span class="co">`# metacell object's metadata`</span>\
\
<span class="va">`pb`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/lapply.html" class="external-link"><code class="sourceCode r R">lapply</code></a></span><span class="op">`(`</span><span class="va">`types`</span>`, `<span class="kw">`function`</span><span class="op">`(`</span><span class="va">`t`</span><span class="op">`)`</span>` `<span class="op">`{`</span>\
`  `<span class="va">`mcs`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">colnames</code></a></span><span class="op">`(`</span><span class="va">`mc`</span><span class="op">`)`</span><span class="op">`[`</span><span class="va">`mc_type`</span>` `<span class="op">`==`</span>` `<span class="va">`t`</span><span class="op">`]`</span>\
`  `<span class="va">`agg`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/rowsum.html" class="external-link"><code class="sourceCode r R">rowsum</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/t.html" class="external-link"><code class="sourceCode r R">t</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/matrix.html" class="external-link"><code class="sourceCode r R">as.matrix</code></a></span><span class="op">`(`</span><span class="va">`counts_mc`</span><span class="op">`[`</span>`, `<span class="va">`mcs`</span>`, drop `<span class="op">`=`</span>` `<span class="cn">`FALSE`</span><span class="op">`]`</span><span class="op">`)`</span><span class="op">`)`</span>`,`\
`                group `<span class="op">`=`</span>` `<span class="va">`mc_donor`</span><span class="op">`[`</span><span class="va">`mcs`</span><span class="op">`]`</span><span class="op">`)`</span>`              `<span class="co">`# donors x genes (counts)`</span>\
`  `<span class="va">`agg`</span>` `<span class="op">`<-`</span>` `<span class="va">`agg`</span><span class="op">`[`</span><span class="va">`donors`</span>`, , drop `<span class="op">`=`</span>` `<span class="cn">`FALSE`</span><span class="op">`]`</span>\
`  `<span class="va">`lib`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/colSums.html" class="external-link"><code class="sourceCode r R">rowSums</code></a></span><span class="op">`(`</span><span class="va">`agg`</span><span class="op">`)`</span>\
`  `<span class="va">`out`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/Log.html" class="external-link"><code class="sourceCode r R">log1p</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/sweep.html" class="external-link"><code class="sourceCode r R">sweep</code></a></span><span class="op">`(`</span><span class="va">`agg`</span>`, `<span class="fl">`1`</span>`, `<span class="va">`lib`</span>`, `<span class="st">`"/"`</span><span class="op">`)`</span>` `<span class="op">`*`</span>` `<span class="fl">`1e4`</span><span class="op">`)`</span>`       `<span class="co">`# donors x genes (logcounts)`</span>\
`  `<span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">rownames</code></a></span><span class="op">`(`</span><span class="va">`out`</span><span class="op">`)`</span>` `<span class="op">`<-`</span>` `<span class="va">`donors`</span>`   `<span class="co">`# same donor id across types: how CMP pairing`</span>\
`  `<span class="va">`out`</span>`                       `<span class="co">`# knows which donors are shared`</span>\
<span class="op">`}`</span><span class="op">`)`</span>\
<span class="fu"><a href="https://rdrr.io/r/base/names.html" class="external-link"><code class="sourceCode r R">names</code></a></span><span class="op">`(`</span><span class="va">`pb`</span><span class="op">`)`</span>` `<span class="op">`<-`</span>` `<span class="va">`types`</span>\
<span class="fu"><a href="https://rdrr.io/r/base/lapply.html" class="external-link"><code class="sourceCode r R">lapply</code></a></span><span class="op">`(`</span><span class="va">`pb`</span>`, `<span class="va">`dim`</span><span class="op">`)`</span>\
\
<span class="va">`sample_metadata`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/data.frame.html" class="external-link"><code class="sourceCode r R">data.frame</code></a></span><span class="op">`(`</span>\
`  observation_id `<span class="op">`=`</span>` `<span class="va">`donors`</span>`,`\
`  planted_activity `<span class="op">`=`</span>` `<span class="va">`a_grad`</span>\
<span class="op">`)`</span>\
<span class="va">`x`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`as_cell_program_data`](../reference/as_cell_program_data.html)</span><span class="op">`(`</span><span class="va">`pb`</span>`, sample_metadata `<span class="op">`=`</span>` `<span class="va">`sample_metadata`</span><span class="op">`)`</span>

</div>

</div>

<div class="section level2">

## Fit and canonicalize programs<a href="#fit-and-canonicalize-programs" class="anchor" aria-label="anchor"></a>

One empirical-Bayes matrix factorization per cell type; the coordinated block’s donor-varying activity and the private blocks’ constant activity should both be recovered as programs.

<div id="cb4" class="sourceCode">

\
<span class="va">`fit`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`fit_celltype_programs`](../reference/fit_celltype_programs.html)</span><span class="op">`(`</span><span class="va">`x`</span>`, loading_prior `<span class="op">`=`</span>` `<span class="st">`"point_laplace"`</span>`,`\
`                             max_factors `<span class="op">`=`</span>` `<span class="fl">`3`</span>`, features `<span class="op">`=`</span>` `<span class="st">`"variable"`</span>`,`\
`                             n_variable_genes `<span class="op">`=`</span>` `<span class="fl">`400`</span>`,`\
`                             feature_scale `<span class="op">`=`</span>` `<span class="st">`"logcounts"`</span><span class="op">`)`</span>\
<span class="va">`fit`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`canonicalize_programs`](../reference/canonicalize_programs.html)</span><span class="op">`(`</span><span class="va">`fit`</span><span class="op">`)`</span>\
<span class="fu">[`program_summary`](../reference/program_summary.html)</span><span class="op">`(`</span><span class="va">`fit`</span><span class="op">`)`</span><span class="op">`[`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link"><code class="sourceCode r R">c</code></a></span><span class="op">`(`</span><span class="st">`"cell_type"`</span>`, `<span class="st">`"n_programs"`</span>`, `<span class="st">`"R2"`</span>`, `<span class="st">`"elbo"`</span><span class="op">`)`</span><span class="op">`]`</span>

</div>

<div id="cb5" class="sourceCode">

\
<span class="fu">[`top_genes`](../reference/top_genes.html)</span><span class="op">`(`</span><span class="va">`fit`</span>`, `<span class="st">`"Tcells_1"`</span>`, n `<span class="op">`=`</span>` `<span class="fl">`6`</span><span class="op">`)`</span>\
<span class="fu">[`top_genes`](../reference/top_genes.html)</span><span class="op">`(`</span><span class="va">`fit`</span>`, `<span class="st">`"Bcells_1"`</span>`, n `<span class="op">`=`</span>` `<span class="fl">`6`</span><span class="op">`)`</span>

</div>

The donor x program score matrix carries the planted gradient: the coordinated programs’ scores track the donor index, the private programs are flat.

<div id="cb6" class="sourceCode">

\
<span class="va">`Z`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`program_scores`](../reference/program_scores.html)</span><span class="op">`(`</span><span class="va">`fit`</span><span class="op">`)`</span>\
<span class="fu">[`plot_program_scores`](../reference/plot_program_scores.html)</span><span class="op">`(`</span><span class="va">`fit`</span>`, annotation `<span class="op">`=`</span>` `<span class="st">`"planted_activity"`</span><span class="op">`)`</span>\
<span class="fu">[`plot_variance_explained`](../reference/plot_variance_explained.html)</span><span class="op">`(`</span><span class="va">`fit`</span><span class="op">`)`</span>

</div>

</div>

<div class="section level2">

## Stability<a href="#stability" class="anchor" aria-label="anchor"></a>

[`assess_program_stability()`](../reference/assess_program_stability.html) refits on donor subsamples and reports how often each program is recovered — the same guarantee as in the pseudobulk workflow, now on metacell-derived inputs.

<div id="cb7" class="sourceCode">

\
<span class="va">`fit`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`assess_program_stability`](../reference/assess_program_stability.html)</span><span class="op">`(`</span><span class="va">`fit`</span>`, n_boot `<span class="op">`=`</span>` `<span class="fl">`4`</span>`, resample_frac `<span class="op">`=`</span>` `<span class="fl">`0.7`</span>`,`\
`                                seed `<span class="op">`=`</span>` `<span class="fl">`7`</span>`, progress `<span class="op">`=`</span>` `<span class="cn">`FALSE`</span><span class="op">`)`</span>\
<span class="va">`fit`</span><span class="op">`$`</span><span class="va">`stability`</span><span class="op">`$`</span><span class="va">`per_program`</span><span class="op">`[`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link"><code class="sourceCode r R">c</code></a></span><span class="op">`(`</span><span class="st">`"program_id"`</span>`, `<span class="st">`"recovery_freq"`</span><span class="op">`)`</span><span class="op">`]`</span>

</div>

</div>

<div class="section level2">

## Coordinated multicellular programs<a href="#coordinated-multicellular-programs" class="anchor" aria-label="anchor"></a>

[`find_coordinated_programs()`](../reference/find_coordinated_programs.html) tests whether programs from *different* cell types co-vary across donors. The planted truth: the two programs of the coordinated block (one per type) share donor activity; the private programs do not.

<div id="cb8" class="sourceCode">

\
<span class="va">`cmp`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`find_coordinated_programs`](../reference/find_coordinated_programs.html)</span><span class="op">`(`</span><span class="va">`fit`</span>`, n_perm `<span class="op">`=`</span>` `<span class="fl">`99`</span>`, seed `<span class="op">`=`</span>` `<span class="fl">`7`</span><span class="op">`)`</span>\
<span class="va">`cmp`</span>\
<span class="fu">[`plot_cmp_network`](../reference/plot_cmp_network.html)</span><span class="op">`(`</span><span class="va">`cmp`</span><span class="op">`)`</span>

</div>

The coordinated-block pair is called; the private programs (constant donor activity, independent noise) are not. [`annotate_cmp_tf()`](../reference/annotate_cmp_tf.html) / [`cmp_triage()`](../reference/cmp_triage.html) add TF-footprint annotation and triage on real data (external CollecTRI download, not run here).

</div>

<div class="section level2">

## Shared vs specific programs<a href="#shared-vs-specific-programs" class="anchor" aria-label="anchor"></a>

[`program_sharing()`](../reference/program_sharing.html) asks whether matched programs share their gene *effects* — here the coordinated block is literally the same genes in both types, so a shared call on the coordinated pair is the planted truth.

<div id="cb9" class="sourceCode">

\
<span class="kw">`if`</span>` `<span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/ns-load.html" class="external-link"><code class="sourceCode r R">requireNamespace</code></a></span><span class="op">`(`</span><span class="st">`"mashr"`</span>`, quietly `<span class="op">`=`</span>` `<span class="cn">`TRUE`</span><span class="op">`)`</span><span class="op">`)`</span>` `<span class="op">`{`</span>\
`  `<span class="va">`shr`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`program_sharing`](../reference/program_sharing.html)</span><span class="op">`(`</span><span class="va">`fit`</span>`, n_perm `<span class="op">`=`</span>` `<span class="fl">`99`</span>`, seed `<span class="op">`=`</span>` `<span class="fl">`7`</span><span class="op">`)`</span>\
`  `<span class="va">`shr`</span><span class="op">`$`</span><span class="va">`matches`</span>\
`  `<span class="fu">[`plot_sharing`](../reference/plot_sharing.html)</span><span class="op">`(`</span><span class="va">`shr`</span><span class="op">`)`</span>\
<span class="op">`}`</span>

</div>

</div>

<div class="section level2">

## Project programs back onto single cells<a href="#project-programs-back-onto-single-cells" class="anchor" aria-label="anchor"></a>

Because the fit declared `feature_scale = "logcounts"`, the fitted programs can be projected onto the original single cells ([`project_programs_to_cells()`](../reference/project_programs_to_cells.html) replays the fit’s preprocessing and infers per-cell activities against the frozen loadings). The per-cell activity of a coordinated program should recover the planted donor gradient — this is the bridge from donor-level discovery to cell-level readout.

<div id="cb10" class="sourceCode">

\
<span class="va">`logcounts_sc`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/matrix.html" class="external-link"><code class="sourceCode r R">as.matrix</code></a></span><span class="op">`(`</span><span class="fu">`SeuratObject`</span><span class="fu">`::`</span><span class="fu"><a href="https://satijalab.github.io/seurat-object/reference/Layers.html" class="external-link"><code class="sourceCode r R">LayerData</code></a></span><span class="op">`(`</span><span class="va">`obj`</span><span class="op">`[[`</span><span class="st">`"RNA"`</span><span class="op">`]`</span><span class="op">`]`</span>`, layer `<span class="op">`=`</span>` `<span class="st">`"data"`</span><span class="op">`)`</span><span class="op">`)`</span>\
<span class="va">`logcounts_sc`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/t.html" class="external-link"><code class="sourceCode r R">t</code></a></span><span class="op">`(`</span><span class="va">`logcounts_sc`</span><span class="op">`)`</span>`                   `<span class="co">`# cells x genes (dense)`</span>\
<span class="va">`proj`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`project_programs_to_cells`](../reference/project_programs_to_cells.html)</span><span class="op">`(`</span><span class="va">`fit`</span>`, expression `<span class="op">`=`</span>` `<span class="va">`logcounts_sc`</span>`,`\
`                                  cell_types `<span class="op">`=`</span>` `<span class="fu"><a href="https://rdrr.io/r/stats/setNames.html" class="external-link"><code class="sourceCode r R">setNames</code></a></span><span class="op">`(`</span><span class="va">`meta`</span><span class="op">`$`</span><span class="va">`celltype`</span>`,`\
`                                                        `<span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">rownames</code></a></span><span class="op">`(`</span><span class="va">`meta`</span><span class="op">`)`</span><span class="op">`)`</span><span class="op">`)`</span>\
<span class="va">`act`</span>` `<span class="op">`<-`</span>` `<span class="va">`proj`</span><span class="op">`$`</span><span class="va">`activity`</span>\
<span class="va">`coord_prog`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/grep.html" class="external-link"><code class="sourceCode r R">grep</code></a></span><span class="op">`(`</span><span class="st">`"^Tcells"`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/colnames.html" class="external-link"><code class="sourceCode r R">colnames</code></a></span><span class="op">`(`</span><span class="va">`act`</span><span class="op">`)`</span>`, value `<span class="op">`=`</span>` `<span class="cn">`TRUE`</span><span class="op">`)`</span><span class="op">`[`</span><span class="fl">`1`</span><span class="op">`]`</span>\
<span class="co">`# Activity is type-specific: cells outside the program's source type get NA,`</span>\
<span class="co">`# so donor means are computed within T cells only.`</span>\
<span class="va">`is_t`</span>` `<span class="op">`<-`</span>` `<span class="va">`meta`</span><span class="op">`$`</span><span class="va">`celltype`</span>` `<span class="op">`==`</span>` `<span class="st">`"Tcells"`</span>\
<span class="va">`donor_mean`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/tapply.html" class="external-link"><code class="sourceCode r R">tapply</code></a></span><span class="op">`(`</span><span class="va">`act`</span><span class="op">`[`</span><span class="va">`is_t`</span>`, `<span class="va">`coord_prog`</span><span class="op">`]`</span>`, `<span class="va">`meta`</span><span class="op">`$`</span><span class="va">`donor`</span><span class="op">`[`</span><span class="va">`is_t`</span><span class="op">`]`</span>`, `<span class="va">`mean`</span><span class="op">`)`</span>\
<span class="fu"><a href="https://rdrr.io/r/base/Round.html" class="external-link"><code class="sourceCode r R">round</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/stats/cor.html" class="external-link"><code class="sourceCode r R">cor</code></a></span><span class="op">`(`</span><span class="va">`donor_mean`</span><span class="op">`[`</span><span class="va">`donors`</span><span class="op">`]`</span>`, `<span class="va">`a_grad`</span>`, method `<span class="op">`=`</span>` `<span class="st">`"spearman"`</span><span class="op">`)`</span>`, `<span class="fl">`3`</span><span class="op">`)`</span>

</div>

The high rank correlation between the projected per-cell activity (averaged within donor) and the planted gradient confirms the metacell -\> donor fit -\> cell projection round trip preserves the biology. Absolute activity levels live on the fit’s scale; the donor ordering is the readout.

</div>

<div class="section level2">

## Session info<a href="#session-info" class="anchor" aria-label="anchor"></a>

</div>
