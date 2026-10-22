<div class="page-header">

# Historical case study: DIVAS COVID-19 multi-omics

<div class="d-none name">

`divas_case_study.Rmd`

</div>

</div>

**Historical, not current acceptance evidence.** The cached analyses below have not been regenerated under the corrected software and source contracts. Prediction and the patient-pairing inferential claims are withdrawn. The source is processed mean expression with unresolved upstream transformation and selection limitations, not verified raw UMI. Reproduction requires approved source/units/rights and independent T1 refits for every retained method.

This article historically reanalysed the DIVAS COVID-19 case study (<a href="https://github.com/ByronSyun/DIVAS_COVID19_CaseStudy" class="external-link uri">https://github.com/ByronSyun/DIVAS_COVID19_CaseStudy</a>) with `cellprograms` and compares the two pipelines’ outputs program by program. DIVAS jointly decomposes 6 data blocks (4 scRNA-seq cell types + proteomics + metabolomics) into components tagged by which blocks they span; the cellprograms workflow fits each block separately, then recovers cross-block structure with coordination (CMPs) and cross-cell-type sharing. All numbers and figures are precomputed; the code chunks show exactly how each result was produced. The full scripts live in `benchmarks/divas_comparison/scripts/` in the package repository.

<div class="section level2">

## Data<a href="#data" class="anchor" aria-label="anchor"></a>

114 samples (85 patients: 29 with paired T1/T2 samples, 56 with a single timepoint), 40 mild / 40 moderate / 34 severe COVID-19, outcome = WHO ordinal score (`wos`). Six blocks: CD4 T, CD8 T, CD14 monocytes, and NK gene blocks (8,634 shared genes each), proteomics (481 proteins), and metabolomics (763 metabolites).

DIVAS was run by the original authors (`DIVASmain(nsim = 500, colCent = FALSE, rowCent = TRUE, seed = 123)`); its precomputed results (Zenodo record 17430294) yield 91 components: 44 individual (1 block), 46 partial (2-5 blocks), 1 joint (all 6 blocks).

</div>

<div class="section level2">

## Fitting at DIVAS parity<a href="#fitting-at-divas-parity" class="anchor" aria-label="anchor"></a>

<div id="cb1" class="sourceCode">

\
<span class="co">`# each block CSV: features x samples -> transpose to samples x features`</span>\
<span class="va">`read_block`</span>` `<span class="op">`<-`</span>` `<span class="kw">`function`</span><span class="op">`(`</span><span class="va">`path`</span><span class="op">`)`</span>` `<span class="op">`{`</span>\
`  `<span class="fu"><a href="https://rdrr.io/r/base/t.html" class="external-link"><code class="sourceCode r R">t</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/matrix.html" class="external-link"><code class="sourceCode r R">as.matrix</code></a></span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/utils/read.table.html" class="external-link"><code class="sourceCode r R">read.csv</code></a></span><span class="op">`(`</span><span class="va">`path`</span>`, row.names `<span class="op">`=`</span>` `<span class="fl">`1`</span>`, check.names `<span class="op">`=`</span>` `<span class="cn">`FALSE`</span><span class="op">`)`</span><span class="op">`)`</span><span class="op">`)`</span>\
<span class="op">`}`</span>\
<span class="va">`blocks`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link"><code class="sourceCode r R">c</code></a></span><span class="op">`(`</span><span class="st">`"CD4"`</span>`, `<span class="st">`"CD8"`</span>`, `<span class="st">`"CD14"`</span>`, `<span class="st">`"NK"`</span>`, `<span class="st">`"proteomics"`</span>`, `<span class="st">`"metabolomics"`</span><span class="op">`)`</span>\
<span class="va">`pseudobulk`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/lapply.html" class="external-link"><code class="sourceCode r R">lapply</code></a></span><span class="op">`(`</span><span class="va">`blocks`</span>`, `<span class="kw">`function`</span><span class="op">`(`</span><span class="va">`b`</span><span class="op">`)`</span>\
`  `<span class="fu">`read_block`</span><span class="op">`(`</span><span class="fu"><a href="https://rdrr.io/r/base/file.path.html" class="external-link"><code class="sourceCode r R">file.path</code></a></span><span class="op">`(`</span><span class="va">`data_dir`</span>`, `<span class="st">`"input"`</span>`, `<span class="fu"><a href="https://rdrr.io/r/base/paste.html" class="external-link"><code class="sourceCode r R">paste0</code></a></span><span class="op">`(`</span><span class="st">`"Combined_"`</span>`, `<span class="va">`b`</span>`, `<span class="st">`".csv"`</span><span class="op">`)`</span><span class="op">`)`</span><span class="op">`)`</span><span class="op">`)`</span>\
<span class="fu"><a href="https://rdrr.io/r/base/names.html" class="external-link"><code class="sourceCode r R">names</code></a></span><span class="op">`(`</span><span class="va">`pseudobulk`</span><span class="op">`)`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link"><code class="sourceCode r R">c</code></a></span><span class="op">`(`</span><span class="st">`"CD4_T"`</span>`, `<span class="st">`"CD8_T"`</span>`, `<span class="st">`"CD14_Monocyte"`</span>`, `<span class="st">`"NK"`</span>`,`\
`                       `<span class="st">`"proteomics"`</span>`, `<span class="st">`"metabolomics"`</span><span class="op">`)`</span>\
<span class="va">`x`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`as_cell_program_data`](../reference/as_cell_program_data.html)</span><span class="op">`(`</span><span class="va">`pseudobulk`</span>`, sample_metadata `<span class="op">`=`</span>` `<span class="va">`sample_metadata`</span><span class="op">`)`</span>\
<span class="va">`fit`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`fit_celltype_programs`](../reference/fit_celltype_programs.html)</span><span class="op">`(`</span><span class="va">`x`</span>`, loading_prior `<span class="op">`=`</span>` `<span class="st">`"point_laplace"`</span>`,`\
`                             max_factors `<span class="op">`=`</span>` `<span class="fl">`30`</span>`, center `<span class="op">`=`</span>` `<span class="cn">`TRUE`</span>`, scale `<span class="op">`=`</span>` `<span class="cn">`FALSE`</span>`,`\
`                             var_type `<span class="op">`=`</span>` `<span class="fl">`1`</span>`, backfit `<span class="op">`=`</span>` `<span class="cn">`TRUE`</span>`, nullcheck `<span class="op">`=`</span>` `<span class="cn">`TRUE`</span>`,`\
`                             seed `<span class="op">`=`</span>` `<span class="fl">`1`</span><span class="op">`)`</span>\
<span class="va">`fit`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`canonicalize_programs`](../reference/canonicalize_programs.html)</span><span class="op">`(`</span><span class="va">`fit`</span><span class="op">`)`</span>

</div>

All 6 blocks retain K = 30 programs (180 scores total).

</div>

<div class="section level2">

## Coordination and sharing<a href="#coordination-and-sharing" class="anchor" aria-label="anchor"></a>

<div id="cb2" class="sourceCode">

\
<span class="va">`cmp`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`find_coordinated_programs`](../reference/find_coordinated_programs.html)</span><span class="op">`(`</span><span class="va">`fit`</span>`, n_perm `<span class="op">`=`</span>` `<span class="fl">`999`</span>`, seed `<span class="op">`=`</span>` `<span class="fl">`1`</span>`, min_donors `<span class="op">`=`</span>` `<span class="fl">`20`</span><span class="op">`)`</span>\
<span class="va">`cmp_within`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`find_coordinated_programs`](../reference/find_coordinated_programs.html)</span><span class="op">`(`</span><span class="va">`fit`</span>`, mode `<span class="op">`=`</span>` `<span class="st">`"within"`</span>`,`\
`                                        context `<span class="op">`=`</span>` `<span class="st">`"severity"`</span>`,`\
`                                        n_perm `<span class="op">`=`</span>` `<span class="fl">`999`</span>`, seed `<span class="op">`=`</span>` `<span class="fl">`1`</span><span class="op">`)`</span>\
<span class="va">`gene_fit`</span>` `<span class="op">`<-`</span>` `<span class="va">`fit`</span>`; `<span class="va">`gene_fit`</span><span class="op">`$`</span><span class="va">`cell_types`</span>` `<span class="op">`<-`</span>` `<span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link"><code class="sourceCode r R">c</code></a></span><span class="op">`(`</span><span class="st">`"CD4_T"`</span>`, `<span class="st">`"CD8_T"`</span>`, `<span class="st">`"CD14_Monocyte"`</span>`, `<span class="st">`"NK"`</span><span class="op">`)`</span>\
<span class="va">`sharing`</span>` `<span class="op">`<-`</span>` `<span class="fu">[`program_sharing`](../reference/program_sharing.html)</span><span class="op">`(`</span><span class="va">`gene_fit`</span>`, n_perm `<span class="op">`=`</span>` `<span class="fl">`199`</span>`, seed `<span class="op">`=`</span>` `<span class="fl">`1`</span>`, z_thresh `<span class="op">`=`</span>` `<span class="fl">`1.5`</span><span class="op">`)`</span>

</div>

- Marginal coordination: 148 significant cross-block edges collapse into 6 coordinated program modules (CMPs); within-severity: 123 edges, 6 CMPs - the structure is not a severity-composition artifact.
- Program sharing: 120 gene-block programs in 31 groups (27 groups span all 4 cell types); all 120 classify as shared (`frac_top_shared >= 0.98`). This is more permissive than DIVAS’s taxonomy (44/46/1 individual/partial/joint): the 4 immune cell types are strongly correlated in this cohort, and the dominant severity axis is present in every cell type, so mash’s lfsr-based classification attributes essentially all structure to shared programs.

</div>

<div class="section level2">

## Concordance with DIVAS<a href="#concordance-with-divas" class="anchor" aria-label="anchor"></a>

<div id="cb3" class="sourceCode">

\
<span class="co">`# score concordance: |Spearman| between each DIVAS component and each`</span>\
<span class="co">`# cellprograms score across the 114 samples`</span>

</div>

23 of 91 DIVAS components have a cellprograms counterpart at \|Spearman\| \> 0.5 (51/180 in the reverse direction); 42 component-program pairs are mutual best matches (19 above 0.5). All 120 gene-block programs are classified **shared** across cell types, and their best-matching DIVAS components span all taxonomy classes (96 partial, 18 individual, 6 joint). The single DIVAS joint component has a strong best match (CD4_5, \|rho\| = 0.78) but is not mutual-best. Score-level agreement does not extend to loadings: median top-100 loading Jaccard for the 40 mutual-best gene-block pairs is 0.16 (IQR 0.14-0.22).

Concordance rises with DIVAS pattern breadth - the median best-match correlation increases from 0.38 (components private to 1 block) to 0.80 (components spanning 5 blocks):

<div id="cb4" class="sourceCode">

\
<span class="fu">`knitr`</span><span class="fu">`::`</span><span class="fu"><a href="https://rdrr.io/pkg/knitr/man/include_graphics.html" class="external-link"><code class="sourceCode r R">include_graphics</code></a></span><span class="op">`(`</span><span class="st">`"figures/fig_breadth_concordance.png"`</span><span class="op">`)`</span>

</div>

<div class="figure">

<img src="figures/fig_breadth_concordance.png" class="r-plt" width="1050" alt="Best-match concordance by DIVAS pattern breadth." />

Best-match concordance by DIVAS pattern breadth.

</div>

The full concordance landscape (DIVAS components x cellprograms scores, row-split by DIVAS taxonomy class):

<div id="cb5" class="sourceCode">

\
<span class="fu">`knitr`</span><span class="fu">`::`</span><span class="fu"><a href="https://rdrr.io/pkg/knitr/man/include_graphics.html" class="external-link"><code class="sourceCode r R">include_graphics</code></a></span><span class="op">`(`</span><span class="st">`"figures/fig_concordance_heatmap.png"`</span><span class="op">`)`</span>

</div>

<div class="figure">

<img src="figures/fig_concordance_heatmap.png" class="r-plt" width="1650" alt="DIVAS component patterns vs cellprograms score concordance." />

DIVAS component patterns vs cellprograms score concordance.

</div>

</div>

<div class="section level2">

## Historical clinical association — prediction withdrawn<a href="#historical-clinical-association-prediction-withdrawn" class="anchor" aria-label="anchor"></a>

The prediction table below is retained solely as an invalid historical result: unsupervised representations and upstream selection were not learned entirely inside patient-held-out folds. These values do not establish predictive performance or superiority. The current benchmark disables that prediction path; no replacement clinical prediction claim is made.

At BH \< 0.05, 44/180 cellprograms scores associate with `wos` versus 8/91 DIVAS components. Cross-validated ridge prediction of `wos`:

| Predictors             | n   | R2   | RMSE |
|------------------------|-----|------|------|
| DIVAS components       | 91  | 0.42 | 1.66 |
| cellprograms scores    | 180 | 0.64 | 1.28 |
| severity-only baseline | 2   | 0.90 | 0.62 |

The historical severity association is not a leakage-free prediction result; the earlier claim that one representation captures more predictive variance is withdrawn. Historical effect-size agreement across the 91 matched pairs is weak overall (Spearman rho = 0.17), but pairs where either method is significant sit mostly above the diagonal - the cellprograms counterpart carries as much or more `wos` signal:

<div id="cb6" class="sourceCode">

\
<span class="fu">`knitr`</span><span class="fu">`::`</span><span class="fu"><a href="https://rdrr.io/pkg/knitr/man/include_graphics.html" class="external-link"><code class="sourceCode r R">include_graphics</code></a></span><span class="op">`(`</span><span class="st">`"figures/fig_effect_agreement.png"`</span><span class="op">`)`</span>

</div>

<div class="figure">

<img src="figures/fig_effect_agreement.png" class="r-plt" width="975" alt="Clinical effect-size agreement for matched pairs." />

Clinical effect-size agreement for matched pairs.

</div>

</div>

<div class="section level2">

## Historical pairing check — inferential conclusion withdrawn<a href="#historical-pairing-check-inferential-conclusion-withdrawn" class="anchor" aria-label="anchor"></a>

The old workflow reported the following check because 29 patients contribute two timepoints. These are historical numbers, not accepted robustness evidence: it checked the headline coordination calls with (i) a T1-only refit (57 samples, no repeated measures): 50/148 full-cohort edges reproduce and both runs find 6 CMPs; and (ii) a patient-blocked permutation null (all samples of a patient permuted together, 999 permutations): 121/148 edges (82%) remain significant at p \< 0.05. **The earlier robustness conclusion is withdrawn.** This was not a valid general patient-cluster permutation design, and a target T1-only refit is not an independent all-method comparison. Even corrected T1 refits cannot undo upstream outcome-balanced/all-visit filtering.

</div>

<div class="section level2">

## Reproducing<a href="#reproducing" class="anchor" aria-label="anchor"></a>

The historical script order was `01_download_data.sh`, `02_fit_cellprograms.R`, `03_cmp_sharing.R`, `04_extract_divas.R`, `05_compare.R`, `06_figures_extra.R`, `07_sensitivity_pairing.R` (all under `benchmarks/divas_comparison/scripts/`). The maintained workflow and its explicit pending/blocking source decisions live in cellprograms-reproducibility. The old script list is not a claim of current end-to-end reproduction.

</div>
