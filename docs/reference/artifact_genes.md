<div class="page-header">

# Expand artifact-gene sets against a gene universe

<div class="d-none name">

`artifact_genes.Rd`

</div>

</div>

<div class="ref-description section level2">

Returns the union of the requested \`artifact_gene_sets\` expanded against \`genes\` (the V(D)J regex is matched; enumerated sets are intersected). Use it to build the \`exclude_genes\` argument of \`fit_celltype_programs()\`.

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
artifact_genes(genes, sets = names(artifact_gene_sets))
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

genes<a href="#arg-genes" class="anchor" aria-label="anchor"></a><br>
Character vector of gene IDs (the universe, e.g. the column names of a pseudobulk matrix).

sets<a href="#arg-sets" class="anchor" aria-label="anchor"></a><br>
Subset of \`names(artifact_gene_sets)\` (default all).

</div>

<div class="section level2">

## Value<a href="#value" class="anchor" aria-label="anchor"></a>

Character vector of gene IDs in \`genes\` belonging to any requested set, in universe order.

</div>

<div class="section level2">

## Examples<a href="#ref-examples" class="anchor" aria-label="anchor"></a>

``` r
universe <- c("IGHV1-2", "IGKC", "GSTM1", "XIST", "MALAT1")
artifact_genes(universe)
artifact_genes(universe, sets = "vdj")
```

</div>
