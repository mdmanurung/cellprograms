<div class="page-header">

# Access metacell results stored by MetacellsByGroupsSEACells

<div class="d-none name">

`GetMetacellSEACells.Rd`

</div>

</div>

<div class="ref-description section level2">

Access metacell results stored by MetacellsByGroupsSEACells

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
GetMetacellSEACells(seurat_obj, wgcna_name = "SEACells")
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

seurat_obj<a href="#arg-seurat-obj" class="anchor" aria-label="anchor"></a><br>
A Seurat object processed by \[MetacellsByGroupsSEACells()\].

wgcna_name<a href="#arg-wgcna-name" class="anchor" aria-label="anchor"></a><br>
Name under which results were stored. Default \`"SEACells"\`.

</div>

<div class="section level2">

## Value<a href="#value" class="anchor" aria-label="anchor"></a>

The stored list: \`metacell_object\`, \`soft_assignments\`, \`qc\`, \`fit_stats\`, \`zerorate_overall\`, \`dubrate_overall\`, \`score_overall\`, \`params\`, \`skipped_groups\`, \`call\`.

</div>
