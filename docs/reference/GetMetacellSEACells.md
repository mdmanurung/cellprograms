# Access metacell results stored by MetacellsByGroupsSEACells

Access metacell results stored by MetacellsByGroupsSEACells

## Usage

``` r
GetMetacellSEACells(seurat_obj, wgcna_name = "SEACells")
```

## Arguments

- seurat_obj:

  A Seurat object processed by \[MetacellsByGroupsSEACells()\].

- wgcna_name:

  Name under which results were stored. Default \`"SEACells"\`.

## Value

The stored list: \`metacell_object\`, \`soft_assignments\`, \`qc\`,
\`fit_stats\`, \`zerorate_overall\`, \`dubrate_overall\`,
\`score_overall\`, \`params\`, \`skipped_groups\`, \`call\`.
