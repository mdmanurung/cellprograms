# Access metacell results stored by MetacellsByGroupsSuperCell

Access metacell results stored by MetacellsByGroupsSuperCell

## Usage

``` r
GetMetacellSuperCell(seurat_obj, wgcna_name = "SuperCell")
```

## Arguments

- seurat_obj:

  A Seurat object processed by \[MetacellsByGroupsSuperCell()\].

- wgcna_name:

  Name under which results were stored. Default \`"SuperCell"\`.

## Value

The stored list: \`metacell_object\`, \`soft_assignments\` (NULL for
this backend), \`qc\`, \`fit_stats\`, \`zerorate_overall\`,
\`dubrate_overall\`, \`score_overall\`, \`params\`, \`skipped_groups\`,
\`call\`.
