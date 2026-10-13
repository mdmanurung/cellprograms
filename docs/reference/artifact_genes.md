# Expand artifact-gene sets against a gene universe

Returns the union of the requested \`artifact_gene_sets\` expanded
against \`genes\` (the V(D)J regex is matched; enumerated sets are
intersected). Use it to build the \`exclude_genes\` argument of
\`fit_celltype_programs()\`.

## Usage

``` r
artifact_genes(genes, sets = names(artifact_gene_sets))
```

## Arguments

- genes:

  Character vector of gene IDs (the universe, e.g. the column names of a
  pseudobulk matrix).

- sets:

  Subset of \`names(artifact_gene_sets)\` (default all).

## Value

Character vector of gene IDs in \`genes\` belonging to any requested
set, in universe order.

## Examples

``` r
universe <- c("IGHV1-2", "IGKC", "GSTM1", "XIST", "MALAT1")
artifact_genes(universe)
#> [1] "IGHV1-2" "IGKC"    "GSTM1"   "XIST"   
artifact_genes(universe, sets = "vdj")
#> [1] "IGHV1-2" "IGKC"   
```
