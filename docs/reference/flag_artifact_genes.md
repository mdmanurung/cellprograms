# Flag programs or CMPs carrying artifact genes

For each CMP (or each program when \`cmp\` is \`NULL\`), takes the
top-\`top_n\` genes by absolute loading of each member program (pooled
and de-duplicated across members for a CMP) and reports, for each
artifact set, the fraction of those genes matching the set and a flag
that is TRUE when at least one gene matches. This mirrors the screening
heuristic used in the project's CMP report; a flag means "inspect before
interpreting", not "discard".

## Usage

``` r
flag_artifact_genes(
  fit,
  cmp = NULL,
  sets = c("vdj", "genotype_marker"),
  top_n = 12L
)
```

## Arguments

- fit:

  A canonicalized \`cell_program_fit\`.

- cmp:

  A \`cmp_result\` or its \`cmps\` table; \`NULL\` (default) flags every
  program of every cell type instead.

- sets:

  Subset of \`names(artifact_gene_sets)\` to screen against (default
  \`c("vdj", "genotype_marker")\`).

- top_n:

  Number of top-\|loading\| genes per program (default 12).

## Value

Data frame. With \`cmp\`: one row per CMP (\`cmp\`, \`n_members\`). With
\`cmp = NULL\`: one row per program (\`cell_type\`, \`program_id\`).
Both have, per set \`\<s\>\`, columns \`\<s\>\_frac\` (fraction of the
top genes matching) and \`\<s\>\_flag\` (TRUE when \>= 1 gene matches).

## Examples

``` r
# Toy fit whose program A_1 is carried by V(D)J / genotype-marker genes
set.seed(1)
g <- c(paste0("G", 1:30), "IGHV1-2", "IGKC", "GSTM1")
W <- matrix(rnorm(length(g) * 2, sd = 0.1), length(g), 2,
            dimnames = list(g, paste0("A_", 1:2)))
W[c("IGHV1-2", "IGKC", "GSTM1"), 1] <- c(5, 4, 3)
fit <- structure(list(loadings = list(A = W), cell_types = "A"),
                 class = c("cell_program_fit", "cell_program_data", "list"))
flag_artifact_genes(fit)
#>   cell_type program_id   vdj_frac vdj_flag genotype_marker_frac
#> 1         A        A_1 0.16666667     TRUE           0.08333333
#> 2         A        A_2 0.08333333     TRUE           0.00000000
#>   genotype_marker_flag
#> 1                 TRUE
#> 2                FALSE
```
