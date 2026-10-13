# Ambient-likeness of programs

For each program, the largest absolute correlation between its gene
loadings and another cell type's expression specificity
(\`profile_other - profile_own\`, log scale). High values indicate a
program that tracks contamination by another cell type's transcripts
(ambient RNA, doublets) rather than biology intrinsic to the cell type.
Genuine biology can also score high (e.g. a monocyte-like program in
DCs), so treat it as a flag.

## Usage

``` r
ambient_program_scores(fit, profiles = NULL, min_genes = 50L)
```

## Arguments

- fit:

  A canonicalized \`cell_program_fit\`.

- profiles:

  Optional named list of per-cell-type mean log-expression vectors named
  by gene; defaults to column means of \`fit\$pseudobulk\`.

- min_genes:

  Minimum shared genes to compare two cell types (default 50).

## Value

Data frame with \`cell_type\`, \`program_id\`, \`ambient_score\` and
\`source\` (the cell type giving the maximum).

## Examples

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- fit_celltype_programs(x, max_factors = 3, features = "variable", n_variable_genes = 150)
fit <- canonicalize_programs(fit)
ambient_program_scores(fit)
#>   cell_type program_id ambient_score source
#> 1         B        B_1    0.49352750     NK
#> 2         B        B_2    0.16956795     NK
#> 3       CD4      CD4_1    0.55106699    CD8
#> 4       CD8      CD8_1    0.09585351    CD4
#> 5        NK       NK_1    0.40634664      B
#> 6        NK       NK_2    0.10544901    CD8
#> 7      Mono     Mono_1    0.65147659    CD8
#> 8      Mono     Mono_2    0.37685525     NK
```
