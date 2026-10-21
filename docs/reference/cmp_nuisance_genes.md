# Default nuisance-gene patterns for CMP triage

Regular expressions for gene classes that produce cross-cell-type donor
co-variation without coordinated biology: immunoglobulin/TCR V(D)J and
constant genes (\`vdj\`), genes with common germline deletion,
copy-number or HLA-haplotype expression variation (\`genotype\`),
sex-chromosome genes (\`sex\`), immediate-early/dissociation stress
genes (\`ieg\`) and platelet or erythrocyte contamination markers
(\`platelet_rbc\`; these cells are usually not fitted, so the ambient
score cannot see them). The \`genotype\` list is a screening heuristic,
not validated against genotype data; extend or replace it for your
cohort.

## Usage

``` r
cmp_nuisance_genes()
```

## Value

Named character vector of regular expressions.

## Examples

``` r
cmp_nuisance_genes()
#> Error in cmp_nuisance_genes(): could not find function "cmp_nuisance_genes"
grep(cmp_nuisance_genes()["vdj"], c("IGHV1-2", "IGKC", "GAPDH"), value = TRUE)
#> Error in cmp_nuisance_genes(): could not find function "cmp_nuisance_genes"
```
