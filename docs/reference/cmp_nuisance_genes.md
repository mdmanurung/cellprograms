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
#>                                                                                     vdj 
#>                                       "^(IG[HKL]V|IG[KL]C|IGH[GAMDE]|JCHAIN|TR[ABGD]V)" 
#>                                                                                genotype 
#>   "^(GSTM1|GSTT1|GSTT2B?|UTS2|ERAP2|CCL3L1|CCL3L3|CCL4L2|HLA-DRB[345]|HLA-DQ[AB][12])$" 
#>                                                                                     sex 
#>                   "^(XIST|TSIX|RPS4Y1|DDX3Y|UTY|KDM5D|EIF1AY|USP9Y|ZFY|TXLNGY|NLGN4Y)$" 
#>                                                                                     ieg 
#>         "^(FOS|FOSB|JUN|JUNB|EGR1|EGR2|EGR3|ATF3|IER2|IER3|NR4A1|DUSP1|PPP1R15A|KLF6)$" 
#>                                                                            platelet_rbc 
#> "^(PF4|PPBP|GP9|GP1BB|TUBB1|ITGA2B|GNG11|CAVIN2|TREML1|PTCRA|HBB|HBA1|HBA2|HBD|ALAS2)$" 
grep(cmp_nuisance_genes()["vdj"], c("IGHV1-2", "IGKC", "GAPDH"), value = TRUE)
#> [1] "IGHV1-2" "IGKC"   
```
