<div class="page-header">

# Default nuisance-gene patterns for CMP triage

<div class="d-none name">

`cmp_nuisance_genes.Rd`

</div>

</div>

<div class="ref-description section level2">

Regular expressions for gene classes that produce cross-cell-type donor co-variation without coordinated biology: immunoglobulin/TCR V(D)J and constant genes (\`vdj\`), genes with common germline deletion, copy-number or HLA-haplotype expression variation (\`genotype\`), sex-chromosome genes (\`sex\`), immediate-early/dissociation stress genes (\`ieg\`) and platelet or erythrocyte contamination markers (\`platelet_rbc\`; these cells are usually not fitted, so the ambient score cannot see them). The \`genotype\` list is a screening heuristic, not validated against genotype data; extend or replace it for your cohort.

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
cmp_nuisance_genes()
```

</div>

<div class="section level2">

## Value<a href="#value" class="anchor" aria-label="anchor"></a>

Named character vector of regular expressions.

</div>
