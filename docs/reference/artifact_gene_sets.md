# Artifact-gene sets

Gene classes that can produce cross-donor co-variation without
coordinated biology, used for flagging called programs and for fit-time
exclusion:

## Usage

``` r
artifact_gene_sets
```

## Format

A named list of character vectors (\`vdj\` is a length-1 regular
expression; \`genotype_marker\` and \`sex\` are gene IDs).

## Details

\* \`vdj\`: immunoglobulin and T-cell-receptor V(D)J and constant genes.
Stored as a single regular expression,
\`^(IG\[HKL\]V\|IG\[KL\]C\|IGH\[GAMDE\]\|JCHAIN\|TR\[ABGD\]V)\`, not an
enumerated gene list, because the full V-gene complement is
genome-annotation-dependent; expand it against a gene universe with
\`.expand_vdj()\` (or use \`artifact_genes()\`, which expands every
set). \* \`genotype_marker\`: genes with common germline deletion,
copy-number or HLA-haplotype expression variation (GSTM1/GSTT1
deletions, CCL3L1/CCL4L1 copy number, HLA-DRB3/4/5 and HLA-DQA/DQB
haplotypes, ...). \* \`sex\`: XIST and the male-specific Y-chromosome
genes RPS4Y1, DDX3Y, UTY, KDM5D, EIF1AY, NLGN4Y, TMSB4Y, USP9Y, ZFY.

Provenance: the \`genotype_marker\` list is a screening heuristic, not
validated against genotype data; extend or replace it for your cohort.
The sets overlap \`cmp_nuisance_genes()\` (used by \`cmp_triage()\`);
unlike that regex vector, \`artifact_gene_sets\` stores explicit gene
IDs where the set is closed, so it can be subtracted from a feature
panel directly.

## Examples

``` r
names(artifact_gene_sets)
#> [1] "vdj"             "genotype_marker" "sex"            
artifact_gene_sets$vdj   # V(D)J genes are stored as a regular expression
#> [1] "^(IG[HKL]V|IG[KL]C|IGH[GAMDE]|JCHAIN|TR[ABGD]V)"
artifact_gene_sets$sex
#>  [1] "XIST"   "RPS4Y1" "DDX3Y"  "UTY"    "KDM5D"  "EIF1AY" "NLGN4Y" "TMSB4Y"
#>  [9] "USP9Y"  "ZFY"   
```
