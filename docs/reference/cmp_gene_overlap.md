# Gene-overlap class of coordinated programs

For each pair of member programs, tests whether their top-\`top_n\`
\|loading\| genes overlap more than expected within the shared gene
universe (hypergeometric, BH across pairs). A CMP is \`"shared-genes"\`
if at least half its pairs overlap, \`"partial-overlap"\` if some do and
\`"distinct-genes"\` otherwise. Distinct-gene CMPs are the case that
gene-matching methods cannot find.

## Usage

``` r
cmp_gene_overlap(fit, cmp, top_n = 100L, q = 0.05)
```

## Arguments

- fit:

  A canonicalized \`cell_program_fit\`.

- cmp:

  A \`cmp_result\` or its \`cmps\` table.

- top_n:

  Top genes per program (default 100).

- q:

  BH threshold (default 0.05).

## Value

Data frame: \`cmp\`, \`n_pairs\`, \`frac_overlapping\`, \`class\`.
