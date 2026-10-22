<div class="page-header">

# Gene-overlap class of coordinated programs

<div class="d-none name">

`cmp_gene_overlap.Rd`

</div>

</div>

<div class="ref-description section level2">

For each pair of member programs, tests whether their top-\`top_n\` \|loading\| genes overlap more than expected within the shared gene universe (hypergeometric, BH across pairs). A CMP is \`"shared-genes"\` if at least half its pairs overlap, \`"partial-overlap"\` if some do and \`"distinct-genes"\` otherwise. Distinct-gene CMPs are the case that gene-matching methods cannot find.

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
cmp_gene_overlap(fit, cmp, top_n = 100L, q = 0.05)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

fit<a href="#arg-fit" class="anchor" aria-label="anchor"></a><br>
A canonicalized \`cell_program_fit\`.

cmp<a href="#arg-cmp" class="anchor" aria-label="anchor"></a><br>
A \`cmp_result\` or its \`cmps\` table.

top_n<a href="#arg-top-n" class="anchor" aria-label="anchor"></a><br>
Top genes per program (default 100).

q<a href="#arg-q" class="anchor" aria-label="anchor"></a><br>
BH threshold (default 0.05).

</div>

<div class="section level2">

## Value<a href="#value" class="anchor" aria-label="anchor"></a>

Data frame: \`cmp\`, \`n_pairs\`, \`frac_overlapping\`, \`class\`.

</div>
