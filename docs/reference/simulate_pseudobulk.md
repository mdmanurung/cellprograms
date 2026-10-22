<div class="page-header">

# Simulate cell-type pseudobulk matrices with planted programs

<div class="d-none name">

`simulate_pseudobulk.Rd`

</div>

</div>

<div class="ref-description section level2">

Simulate cell-type pseudobulk matrices with planted programs

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
simulate_pseudobulk(
  scenario = "mixed",
  n_samples = 150,
  n_genes = 1500,
  program_size = 50,
  signal_sd = 1.5,
  noise_sd = 1,
  cell_types = c("B", "CD4", "CD8", "NK", "Mono"),
  seed = 1
)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

scenario<a href="#arg-scenario" class="anchor" aria-label="anchor"></a><br>
One of "private" (Scenario 1), "mixed" (Scenario 4, principal), "same_genes_independent" (Scenario 7 negative control), or "global" (Scenario 2).

n_samples<a href="#arg-n-samples" class="anchor" aria-label="anchor"></a><br>
Number of biological observations.

n_genes<a href="#arg-n-genes" class="anchor" aria-label="anchor"></a><br>
Total gene pool size.

program_size<a href="#arg-program-size" class="anchor" aria-label="anchor"></a><br>
Number of active genes per program.

signal_sd<a href="#arg-signal-sd" class="anchor" aria-label="anchor"></a><br>
Program activity standard deviation.

noise_sd<a href="#arg-noise-sd" class="anchor" aria-label="anchor"></a><br>
Residual noise standard deviation.

cell_types<a href="#arg-cell-types" class="anchor" aria-label="anchor"></a><br>
Character vector of cell-type names.

seed<a href="#arg-seed" class="anchor" aria-label="anchor"></a><br>
Random seed.

</div>

<div class="section level2">

## Value<a href="#value" class="anchor" aria-label="anchor"></a>

A list with \`pseudobulk\` (named list of observations x genes matrices), \`truth\` (per-cell-type list of true score vectors, loading indices, and program labels), and \`scenario\`.

</div>

<div class="section level2">

## Examples<a href="#ref-examples" class="anchor" aria-label="anchor"></a>

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
names(sim$pseudobulk)
dim(sim$pseudobulk$B)
```

</div>
