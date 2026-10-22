<div class="page-header">

# Plot the coordinated-program network

<div class="d-none name">

`plot_cmp_network.Rd`

</div>

</div>

<div class="ref-description section level2">

Network of significant cross-cell-type program edges. Node size = node statistic A_v, node color = CMP membership, edge width = \|T\|, edge color = edge sign. Programs flagged as ambient get hollow nodes.

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
plot_cmp_network(cmp, max_edges = NULL, seed = 1)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

cmp<a href="#arg-cmp" class="anchor" aria-label="anchor"></a><br>
Output of \`find_coordinated_programs\`.

max_edges<a href="#arg-max-edges" class="anchor" aria-label="anchor"></a><br>
Draw only the strongest edges (default: all).

seed<a href="#arg-seed" class="anchor" aria-label="anchor"></a><br>
Layout seed (the spring layout is stochastic).

</div>

<div class="section level2">

## Value<a href="#value" class="anchor" aria-label="anchor"></a>

A ggplot object.

</div>

<div class="section level2">

## Examples<a href="#ref-examples" class="anchor" aria-label="anchor"></a>

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 3, seed = 1))
cmp <- find_coordinated_programs(fit, n_perm = 200, seed = 1)
plot_cmp_network(cmp)
```

</div>
