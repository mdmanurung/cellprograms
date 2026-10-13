# Plot the coordinated-program network

Network of significant cross-cell-type program edges. Node size = node
statistic A_v, node color = CMP membership, edge width = \|T\|, edge
color = edge sign. Programs flagged as ambient get hollow nodes.

## Usage

``` r
plot_cmp_network(cmp, max_edges = NULL, seed = 1)
```

## Arguments

- cmp:

  Output of \`find_coordinated_programs\`.

- max_edges:

  Draw only the strongest edges (default: all).

- seed:

  Layout seed (the spring layout is stochastic).

## Value

A ggplot object.

## Examples

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 3, seed = 1))
cmp <- find_coordinated_programs(fit, n_perm = 200, seed = 1)
plot_cmp_network(cmp)
```
