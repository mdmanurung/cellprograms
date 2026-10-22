<div class="page-header">

# Themes of coordinated programs

<div class="d-none name">

`annotate_cmp_themes.Rd`

</div>

</div>

<div class="ref-description section level2">

Univariate linear model (decoupler-style ULM) of each member program's absolute loadings on gene-set membership, for every set with at least \`min_size\` genes in the program's gene universe; BH across sets within a program. Absolute loadings are used because program signs are arbitrary up to canonicalization and a set can be enriched at either pole.

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
annotate_cmp_themes(fit, cmp, sets, min_size = 10L, top = 3L, q = 0.05)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

fit<a href="#arg-fit" class="anchor" aria-label="anchor"></a><br>
A canonicalized \`cell_program_fit\`.

cmp<a href="#arg-cmp" class="anchor" aria-label="anchor"></a><br>
A \`cmp_result\` or its \`cmps\` table.

sets<a href="#arg-sets" class="anchor" aria-label="anchor"></a><br>
Named list of gene vectors (e.g. MSigDB Hallmark).

min_size<a href="#arg-min-size" class="anchor" aria-label="anchor"></a><br>
Minimum set size within the universe (default 10).

top<a href="#arg-top" class="anchor" aria-label="anchor"></a><br>
Top sets per program (default 3).

q<a href="#arg-q" class="anchor" aria-label="anchor"></a><br>
BH threshold (default 0.05).

</div>

<div class="section level2">

## Value<a href="#value" class="anchor" aria-label="anchor"></a>

Data frame: \`cmp\`, \`cell_type\`, \`program_id\`, \`set\`, \`t\`, \`padj\`.

</div>

<div class="section level2">

## Examples<a href="#ref-examples" class="anchor" aria-label="anchor"></a>

``` r
# Toy canonicalized fit: program 1 of cell types A and B shares donor
# activity z (a planted coordinated program) and top genes G1-G15.
set.seed(1)
ids <- paste0("d", seq_len(40))
genes <- paste0("G", seq_len(100))
z <- rnorm(40)
mk <- function(ct) {
  S <- matrix(rnorm(40 * 2), 40, 2, dimnames = list(ids, paste0(ct, "_", 1:2)))
  S[, 1] <- z + rnorm(40, sd = 0.3)
  W <- matrix(rnorm(length(genes) * 2, sd = 0.1), length(genes), 2,
              dimnames = list(genes, colnames(S)))
  W[1:15, 1] <- W[1:15, 1] + 1
  list(S = S, W = W)
}
ta <- mk("A"); tb <- mk("B")
fit <- structure(list(
  scores = list(A = ta$S, B = tb$S), loadings = list(A = ta$W, B = tb$W),
  cell_types = c("A", "B"), observation_ids = ids,
  sample_metadata = data.frame(observation_id = ids,
                               group = rep(c("ctrl", "stim"), each = 20))),
  class = c("cell_program_fit", "cell_program_data", "list"))
res <- find_coordinated_programs(fit, n_perm = 99, ambient = "none")
sets <- list(planted = paste0("G", 1:15), other = paste0("G", 50:65))
annotate_cmp_themes(fit, res, sets, min_size = 5)
```

</div>
