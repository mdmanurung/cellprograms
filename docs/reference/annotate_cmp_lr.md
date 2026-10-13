# Ligand-receptor compatibility of coordinated programs

For member programs a and b of a CMP in different cell types and a
ligand-receptor pair (L, R), scores \`w_a\[L\] \* w_b\[R\]\`, where
\`w\` is the absolute gene loading times (1 - lfsr), scaled to a maximum
of 1 within each program. Significance uses a degree-preserving resource
permutation (ligand labels shuffled among ligands, receptor labels among
receptors) with the maximum over the CMP as the null statistic.
Optionally, the donor-level co-variation of the ligand in the sender and
the receptor in the receiver is reported from the pseudobulk (a
MOFAtalk-style check).

## Usage

``` r
annotate_cmp_lr(
  fit,
  cmp,
  lr,
  n_perm = 499L,
  alpha = 0.05,
  top = 10L,
  control = c("fwer", "fdr"),
  donor_check = TRUE,
  seed = 1L
)
```

## Arguments

- fit:

  A canonicalized \`cell_program_fit\`.

- cmp:

  A \`cmp_result\` or its \`cmps\` table.

- lr:

  Data frame with \`ligand\` and \`receptor\` columns (gene IDs matching
  the fit), e.g. a LIANA consensus resource.

- n_perm:

  Permutations (default 499).

- alpha:

  Level for reported pairs (default 0.05).

- top:

  Maximum pairs reported per CMP (default 10).

- control:

  \`"fwer"\` (default): each pair is compared with the permutation
  maximum over all pairs of the CMP (family-wise). \`"fdr"\`: each pair
  is compared with the pooled permutation null of all pair scores, then
  Benjamini-Hochberg over all tested (pair, sender, receiver)
  combinations; \`p\` is then the BH-adjusted value.

- donor_check:

  Add the donor-level ligand/receptor correlation (default TRUE).

- seed:

  Random seed.

## Value

Data frame: \`cmp\`, \`sender\`, \`receiver\`, \`ligand\`, \`receptor\`,
\`score\`, \`p\` and, with \`donor_check\`, \`donor_cor\`, \`donor_p\`
(Spearman test of ligand in sender vs receptor in receiver across shared
donors) and \`n_donors\`.

## Details

A hit means the coordinated programs are \*compatible\* with that
interaction: it is not evidence that the interaction causes the
coordination. Programs rich in secreted or surface genes (e.g.
inflammatory programs) yield hits more readily.

## Examples

``` r
# Toy canonicalized fit: program 1 of cell types A and B shares donor
# activity z (a planted coordinated program) and top genes G1-G15.
set.seed(1)
ids <- paste0("d", seq_len(40))
genes <- paste0("G", seq_len(600))
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
# plant a ligand at the top of A_1 and its receptor at the top of B_1
fit$loadings$A["G150", "A_1"] <- 5
fit$loadings$B["G151", "B_1"] <- 5
set.seed(2)
lr <- data.frame(ligand = c("G150", sample(paste0("G", 60:349), 300, replace = TRUE)),
                 receptor = c("G151", sample(paste0("G", 352:600), 300, replace = TRUE)))
annotate_cmp_lr(fit, res, lr, n_perm = 99)
#>   cmp sender receiver ligand receptor score    p
#> 1   1  A|A_1    B|B_1   G150     G151     1 0.02
```
