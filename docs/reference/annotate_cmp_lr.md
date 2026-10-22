<div class="page-header">

# Ligand-receptor compatibility of coordinated programs

<div class="d-none name">

`annotate_cmp_lr.Rd`

</div>

</div>

<div class="ref-description section level2">

For member programs a and b of a CMP in different cell types and a ligand-receptor pair (L, R), scores \`w_a\[L\] \* w_b\[R\]\`, where \`w\` is the absolute gene loading times (1 - lfsr), scaled to a maximum of 1 within each program. Significance uses a degree-preserving resource permutation (ligand labels shuffled among ligands, receptor labels among receptors) with the maximum over the CMP as the null statistic. Optionally, the donor-level co-variation of the ligand in the sender and the receptor in the receiver is reported from the pseudobulk (a MOFAtalk-style check).

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

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

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

fit<a href="#arg-fit" class="anchor" aria-label="anchor"></a><br>
A canonicalized \`cell_program_fit\`.

cmp<a href="#arg-cmp" class="anchor" aria-label="anchor"></a><br>
A \`cmp_result\` or its \`cmps\` table.

lr<a href="#arg-lr" class="anchor" aria-label="anchor"></a><br>
Data frame with \`ligand\` and \`receptor\` columns (gene IDs matching the fit), e.g. a LIANA consensus resource.

n_perm<a href="#arg-n-perm" class="anchor" aria-label="anchor"></a><br>
Permutations (default 499).

alpha<a href="#arg-alpha" class="anchor" aria-label="anchor"></a><br>
Level for reported pairs (default 0.05).

top<a href="#arg-top" class="anchor" aria-label="anchor"></a><br>
Maximum pairs reported per CMP (default 10).

control<a href="#arg-control" class="anchor" aria-label="anchor"></a><br>
\`"fwer"\` (default): each pair is compared with the permutation maximum over all pairs of the CMP (family-wise). \`"fdr"\`: each pair is compared with the pooled permutation null of all pair scores, then Benjamini-Hochberg over all tested (pair, sender, receiver) combinations; \`p\` is then the BH-adjusted value.

donor_check<a href="#arg-donor-check" class="anchor" aria-label="anchor"></a><br>
Add the donor-level ligand/receptor correlation (default TRUE).

seed<a href="#arg-seed" class="anchor" aria-label="anchor"></a><br>
Random seed.

</div>

<div class="section level2">

## Value<a href="#value" class="anchor" aria-label="anchor"></a>

Data frame: \`cmp\`, \`sender\`, \`receiver\`, \`ligand\`, \`receptor\`, \`score\`, \`p\` and, with \`donor_check\`, \`donor_cor\`, \`donor_p\` (Spearman test of ligand in sender vs receptor in receiver across shared donors) and \`n_donors\`.

</div>

<div class="section level2">

## Details<a href="#details" class="anchor" aria-label="anchor"></a>

A hit means the coordinated programs are \*compatible\* with that interaction: it is not evidence that the interaction causes the coordination. Programs rich in secreted or surface genes (e.g. inflammatory programs) yield hits more readily.

</div>

<div class="section level2">

## Examples<a href="#ref-examples" class="anchor" aria-label="anchor"></a>

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
```

</div>
