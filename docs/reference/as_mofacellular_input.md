# Export pseudobulk profiles to MOFAcellulaR input format

Builds the \`counts\` matrix (genes in rows, pseudobulk profiles in
columns, column names \`\<cell_type\>\_\<observation_id\>\`) and the
accompanying \`coldata\` data.frame (\`donor_id\`, \`cell_type\`,
\`cell_counts\`) expected by \`MOFAcellulaR::create_init_exp()\`. Works
on both \`cell_program_data\` and \`cell_program_fit\` objects (the
pseudobulk matrices are carried through fitting), guaranteeing that
local EBMF and joint MOFA2 analyses see identical inputs.

## Usage

``` r
as_mofacellular_input(x, cell_counts = NULL)
```

## Arguments

- x:

  A \`cell_program_data\` or \`cell_program_fit\` object.

- cell_counts:

  Optional data.frame with columns \`observation_id\`, \`cell_type\`,
  and \`cell_counts\` (number of cells behind each pseudobulk profile).
  If NULL (default), \`cell_counts\` is filled with NA; supply it if you
  intend to use MOFAcellulaR's cell-count-based QC
  (\`filt_profiles(ncells = ...)\`).

## Value

A list with elements \`counts\` (genes x profiles matrix) and
\`coldata\` (data.frame with rownames matching \`colnames(counts)\`),
ready for \`MOFAcellulaR::create_init_exp(counts = ..., coldata =
...)\`.

## Details

MOFAcellulaR's TMM normalization (\`tmm_trns\`) expects \*raw-count\*
pseudobulk; pass raw-count pseudobulk in \`x\` if you plan to use it.
Log-normalized pseudobulk is fine for the factor analysis itself but not
for TMM.

## Examples

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
mofa_in <- as_mofacellular_input(x)
dim(mofa_in$counts)
#> [1] 300 120
head(mofa_in$coldata)
#>        donor_id cell_type cell_counts
#> B_obs1     obs1         B          NA
#> B_obs2     obs2         B          NA
#> B_obs3     obs3         B          NA
#> B_obs4     obs4         B          NA
#> B_obs5     obs5         B          NA
#> B_obs6     obs6         B          NA
```
