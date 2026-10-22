<div class="page-header">

# Export pseudobulk profiles to MOFAcellulaR input format

<div class="d-none name">

`as_mofacellular_input.Rd`

</div>

</div>

<div class="ref-description section level2">

Builds the \`counts\` matrix (genes in rows, pseudobulk profiles in columns, column names \`\<cell_type\>\_\<observation_id\>\`) and the accompanying \`coldata\` data.frame (\`donor_id\`, \`cell_type\`, \`cell_counts\`) expected by \`MOFAcellulaR::create_init_exp()\`. Works on both \`cell_program_data\` and \`cell_program_fit\` objects (the pseudobulk matrices are carried through fitting), guaranteeing that local EBMF and joint MOFA2 analyses see identical inputs.

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
as_mofacellular_input(x, cell_counts = NULL)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

x<a href="#arg-x" class="anchor" aria-label="anchor"></a><br>
A \`cell_program_data\` or \`cell_program_fit\` object.

cell_counts<a href="#arg-cell-counts" class="anchor" aria-label="anchor"></a><br>
Optional data.frame with columns \`observation_id\`, \`cell_type\`, and \`cell_counts\` (number of cells behind each pseudobulk profile). If NULL (default), \`cell_counts\` is filled with NA; supply it if you intend to use MOFAcellulaR's cell-count-based QC (\`filt_profiles(ncells = ...)\`).

</div>

<div class="section level2">

## Value<a href="#value" class="anchor" aria-label="anchor"></a>

A list with elements \`counts\` (genes x profiles matrix) and \`coldata\` (data.frame with rownames matching \`colnames(counts)\`), ready for \`MOFAcellulaR::create_init_exp(counts = ..., coldata = ...)\`.

</div>

<div class="section level2">

## Details<a href="#details" class="anchor" aria-label="anchor"></a>

MOFAcellulaR's TMM normalization (\`tmm_trns\`) expects \*raw-count\* pseudobulk; pass raw-count pseudobulk in \`x\` if you plan to use it. Log-normalized pseudobulk is fine for the factor analysis itself but not for TMM.

</div>

<div class="section level2">

## Examples<a href="#ref-examples" class="anchor" aria-label="anchor"></a>

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
mofa_in <- as_mofacellular_input(x)
dim(mofa_in$counts)
head(mofa_in$coldata)
```

</div>
