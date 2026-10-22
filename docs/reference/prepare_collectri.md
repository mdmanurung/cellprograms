<div class="page-header">

# Format an OmniPath CollecTRI table as a signed TF-target network

<div class="d-none name">

`prepare_collectri.Rd`

</div>

</div>

<div class="ref-description section level2">

Follows the decoupler convention: one row per TF-target pair, mode of regulation (\`mor\`) -1 when the consensus is inhibition only and +1 otherwise (stimulation, both, or unknown); NF-kB and AP-1 complex members are merged into \`NFKB\` and \`AP1\`.

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
prepare_collectri(raw)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

raw<a href="#arg-raw" class="anchor" aria-label="anchor"></a><br>
Data frame from the OmniPath \`interactions\` endpoint with \`datasets=collectri&genesymbols=yes\` (columns \`source\`, \`source_genesymbol\`, \`target_genesymbol\`, \`consensus_stimulation\`, \`consensus_inhibition\`).

</div>

<div class="section level2">

## Value<a href="#value" class="anchor" aria-label="anchor"></a>

Data frame with \`source\`, \`target\`, \`mor\`.

</div>

<div class="section level2">

## Examples<a href="#ref-examples" class="anchor" aria-label="anchor"></a>

``` r
# `raw` mimics the OmniPath interactions endpoint (datasets=collectri)
raw <- data.frame(
  source = c("P1", "COMPLEX:P2_P3", "P4"),
  source_genesymbol = c("STAT1", "JUN_FOS", "MYC"),
  target_genesymbol = c("IRF1", "MMP1", "CDKN1A"),
  consensus_stimulation = c(TRUE, TRUE, FALSE),
  consensus_inhibition = c(FALSE, FALSE, TRUE)
)
prepare_collectri(raw)
```

</div>
