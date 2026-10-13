# Format an OmniPath CollecTRI table as a signed TF-target network

Follows the decoupler convention: one row per TF-target pair, mode of
regulation (\`mor\`) -1 when the consensus is inhibition only and +1
otherwise (stimulation, both, or unknown); NF-kB and AP-1 complex
members are merged into \`NFKB\` and \`AP1\`.

## Usage

``` r
prepare_collectri(raw)
```

## Arguments

- raw:

  Data frame from the OmniPath \`interactions\` endpoint with
  \`datasets=collectri&genesymbols=yes\` (columns \`source\`,
  \`source_genesymbol\`, \`target_genesymbol\`,
  \`consensus_stimulation\`, \`consensus_inhibition\`).

## Value

Data frame with \`source\`, \`target\`, \`mor\`.

## Examples

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
#>   source target mor
#> 1  STAT1   IRF1   1
#> 2    AP1   MMP1   1
#> 3    MYC CDKN1A  -1
```
