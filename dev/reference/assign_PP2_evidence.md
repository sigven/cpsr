# Function that assigns ACMG evidence indicators for PP2

Function that assigns ACMG evidence indicators for PP2

## Usage

``` r
assign_PP2_evidence(
  var_calls = NULL,
  max_benign_missense_frac = 0.1,
  max_truncating_frac = 0.5
)
```

## Arguments

- var_calls:

  variants in cancer predisposition genes

- max_benign_missense_frac:

  maximum tolerated fraction of benign missense variants for gene

- max_truncating_frac:

  maximum tolerated fraction of pathogenic truncating variants for gene

## Value

data.frame with ACMG evidence indicators
