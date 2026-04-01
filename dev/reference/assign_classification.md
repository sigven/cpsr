# Function that assigns final pathogenicity classification (B, LB, VUS, P, LP) based on accumulated scores from different ACMG criteria and pre-defined cutoffs (calibrated against ClinVar)

Function that assigns final pathogenicity classification (B, LB, VUS, P,
LP) based on accumulated scores from different ACMG criteria and
pre-defined cutoffs (calibrated against ClinVar)

## Usage

``` r
assign_classification(var_calls)
```

## Arguments

- var_calls:

  data frame with variant calls in predisposition genes

## Value

var_calls data frame with pathogenicity classification appended
