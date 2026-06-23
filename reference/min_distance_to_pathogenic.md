# Function that calculates the closest distance to a known pathogenic variant in the same gene

Function that calculates the closest distance to a known pathogenic
variant in the same gene

## Usage

``` r
min_distance_to_pathogenic(var_calls = NULL, ref_data = NULL)
```

## Arguments

- var_calls:

  variants in cancer predisposition genes

- ref_data:

  reference data object with ClinVar annotations

## Value

data.frame with closest distance to known pathogenic variant in same
gene
