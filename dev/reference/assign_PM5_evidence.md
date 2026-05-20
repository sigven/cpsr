# Function that assigns ACMG evidence indicators for PM5

Function that assigns ACMG evidence indicators for PM5

## Usage

``` r
assign_PM5_evidence(var_calls = NULL, pathogenic_codons = NULL)
```

## Arguments

- var_calls:

  variants in cancer predisposition genes

- pathogenic_codons:

  data.frame with codons where pathogenic variants have been observed
  (from ClinVar)

## Value

data.frame with ACMG evidence indicators
