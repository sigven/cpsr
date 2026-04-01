# Function that appends mechanism-of-inheritance (MOI) and mechanism of disease (MOD) annotations to cancer predisposition genes (cpg's), as well as estimated fractions of truncation vs. non-truncating variants etc per predisposition gene (from ClinVar)

Function that appends mechanism-of-inheritance (MOI) and mechanism of
disease (MOD) annotations to cancer predisposition genes (cpg's), as
well as estimated fractions of truncation vs. non-truncating variants
etc per predisposition gene (from ClinVar)

## Usage

``` r
append_cpg_properties(cpg_calls, ref_data = NULL)
```

## Arguments

- cpg_calls:

  data frame with variant calls in predisposition genes

- ref_data:

  PCGR/CPSR reference data object

## Value

cpg_calls data frame with cancer predisposition gene properties appended
(mechanism of disease, inheritance patterns etc)
