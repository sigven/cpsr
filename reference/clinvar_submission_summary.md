# Summarise ClinVar germline submissions per classification

Converts ClinVar contributing classifications (e.g.
'Likely_pathogenic:3,Uncertain_significance:1') into a short label
('LP×3, VUS×1').

## Usage

``` r
clinvar_submission_summary(x)
```

## Arguments

- x:

  character vector with ClinVar contributing classifications

## Value

character vector (NA where no submissions are available)
