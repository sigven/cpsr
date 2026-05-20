# Prepare unified variant dataset for single-table display

Combines ClinVar and CPSR-classified variants across ALL significance
classes into a single data frame. Includes a variant_key for crosstalk
linking.

## Usage

``` r
prepare_unified_all_variants(cps_report = NULL, max_rows = 1000)
```

## Arguments

- cps_report:

  CPSR report object

- max_rows:

  Integer. Max rows per source

## Value

Data frame with all variants from both sources, keyed for crosstalk
