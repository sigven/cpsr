# Create crosstalk-linked unified filters for all variants

Shared filters across ClinVar and CPSR variants that control the single
table.

## Usage

``` r
create_unified_variant_filters(
  shared_data = NULL,
  cps_report = NULL,
  coding_status = "coding"
)
```

## Arguments

- shared_data:

  Crosstalk SharedData object

- cps_report:

  CPSR report object

- coding_status:

  Character. "coding" or "noncoding". Determines which filters to
  include.

## Value

bscols HTML widget with filters
