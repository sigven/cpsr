# Create the complete unified variant reactable

Simple, fast reactable showing primary columns only. All other columns
accessible via nested details.

## Usage

``` r
create_unified_variant_reactable(
  shared_data = NULL,
  primary_cols = c("SYMBOL", "gnomADg_AF", "CONSEQUENCE", "ALTERATION", "CLASSIFICATION",
    "GENOTYPE"),
  color_palette = NULL,
  flag_low_depth = FALSE,
  dp_threshold = 10
)
```

## Arguments

- shared_data:

  Crosstalk SharedData object with unified variants

- primary_cols:

  Character vector of primary column names to display

- color_palette:

  CPSR color_palette object

- flag_low_depth:

  Logical. If TRUE, flag variants with a control-sample sequencing depth
  (DP_CONTROL) below `dp_threshold` with a marker on the genotype cell.
  Should only be set to TRUE when `DP_CONTROL` carries real depth values
  for at least some variants (i.e. not the "-1" sentinel used when CPSR
  is run without a matched control sample).

- dp_threshold:

  Integer. Depth threshold below which a variant is considered
  low-depth. Default: 10.

## Value

reactable widget
