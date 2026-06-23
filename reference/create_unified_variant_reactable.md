# Create the complete unified variant reactable

Simple, fast reactable showing primary columns only. All other columns
accessible via nested details.

## Usage

``` r
create_unified_variant_reactable(
  shared_data = NULL,
  primary_cols = c("SYMBOL", "gnomADg_AF", "CONSEQUENCE", "ALTERATION", "CLASSIFICATION",
    "GENOTYPE"),
  color_palette = NULL
)
```

## Arguments

- shared_data:

  Crosstalk SharedData object with unified variants

- primary_cols:

  Character vector of primary column names to display

- color_palette:

  CPSR color_palette object

## Value

reactable widget
