# Create a ClinVar variant reactable with expandable detail rows

Displays ClinVar-classified variants with primary columns in the main
row and all remaining columns accessible via an expandable details row.
Used for secondary findings and pharmacogenomic tables.

## Usage

``` r
create_clinvar_reactable(
  data = NULL,
  primary_cols = c("SYMBOL", "ALTERATION", "CLINVAR_CLASSIFICATION", "CLINVAR_PHENOTYPE",
    "GENOTYPE", "CONSEQUENCE"),
  color_palette = NULL
)
```

## Arguments

- data:

  Data frame of variants to display

- primary_cols:

  Character vector of primary column names to display

- color_palette:

  CPSR color_palette object

## Value

reactable widget
