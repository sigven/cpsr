# Variant table (classified variants or secondary findings) for the PDF report, with column names as in the HTML report

Variant table (classified variants or secondary findings) for the PDF
report, with column names as in the HTML report

## Usage

``` r
variant_table_pdf(
  variants,
  classification_col = "CLASSIFICATION",
  source_col = "ASSERTION_AUTHORITY"
)
```

## Arguments

- variants:

  data frame with variants (variant_display)

- classification_col:

  column with the variant classification

- source_col:

  column with the assertion authority (ClinVar/CPSR), or NULL

## Value

tinytable object
