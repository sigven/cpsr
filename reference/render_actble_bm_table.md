# Build biomarker reactable with category-aware styling

Build biomarker reactable with category-aware styling

## Usage

``` r
render_actble_bm_table(
  rctbl_recs = NULL,
  variant_category = "snv_indel",
  color_palette = NULL
)
```

## Arguments

- rctbl_recs:

  List with \$main and \$nested data frames containing prepared
  biomarker table records, for example as returned by
  [`prep_biomarker_tbl()`](https://sigven.github.io/cpsr/reference/prep_biomarker_tbl.md).

- variant_category:

  One of "snv_indel", "cnv", "fusion"

- color_palette:

  color palette for therapeutic biomarkers

## Value

A reactable object with the biomarker table
