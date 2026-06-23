# Cell renderer factory for biomarker clinical significance Returns a function that renders biomarker clinical significance values as colored pills based on the provided color palette.

Cell renderer factory for biomarker clinical significance Returns a
function that renders biomarker clinical significance values as colored
pills based on the provided color palette.

## Usage

``` r
rt_cell_bm_significance(color_palette)
```

## Arguments

- color_palette:

  CPSR color palette object containing biomarker_types styling
  information

## Value

A function that takes a biomarker clinical significance value and
returns an HTML div element containing styled pills for each
significance category
