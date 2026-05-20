# Cell renderer factory for classification column Returns a function that renders classification values as colored pills based on the provided color palette.

Cell renderer factory for classification column Returns a function that
renders classification values as colored pills based on the provided
color palette.

## Usage

``` r
rt_cell_classification(color_palette)
```

## Arguments

- color_palette:

  CPSR color palette object containing pathogenicity styling information

## Value

A function that takes a classification value and returns an HTML span
element with appropriate styling
