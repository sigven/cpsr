# Variant classification as colored pills in a PDF report table - pathogenicity colors as in the HTML report

Variant classification as colored pills in a PDF report table -
pathogenicity colors as in the HTML report

## Usage

``` r
classification_pill_pdf(classifications)
```

## Arguments

- classifications:

  classification per table row (e.g. 'Pathogenic', 'Likely Pathogenic',
  'VUS')

## Value

character vector with Typst markup (see pill_pdf); unknown categories in
grey
