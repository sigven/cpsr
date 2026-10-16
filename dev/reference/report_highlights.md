# Key findings of a CPSR report (highlights) as markdown

The most critical findings come first: pathogenic and likely pathogenic
variants in the cancer predisposition genes of the virtual panel. These
are followed by supporting findings, each included only when assessed:
an overview of variant classifications (including VUS), secondary
findings (ACMG SF), genomic biomarkers (evidence levels A/B), and
pharmacogenomic phenotypes (CPIC).

## Usage

``` r
report_highlights(
  report,
  max_show = 5,
  pill = FALSE,
  vus_coding_only = FALSE,
  links = FALSE
)
```

## Arguments

- report:

  CPSR report object

- max_show:

  maximum number of items listed per finding

- pill:

  label 'Pathogenic' / 'Likely Pathogenic' as colored pills (white bold
  text on pathogenicity color, as in the Clinical significance column):
  TRUE gives raw Typst (PDF report), "html" gives an HTML span (HTML
  report), FALSE (default) gives bold text

- vus_coding_only:

  count protein-coding VUS only (PDF report - more condensed; non-coding
  VUS are noted separately)

- links:

  finding titles as links to the corresponding report sections (PDF
  report - section labels as in the quarto/pdf templates, e.g.
  'cpsr-plp-variants'); only for sections that are included in the
  report

## Value

list with 'lead' (markdown on pathogenic/likely pathogenic variants) and
'supporting' (markdown bullet list, "" if none)
