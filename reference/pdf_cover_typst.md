# Cover page of the PDF report: sample metadata (cpsr-cover(), defined in typst/\_preamble.typ) and key findings (report_highlights()), followed by the table of contents

Cover page of the PDF report: sample metadata (cpsr-cover(), defined in
typst/\_preamble.typ) and key findings (report_highlights()), followed
by the table of contents

## Usage

``` r
pdf_cover_typst(report)
```

## Arguments

- report:

  CPSR report object

## Value

raw Typst blocks around markdown (character), to be printed with cat()
in a chunk with output 'asis'
