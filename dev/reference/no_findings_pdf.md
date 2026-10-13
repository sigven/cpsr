# Message box for a PDF report section without findings

Message box for a PDF report section without findings

## Usage

``` r
no_findings_pdf(text, title = "No findings")
```

## Arguments

- text:

  message (plain text; Typst markup characters are escaped)

- title:

  short title shown in bold before the message

## Value

raw Typst block (character), to be printed with cat() in a chunk with
output 'asis'
