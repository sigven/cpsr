# Function that assigns ACMG evidence indicators for PS1 ('ACMG_PS1') - Same amino acid change as a previously established pathogenic variant regardless of nucleotide change

Function that assigns ACMG evidence indicators for PS1 ('ACMG_PS1') -
Same amino acid change as a previously established pathogenic variant
regardless of nucleotide change

## Usage

``` r
assign_PS1_evidence(
  var_calls = NULL,
  pathogenic_codons = NULL,
  pathogenic_nucleotides = NULL
)
```

## Arguments

- var_calls:

  variants in cancer predisposition genes

- pathogenic_codons:

  data.frame with known pathogenic codons pathogenic variants have been
  observed (from ClinVar)

- pathogenic_nucleotides:

  data.frame with pathogenic nucleotides (non-coding, e.g. splice)

## Value

data.frame with ACMG evidence indicators
