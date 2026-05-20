# Function that extracts known benign peptide changes from ClinVar (\>= 2 gold stars), as well as known pathogenic nucleotide and codon sites

Function that extracts known benign peptide changes from ClinVar (\>= 2
gold stars), as well as known pathogenic nucleotide and codon sites

## Usage

``` r
known_path_benign_sites(ref_data = NULL, min_gold_stars = 2)
```

## Arguments

- ref_data:

  reference data object with ClinVar annotations

- min_gold_stars:

  minimum number of gold stars for considering a variant as known
  benign/pathogenic

## Value

list with data.frames for benign peptide changes, pathogenic nucleotide
sites, and pathogenic codon sites
