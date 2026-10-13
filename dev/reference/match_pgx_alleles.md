# Match sample variants against CPIC allele definitions

Sample variants (CHROM/POS/REF/ALT) are matched against CPIC allele
definitions (pgx_allele, reference data bundle). For each matched
variant, all allele definitions that include the variant are returned,
together with the number of allele-defining variants that are assessable
from VCF (N_ASSESSABLE) and present in the sample (N_PRESENT).

## Usage

``` r
match_pgx_alleles(calls, pgx_allele)
```

## Arguments

- calls:

  data frame with sample variants

- pgx_allele:

  data frame with CPIC allele definitions

## Value

data frame with one row per sample variant and allele definition
