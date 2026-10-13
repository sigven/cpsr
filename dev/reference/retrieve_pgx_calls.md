# Function that retrieves pharmacogenomic variants, i.e. variants that define CPIC alleles (DPYD, TPMT, NUDT15, G6PD) with altered function

Sample variants are matched (CHROM/POS/REF/ALT) against CPIC allele
definitions in the reference data bundle (pgx_allele). Each matched
variant is annotated with the CPIC allele(s) it defines, together with
allele function and activity value. Alleles defined by multiple variants
(e.g. TPMT \*3A, DPYD HapB3) are flagged as partial when not all
defining variants are found in the sample (phase is not assessed).

## Usage

``` r
retrieve_pgx_calls(calls, ref_data = NULL, include_normal_function = FALSE)
```

## Arguments

- calls:

  data frame with all calls found

- ref_data:

  object with PCGR/CPSR reference data

- include_normal_function:

  logical, report also variants that only define normal function alleles
