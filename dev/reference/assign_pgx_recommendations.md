# Map inferred pharmacogene phenotypes to CPIC prescribing recommendations

Recommendations are retrieved (pgx_guideline, reference data bundle) for
all drugs whose CPIC lookup key involves at least one gene with an
inferred phenotype (assign_pgx_phenotypes). For multi-gene
recommendations (e.g. thiopurines - TPMT and NUDT15), genes without
detected altered-function alleles are assumed to be normal metabolizers,
which is stated in PGX_NOTE.

## Usage

``` r
assign_pgx_recommendations(gene_phenotypes, ref_data = NULL)
```

## Arguments

- gene_phenotypes:

  data frame with inferred gene phenotypes

- ref_data:

  object with PCGR/CPSR reference data

## Value

data frame with one row per drug (and CPIC population)
