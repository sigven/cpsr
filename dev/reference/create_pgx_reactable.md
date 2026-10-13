# Function that creates a reactable of CPIC pharmacogenomic phenotypes or prescribing recommendations, for display in the germline report

The phenotype table has one row per gene (detected variants, CPIC
diplotype, phenotype); each row can be expanded to show the
variant-level details (CPIC allele/function, consequence, ClinVar,
gnomAD) and the CPIC consultation text. The recommendation table has one
row per drug and phenotype; each row can be expanded to show the
implications, CPIC comments and guideline. Phenotypes are shown as blue
badges, darker for higher impact.

## Usage

``` r
create_pgx_reactable(
  data = NULL,
  type = "phenotype",
  color_palette = NULL,
  variants = NULL
)
```

## Arguments

- data:

  data frame with gene phenotypes (assign_pgx_phenotypes) or drug
  recommendations (assign_pgx_recommendations)

- type:

  'phenotype' or 'recommendation'

- color_palette:

  CPSR color_palette object

- variants:

  data frame with pharmacogenomic variants (display version of
  retrieve_pgx_calls output), used for row details of the phenotype
  table
