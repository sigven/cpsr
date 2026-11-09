# Impact level of a CPIC phenotype (for color coding)

Impact level of a CPIC phenotype (for color coding)

## Usage

``` r
pgx_phenotype_level(phenotype)
```

## Arguments

- phenotype:

  CPIC phenotype, e.g. 'Poor Metabolizer'. Multi-gene phenotypes ('TPMT:
  Normal Metabolizer; NUDT15: Poor Metabolizer') are assigned the
  highest impact level among the genes

## Value

integer, 1 (normal) to 5 (highest impact), or 0 when the phenotype is
indeterminate/not determined
