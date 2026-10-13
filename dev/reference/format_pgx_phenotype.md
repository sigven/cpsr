# Readable phenotype of a CPIC recommendation

Combines gene-keyed CPIC phenotype and activity score strings, e.g.
'DPYD:Intermediate Metabolizer' and 'DPYD:1.0' into 'DPYD: Intermediate
Metabolizer (activity score 1.0)'; multi-gene phenotypes are separated
by '; '.

## Usage

``` r
format_pgx_phenotype(phenotype, activity_score = NA)
```

## Arguments

- phenotype:

  gene-keyed CPIC phenotype(s)

- activity_score:

  gene-keyed CPIC activity score(s), or NA
