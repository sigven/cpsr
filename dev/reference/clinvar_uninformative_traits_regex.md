# Regular expression matching uninformative ClinVar trait names

Trait names that carry no information about the phenotype (e.g. 'not
provided'). Records with only such traits are given an unknown (NA)
cancer phenotype status, rather than a non-cancer status.

## Usage

``` r
clinvar_uninformative_traits_regex
```
