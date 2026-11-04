# Regular expression matching cancer-related ClinVar trait names

Applied (case-insensitive, to lower-cased trait names) when determining
whether a ClinVar record is attributed to a cancer-related phenotype.
Moved from pcgrr
([`pcgrr::cancer_phenotypes_regex`](https://sigven.github.io/pcgr/reference/cancer_phenotypes_regex.html)),
where it was used by CPSR only.

## Usage

``` r
cancer_phenotypes_regex
```
