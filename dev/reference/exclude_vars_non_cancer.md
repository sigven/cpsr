# Function that excludes variants attributed to non-cancer related phenotypes in ClinVar (option –clinvar_cancer_only)

Only variants with a ClinVar record where all informative traits are
non-cancer related (CLINVAR_PHENOTYPE_CANCER == 0) are excluded.
Variants with cancer-related traits (1), uninformative/unmapped traits
only (NA), or no ClinVar record are retained.

## Usage

``` r
exclude_vars_non_cancer(var_calls = NULL, conf = NULL)
```

## Arguments

- var_calls:

  variants in cancer predisposition genes

- conf:

  CPSR configuration object with run settings
