# Get documentation string for the ClinVar trust levels

Table of the levels of the 'clinvar_trust_level' setting, i.e. when the
CPSR classification takes precedence over an existing ClinVar
classification (as implemented in assign_classification_authority())

## Usage

``` r
clinvar_trust_levels_doc_note(trust_level = NULL)
```

## Arguments

- trust_level:

  trust level used for the report (0-4), marked in the table; NULL for
  none

## Value

A documentation string (markdown)
