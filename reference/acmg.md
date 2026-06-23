# Scores and documentation of ACMG evidence criteria used for variant classification in CPSR

Scores and documentation of ACMG evidence criteria used for variant
classification in CPSR

## Usage

``` r
acmg
```

## Format

A list object with three elements: 'score2tier', 'evidence_codes',
'pathogenic_range_gnomad'

**score2tier** - A data frame with 5 rows and two columns that indicate
current score thresholds for variant classification in CPSR:

- *CPSR_CLASSIFICATION* - variant classification level (P, LP, VUS etc)

- *CPSR_PATHOGENICITY_SCORE* - indication of CPSR "score bucket" for a
  given classification (HTML string)

\#' @format **evidence_codes** - A data frame with 34 rows and 7 columns
that document all ACMG evidence criteria that are used for variant
classification in CPSR:

- *cpsr_evidence_code* - code for evidence criterion ('ACMG_BA1_AD' etc)

- *category* - type of evidence feature
  ('clinpop','funcvarpop','funcvar','funccomp')

- *pathogenicity_pole* - whether the given evidence support a benign
  variant character ('B'), or pathogenic character ('P')

- *category_long* - long version of 'category' column

- *description* - Verbose description for the given evidence criterion

- *sherloc_code* - Corresponding code identifier in SherLoc (Nykamp et
  al., GiM, 2017)

- *path_score* - Score associated with the given evidence criterion
  (negative or positive)
