# Function that assigns a summary string of all ACMG evidence codes met for a given variant (e.g. "PVS1\|PS1\|PM2_supporting\|PP3")

Function that assigns a summary string of all ACMG evidence codes met
for a given variant (e.g. "PVS1\|PS1\|PM2_supporting\|PP3")

## Usage

``` r
assign_acmg_consensus(var_calls = NULL)
```

## Arguments

- var_calls:

  data frame with variant calls in predisposition genes

## Value

var_calls data frame with ACMG_CODE column appended
