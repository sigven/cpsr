# Short display name of the virtual gene panel of a CPSR run

Curated short name (`panel_short_names`) for single panels, a compact
label for multiple panels and custom gene lists, and otherwise the panel
name without the '(GEP)' suffix, truncated to `max_chars` characters.

## Usage

``` r
panel_display_name(conf = NULL, max_chars = 42)
```

## Arguments

- conf:

  CPSR configuration object with run settings

- max_chars:

  maximum number of characters

## Value

list with 'short' (display name) and 'full' (full panel name)
