# Use Bootstrap icons (as inline SVG)

Copied from https://github.com/rstudio/bsicons.

## Usage

``` r
bs_icon2(
  name,
  size = "1em",
  class = NULL,
  title = NULL,
  a11y = c("auto", "deco", "sem", "none"),
  ...
)
```

## Arguments

- name:

  The name of the Bootstrap icon. Whitespace is replaced with `-` (that
  way, `"arrow up"` can be used to refer to the "actual name" of
  `"arrow-up"`). For a searchable list of names, see
  <https://icons.getbootstrap.com/>

- size:

  Any valid CSS unit defining both the height and width of the icon.

- class:

  Additional CSS classes to add to the `<svg>` element. Consider
  providing Bootstrap 5+ utility classes (e.g., `text-success`) here to
  stylize the icon (but also note that those utility classes will only
  work when Bootstrap 5+ is on the page).

- title:

  If provided (highly recommended), `a11y` defaults to `"sem"`, meaning
  the title is used for on-hover text and screen reader announcements.

- a11y:

  Cases that distinguish the role of the icon and inform which
  accessibility attributes to be used. Icons can either be `"deco"`
  (decorative, the default case), `"sem"` (semantic), `"none"` (no
  accessibility features). The default, `"auto"`, resolves to `"sem"` if
  a `title` is provided (and `"deco"` otherwise).

- ...:

  additional CSS properties (e.g., `margin`, `position`, etc.) placed on
  the `<svg>` tag.

## Value

An
[`htmltools::HTML()`](https://rstudio.github.io/htmltools/reference/HTML.html)
string containing the SVG icon.
