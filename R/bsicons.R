# Copied from https://github.com/rstudio/bsicons/blob/f7dca31/R/icons.R
icon_info_mini <- list(
  name = c(
    "bullseye",
    "check-square"
  ),
  contents = c(
    "<path d=\"M8 15A7 7 0 1 1 8 1a7 7 0 0 1 0 14zm0 1A8 8 0 1 0 8 0a8 8 0 0 0 0 16z\"></path>\n<path d=\"M8 13A5 5 0 1 1 8 3a5 5 0 0 1 0 10zm0 1A6 6 0 1 0 8 2a6 6 0 0 0 0 12z\"></path>\n<path d=\"M8 11a3 3 0 1 1 0-6 3 3 0 0 1 0 6zm0 1a4 4 0 1 0 0-8 4 4 0 0 0 0 8z\"></path>\n<path d=\"M9.5 8a1.5 1.5 0 1 1-3 0 1.5 1.5 0 0 1 3 0z\"></path>",
    "<path d=\"M14 1a1 1 0 0 1 1 1v12a1 1 0 0 1-1 1H2a1 1 0 0 1-1-1V2a1 1 0 0 1 1-1h12zM2 0a2 2 0 0 0-2 2v12a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V2a2 2 0 0 0-2-2H2z\"></path>\n<path d=\"M10.97 4.97a.75.75 0 0 1 1.071 1.05l-3.992 4.99a.75.75 0 0 1-1.08.02L4.324 8.384a.75.75 0 1 1 1.06-1.06l2.094 2.093 3.473-4.425a.235.235 0 0 1 .02-.022z\"></path>"
  )
)

#' Use Bootstrap icons (as inline SVG)
#'
#' Copied from https://github.com/rstudio/bsicons.
#' @param name The name of the Bootstrap icon. Whitespace is replaced with `-`
#'   (that way, `"arrow up"` can be used to refer to the "actual name" of
#'   `"arrow-up"`). For a searchable list of names, see <https://icons.getbootstrap.com/>
#' @param size Any valid CSS unit defining both the height and width of the
#'   icon.
#' @param class Additional CSS classes to add to the `<svg>` element. Consider
#'   providing Bootstrap 5+ utility classes (e.g., `text-success`) here to
#'   stylize the icon (but also note that those utility classes will only work
#'   when Bootstrap 5+ is on the page).
#' @param title If provided (highly recommended), `a11y` defaults to `"sem"`,
#'   meaning the title is used for on-hover text and screen reader
#'   announcements.
#' @param a11y Cases that distinguish the role of the icon and inform which
#'   accessibility attributes to be used. Icons can either be `"deco"`
#'   (decorative, the default case), `"sem"` (semantic), `"none"` (no
#'   accessibility features). The default, `"auto"`, resolves to `"sem"` if a
#'   `title` is provided (and `"deco"` otherwise).
#' @param ... additional CSS properties (e.g., `margin`, `position`, etc.)
#'   placed on the `<svg>` tag.
#'
#' @return An [htmltools::HTML()] string containing the SVG icon.
#' @export
bs_icon2 <- function(name, size = "1em", class = NULL, title = NULL,
                     a11y = c("auto", "deco", "sem", "none"), ...) {
  if (length(name) != 1) {
    rlang::abort("The number of icons specified in `name` must be 1.")
  }

  "%||%" <- function(x, y) {
    if (is.null(x)) y else x
  }

  name <- sub("\\s+", "-", tolower(name))
  idx <- match(name, tolower(icon_info_mini$name))
  stopifnot(!is.na(idx))
  svg_children <- icon_info_mini$contents[idx]
  size <- htmltools::validateCssUnit(size)
  style_attr <- paste0(
    "height:", size, ";",
    "width:", size, ";",
    "fill:currentColor;",
    # Better default vertical positioning of icons in a inline context (inspired by fontawesome::fa())
    "vertical-align:-0.125em;",
    htmltools::css(...)
  )

  # Generate accessibility attributes if either of
  # the "deco" or "sem" cases are chosen
  a11y <- rlang::arg_match(a11y)
  a11y_attrs <- ""

  if (a11y == "auto") {
    a11y <- if (is.null(title)) "deco" else "sem"
  }

  if (a11y == "deco") {
    a11y_attrs <- 'aria-hidden="true" role="img" '
  } else if (a11y == "sem") {
    title <- title %||% name
    a11y_attrs <- sprintf(
      'aria-label="%s" role="img" ',
      htmltools::htmlEscape(title, attribute = TRUE)
    )
  }

  res <- sprintf(
    '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 16 16" class="bi bi-%s %s" style="%s" %s>%s%s</svg>',
    name,
    paste(class, collapse = " "),
    style_attr,
    a11y_attrs,
    if (is.null(title)) "" else paste0("<title>", htmltools::htmlEscape(title), "</title>"),
    svg_children
  )

  htmltools::browsable(htmltools::HTML(res))
}
