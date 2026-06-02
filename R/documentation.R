#' Get documentation string for the CPSR report introduction
#'
#' @return A documentation string
#' @export
cpsr_intro_doc_note <- function() {

  doc_md_file <- system.file(
    "templates", "doc_notes_md", "cpsr_intro.md", package = "cpsr")
  template <- paste0(readLines(doc_md_file, warn = FALSE), collapse = "\n")

  return(glue::glue(template))
}

#' Get documentation string for variant classification synopsis
#'
#' @param quarto Logical; if TRUE, use Quarto cross-reference syntax for
#'   internal links (@nte-*). Set to FALSE when rendering in plain R Markdown
#'   contexts (e.g. vignettes) where those references are not resolved.
#'
#' @return A documentation string
#' @export
variant_classification_synopsis_doc_note <- function(quarto = TRUE) {

  doc_md_file <- system.file(
    "templates", "doc_notes_md", "variant_classification_synopsis.md",
    package = "cpsr")
  template <- paste0(readLines(doc_md_file, warn = FALSE), collapse = "\n")

  if (quarto) {
    criteria_ref  <- "@nte-table-criteria"
    threshold_ref <- "@nte-threshold-calibration"
  } else {
    criteria_ref  <- "the criteria table below"
    threshold_ref <- "the calibration section below"
  }

  return(glue::glue(
    template,
    criteria_ref  = criteria_ref,
    threshold_ref = threshold_ref
  ))
}

#' Get documentation string for classification threshold calibration (intro)
#'
#' Returns the prose that precedes the calibration plot and score-tier table
#' in the threshold calibration section.
#'
#' @return A documentation string
#' @export
classification_thresholds_intro_doc_note <- function() {

  doc_md_file <- system.file(
    "templates", "doc_notes_md", "classification_thresholds_intro.md",
    package = "cpsr")
  template <- paste0(readLines(doc_md_file, warn = FALSE), collapse = "\n")

  return(glue::glue(template))
}

#' Get documentation string for classification threshold calibration (outro)
#'
#' Returns the prose that follows the score-tier table in the threshold
#' calibration section (score-weight downgrading rules and caveats).
#'
#' @return A documentation string
#' @export
classification_thresholds_outro_doc_note <- function() {

  doc_md_file <- system.file(
    "templates", "doc_notes_md", "classification_thresholds_outro.md",
    package = "cpsr")
  template <- paste0(readLines(doc_md_file, warn = FALSE), collapse = "\n")

  return(glue::glue(template))
}

#' Get documentation string for gnomAD usage in germline variant classification
#'
#' @return A documentation string
#' @export
gnomad_germline_doc_note <- function() {

  doc_md_file <- system.file(
    "templates", "doc_notes_md", "gnomad_germline.md", package = "cpsr")
  template <- paste0(readLines(doc_md_file, warn = FALSE), collapse = "\n")

  return(glue::glue(template))
}
