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
variant_classification_synopsis_doc_note <- function(
    quarto = TRUE, trust_level = NULL) {

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
    threshold_ref = threshold_ref,
    trust_levels  = clinvar_trust_levels_doc_note(trust_level)
  ))
}

#' Get documentation string for the ClinVar trust levels
#'
#' Table of the levels of the 'clinvar_trust_level' setting, i.e. when the
#' CPSR classification takes precedence over an existing ClinVar
#' classification (as implemented in assign_classification_authority())
#'
#' @param trust_level trust level used for the report (0-4), marked in the
#' table; NULL for none
#'
#' @return A documentation string (markdown)
#' @export
clinvar_trust_levels_doc_note <- function(trust_level = NULL) {
  levels <- c(
    "0" = "Conflicting ClinVar interpretations only - ClinVar is trusted otherwise (default)",
    "1" = "As level 0, and ClinVar records with zero gold stars (no assertion criteria provided)",
    "2" = "As level 1, and ClinVar records with one gold star (e.g. criteria provided, single submitter)",
    "3" = "As level 2, and ClinVar records with non-cancer phenotypes only (regardless of review status)",
    "4" = "All ClinVar records - CPSR always classifies, ClinVar never takes precedence")
  level_col <- names(levels)
  if (!is.null(trust_level) && as.character(trust_level) %in% level_col) {
    i <- match(as.character(trust_level), level_col)
    level_col[i] <- paste0(level_col[i], " (this report)")
  }
  paste0(
    "The *clinvar_trust_level* setting (`--clinvar_trust_level`) determines ",
    "for which ClinVar records the CPSR classification takes precedence. ",
    "At all levels, variants absent from ClinVar are classified by CPSR.\n\n",
    "| Trust level | CPSR classification takes precedence over |\n",
    ## column widths (pandoc): ratio of the dashes
    "|:-----|:-----------------------------|\n",
    paste0("| **", level_col, "** | ", levels, " |", collapse = "\n"),
    "\n")
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
