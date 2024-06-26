
#' Scores and documentation of ACMG evidence criteria used for variant classification
#' in CPSR
#'
#' @format A list object with three elements: 'score2tier', 'evidence_codes',
#' 'pathogenic_range_gnomad'
#'
#' @format \bold{score2tier} - A data frame with 5 rows and two columns that indicate
#' current score thresholds for variant classification in CPSR:
#' \itemize{
#'   \item \emph{CPSR_CLASSIFICATION} - variant classification level (P, LP, VUS etc)
#'   \item \emph{CPSR_PATHOGENICITY_SCORE} - indication of CPSR "score bucket" for a given
#'   classification (HTML string)
#' }
#'
#' #' @format \bold{evidence_codes} - A data frame with 34 rows and 7 columns that document
#' all ACMG evidence criteria that are used for variant classification in CPSR:
#' \itemize{
#'   \item \emph{cpsr_evidence_code} - code for evidence criterion ('ACMG_BA1_AD' etc)
#'   \item \emph{category} - type of evidence feature ('clinpop','funcvarpop','funcvar','funccomp')
#'   \item \emph{pathogenicity_pole} - whether the given evidence support a benign variant character ('B'), or
#'   pathogenic character ('P')
#'   \item \emph{category_long} - long version of 'category' column
#'   \item \emph{description} - Verbose description for the given evidence criterion
#'   \item \emph{sherloc_code} - Corresponding code identifier in SherLoc (Nykamp et al., GiM, 2017)
#'   \item \emph{path_score} - Score associated with the given evidence criterion (negative or positive)
#' }
#'
"acmg"

#' Format of CPSR output data frames (HTML, TSV)
#'
#' @format A list object with four elements: 'html_tier', 'html_sf',
#' 'tsv','html_gwas'
"col_format_output"
