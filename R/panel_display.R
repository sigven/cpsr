#' Short display names of the CPSR virtual gene panels
#'
#' Named character vector (names: full panel names as listed in the
#' reference data bundle, values: short names of max. 42 characters) used
#' when space is limited in the reports (e.g. the virtual gene panel card at
#' the top of the HTML report). Full panel names are used elsewhere.
#'
#' @keywords internal
panel_short_names <- c(
  "CPSR exploratory cancer predisposition panel (PanelApp genes / TCGA's germline study / Cancer Gene Census / Other)" =
    "Exploratory - all predisposition genes",
  "Adult solid tumours cancer susceptibility (GEP)" =
    "Adult solid tumours susceptibility",
  "Adult solid tumours for rare disease (GEP)" =
    "Adult solid tumours (rare disease)",
  "Bladder cancer pertinent cancer susceptibility (GEP)" =
    "Bladder cancer susceptibility",
  "Brain cancer pertinent cancer susceptibility (GEP)" =
    "Brain cancer susceptibility",
  "Breast cancer pertinent cancer susceptibility (GEP)" =
    "Breast cancer susceptibility",
  "Childhood solid tumours cancer susceptibility (GEP)" =
    "Childhood solid tumours susceptibility",
  "Colorectal cancer pertinent cancer susceptibility (GEP)" =
    "Colorectal cancer susceptibility",
  "Endometrial cancer pertinent cancer susceptibility (GEP)" =
    "Endometrial cancer susceptibility",
  "Familial Tumours Syndromes of the central & peripheral Nervous system (GEP)" =
    "Familial nervous system tumour syndromes",
  "Familial breast cancer (GEP)" =
    "Familial breast cancer",
  "Familial melanoma (GEP)" =
    "Familial melanoma",
  "Familial prostate cancer (GEP)" =
    "Familial prostate cancer",
  "Familial rhabdomyosarcoma (GEP)" =
    "Familial rhabdomyosarcoma",
  "GI tract tumours (GEP)" =
    "GI tract tumours",
  "Genodermatoses with malignancies (GEP)" =
    "Genodermatoses with malignancies",
  "Haematological malignancies cancer susceptibility (GEP)" =
    "Haematological malignancies susceptibility",
  "Haematological malignancies for rare disease (GEP)" =
    "Haematological malignancies (rare disease)",
  "Head and neck cancer pertinent cancer susceptibility (GEP)" =
    "Head and neck cancer susceptibility",
  "Inherited MMR deficiency (Lynch Syndrome) - GEP" =
    "Lynch syndrome (inherited MMR deficiency)",
  "Inherited non-medullary thyroid cancer (GEP)" =
    "Inherited non-medullary thyroid cancer",
  "Inherited ovarian cancer (without breast cancer) (GEP)" =
    "Inherited ovarian cancer (without breast)",
  "Inherited pancreatic cancer (GEP)" =
    "Inherited pancreatic cancer",
  "Inherited polyposis and early onset colorectal cancer (GEP)" =
    "Polyposis / early-onset colorectal cancer",
  "Inherited predisposition to acute myeloid leukaemia (AML) (GEP)" =
    "Inherited predisposition to AML",
  "Inherited susceptibility to acute lymphoblastoid leukaemia (ALL) (GEP)" =
    "Inherited susceptibility to ALL",
  "Inherited predisposition to GIST (GEP)" =
    "Inherited predisposition to GIST",
  "Inherited renal cancer (GEP)" =
    "Inherited renal cancer",
  "Inherited phaeochromocytoma and paraganglioma (GEP)" =
    "Phaeochromocytoma and paraganglioma",
  "Melanoma pertinent cancer susceptibility (GEP)" =
    "Melanoma susceptibility",
  "Multiple endocrine tumours (GEP)" =
    "Multiple endocrine tumours",
  "Multiple monogenic benign skin tumours (GEP)" =
    "Monogenic benign skin tumours",
  "Neuroendocrine cancer pertinent cancer susceptibility (GEP)" =
    "Neuroendocrine cancer susceptibility",
  "Neurofibromatosis Type 1 (GEP)" =
    "Neurofibromatosis type 1",
  "Ovarian cancer pertinent cancer susceptibility (GEP)" =
    "Ovarian cancer susceptibility",
  "Parathyroid Cancer (GEP)" =
    "Parathyroid cancer",
  "Prostate cancer pertinent cancer susceptibility (GEP)" =
    "Prostate cancer susceptibility",
  "Renal cancer pertinent cancer susceptibility (GEP)" =
    "Renal cancer susceptibility",
  "Sarcoma cancer susceptibility (GEP)" =
    "Sarcoma cancer susceptibility",
  "Sarcoma susceptibility (GEP)" =
    "Sarcoma susceptibility",
  "Thyroid cancer pertinent cancer susceptibility (GEP)" =
    "Thyroid cancer susceptibility",
  "Tumour predisposition - childhood onset (GEP)" =
    "Childhood-onset tumour predisposition",
  "Upper gastrointestinal cancer pertinent cancer susceptibility (GEP)" =
    "Upper GI cancer susceptibility",
  "DNA repair genes pertinent cancer susceptibility (GEP)" =
    "DNA repair cancer susceptibility"
)

#' Shorten a text to a maximum number of characters
#'
#' Cut at the last word boundary within the limit, followed by an ellipsis.
#'
#' @param x character string
#' @param max_chars maximum number of characters (including the ellipsis)
#'
#' @keywords internal
truncate_label <- function(x, max_chars = 42) {
  if (is.na(x) || nchar(x) <= max_chars) {
    return(x)
  }
  cut <- substr(x, 1, max_chars - 1)
  ## last word boundary, if not too far back
  ws <- regexpr("[[:space:]][^[:space:]]*$", cut)
  if (ws > max_chars / 2) {
    cut <- substr(cut, 1, ws - 1)
  }
  paste0(sub("[[:space:][:punct:]]+$", "", cut), "…")
}

#' Short display name of the virtual gene panel of a CPSR run
#'
#' Curated short name (\code{panel_short_names}) for single panels, a
#' compact label for multiple panels and custom gene lists, and otherwise
#' the panel name without the '(GEP)' suffix, truncated to
#' \code{max_chars} characters.
#'
#' @param conf CPSR configuration object with run settings
#' @param max_chars maximum number of characters
#'
#' @return list with 'short' (display name) and 'full' (full panel name)
#' @keywords internal
panel_display_name <- function(conf = NULL, max_chars = 42) {
  gp <- conf[["gene_panel"]]
  panel_id <- as.character(gp[["panel_id"]])
  if (length(panel_id) != 1 || is.na(panel_id)) panel_id <- ""

  is_set <- function(x) {
    !is.null(x) && length(x) == 1 && !is.na(x) && !x %in% c("", "None")
  }
  full <- if (is_set(gp[["description"]])) gp[["description"]] else ""
  if (!panel_id %in% c("-1", "0") && !grepl(",", panel_id) &&
      is_set(gp[["description_trait"]])) {
    full <- gp[["description_trait"]]
  }

  short <- if (panel_id == "-1") {
    "Custom gene list (from panel 0)"
  } else if (grepl(",", panel_id)) {
    ids <- trimws(strsplit(panel_id, ",")[[1]])
    lbl <- paste0(length(ids), " PanelApp panels (",
                  paste(ids, collapse = ", "), ")")
    ## panel identifiers only if they fit
    if (nchar(lbl) > max_chars) paste0(length(ids), " PanelApp panels") else lbl
  } else if (full %in% names(panel_short_names)) {
    panel_short_names[[full]]
  } else if (panel_id == "0") {
    panel_short_names[[1]]
  } else {
    stringr::str_trim(stringr::str_remove(full, "\\s*(\\(GEP\\)|- GEP)$"))
  }

  list(short = truncate_label(short, max_chars), full = full)
}
