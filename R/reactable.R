# Cell renderer factories

#' Cell renderer factory for classification column
#' Returns a function that renders classification values as
#' colored pills based on the provided color palette.
#'
#' @param color_palette CPSR color palette object containing
#' pathogenicity styling information
#' @return A function that takes a classification value and returns
#' an HTML span element with appropriate styling
#'
#' @export
#'
rt_cell_classification <- function(color_palette) {
  function(value) {
    if (is.na(value)) return("-")
    idx <- match(value, color_palette$pathogenicity$levels)
    bg <- if (!is.na(idx)) color_palette$pathogenicity$values[idx] else "#999"
    htmltools::span(
      style = list(
        background = bg, color = "white",
        padding = "6px 10px", borderRadius = "4px",
        fontWeight = "bold", display = "inline-block",
        whiteSpace = "nowrap", fontSize = "0.92em"
      ),
      value
    )
  }
}

#' Cell renderer factory for classification rank column
#' Maps a numeric rank (1=Benign ... 5=Pathogenic) to a labelled
#' colored pill, enabling correct sort order via JavaScript.
#'
#' @param color_palette CPSR color palette object
#' @return A function(value) returning an HTML span
#'
#' @export
#'
rt_cell_classification_rank <- function(color_palette) {
  labels <- c("Benign", "Likely Benign", "VUS", "Likely Pathogenic", "Pathogenic")
  function(value) {
    if (is.na(value)) return("-")
    label <- labels[as.integer(value)]
    if (is.na(label)) return("-")
    idx <- match(label, color_palette$pathogenicity$levels)
    bg <- if (!is.na(idx)) color_palette$pathogenicity$values[idx] else "#999"
    htmltools::span(
      style = list(
        background = bg, color = "white",
        padding = "6px 10px", borderRadius = "4px",
        fontWeight = "bold", display = "inline-block",
        whiteSpace = "nowrap", fontSize = "0.92em"
      ),
      label
    )
  }
}

#' Display label of a CPSR genotype code
#'
#' @param genotype genotype code(s) ('het', 'hom_alt', 'hom_ref',
#' 'undefined')
#'
#' @return readable genotype(s), e.g. 'heterozygous', 'homozygous'
#' @keywords internal
genotype_label <- function(genotype) {
  dplyr::case_when(
    genotype == "het" ~ "heterozygous",
    genotype == "hom_alt" ~ "homozygous",
    genotype == "hom_ref" ~ "homozygous (ref)",
    TRUE ~ as.character(genotype))
}

#' Cell renderer factory for genotype column
#' Returns a function that renders genotype values as colored
#' pills based on the provided color palette.
#'
#' @param color_palette CPSR color palette object containing genotype styling information
#' @return A function that takes a genotype value and returns an HTML span element with appropriate styling
#'
#' @export
#'
rt_cell_genotype <- function(color_palette) {
  function(value) {
    if (is.na(value)) return("-")
    idx <- match(value, color_palette$genotypes$levels)
    if (is.na(idx)) return(value)
    value <- genotype_label(value)
    htmltools::span(
      style = list(
        background = color_palette$genotypes$bgcolor_values[idx],
        color      = color_palette$genotypes$color_values[idx],
        padding = "4px 8px", borderRadius = "3px",
        fontWeight = "bold", display = "inline-block",
        fontSize = "0.92em"
      ),
      value
    )
  }
}

#' Cell renderer factory for biomarker clinical significance
#' Returns a function that renders biomarker clinical significance
#' values as colored pills based on the provided color palette.
#'
#' @param color_palette CPSR color palette object containing
#' biomarker_types styling information
#' @return A function that takes a biomarker clinical significance value and returns an HTML div element containing styled pills for each significance category
#'
#' @export
#'
rt_cell_bm_significance <- function(color_palette) {
  function(value) {
    if (is.na(value) || value == "") return("-")
    items <- trimws(strsplit(value, ",")[[1]])
    pills <- lapply(items, function(item) {
      idx <- match(item, color_palette$biomarker_types$levels)
      bg <- if (!is.na(idx)) color_palette$biomarker_types$values[idx] else "#999"
      htmltools::span(
        style = list(
          background = bg, color = "white",
          padding = "3px 8px", borderRadius = "4px",
          fontWeight = "bold", display = "inline-block",
          whiteSpace = "nowrap", fontSize = "0.92em",
          marginRight = "4px", marginBottom = "3px"
        ),
        item
      )
    })
    htmltools::div(
      style = list(display = "flex", flexWrap = "wrap", gap = "3px"),
      pills
    )
  }
}


#' Prepare unified variant dataset for single-table display
#'
#' Combines ClinVar and CPSR-classified variants across ALL significance classes
#' into a single data frame. Includes a variant_key for crosstalk linking.
#'
#' @param cps_report CPSR report object
#' @param max_rows Integer. Max rows per source
#'
#' @return Data frame with all variants from both sources, keyed for crosstalk
#'
#' @export
#'
prepare_unified_all_variants <- function(
    cps_report = NULL,
    max_rows = 1000) {

  assertthat::assert_that(
    !is.null(cps_report), msg = "cps_report is NULL")

  all_variants <-
    cps_report[["content"]][["snv_indel"]]$callset$variant_display$cpg_non_sf

  required_cols <- c(
    "ASSERTION_AUTHORITY",
    "CLASSIFICATION",
    "CLINVAR_GOLD_STARS",
    "CPSR_PATHOGENICITY_SCORE",
    "SYMBOL", "GENOMIC_CHANGE",
    "ALTERATION",
    "GERP_SCORE", "HGVSc",
    "HGVSp", "HGVSc_RefSeq",
    "CLINVAR_PHENOTYPE",
    "CLINVAR", "PROTEIN_DOMAIN",
    "GENENAME", "CODING_STATUS",
    "ACMG_CODE","PROTEIN_CHANGE"
  )

  if (is.null(all_variants) ||
      !all(required_cols %in% colnames(all_variants))) {
    return(
      data.frame(matrix(ncol = length(required_cols), nrow = 0)) |>
        stats::setNames(required_cols)
    )
  }

  all_variants <- all_variants |>
    dplyr::mutate(
      ACMG_CODE = dplyr::if_else(
        is.na(.data$ACMG_CODE) | .data$ACMG_CODE == "",
        "\u2014",
        .data$ACMG_CODE
      )
    ) |>
    dplyr::mutate(
      SEARCH_INDEX = stringr::str_replace_all(
        paste(
          pcgrr::strip_html(.data$PROTEIN_DOMAIN),
          pcgrr::strip_html(.data$GENENAME),
          .data$PROTEIN_CHANGE,
          .data$HGVSc,
          .data$HGVSp,
          .data$HGVSc_RefSeq,
          .data$CLINVAR_PHENOTYPE,
          pcgrr::strip_html(.data$CLINVAR),
          .data$GENOMIC_CHANGE,
          .data$ACMG_CODE,
          sep = " "
        ),
        "NA ",""
      )
    )


  # Split by source and arrange by quality metrics
  clinvar_data <- all_variants |>
    dplyr::filter(.data$ASSERTION_AUTHORITY == "ClinVar") |>
    dplyr::arrange(
      dplyr::desc(.data$CLASSIFICATION),
      dplyr::desc(.data$CLINVAR_GOLD_STARS),
      .data$SYMBOL
    ) |>
    utils::head(max_rows)

  cpsr_data <- all_variants |>
    dplyr::filter(.data$ASSERTION_AUTHORITY == "CPSR") |>
    dplyr::arrange(
      dplyr::desc(.data$CLASSIFICATION),
      dplyr::desc(.data$CPSR_PATHOGENICITY_SCORE),
      .data$SYMBOL
    ) |>
    utils::head(max_rows)

  # Combine and create unique key
  unified_data <-
    dplyr::bind_rows(
      clinvar_data,
      cpsr_data) |>
    dplyr::arrange(
      dplyr::desc(.data$CLASSIFICATION),
      .data$SYMBOL
    ) |>
    dplyr::mutate(
      VARKEY_CLASSIFICATION = paste(
        .data$SYMBOL,
        .data$GENOMIC_CHANGE,
        .data$ALTERATION,
        .data$ASSERTION_AUTHORITY, sep = "|"),
      .after = "SYMBOL"
    ) |>
    dplyr::mutate(
      dplyr::across(
        c("GERP_SCORE"),
        \(x) round(x, 3)
      )
    ) |>
    dplyr::mutate(
      CLASSIFICATION = factor(
        .data$CLASSIFICATION,
        levels = c(
          "Benign", "Likely Benign", "VUS",
          "Likely Pathogenic", "Pathogenic"),
        ordered = TRUE
      )
    ) |>
    dplyr::mutate(
      CLASSIFICATION_RANK = dplyr::case_when(
        as.character(.data$CLASSIFICATION) == "Pathogenic"        ~ 5L,
        as.character(.data$CLASSIFICATION) == "Likely Pathogenic" ~ 4L,
        as.character(.data$CLASSIFICATION) == "VUS"               ~ 3L,
        as.character(.data$CLASSIFICATION) == "Likely Benign"     ~ 2L,
        as.character(.data$CLASSIFICATION) == "Benign"            ~ 1L,
        TRUE ~ NA_integer_
      )
    )

  return(unified_data)
}


#' Create crosstalk-linked unified filters for all variants
#'
#' Shared filters across ClinVar and CPSR variants that control the single table.
#'
#' @param shared_data Crosstalk SharedData object
#' @param cps_report CPSR report object
#' @param coding_status Character. "coding" or "noncoding".
#' Determines which filters to include.
#'
#' @return bscols HTML widget with filters
#'
#' @export
#'
create_unified_variant_filters <- function(
    shared_data = NULL,
    cps_report = NULL,
    coding_status = "coding") {

  assertthat::assert_that(
    !is.null(shared_data), msg = "shared_data is NULL")
  assertthat::assert_that(
    !is.null(cps_report), msg = "cps_report is NULL")

  if (nrow(shared_data$data()) == 0) {
    return(htmltools::div())
  }

  crosstalk::bscols(
    list(
      crosstalk::filter_select(
        "CLASSIFICATION",
        "Clinical significance",
        shared_data, ~CLASSIFICATION
      ),
      crosstalk::filter_select(
        "ASSERTION_AUTHORITY",
        "Assertion authority",
        shared_data,
        ~ASSERTION_AUTHORITY
      )
    ),
    list(
      crosstalk::filter_select(
        "ACMG_CODE",
        "ACMG/AMP classification criteria",
        shared_data,
        ~ACMG_CODE
      ),
      if(coding_status == "noncoding"){
        crosstalk::filter_select(
          "TF_BINDING_SITE_VARIANT",
          "TF binding site alteration",
          shared_data, ~TF_BINDING_SITE_VARIANT)
      }
    #)
      # if(cps_report$settings$conf$sample_properties$dp_detected == 1) {
      #   crosstalk::filter_slider(
      #     "DP_CONTROL",
      #     "Sequencing depth",
      #     shared_data,
      #     ~DP_CONTROL
      #   )
      # }

      # if(coding_status == "noncoding"){
      #   crosstalk::filter_slider(
      #     "GERP_SCORE",
      #     "GERP conservation score",
      #     shared_data,
      #     ~GERP_SCORE
      #   )
      # }
    )
  )

}
#'
#' Applies consistent styling using CPSR color palette.
#'
#' @param color_palette CPSR color_palette object
#' @param header_color Optional hex color for the header background.
#'   Defaults to color_palette$report.
#'
#' @return reactableTheme object
#'
#' @export
#'
create_variant_table_theme <- function(
    color_palette = NULL,
    header_color = NULL) {

  assertthat::assert_that(
    !is.null(color_palette),
    msg = "color_palette is NULL")

  hdr_bg <- if (!is.null(header_color)) header_color else color_palette$report

  reactable::reactableTheme(
    style = list(
      fontFamily = "inherit",
      fontSize = "0.99em"
    ),
    headerStyle = list(
      background = hdr_bg,
      color = "white",
      fontFamily = "inherit",
      fontWeight = "bold",
      fontSize = "0.99em",
      borderRight = "1px solid rgba(255,255,255,0.3)",
      padding = "10px 8px",
      display = "flex",
      alignItems = "center"
    ),
    cellStyle = list(
      borderRight = "1px solid #e8e8e8",
      padding = "8px 10px",
      display = "flex",
      alignItems = "center"
    ),
    rowStyle = list(
      borderBottom = "1px solid #f0f0f0"
    ),
    stripedColor = "#fafafa",
    highlightColor = "#f0f7ff",
    borderColor = "#e0e0e0",
    searchInputStyle = list(
      borderColor = color_palette$report,
      fontSize = "0.99em"
    )
  )
}


#' Create the complete unified variant reactable
#'
#' Simple, fast reactable showing primary columns only.
#' All other columns accessible via nested details.
#'
#' @param shared_data Crosstalk SharedData object with unified variants
#' @param primary_cols Character vector of primary column names to display
#' @param color_palette CPSR color_palette object
#' @param flag_low_depth Logical. If TRUE, flag variants with a control-sample
#' sequencing depth (DP_CONTROL) below \code{dp_threshold} with a marker on the
#' genotype cell. Should only be set to TRUE when \code{DP_CONTROL} carries
#' real depth values for at least some variants (i.e. not the "-1" sentinel
#' used when CPSR is run without a matched control sample).
#' @param dp_threshold Integer. Depth threshold below which a variant is
#' considered low-depth. Default: 10.
#'
#' @return reactable widget
#'
#' @export
#'
create_unified_variant_reactable <- function(
    shared_data = NULL,
    primary_cols =
      c("SYMBOL",
        "gnomADg_AF",
        "CONSEQUENCE",
        "ALTERATION",
        "CLASSIFICATION",
        "GENOTYPE"),
    color_palette = NULL,
    flag_low_depth = FALSE,
    dp_threshold = 10) {

  assertthat::assert_that(
    !is.null(shared_data), msg = "data is NULL")
  assertthat::assert_that(
    !is.null(color_palette), msg = "color_palette is NULL")

  # # Extract actual data if SharedData object is passed
  # if (inherits(data, "SharedData")) {
  #   actual_data <- data$origData()
  # } else {
  #   actual_data <- data
  # }

  # Build column definitions for primary columns only
  col_defs <- list(
    VARKEY_CLASSIFICATION =
      reactable::colDef(show = FALSE),
    ASSERTION_AUTHORITY =
      reactable::colDef(show = FALSE),
    SEARCH_INDEX = reactable::colDef(
      show = FALSE,
      searchable = TRUE
    ),

    SYMBOL = reactable::colDef(
      name = "Gene",
      minWidth = 90,
      sticky = "left",
      style = list(fontWeight = "bold")
    ),
    gnomADg_AF = reactable::colDef(
      name = "gnomADg AF",
      minWidth = 100,
      sticky = "left",
      cell = function(value) {
        if (is.na(value)) return("-")
        if (value == 0) return("0")
        formatC(value, format = "e", digits = 2)
      }
    ),
    CONSEQUENCE = reactable::colDef(
      name = "Consequence",
      minWidth = 140
    ),
    ALTERATION = reactable::colDef(
      name = "Alteration",
      minWidth = 140,
      cell = function(value, index) {
        # Get the variant_source value for this row
        source <- shared_data$origData()[index, "ASSERTION_AUTHORITY"]

        # Choose border color based on source
        border_color <- if (source == "ClinVar") {
          "#0277bd"  # Blue for ClinVar
        } else {
          "#e65100"  # Orange for CPSR
        }

        htmltools::span(
          style = list(
            fontSize = "0.98em",
            borderLeft = paste0("6px solid ", border_color),
            paddingLeft = "8px",
            display = "inline-block"
          ),
          value
        )
      }
    ),
    CLASSIFICATION = reactable::colDef(show = FALSE),
    CLASSIFICATION_RANK = reactable::colDef(
      name = "Clinical significance",
      minWidth = 140,
      align = "center",
      defaultSortOrder = "desc",
      cell = rt_cell_classification_rank(color_palette)
    ),
    GENOTYPE = reactable::colDef(
      name = "Genotype",
      minWidth = 120,
      align = "center",
      cell = function(value, index) {
        if (is.na(value)) return("-")
        gidx <- match(value, color_palette$genotypes$levels)
        style <- if (!is.na(gidx)) {
          list(
            background = color_palette$genotypes$bgcolor_values[gidx],
            color      = color_palette$genotypes$color_values[gidx],
            padding = "4px 8px", borderRadius = "3px",
            fontWeight = "bold", display = "inline-block",
            fontSize = "0.92em"
          )
        } else {
          list()
        }
        if (isTRUE(flag_low_depth)) {
          dp <- shared_data$origData()[index, "DP_CONTROL"]
          if (!is.na(dp) && dp >= 0 && dp < dp_threshold) {
            style$borderLeft <- "4px solid #E69F00"
            style$paddingLeft <- "6px"
          }
        }
        if (length(style) == 0) return(value)
        htmltools::span(style = style, genotype_label(value))
      }
    )
  )

  # Hide all remaining non-primary columns
  for (col in colnames(shared_data$origData())) {
    if (!col %in% names(col_defs)) {
      col_defs[[col]] <- reactable::colDef(show = FALSE)
    }
  }

  theme <- create_variant_table_theme(
    color_palette = color_palette)

  ## make primary_cols2 (not including ASSERTION_AUTHORITY)
  ## - i want this shown in the row details
  primary_cols2 <- setdiff(primary_cols, "ASSERTION_AUTHORITY")

  reactable::reactable(
    shared_data,
    columns = col_defs,
    defaultColDef = reactable::colDef(html = TRUE),
    details = pcgrr::build_rt_row_details(
      primary_cols = c(primary_cols2,
                       "SEARCH_INDEX",
                       "CLASSIFICATION_RANK"),
      font_size = "0.97em"),
    searchable = TRUE,
    filterable = TRUE,
    highlight = TRUE,
    striped = TRUE,
    compact = TRUE,
    wrap = TRUE,
    defaultPageSize = 10,
    theme = theme,
    language = reactable::reactableLang(
      searchPlaceholder = "Search variants..."
    )
  )
}


#' Create a ClinVar variant reactable with expandable detail rows
#'
#' Displays ClinVar-classified variants with primary columns in the main row
#' and all remaining columns accessible via an expandable details row.
#' Used for secondary findings and pharmacogenomic tables.
#'
#' @param data Data frame of variants to display
#' @param primary_cols Character vector of primary column names to display
#' @param color_palette CPSR color_palette object
#'
#' @return reactable widget
#'
#' @export
#'
create_clinvar_reactable <- function(
    data = NULL,
    primary_cols = c(
      "SYMBOL",
      "ALTERATION",
      "CLINVAR_CLASSIFICATION",
      "CLINVAR_PHENOTYPE",
      "GENOTYPE",
      "CONSEQUENCE"),
    color_palette = NULL) {

  assertthat::assert_that(!is.null(data), msg = "data is NULL")
  assertthat::assert_that(!is.null(color_palette), msg = "color_palette is NULL")

  col_defs <- list(
    SYMBOL = reactable::colDef(
      name = "Gene",
      minWidth = 90,
      sticky = "left",
      style = list(fontWeight = "bold")
    ),
    ALTERATION = reactable::colDef(
      name = "Alteration",
      minWidth = 150
    ),
    CLINVAR_CLASSIFICATION = reactable::colDef(
      name = "Clinical significance",
      minWidth = 160,
      align = "center",
      cell = rt_cell_classification(color_palette)
    ),
    CLINVAR_PHENOTYPE = reactable::colDef(
      name = "ClinVar phenotype",
      minWidth = 220
    ),
    GENOTYPE = reactable::colDef(
      name = "Genotype",
      minWidth = 120,
      align = "center",
      cell = rt_cell_genotype(color_palette)
    ),
    CONSEQUENCE = reactable::colDef(
      name = "Consequence",
      minWidth = 130
    ),
    PGX_ALLELE = reactable::colDef(
      name = "CPIC allele",
      minWidth = 180
    ),
    PGX_ALLELE_FUNCTION = reactable::colDef(
      name = "Allele function",
      minWidth = 140
    ),
    PGX_ACTIVITY_VALUE = reactable::colDef(
      name = "Activity value",
      minWidth = 90,
      align = "center"
    )
  )
  col_defs <- col_defs[names(col_defs) %in% colnames(data)]
  ## defined columns that are not primary are shown in row details
  for (col in setdiff(names(col_defs), primary_cols)) {
    col_defs[[col]]$show <- FALSE
  }

  for (col in colnames(data)) {
    if (!col %in% names(col_defs)) {
      col_defs[[col]] <- reactable::colDef(show = FALSE)
    }
  }

  theme <- create_variant_table_theme(
    color_palette = color_palette,
    header_color = "#2c313c")

  reactable::reactable(
    data,
    columns = col_defs,
    defaultColDef = reactable::colDef(html = TRUE),
    details = pcgrr::build_rt_row_details(
      primary_cols = primary_cols,
      font_size = "0.97em"),
    searchable = TRUE,
    filterable = TRUE,
    highlight = TRUE,
    striped = TRUE,
    compact = TRUE,
    wrap = TRUE,
    defaultPageSize = 10,
    theme = theme,
    language = reactable::reactableLang(
      searchPlaceholder = "Search variants..."
    )
  )
}


#' Function that gathers data table on biomarker variants
#' for display in germline report
#'
#' @param rep report object
#' @param variant_category variant category
#' @export
#'
prep_biomarker_tbl <- function(
    rep = NULL,
    variant_category = "snv_indel") {

  if (is.null(rep)) {
    pcgrr::log4r_fatal("report object is NULL")
  }
  if (!variant_category %in% names(rep$content)) {
    pcgrr::log4r_fatal(paste0(
      "rep$content object does not contain '", variant_category,"'"))
  }

  ## check variant_category is valid
  if (!variant_category %in% c("snv_indel", "cna")) {
    pcgrr::log4r_fatal(
      "variant_category must be one of 'snv_indel', 'cna'")
  }

  if (!"callset" %in% names(rep$content[[variant_category]])) {
    pcgrr::log4r_fatal("rep$content$variant_category object does not contain 'callset'")
  }

  callset <- rep$content[[variant_category]]$callset

  if (is.null(callset$variant_display$biomarker) ||
      NROW(callset$variant_display$biomarker) == 0) {
    return(list(main = data.frame(), nested = data.frame()))
  }

  vars <- callset$variant_display$biomarker |>
    dplyr::filter(.data$BM_EVIDENCE_DIRECTION == "Supports") |>
    dplyr::filter(.data$BM_CLINICAL_SIGNIFICANCE %in%
                    c("Sensitivity/Response",
                      "Poor Outcome",
                      "Better Outcome",
                      "Resistance/Non-response",
                      "Predisposition",
                      "Positive",
                      "Toxicity",
                      "Adverse Response"))

  if (NROW(vars) == 0) {
    pcgrr::log4r_info(
      paste0("No bioimarker variants found."))
    return(list(main = data.frame(), nested = data.frame()))
  }

  rctbl_recs <- list()
  rctbl_recs[['main']] <- data.frame()
  rctbl_recs[['nested']] <- data.frame()

  if (NROW(vars) > 0) {
    ## debug: print colnames of vars and eitems

    vars <- vars |>
      dplyr::select(
        dplyr::any_of(
          c("VAR_ID",
            "SAMPLE_ALTERATION",
            "VARIANT_CLASS",
            "CONSEQUENCE",
            "GENOTYPE",
            "ENTREZGENE",
            "ASSERTION_AUTHORITY",
            "CLASSIFICATION",
            "BM_SOURCE_DB",
            "BM_REFERENCE",
            "BM_RATING",
            "BM_MOLECULAR_PROFILE",
            "BM_CANCER_TYPE",
            "BM_EVIDENCE_DESCRIPTION",
            "BM_EVIDENCE_TYPE",
            "BM_EVIDENCE_LEVEL",
            "BM_THERAPEUTIC_CONTEXT",
            "BM_CLINICAL_SIGNIFICANCE",
            "BM_PRIMARY_SITE",
            "BM_MAPPING_CONFIDENCE")
        )
      ) |>
      dplyr::mutate(
        BM_CLINICAL_SIGNIFICANCE = dplyr::if_else(
          .data$BM_CLINICAL_SIGNIFICANCE == "Positive",
          "Diagnostic",
          .data$BM_CLINICAL_SIGNIFICANCE
        )
      ) |>
      dplyr::mutate(
        CLASSIFICATION_RANK = dplyr::case_when(
          .data$CLASSIFICATION == "Pathogenic"        ~ 5L,
          .data$CLASSIFICATION == "Likely Pathogenic" ~ 4L,
          .data$CLASSIFICATION == "VUS"               ~ 3L,
          .data$CLASSIFICATION == "Likely Benign"     ~ 2L,
          .data$CLASSIFICATION == "Benign"            ~ 1L,
          TRUE ~ NA_integer_
        )
      )

    if (NROW(vars) > 0) {


      ## get all clinical significance categories for each
      ## variant

      biomarker_clinical_significance <-
        vars |>
        dplyr::select(
          c("VAR_ID",
            "ENTREZGENE",
            "BM_CLINICAL_SIGNIFICANCE")
        ) |>
        dplyr::group_by(
          .data$ENTREZGENE,
          .data$VAR_ID,
        ) |>
        dplyr::summarise(
          CLINICAL_SIGNIFICANCE = paste(
            unique(sort(.data$BM_CLINICAL_SIGNIFICANCE)),
            collapse = ", "),
          .groups = "drop"
        ) |>
        dplyr::distinct()

      ## get the highest confidence level and resolution for
      ## each variant, based on all evidence items associated
      ## with the variant
      biomarker_top_resolution <-
        vars |>
        dplyr::select(
          c("VAR_ID",
            "ENTREZGENE",
            "BM_MAPPING_CONFIDENCE")
        ) |>
        dplyr::group_by(
          .data$ENTREZGENE,
          .data$VAR_ID,
        ) |>
        dplyr::summarise(
          BM_MAPPING_CONFIDENCE = paste(
            unique(sort(.data$BM_MAPPING_CONFIDENCE)),
            collapse = ","),
          .groups = "drop"
        ) |>
        dplyr::mutate(
          BM_TOP_MAPPING_CONFIDENCE = dplyr::case_when(
            stringr::str_detect(
              .data$BM_MAPPING_CONFIDENCE, "high") ~ "high",
            stringr::str_detect(
              .data$BM_MAPPING_CONFIDENCE, "medium") ~ "medium",
            TRUE ~ "low"
          )
        ) |>
        dplyr::select(
          -c("BM_MAPPING_CONFIDENCE")) |>
        dplyr::distinct()

      ## across all evidence items, get the unique sources
      ## supporting biomarker evidence for each variant
      biomarker_source_support <-
        vars |>
        dplyr::group_by(
          .data$VAR_ID, .data$ENTREZGENE
        ) |>
        dplyr::summarise(
          BM_SOURCES = paste(
            unique(sort(.data$BM_SOURCE_DB)),
            collapse = "|"),
          .groups = "drop") |>
        dplyr::distinct()

      ## search index: concatenate nested text fields per variant so that
      ## the main-table search box can match content from evidence items
      biomarker_search_index <-
        vars |>
        dplyr::group_by(.data$VAR_ID, .data$ENTREZGENE) |>
        dplyr::summarise(
          SEARCH_INDEX = stringr::str_squish(paste(
            c(unique(stats::na.omit(.data$BM_EVIDENCE_DESCRIPTION)),
              unique(stats::na.omit(.data$BM_CANCER_TYPE)),
              unique(stats::na.omit(.data$BM_THERAPEUTIC_CONTEXT))),
            collapse = " "
          )),
          .groups = "drop"
        )

      ## for the main report table, we want to aggregate evidence items
      ## for each variant, and show the unique therapeutic contexts
      ## and primary sites associated with the variant, as well as the
      ## highest mapping confidence and resolution across all evidence items.
      ## For the nested table, we want to show all evidence items for each variant.
      ##
      rctbl_recs[['main']] <-
        vars |>
        dplyr::left_join(
          biomarker_top_resolution,
          by = c("VAR_ID","ENTREZGENE")
        ) |>
        dplyr::left_join(
          biomarker_source_support,
          by = c("VAR_ID","ENTREZGENE")
        ) |>
        dplyr::left_join(
          biomarker_clinical_significance,
          by = c("VAR_ID","ENTREZGENE")
        ) |>
        dplyr::arrange(
          .data$BM_TOP_MAPPING_CONFIDENCE,
          dplyr::desc(.data$BM_RATING)
        ) |>
        dplyr::select(
          c("VAR_ID",
            "ENTREZGENE",
            "SAMPLE_ALTERATION",
            "BM_SOURCES",
            "CLINICAL_SIGNIFICANCE",
            "ASSERTION_AUTHORITY",
            "CLASSIFICATION",
            "CLASSIFICATION_RANK",
            "GENOTYPE",
            "BM_TOP_MAPPING_CONFIDENCE")
        ) |>
        dplyr::distinct() |>
        dplyr::left_join(
          biomarker_search_index,
          by = c("VAR_ID", "ENTREZGENE")
        )

      rctbl_recs[['nested']] <- vars |>
        dplyr::select(
          c("VAR_ID",
            "ENTREZGENE",
            "BM_CLINICAL_SIGNIFICANCE",
            "BM_THERAPEUTIC_CONTEXT",
            "BM_REFERENCE",
            "BM_MOLECULAR_PROFILE",
            "BM_EVIDENCE_LEVEL",
            "BM_SOURCE_DB",
            "BM_CANCER_TYPE",
            "BM_EVIDENCE_DESCRIPTION")
        ) |>
        dplyr::arrange(
          .data$BM_EVIDENCE_LEVEL
        ) |>
        dplyr::mutate(
          BM_THERAPEUTIC_CONTEXT = stringr::str_replace_all(
            .data$BM_THERAPEUTIC_CONTEXT, ",", ", "
          )
        ) |>
        dplyr::distinct()
    }
  }

  return(rctbl_recs)

}



#' Build biomarker reactable with category-aware styling
#'
#' @param rctbl_recs List with $main and $nested data frames
#'   containing prepared biomarker table records, for example as
#'   returned by [prep_biomarker_tbl()].
#' @param variant_category One of "snv_indel", "cnv", "fusion"
#' @param color_palette color palette for therapeutic biomarkers
#'
#' @return A reactable object with the biomarker table
#' @export
#'
render_actble_bm_table <- function(
    rctbl_recs = NULL,
    variant_category = "snv_indel",
    color_palette = NULL) {

  ## check that rctbl_recs is
  ## 1. non-null
  ## 2. is a list object that contains two elements
  ## 3. both elements are data frames
  ## 4. main data frame contains required columns
  ##.   (pending upon variant_category)
  if (is.null(rctbl_recs) ||
      !is.list(rctbl_recs) ||
      !all(c("main", "nested") %in% names(rctbl_recs)) ||
      !is.data.frame(rctbl_recs$main) ||
      !is.data.frame(rctbl_recs$nested)) {
    pcgrr::log4r_fatal(
      "rctbl_recs must be a list with 'main' and 'nested' data frames")
  }


  if (NROW(rctbl_recs$main) == 0) {
    return(htmltools::div(
      style = "color:#666; font-style:italic; padding:8px;",
      "No biomarker evidence found."
    ))
  }

  required_cols <-
    c("VAR_ID",
      "ENTREZGENE",
      "BM_SOURCES",
      "GENOTYPE",
      "BM_TOP_MAPPING_CONFIDENCE",
      "CLINICAL_SIGNIFICANCE",
      "CLASSIFICATION",
      "ASSERTION_AUTHORITY")


  assertable::assert_colnames(
    rctbl_recs$main,
    required_cols,
    only_colnames = FALSE,
    quiet = TRUE
  )

  assertable::assert_colnames(
    rctbl_recs$nested,
    c("VAR_ID",
      "ENTREZGENE",
      "BM_MOLECULAR_PROFILE",
      "BM_REFERENCE",
      "BM_CLINICAL_SIGNIFICANCE",
      "BM_SOURCE_DB",
      "BM_CANCER_TYPE",
      "BM_THERAPEUTIC_CONTEXT",
      "BM_EVIDENCE_LEVEL",
      "BM_EVIDENCE_DESCRIPTION"),
    only_colnames = FALSE,
    quiet = TRUE
  )


  theme <- create_variant_table_theme(
    color_palette = color_palette,
    header_color = "#2c313c"
  )

  if(variant_category == "snv_indel"){
    main_cols = list(
      SAMPLE_ALTERATION = reactable::colDef(
        name = "Alteration",
        cell = pcgrr::render_alteration_cell(
          rctbl_recs$main),
        minWidth = 140,       # icon + monospace text like "BCR::ABL1 fusion"
        maxWidth = 200
      ),
      BM_SOURCES = reactable::colDef(
        name = "Source",
        cell = pcgrr::render_source_logos(),
        align = "center",
        minWidth = 70
      ),
      CLINICAL_SIGNIFICANCE = reactable::colDef(
        name = "Biomarker relevance",
        minWidth = 200,
        cell = rt_cell_bm_significance(color_palette)
      ),
      GENOTYPE = reactable::colDef(
        name = "Genotype",
        minWidth = 120,
        align = "center",
        cell = rt_cell_genotype(color_palette)
      ),
      CLASSIFICATION = reactable::colDef(show = FALSE),
      CLASSIFICATION_RANK = reactable::colDef(
        name = "Clinical significance",
        minWidth = 140,
        align = "center",
        defaultSortOrder = "desc",
        cell = rt_cell_classification_rank(color_palette)
      ),
      ASSERTION_AUTHORITY = reactable::colDef(show = FALSE),
      VAR_ID = reactable::colDef(show = FALSE),
      ENTREZGENE = reactable::colDef(show = FALSE),
      BM_TOP_MAPPING_CONFIDENCE = reactable::colDef(show = FALSE),
      SEARCH_INDEX = reactable::colDef(show = FALSE, searchable = TRUE)
    )
  }

  result_table <- reactable::reactable(
    rctbl_recs$main,
    columns = main_cols,
    details = function(index) {
      eg <- rctbl_recs$main$ENTREZGENE[index]
      vid <- rctbl_recs$main$VAR_ID[index]
      nested <- rctbl_recs$nested[
        rctbl_recs$nested$ENTREZGENE == eg &
          rctbl_recs$nested$VAR_ID == vid,
      ]
      if (nrow(nested) == 0) return(NULL)
      display_cols <- setdiff(
        names(nested),
        c("VAR_ID",
          "ENTREZGENE")
      )

      htmltools::div(
        style = "padding: 4px 40px 12px 40px; background: #f9f9f9;",
        reactable::reactable(
          nested[, display_cols],
          columns = c(
            list(
              BM_MOLECULAR_PROFILE = reactable::colDef(
                name = "Molecular Profile",
                html = TRUE,
                minWidth = 140
              ),
              BM_REFERENCE = reactable::colDef(
                name = "Reference",
                html = TRUE,
                minWidth = 180
              ),
              BM_SOURCE_DB = reactable::colDef(
                name = "Source",
                maxWidth = 100
              ),
              BM_CANCER_TYPE = reactable::colDef(
                name = "Cancer Type",
                minWidth = 130
              ),
              BM_CLINICAL_SIGNIFICANCE = reactable::colDef(
                name = "Clinical Significance",
                minWidth = 170,
              ),
              BM_THERAPEUTIC_CONTEXT = reactable::colDef(
                name = "Therapy Match",
                minWidth = 120
              ),
              BM_EVIDENCE_LEVEL = reactable::colDef(
                name = "Evidence level",
                cell = pcgrr::render_evidence_level_cell(color_palette),
                align = "center",
                minWidth = 100
              ),
              BM_EVIDENCE_DESCRIPTION = reactable::colDef(
                name = "Evidence Description",
                cell = pcgrr::render_evidence_desc_cell(),
                minWidth = 350
              )
            )
          ),
          outlined = TRUE,
          compact = TRUE,
          highlight = TRUE,
          wrap = TRUE,
          defaultPageSize = 5,
          theme = reactable::reactableTheme(
            backgroundColor = "#f9f9f9"
          )
        )
      )
    },
    searchable = TRUE,
    striped = TRUE,
    highlight = TRUE,
    compact = TRUE,
    filterable = TRUE,
    defaultPageSize = 5,
    theme = theme
  )


  return(result_table)
}



#' Impact level of a CPIC phenotype (for color coding)
#'
#' @param phenotype CPIC phenotype, e.g. 'Poor Metabolizer'. Multi-gene
#' phenotypes ('TPMT: Normal Metabolizer; NUDT15: Poor Metabolizer') are
#' assigned the highest impact level among the genes
#'
#' @return integer, 1 (normal) to 5 (highest impact), or 0 when the
#' phenotype is indeterminate/not determined
#' @keywords internal
pgx_phenotype_level <- function(phenotype) {
  vapply(phenotype, function(p) {
    if (is.na(p) || p == "") return(0L)
    parts <- trimws(sub("^[A-Z0-9]+:", "", trimws(strsplit(p, ";")[[1]])))
    parts <- trimws(sub("\\(activity score.*\\)$", "", parts))
    levels <- vapply(parts, function(x) {
      x <- tolower(x)
      if (grepl("cnsha", x)) return(5L)
      if (grepl("^poor|^deficient", x)) return(4L)
      if (grepl("^intermediate|^variable", x)) return(3L)
      if (grepl("^possible intermediate", x)) return(2L)
      if (grepl("^normal", x)) return(1L)
      0L
    }, integer(1))
    max(levels)
  }, integer(1), USE.NAMES = FALSE)
}

#' Color of a CPIC phenotype - blues, darker = more impact (HTML and PDF)
#'
#' @param phenotype CPIC phenotype(s)
#'
#' @return hex color(s); grey for indeterminate/not determined
#' @keywords internal
pgx_phenotype_color <- function(phenotype) {
  pgx_phenotype_colors <- c(
    "#9e9e9e",  # indeterminate / not determined
    "#6baed6",  # normal
    "#4292c6",  # possible intermediate
    "#2171b5",  # intermediate / variable
    "#08519c",  # poor / deficient
    "#08306b")  # deficient with CNSHA
  pgx_phenotype_colors[pgx_phenotype_level(phenotype) + 1]
}

#' Cell renderer for CPIC phenotypes - blue badge, darker = more impact
#'
#' @keywords internal
rt_cell_pgx_phenotype <- function(value) {
  if (is.na(value) || value == "") return("-")
  htmltools::span(
    style = list(
      background = pgx_phenotype_color(value), color = "#ffffff",
      padding = "3px 8px", borderRadius = "3px", fontWeight = "bold",
      display = "inline-block", fontSize = "0.92em", lineHeight = "1.4"),
    value)
}

#' Cell renderer for the detected variants of a gene ('; '-separated),
#' one variant per line (aligned with rt_cell_pgx_genotype)
#'
#' @keywords internal
rt_cell_pgx_variants <- function(value) {
  if (is.na(value) || value == "") return("-")
  htmltools::div(unname(lapply(strsplit(value, "; ")[[1]], function(v)
    htmltools::div(style = list(padding = "3px 0", lineHeight = "1.4"), v))))
}

#' Cell renderer factory for the genotypes of the detected variants of a
#' gene ('; '-separated: heterozygous, homozygous, hemizygous) - pills
#' styled as the genotypes of the variant classification table
#'
#' @param color_palette CPSR color palette object
#' @keywords internal
rt_cell_pgx_genotype <- function(color_palette) {
  function(value) {
    if (is.na(value) || value == "") return("-")
    pills <- lapply(strsplit(value, "; ")[[1]], function(gt) {
      level <- dplyr::case_when(
        gt == "heterozygous" ~ "het",
        gt %in% c("homozygous", "hemizygous") ~ "hom_alt",
        TRUE ~ "undefined")
      idx <- match(level, color_palette$genotypes$levels)
      htmltools::div(
        style = list(padding = "1px 0"),
        htmltools::span(
          style = list(
            background = color_palette$genotypes$bgcolor_values[idx],
            color = color_palette$genotypes$color_values[idx],
            padding = "2px 8px", borderRadius = "3px",
            fontWeight = "bold", display = "inline-block",
            fontSize = "0.92em", lineHeight = "1.4"),
          gt))
    })
    htmltools::div(unname(pills))
  }
}

#' Convert an HTML link string ("<a href='...'>text</a>") into an
#' htmltools tag (raw HTML is not rendered in reactable row details)
#'
#' @keywords internal
html_link_to_tag <- function(value) {
  if (is.na(value) || !grepl("<a ", value, fixed = TRUE)) return(value)
  links <- stringr::str_match_all(
    value, "<a [^>]*href=['\"]([^'\"]+)['\"][^>]*>([^<]*)</a>")[[1]]
  if (NROW(links) == 0) return(value)
  tags <- lapply(seq_len(NROW(links)), function(i)
    htmltools::tags$a(href = links[i, 2], target = "_blank", links[i, 3]))
  ## multiple links separated by ', '
  unname(unlist(lapply(seq_along(tags), function(i)
    if (i == 1) list(tags[[i]]) else list(", ", tags[[i]])),
    recursive = FALSE))
}

#' Function that creates a reactable of CPIC pharmacogenomic phenotypes
#' or prescribing recommendations, for display in the germline report
#'
#' The phenotype table has one row per gene (detected variants, CPIC
#' diplotype, phenotype); each row can be expanded to show the
#' variant-level details (CPIC allele/function, consequence, ClinVar,
#' gnomAD) and the CPIC consultation text. The recommendation table has one
#' row per drug and phenotype; each row can be expanded to show the
#' implications, CPIC comments and guideline. Phenotypes are shown as
#' blue badges, darker for higher impact.
#'
#' @param data data frame with gene phenotypes (assign_pgx_phenotypes)
#' or drug recommendations (assign_pgx_recommendations)
#' @param type 'phenotype' or 'recommendation'
#' @param color_palette CPSR color_palette object
#' @param variants data frame with pharmacogenomic variants (display
#' version of retrieve_pgx_calls output), used for row details of the
#' phenotype table
#'
#' @export
create_pgx_reactable <- function(
    data = NULL,
    type = "phenotype",
    color_palette = NULL,
    variants = NULL) {

  assertthat::assert_that(!is.null(data), msg = "data is NULL")
  assertthat::assert_that(type %in% c("phenotype", "recommendation"))
  assertthat::assert_that(!is.null(color_palette), msg = "color_palette is NULL")

  if (type == "phenotype") {
    col_defs <- list(
      SYMBOL = reactable::colDef(
        name = "Gene", minWidth = 70, sticky = "left",
        style = list(fontWeight = "bold")),
      PGX_VARIANTS = reactable::colDef(
        name = "Alteration", minWidth = 110, html = FALSE,
        cell = rt_cell_pgx_variants),
      PGX_DIPLOTYPE = reactable::colDef(
        name = "Diplotype", minWidth = 150),
      PGX_GENOTYPES = reactable::colDef(
        name = "Genotype", minWidth = 115, html = FALSE,
        cell = rt_cell_pgx_genotype(color_palette)),
      PGX_PHENOTYPE = reactable::colDef(
        name = "Phenotype", minWidth = 150, html = FALSE,
        cell = rt_cell_pgx_phenotype),
      PGX_ACTIVITY_SCORE = reactable::colDef(
        name = "Activity score", minWidth = 95, align = "center",
        cell = function(value) {
          if (is.na(value) || value %in% c("n/a", "")) "-" else value
        }),
      ## interpretation notes are shown in the row details (keeps the
      ## table within the page width)
      PGX_PHENOTYPE_NOTE = reactable::colDef(show = FALSE)
    )
    details <- pgx_phenotype_row_details(data, variants)
  } else {
    data <- data |>
      dplyr::mutate(
        DRUG_NAME = dplyr::if_else(
          !is.na(.data$GUIDELINE_URL),
          paste0("<a href='", .data$GUIDELINE_URL, "' target='_blank'>",
                 .data$DRUG_NAME, "</a>"),
          .data$DRUG_NAME),
        PHENOTYPE_DISPLAY = dplyr::coalesce(
          format_pgx_phenotype(.data$PHENOTYPE, .data$ACTIVITY_SCORE),
          stringr::str_replace_all(.data$LOOKUP_KEY, c(":" = ": ", ";" = "; "))),
        RECOMMENDATION = dplyr::coalesce(
          .data$RECOMMENDATION,
          paste0("<i>", .data$PGX_NOTE, "</i>")),
        GENES = stringr::str_replace_all(.data$GENES, "\\|", ", "),
        IMPLICATIONS = format_pgx_gene_text(.data$IMPLICATIONS))
    col_defs <- list(
      DRUG_NAME = reactable::colDef(
        name = "Drug", minWidth = 100, sticky = "left", html = TRUE,
        style = list(fontWeight = "bold")),
      PHENOTYPE_DISPLAY = reactable::colDef(
        name = "Phenotype", minWidth = 170, html = FALSE,
        cell = rt_cell_pgx_phenotype),
      RECOMMENDATION = reactable::colDef(
        name = "CPIC recommendation", minWidth = 225),
      CLASSIFICATION = reactable::colDef(
        name = "Strength", minWidth = 105, align = "center",
        cell = function(value) {
          if (is.na(value) || value %in% c("", ".")) "-" else value
        })
    )
    details <- pgx_recommendation_row_details(data)
  }

  col_defs <- col_defs[names(col_defs) %in% colnames(data)]
  table_data <- dplyr::select(data, dplyr::all_of(names(col_defs)))

  reactable::reactable(
    table_data,
    columns = col_defs,
    defaultColDef = reactable::colDef(html = TRUE),
    ## details rendered from R (htmltools tags) - must not inherit
    ## html = TRUE from defaultColDef (shown as '[object Object]')
    details = reactable::colDef(details = details, html = FALSE, width = 45),
    searchable = FALSE,
    highlight = TRUE,
    striped = TRUE,
    compact = TRUE,
    wrap = TRUE,
    defaultPageSize = 10,
    theme = create_variant_table_theme(
      color_palette = color_palette,
      header_color = "#2c313c")
  )
}

#' Row details (variant level) for the CPIC phenotype table
#'
#' @param phenotypes data frame with gene phenotypes
#' @param variants data frame with pharmacogenomic variants
#'
#' @return function(index) for reactable details
#' @keywords internal
pgx_phenotype_row_details <- function(phenotypes, variants = NULL) {

  ## genotype and activity (score) are shown in the main row
  detail_fields <- c(
    "ALTERATION" = "Alteration",
    "PGX_ALLELE" = "CPIC allele(s)",
    "PGX_ALLELE_FUNCTION" = "CPIC allele function",
    "CONSEQUENCE" = "Consequence",
    "CLINVAR_CLASSIFICATION" = "ClinVar",
    "DBSNP_RSID" = "dbSNP",
    "gnomADe_AF" = "gnomAD AF (exomes)",
    "gnomADg_AF" = "gnomAD AF (genomes)")

  function(index) {
    gene <- phenotypes$SYMBOL[index]
    blocks <- list()
    if (NROW(variants) > 0 && "SYMBOL" %in% colnames(variants)) {
      gene_vars <- variants[variants$SYMBOL == gene, , drop = FALSE]
      fields <- detail_fields[names(detail_fields) %in% colnames(gene_vars)]
      ## alteration only needed to tell multiple variants apart
      if (NROW(gene_vars) == 1) {
        fields <- fields[names(fields) != "ALTERATION"]
      }
      if (NROW(gene_vars) > 0 && length(fields) > 0) {
        ## unnamed lists - named arguments become HTML attributes
        header <- htmltools::tags$tr(unname(lapply(fields, function(f)
          htmltools::tags$th(
            f, style = "text-align:left; padding:3px 12px 3px 0;"))))
        rows <- unname(lapply(seq_len(NROW(gene_vars)), function(i) {
          htmltools::tags$tr(unname(lapply(names(fields), function(f) {
            ## single value - gene_vars may be a tibble (gene_vars[i, f]
            ## would then be a 1x1 tibble)
            value <- as.character(gene_vars[[f]][i])
            htmltools::tags$td(
              if (is.na(value) || value == "") "-" else html_link_to_tag(value),
              style = "padding:3px 12px 3px 0; vertical-align:top;")
          })))
        }))
        blocks[[length(blocks) + 1]] <- htmltools::tags$table(
          style = "font-size:0.92em; margin-bottom:8px;",
          htmltools::tags$thead(header), htmltools::tags$tbody(rows))
      }
    }
    if ("PGX_PHENOTYPE_NOTE" %in% colnames(phenotypes) &&
        !is.na(phenotypes$PGX_PHENOTYPE_NOTE[index])) {
      blocks[[length(blocks) + 1]] <- htmltools::div(
        style = "font-size:0.92em; margin-bottom:8px;",
        htmltools::tags$b("Note: "),
        htmltools::tags$i(phenotypes$PGX_PHENOTYPE_NOTE[index]))
    }
    if ("CONSULTATION_TEXT" %in% colnames(phenotypes) &&
        !is.na(phenotypes$CONSULTATION_TEXT[index])) {
      blocks[[length(blocks) + 1]] <- htmltools::div(
        style = "font-size:0.92em; color:#444;",
        htmltools::tags$b("CPIC consultation text: "),
        phenotypes$CONSULTATION_TEXT[index])
    }
    if (length(blocks) == 0) return(NULL)
    htmltools::div(style = "padding:10px 16px;", blocks)
  }
}

#' Row details for the CPIC prescribing recommendation table
#'
#' @param recommendations data frame with drug recommendations
#'
#' @return function(index) for reactable details
#' @keywords internal
pgx_recommendation_row_details <- function(recommendations) {

  detail_fields <- c(
    "IMPLICATIONS" = "Implications",
    "COMMENTS" = "CPIC comments",
    "PGX_NOTE" = "Note",
    "CPIC_LEVEL" = "CPIC level (gene-drug pair)",
    "GUIDELINE_NAME" = "CPIC guideline")

  function(index) {
    items <- list()
    for (f in names(detail_fields)) {
      if (!f %in% colnames(recommendations)) next
      value <- recommendations[[f]][index]
      if (is.na(value) || value %in% c("", "n/a")) next
      if (f == "GUIDELINE_NAME" && "GUIDELINE_URL" %in% colnames(recommendations) &&
          !is.na(recommendations$GUIDELINE_URL[index])) {
        value <- htmltools::tags$a(
          href = recommendations$GUIDELINE_URL[index], target = "_blank", value)
      }
      if (f == "CPIC_LEVEL") {
        value <- stringr::str_replace_all(value, c(":" = ": ", ";" = "; "))
      }
      items[[length(items) + 1]] <- htmltools::div(
        style = "margin-bottom:6px;",
        htmltools::tags$b(paste0(detail_fields[[f]], ": ")), value)
    }
    if (length(items) == 0) return(NULL)
    htmltools::div(style = "padding:10px 16px; font-size:0.92em;", items)
  }
}
