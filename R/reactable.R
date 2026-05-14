# ── Cell renderer factories ──────────────────────────────────────────────────

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
#' Maps a numeric rank (1=Benign … 5=Pathogenic) to a labelled
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

# ─────────────────────────────────────────────────────────────────────────────

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
    "ASSERTION_AUTHORITY", "CLASSIFICATION", "CLINVAR_GOLD_STARS",
    "CPSR_PATHOGENICITY_SCORE", "SYMBOL", "GENOMIC_CHANGE", "ALTERATION",
    "GERP_SCORE", "HGVSc", "HGVSp", "HGVSc_RefSeq", "CLINVAR_PHENOTYPE",
    "CLINVAR", "PROTEIN_DOMAIN", "GENENAME", "CODING_STATUS"
  )

  if (is.null(all_variants) ||
      !all(required_cols %in% colnames(all_variants))) {
    return(
      data.frame(matrix(ncol = length(required_cols), nrow = 0)) |>
        setNames(required_cols)
    )
  }

  all_variants <- all_variants |>
    dplyr::mutate(
      ACMG_CODE = dplyr::if_else(
        is.na(.data$ACMG_CODE) | .data$ACMG_CODE == "",
        "—",
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
          GENOMIC_CHANGE,
          .data$ACMG_CODE,
          sep = " "
        ),
        "NA ",""
      )
    )


  # Split by source and arrange by quality metrics
  clinvar_data <- all_variants |>
    dplyr::filter(ASSERTION_AUTHORITY == "ClinVar") |>
    dplyr::arrange(
      dplyr::desc(CLASSIFICATION),
      dplyr::desc(CLINVAR_GOLD_STARS),
      SYMBOL
    ) |>
    utils::head(max_rows)

  cpsr_data <- all_variants |>
    dplyr::filter(ASSERTION_AUTHORITY == "CPSR") |>
    dplyr::arrange(
      dplyr::desc(CLASSIFICATION),
      dplyr::desc(CPSR_PATHOGENICITY_SCORE),
      SYMBOL
    ) |>
    utils::head(max_rows)

  # Combine and create unique key
  unified_data <-
    dplyr::bind_rows(
      clinvar_data,
      cpsr_data) |>
    dplyr::arrange(
      dplyr::desc(.data$CLASSIFICATION),
      SYMBOL
    ) |>
    dplyr::mutate(
      VARKEY_CLASSIFICATION = paste(
        SYMBOL,
        GENOMIC_CHANGE,
        ALTERATION,
        ASSERTION_AUTHORITY, sep = "|"),
      .after = SYMBOL
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
    color_palette = NULL) {

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
      minWidth = 120,
      sticky = "left",
      cell = function(value) {
        if (is.na(value)) return("-")
        if (value == 0) return("0")
        formatC(value, format = "e", digits = 2)
      }
    ),
    CONSEQUENCE = reactable::colDef(
      name = "Consequence",
      minWidth = 120
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
      minWidth = 90,
      align = "center",
      cell = rt_cell_genotype(color_palette)
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
      minWidth = 90,
      align = "center",
      cell = rt_cell_genotype(color_palette)
    ),
    CONSEQUENCE = reactable::colDef(
      name = "Consequence",
      minWidth = 130
    )
  )

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
      paste0("No tier ", paste(tier, collapse = "/"), " variants found."))
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
          VAR_ID, ENTREZGENE
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
            c(unique(na.omit(.data$BM_EVIDENCE_DESCRIPTION)),
              unique(na.omit(.data$BM_CANCER_TYPE)),
              unique(na.omit(.data$BM_THERAPEUTIC_CONTEXT))),
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
#' Combines tier 1 and tier 2 records in one table.
#' Header uses tier 1 color; THERAPY_MATCH cell background
#' reflects the row's tier (1 or 2).
#' @param rctbl_recs List with $main and $nested data frames.
#'   $main must contain ACTIONABILITY_TIER with values 1 and 2.
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
        minWidth = 90,
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
                cell = pcgrr::render_evidence_level_cell(),
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


