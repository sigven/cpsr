#' Function that assigns final pathogenicity classification (B, LB, VUS, P, LP)
#' based on accumulated scores from different ACMG criteria and pre-defined
#' cutoffs (calibrated against ClinVar)
#'
#' @param var_calls data frame with variant calls in predisposition genes
#'
#' @return var_calls data frame with pathogenicity classification appended
#'
#' @export
assign_classification <- function(var_calls) {


  evidence_codes <- cpsr::acmg[["evidence_codes"]]

  pcgrr::log4r_info(paste0(
    "Assigning five-tier classifications (P, LP, VUS, LB, B) based on ",
    "aggregated ACMG points"))

  path_cols <- c(
    "CPSR_CLASSIFICATION",
    "CPSR_CLASSIFICATION_DOC",
    "CPSR_CLASSIFICATION_CODE",
    "cpsr_score_pathogenic",
    "cpsr_score_benign"
  )
  var_calls <- var_calls[, !(colnames(var_calls) %in% path_cols)]

  var_calls$CPSR_CLASSIFICATION <- "VUS"
  var_calls$CPSR_CLASSIFICATION_DOC <- ""
  var_calls$CPSR_CLASSIFICATION_CODE <- ""
  var_calls$cpsr_score_pathogenic <- 0
  var_calls$cpsr_score_benign <- 0

  i <- 1
  while (i <= nrow(evidence_codes)) {
    category <- evidence_codes[i, ]$category
    pole <- evidence_codes[i, ]$pathogenicity_pole
    description <- evidence_codes[i, ]$description
    cpsr_evidence_code <- evidence_codes[i, ]$cpsr_evidence_code
    score <- evidence_codes[i, ]$path_score
    if (cpsr_evidence_code %in% colnames(var_calls)) {
      var_calls <- var_calls |>
        dplyr::mutate(
          cpsr_score_benign = .data$cpsr_score_benign +
            dplyr::if_else(
              pole == "B" & !!rlang::sym(cpsr_evidence_code) == T,
              score, 0
            )
        ) |>
        dplyr::mutate(
          cpsr_score_pathogenic = .data$cpsr_score_pathogenic +
            dplyr::if_else(
              pole == "P" & !!rlang::sym(cpsr_evidence_code) == T,
              score, 0
            )
        ) |>
        dplyr::mutate(
          CPSR_CLASSIFICATION_DOC =
            paste0(
              .data$CPSR_CLASSIFICATION_DOC,
              dplyr::if_else(
                !!rlang::sym(cpsr_evidence_code) == T,
                paste0("- <i>", cpsr_evidence_code,
                       "</i>: ", description, " (<b>", score, "</b>)"), ""
              ),
              sep = "<br>"
            )
        ) |>
        dplyr::mutate(
          CPSR_CLASSIFICATION_CODE =
            paste0(
              .data$CPSR_CLASSIFICATION_CODE,
              dplyr::if_else(
                !!rlang::sym(cpsr_evidence_code) == T,
                cpsr_evidence_code, ""
              ),
              sep = "|"
            )
        )
    }
    i <- i + 1
  }

  p_lower_limit <- cpsr::acmg[['score_thresholds']][['p_lower']]
  lp_upper_limit <- cpsr::acmg[['score_thresholds']][['lp_upper']]
  lp_lower_limit <- cpsr::acmg[['score_thresholds']][['lp_lower']]
  vus_upper_limit <- cpsr::acmg[['score_thresholds']][['vus_upper']]
  vus_lower_limit <- cpsr::acmg[['score_thresholds']][['vus_lower']]
  lb_upper_limit <- cpsr::acmg[['score_thresholds']][['lb_upper']]
  lb_lower_limit <- cpsr::acmg[['score_thresholds']][['lb_lower']]
  b_upper_limit <- cpsr::acmg[['score_thresholds']][['b_upper']]

  var_calls <- var_calls |>
    dplyr::mutate(
      CPSR_CLASSIFICATION_CODE =
        stringr::str_replace_all(
          stringr::str_replace_all(
            .data$CPSR_CLASSIFICATION_CODE,
            "(\\|{2,})", "|"
          ),
          "(^\\|)|(\\|$)", ""
        )
    ) |>
    dplyr::mutate(
      CPSR_CLASSIFICATION_DOC =
        stringr::str_replace_all(
          stringr::str_replace_all(
            .data$CPSR_CLASSIFICATION_DOC,
            "(<br>){2,}", "<br>"
          ), "(^(<br>))|((<br>)$)", ""
        )
    ) |>
    ## Adjust scores in cases where critera are acting as a
    ## prerequisite for other criteria
    dplyr::mutate(
      cpsr_score_pathogenic =
        dplyr::if_else(
          stringr::str_detect(
            .data$CPSR_CLASSIFICATION_CODE, "ACMG_PVS") &
            stringr::str_detect(
              .data$CPSR_CLASSIFICATION_CODE, "ACMG_PM2_2"),
          .data$cpsr_score_pathogenic - 1,
          .data$cpsr_score_pathogenic
        )
    ) |>
    dplyr::mutate(
      cpsr_score_pathogenic =
        dplyr::if_else(
          stringr::str_detect(
            .data$CPSR_CLASSIFICATION_CODE, "ACMG_PVS") &
            stringr::str_detect(
              .data$CPSR_CLASSIFICATION_CODE, "ACMG_PM2_1"),
          .data$cpsr_score_pathogenic - 0.5,
          .data$cpsr_score_pathogenic
        )
    ) |>
    dplyr::mutate(
      cpsr_score_pathogenic =
        dplyr::if_else(
          stringr::str_detect(
            .data$CPSR_CLASSIFICATION_CODE, "ACMG_PVS1_10") &
            stringr::str_detect(
              .data$CPSR_CLASSIFICATION_CODE, "ACMG_PP3"),
          .data$cpsr_score_pathogenic - 0.5,
          .data$cpsr_score_pathogenic
        )
    ) |>
    ## Add scores accumulated with benign criteria and pathogenic criteria
    dplyr::mutate(
      CPSR_PATHOGENICITY_SCORE =
        dplyr::if_else(
          .data$cpsr_score_benign == 0,
          .data$cpsr_score_pathogenic,
          .data$cpsr_score_benign
        )
    ) |>
    dplyr::mutate(
      CPSR_PATHOGENICITY_SCORE =
        dplyr::if_else(
          .data$cpsr_score_benign < 0 &
            .data$cpsr_score_pathogenic > 0,
          .data$cpsr_score_benign + .data$cpsr_score_pathogenic,
          .data$CPSR_PATHOGENICITY_SCORE
        )
    ) |>
    dplyr::mutate(
      CPSR_CLASSIFICATION =
        dplyr::case_when(
          .data$CPSR_PATHOGENICITY_SCORE <= lb_upper_limit &
            .data$CPSR_PATHOGENICITY_SCORE >= lb_lower_limit ~ "Likely_Benign",
          .data$CPSR_PATHOGENICITY_SCORE <= b_upper_limit ~ "Benign",
          .data$CPSR_PATHOGENICITY_SCORE <= vus_upper_limit &
            .data$CPSR_PATHOGENICITY_SCORE >= vus_lower_limit ~ "VUS",
          .data$CPSR_PATHOGENICITY_SCORE >= p_lower_limit ~ "Pathogenic",
          .data$CPSR_PATHOGENICITY_SCORE >= lp_lower_limit &
            .data$CPSR_PATHOGENICITY_SCORE <= lp_upper_limit ~ "Likely_Pathogenic",
          TRUE ~ as.character("VUS")
        )
    ) |>
    dplyr::select(-c(.data$cpsr_score_benign,
                     .data$cpsr_score_pathogenic))

  return(var_calls)
}


#' Function that assigns variant pathogenicity evidence based on ACMG guidelines
#'
#' @param var_calls sample calls with dbnsfp annotations
#' @param settings cpsr settings object
#' @param ref_data pcgr data object
#'
#' @return calls
#'
#' @export
assign_acmg_evidence <- function(var_calls, settings, ref_data) {

  invisible(assertthat::assert_that(!is.null(var_calls)))
  invisible(assertthat::assert_that(!is.null(settings)))
  invisible(assertthat::assert_that(!is.null(settings$conf)))
  invisible(assertthat::assert_that(
    !is.null(settings$conf$variant_classification)))
  invisible(assertthat::assert_that(!is.null(ref_data)))
  invisible(assertthat::assert_that(is.data.frame(var_calls)))

  pcgrr::log4r_info(
    "Assigning variant classification codes according to refined ACMG criteria")
  classification_settings <-
    settings$conf$variant_classification

  acmg_ev_codes <-
    c(
      "ACMG_PM1_MOD",
      "ACMG_PM1_SUPP",
      "ACMG_PM2_SUPP",
      "ACMG_BA1",
      "ACMG_BS1",
      "ACMG_BS1_SUPP",
      "ACMG_PVS1",
      "ACMG_PVS1_STR",
      "ACMG_PVS1_MOD",
      "ACMG_PS1",
      "ACMG_PP3",
      "ACMG_PS1",
      "ACMG_PM5",
      "ACMG_PP3",
      "ACMG_BP4",
      "ACMG_BP7",
    )

  var_calls <- var_calls[, !(colnames(var_calls) %in% acmg_ev_codes)]

  if("CONSEQUENCE" %in% colnames(var_calls) &
     "HGVSP" %in% colnames(var_calls) &
     "HGVSp" %in% colnames(var_calls) &
     "ENTREZGENE" %in% colnames(var_calls)){
    var_calls <- var_calls |>
      dplyr::mutate(
        CODON = dplyr::if_else(
          !is.na(.data$CONSEQUENCE) &
            stringr::str_detect(
              .data$CONSEQUENCE,
              "^missense_variant"
            ) &
            !is.na(.data$ENTREZGENE) &
            !is.na(.data$HGVSP),
          stringr::str_match(
            .data$HGVSP,
            "p\\.[A-Z]{1}[0-9]{1,}"
          )[,1],
          as.character(NA)
        )
      ) |>
      tidyr::separate(
        "HGVSp",c("tmpENSP","HGVSp_long"),sep=":", remove = F) |>
      dplyr::rename("HGVSP_query" = "HGVSp_long") |>
      dplyr::select(-c("tmpENSP"))
  }

  known_sites <- cpsr::known_path_benign_sites(
    ref_data = ref_data
  )

  var_calls3 <- var_calls |>
    cpsr::assign_BA1_evidence() |>
    cpsr::assign_BS1_evidence() |>
    cpsr::assign_PM2_evidence() |>
    cpsr::assign_PM1_evidence() |>
    cpsr::assign_PVS1_evidence() |>
    cpsr::assign_PP3_evidence() |>
    cpsr::assign_BP4_evidence() |>
    cpsr::assign_PM5_evidence(
      pathogenic_codons = known_sites[['pathogenic_codon']]) |>
    cpsr::assign_PS1_evidence(
      pathogenic_codons = known_sites[['pathogenic_codon']],
      pathogenic_nucleotides =
        known_sites[['pathogenic_nucleotide']]
    )

  # Assign logical ACMG evidence indicators
  # # TODO - BA1 -  exceptions for high population germline frequency
  #  (gnomAD) - HFE/SERPINA1

  ## Assign logical ACMG evidence indicator
  # PM4 - Protein length changes (in non-repetitive regions) due to
  # inframe indels or nonstop variant of genes that harbor variants with
  # a dominant mode of inheritance
  #
  if ("RMSK_HIT" %in% colnames(var_calls)) {
    var_calls <- var_calls |>
      dplyr::mutate(
        ACMG_PM4 =
          dplyr::if_else(
            stringr::str_detect(
              .data$CONSEQUENCE,
              "stop_lost|inframe_deletion|inframe_insertion") &
              is.na(.data$RMSK_HIT) &
              .data$CPG_MOI == "AD",
            TRUE, FALSE, FALSE
          )
      )
  }

  ## Assign logical ACMG evidence indicator
  # ACMG_PP2 - Missense variant in a gene that has a relatively low rate
  # of benign missense variation and where missense variants are a
  # common mechanism of disease
  var_calls <- var_calls |>
    dplyr::mutate(
      ACMG_PP2 =
        dplyr::if_else(
          (is.na(.data$BENIGN_MISSENSE_FRAC) |
             .data$BENIGN_MISSENSE_FRAC <= 0.1) &
            (is.na(.data$PATH_TRUNC_FRAC) |
               .data$PATH_TRUNC_FRAC < 0.5) &
            stringr::str_detect(
              .data$CONSEQUENCE, "^missense_variant"),
          TRUE, FALSE, FALSE
        )
    )

  ## Assign a logical ACMG evidence indicator
  # ACMG_BP1 - Missense variant in a gene for which primarily truncating
  # variants (> 90%, as given in Maxwell et al.) are known to cause disease
  var_calls <- var_calls |>
    dplyr::mutate(
      ACMG_BP1 =
        dplyr::if_else(
          .data$PATH_TRUNC_FRAC > 0.90 &
            stringr::str_detect(
              .data$CONSEQUENCE, "^missense_variant"),
          TRUE, FALSE, FALSE
        )
    )


  # var_calls <- var_calls |>
  #   pcgrr::remove_cols_from_df(
  #     cnames = c(
  #       "PATHOGENIC_CODON",
  #       "BENIGN_CODON",
  #       "NUC_SITE",
  #       "ESSENTIAL_NUCLEOTIDE",
  #       "PATHOGENIC_PEPTIDE_CHANGE",
  #       "BENIGN_PEPTIDE_CHANGE",
  #       "VAR_ID_PATH_CHANGE",
  #       "VAR_ID_BENIGN_CHANGE",
  #       "VAR_ID_PATH_CODON",
  #       "VAR_ID_BENIGN_CODON",
  #       "CODON",
  #       "gad_af",
  #       "hotspot_region",
  #       "hotspot_entrezgene",
  #       "hotspot_symbol",
  #       "hotspot_codon",
  #       "hotspot_aa",
  #       "hotspot_pvalue"
  #     )
  #   ) |>
  #   dplyr::distinct()

  return(var_calls)
}

#' Function that assigns a summary string of
#' all ACMG evidence codes met for a given variant
#' (e.g. "PVS1|PS1|PM2_supporting|PP3")
#'
#' @param var_calls data frame with variant calls in predisposition genes
#'
#' @return var_calls data frame with ACMG_CONCENSUS column appended
#' @export
#'
assign_acmg_concensus <- function(
    var_calls = NULL){

  if(is.data.frame(var_calls) &
     all(c("ACMG_PVS1",
           "ACMG_PVS1_STR",
           "ACMG_PVS1_MOD",
           "ACMG_PS1",
           "ACMG_PM1",
           "ACMG_PM1_SUPP",
           "ACMG_PM2_SUPP",
           "ACMG_PM5",
           "ACMG_PP3",
           "ACMG_BS1",
           "ACMG_BS1_SUPP",
           "ACMG_BP4",
           "ACMG_BP7") %in%
         colnames(var_calls))){

    var_calls$ACMG_CONCENSUS <- ""

    var_calls <- var_calls |>
      dplyr::mutate(ACMG_CONCENSUS = dplyr::if_else(
        .data$ACMG_PVS1 == TRUE,
        paste0(.data$ACMG_CONCENSUS, "PVS1|"),
        as.character(.data$ACMG_CONCENSUS)
      )) |>
      dplyr::mutate(ACMG_CONCENSUS = dplyr::if_else(
        .data$ACMG_PVS1_STR == TRUE,
        paste0(.data$ACMG_CONCENSUS, "PVS1_strong|"),
        as.character(.data$ACMG_CONCENSUS)
      )) |>
      dplyr::mutate(ACMG_CONCENSUS = dplyr::if_else(
        .data$ACMG_PVS1_MOD == TRUE,
        paste0(.data$ACMG_CONCENSUS, "PVS1_moderate|"),
        as.character(.data$ACMG_CONCENSUS)
      )) |>
      dplyr::mutate(ACMG_CONCENSUS = dplyr::if_else(
        .data$ACMG_PS1 == TRUE,
        paste0(.data$ACMG_CONCENSUS, "PS1|"),
        as.character(.data$ACMG_CONCENSUS)
      )) |>
      dplyr::mutate(ACMG_CONCENSUS = dplyr::if_else(
        .data$ACMG_PM1 == TRUE,
        paste0(.data$ACMG_CONCENSUS, "PM1|"),
        as.character(.data$ACMG_CONCENSUS)
      )) |>
      dplyr::mutate(ACMG_CONCENSUS = dplyr::if_else(
        .data$ACMG_PM1_SUPP == TRUE,
        paste0(.data$ACMG_CONCENSUS, "PM1_supporting|"),
        as.character(.data$ACMG_CONCENSUS)
      )) |>
      dplyr::mutate(ACMG_CONCENSUS = dplyr::if_else(
        .data$ACMG_PM2_SUPP == TRUE,
        paste0(.data$ACMG_CONCENSUS, "PM2_supporting|"),
        as.character(.data$ACMG_CONCENSUS)
      )) |>
      dplyr::mutate(ACMG_CONCENSUS = dplyr::if_else(
        .data$ACMG_PM5 == TRUE,
        paste0(.data$ACMG_CONCENSUS, "PM5|"),
        as.character(.data$ACMG_CONCENSUS)
      )) |>
      dplyr::mutate(ACMG_CONCENSUS = dplyr::if_else(
        .data$ACMG_PP3 == TRUE,
        paste0(.data$ACMG_CONCENSUS, "PP3|"),
        as.character(.data$ACMG_CONCENSUS)
      )) |>
      dplyr::mutate(ACMG_CONCENSUS = dplyr::if_else(
        .data$ACMG_BA1 == TRUE,
        paste0(.data$ACMG_CONCENSUS, "BA1|"),
        as.character(.data$ACMG_CONCENSUS)
      )) |>
      dplyr::mutate(ACMG_CONCENSUS = dplyr::if_else(
        .data$ACMG_BS1 == TRUE,
        paste0(.data$ACMG_CONCENSUS, "BS1|"),
        as.character(.data$ACMG_CONCENSUS)
      )) |>
      dplyr::mutate(ACMG_CONCENSUS = dplyr::if_else(
        .data$ACMG_BS1_SUPP == TRUE,
        paste0(.data$ACMG_CONCENSUS, "BS1_supporting|"),
        as.character(.data$ACMG_CONCENSUS)
      )) |>
      dplyr::mutate(ACMG_CONCENSUS = dplyr::if_else(
        .data$ACMG_BP4 == TRUE,
        paste0(.data$ACMG_CONCENSUS, "BP4|"),
        as.character(.data$ACMG_CONCENSUS)
      )) |>
      dplyr::mutate(ACMG_CONCENSUS = dplyr::if_else(
        .data$ACMG_BP7 == TRUE,
        paste0(.data$ACMG_CONCENSUS, "BP7|"),
        as.character(.data$ACMG_CONCENSUS)
      )) |>
      dplyr::mutate(ACMG_CONCENSUS = stringr::str_replace_all(
        .data$ACMG_CONCENSUS,
        "(\\|$)", ""
      ))

    var_calls <- var_calls |>
      dplyr::mutate(CPSR_CLASSIFICATION_SCORE = 0) |>
      dplyr::mutate(CPSR_CLASSIFICATION_SCORE = dplyr::if_else(
        .data$ACMG_PVS1 == TRUE,
        .data$CPSR_CLASSIFICATION_SCORE + 8,
        .data$CPSR_CLASSIFICATION_SCORE
      )) |>
      dplyr::mutate(CPSR_CLASSIFICATION_SCORE = dplyr::if_else(
        .data$ACMG_PVS1_STR == TRUE,
        .data$CPSR_CLASSIFICATION_SCORE + 4,
        .data$CPSR_CLASSIFICATION_SCORE
      )) |>
      dplyr::mutate(CPSR_CLASSIFICATION_SCORE = dplyr::if_else(
        .data$ACMG_PVS1_MOD == TRUE,
        .data$CPSR_CLASSIFICATION_SCORE + 2,
        .data$CPSR_CLASSIFICATION_SCORE
      )) |>
      dplyr::mutate(CPSR_CLASSIFICATION_SCORE = dplyr::if_else(
        .data$ACMG_PS1 == TRUE,
        .data$CPSR_CLASSIFICATION_SCORE + 4,
        .data$CPSR_CLASSIFICATION_SCORE
      )) |>
      dplyr::mutate(CPSR_CLASSIFICATION_SCORE = dplyr::if_else(
        .data$ACMG_PM1 == TRUE,
        .data$CPSR_CLASSIFICATION_SCORE + 2,
        .data$CPSR_CLASSIFICATION_SCORE
      )) |>
      dplyr::mutate(CPSR_CLASSIFICATION_SCORE = dplyr::if_else(
        .data$ACMG_PM1_SUPP == TRUE,
        .data$CPSR_CLASSIFICATION_SCORE + 1,
        .data$CPSR_CLASSIFICATION_SCORE
      )) |>
      dplyr::mutate(CPSR_CLASSIFICATION_SCORE = dplyr::if_else(
        .data$ACMG_PM2_SUPP == TRUE,
        .data$CPSR_CLASSIFICATION_SCORE + 0,
        .data$CPSR_CLASSIFICATION_SCORE
      )) |>
      dplyr::mutate(CPSR_CLASSIFICATION_SCORE = dplyr::if_else(
        .data$ACMG_PM5 == TRUE,
        .data$CPSR_CLASSIFICATION_SCORE + 2,
        .data$CPSR_CLASSIFICATION_SCORE
      )) |>
      dplyr::mutate(CPSR_CLASSIFICATION_SCORE = dplyr::if_else(
        .data$ACMG_PP3 == TRUE,
        .data$CPSR_CLASSIFICATION_SCORE + 1,
        .data$CPSR_CLASSIFICATION_SCORE
      )) |>
      dplyr::mutate(CPSR_CLASSIFICATION_SCORE = dplyr::if_else(
        .data$ACMG_BA1 == TRUE,
        CPSR_CLASSIFICATION_SCORE - 8,
        CPSR_CLASSIFICATION_SCORE
      )) |>
      dplyr::mutate(CPSR_CLASSIFICATION_SCORE = dplyr::if_else(
        .data$ACMG_BS1 == TRUE,
        CPSR_CLASSIFICATION_SCORE - 4,
        CPSR_CLASSIFICATION_SCORE
      )) |>
      dplyr::mutate(CPSR_CLASSIFICATION_SCORE = dplyr::if_else(
        .data$ACMG_BS1_SUPP == TRUE,
        CPSR_CLASSIFICATION_SCORE - 1,
        CPSR_CLASSIFICATION_SCORE
      )) |>
      dplyr::mutate(CPSR_CLASSIFICATION_SCORE = dplyr::if_else(
        .data$ACMG_BP4 == TRUE,
        CPSR_CLASSIFICATION_SCORE - 1,
        CPSR_CLASSIFICATION_SCORE
      )) |>
      dplyr::mutate(CPSR_CLASSIFICATION_SCORE = dplyr::if_else(
        .data$ACMG_BP7 == TRUE,
        CPSR_CLASSIFICATION_SCORE - 1,
        CPSR_CLASSIFICATION_SCORE
      ))
  }
  return(var_calls)

}


#' Function that combines classifications of novel and
#' pre-classified variants
#'
#' @param var_calls variants in cancer predisposition genes
#' @param conf CPSR configuration object with run settings
#'
#' @export
combine_novel_and_preclassified <-
  function(var_calls,
           conf = NULL){
           #col_format_output = NULL) {

    #if("variant" %in% names(callset)){
    #var_df <-
    #callset$variant

    assertable::assert_colnames(
      var_calls,
      c("CPSR_CLASSIFICATION",
        "CPSR_PATHOGENICITY_SCORE",
        "CPSR_CLASSIFICATION_CODE",
        "CPSR_CLASSIFICATION_DOC",
        "gnomADe_AF",
        "CANCER_PHENOTYPE",
        "CLINVAR_CLASSIFICATION"),
      only_colnames = F,quiet = T
    )

    pcgrr::log4r_info(
      "Combining pre-classified (ClinVar) and novel variants")

    #snv_indel_report <- pcgrr::init_germline_content()

    var_calls <- var_calls |>
      dplyr::mutate(CPSR_CLASSIFICATION_SOURCE = dplyr::case_when(
        !is.na(.data$CLINVAR_CLASSIFICATION) ~ "ClinVar",
        is.na(.data$CLINVAR_CLASSIFICATION) &
          !is.na(CPSR_CLASSIFICATION) ~ "CPSR_ACMG",
        TRUE ~ as.character(NA)
      ))


    ## set FINAL_CLASSIFICATION col
    var_calls <-
      var_calls |>
      dplyr::mutate(
        FINAL_CLASSIFICATION = dplyr::case_when(
          !is.na(.data$CLINVAR_CLASSIFICATION) ~
            as.character(.data$CLINVAR_CLASSIFICATION),
          is.na(.data$CLINVAR_CLASSIFICATION) ~
            as.character(.data$CPSR_CLASSIFICATION),
          TRUE ~ as.character(NA)
        )
      ) |>
      dplyr::mutate(FINAL_CLASSIFICATION = dplyr::if_else(
        FINAL_CLASSIFICATION != "Benign" &
          FINAL_CLASSIFICATION != "Likely_Benign" &
          FINAL_CLASSIFICATION != "VUS" &
          FINAL_CLASSIFICATION != "Likely_Pathogenic" &
          FINAL_CLASSIFICATION != "Pathogenic",
        "VUS",
        as.character(FINAL_CLASSIFICATION)
      )) |>
      dplyr::mutate(FINAL_CLASSIFICATION = dplyr::if_else(
        is.na(FINAL_CLASSIFICATION),"VUS",
        as.character(FINAL_CLASSIFICATION)
      )) |>
      dplyr::mutate(FINAL_CLASSIFICATION = factor(
        .data$FINAL_CLASSIFICATION,
        levels = c("Benign","Likely_Benign","VUS",
                   "Likely_Pathogenic","Pathogenic")
      )) |>
      dplyr::arrange(
        dplyr::desc(
          .data$FINAL_CLASSIFICATION),
        dplyr::desc(.data$CANCER_PHENOTYPE),
        dplyr::desc(.data$CPSR_PATHOGENICITY_SCORE))


    ## If not 'classify_all' is turned on,
    ## remove CPSR classifications for existing
    ## ClinVar classifications
    if (conf[["variant_classification"]][["classify_all"]] == 0) {
      var_calls <-
        var_calls |>
        dplyr::mutate(
          CPSR_CLASSIFICATION =
            dplyr::if_else(
              !is.na(.data$CLINVAR_CLASSIFICATION),
              "",
              as.character(.data$CPSR_CLASSIFICATION)
            )
        ) |>
        dplyr::mutate(
          CPSR_CLASSIFICATION_DOC =
            dplyr::if_else(
              !is.na(.data$CLINVAR_CLASSIFICATION),
              "",
              as.character(.data$CPSR_CLASSIFICATION_DOC)
            )
        ) |>
        dplyr::mutate(
          CPSR_CLASSIFICATION_CODE =
            dplyr::if_else(
              !is.na(.data$CLINVAR_CLASSIFICATION),
              "",
              as.character(.data$CPSR_CLASSIFICATION_CODE)
            )
        ) |>
        dplyr::mutate(
          CPSR_PATHOGENICITY_SCORE =
            dplyr::if_else(
              !is.na(.data$CLINVAR_CLASSIFICATION),
              as.numeric(NA),
              as.numeric(.data$CPSR_PATHOGENICITY_SCORE)
            )
        )
    }

    var_calls$CPSR_CLASSIFICATION_DOC <-
      stringr::str_replace_all(
        var_calls$CPSR_CLASSIFICATION_DOC,
        "<br>-", ","
      )
    var_calls$CPSR_CLASSIFICATION_DOC <-
      stringr::str_replace_all(
        var_calls$CPSR_CLASSIFICATION_DOC,
        "^, ", ""
      )

    for (col in colnames(var_calls)) {
      if (nrow(var_calls[!is.na(
        var_calls[, col]
      ) &
      var_calls[, col] == "", ]) > 0) {
        var_calls[!is.na(
          var_calls[, col]
        ) &
          var_calls[, col] == "", col] <- NA
      }
    }


    return(var_calls)
  }

#' Function that assigns ACMG evidence indicators for
#' BA1
#' @param var_calls variants in cancer predisposition genes
#' @return data.frame with ACMG evidence indicators
#'
#' @export
#'
assign_BA1_evidence <- function(
    var_calls = NULL){

  var_calls$ACMG_BA1 <- FALSE

  ## Assign logical ACMG evidence indicator
  ## BA1 -
  if(is.data.frame(var_calls) &
     "gnomAD_NC_FAF_GRPMAX" %in% colnames(var_calls) &
     "gnomADj_FAF_GRPMAX" %in% colnames(var_calls) &
     "CPG_MOI" %in% colnames(var_calls)){

    var_calls <- var_calls |>
      dplyr::mutate(
        ACMG_BA1 = dplyr::if_else(
          ((.data$gnomADj_FAF_GRPMAX >= 0.001 &
              !(.data$SYMBOL %in% c("BRCA1", "BRCA2", "PALB2","ATM","PTEN",
                                  "MSH2", "MLH1", "MSH6", "PMS2","APC","TP53")) &
              .data$CPG_MOI == "AD") |
             (.data$gnomADj_FAF_GRPMAX >= 0.01 &
                !(.data$SYMBOL %in% c("BRCA1", "BRCA2", "PALB2","ATM","PTEN",
                                      "MSH2", "MLH1", "MSH6", "PMS2","APC","TP53")) &
             (.data$CPG_MOI == "AR" |
                is.na(.data$CPG_MOI)))),
          TRUE, FALSE, FALSE
        )
      ) |>

      ## gene-specific BA1 thresholds for core cancer predisposition genes
      ## - as found in ClinGen VCEPs
      dplyr::mutate(
        ACMG_BA1 = dplyr::case_when(
          .data$SYMBOL == "TP53" & .data$gnomADj_FAF_GRPMAX >= 0.001 ~ TRUE,
          .data$SYMBOL == "BRCA1" & .data$gnomAD_NC_FAF_GRPMAX > 0.001 ~ TRUE,
          .data$SYMBOL == "BRCA2" & .data$gnomAD_NC_FAF_GRPMAX > 0.001 ~ TRUE,
          .data$SYMBOL == "PALB2" & .data$gnomAD_NC_FAF_GRPMAX > 0.001 ~ TRUE,
          .data$SYMBOL == "ATM" & .data$gnomADj_FAF_GRPMAX > 0.005 ~ TRUE,
          .data$SYMBOL == "APC" & .data$gnomAD_NC_FAF_GRPMAX >= 0.001 ~ TRUE,
          .data$SYMBOL == "MLH1" & .data$gnomADj_FAF_GRPMAX >= 0.001 ~ TRUE,
          .data$SYMBOL == "MSH2" & .data$gnomADj_FAF_GRPMAX >= 0.001 ~ TRUE,
          .data$SYMBOL == "PTEN" & .data$gnomADj_FAF_GRPMAX > 0.00056 ~ TRUE,
          .data$SYMBOL == "MSH6" & .data$gnomADj_FAF_GRPMAX >= 0.0022 ~ TRUE,
          .data$SYMBOL == "PMS2" & .data$gnomADj_FAF_GRPMAX >= 0.0028 ~ TRUE,
          TRUE ~ as.logical(.data$ACMG_BA1)
        ))
  }
  return(var_calls)
}

#' Function that assigns ACMG evidence indicators for
#' BS1_ST and BS1_SUP
#'
#' @param var_calls variants in cancer predisposition genes
#' @return data.frame with ACMG evidence indicators
#'
#' @export
assign_BS1_evidence <- function(
    var_calls = NULL){

  var_calls$ACMG_BS1 <- FALSE
  var_calls$ACMG_BS1_SUPP <- FALSE

  ## Assign logical ACMG evidence indicators
  ## BS1_strong
  ## BS1_supporting
  if(is.data.frame(var_calls) &
     "gnomAD_NC_FAF_GRPMAX" %in% colnames(var_calls) &
     "gnomADj_FAF_GRPMAX" %in% colnames(var_calls) &
     "ACMG_BA1" %in% colnames(var_calls) &
     "CPG_MOI" %in% colnames(var_calls)){

    var_calls <- var_calls |>
      dplyr::mutate(
        ACMG_BS1 = dplyr::if_else(
          ((.data$gnomADj_FAF_GRPMAX >= 0.0001 &
              .data$CPG_MOI == "AD" &
              !(.data$SYMBOL %in%
                  c("BRCA1", "BRCA2", "PALB2",
                    "ATM","PTEN", "MSH2", "MLH1",
                    "MSH6", "PMS2","APC", "TP53")) &
              .data$ACMG_BA1 == FALSE) |
             (.data$gnomADj_FAF_GRPMAX >= 0.001 &
             .data$ACMG_BA1 == FALSE &
               !(.data$SYMBOL %in%
                   c("BRCA1", "BRCA2", "PALB2",
                     "ATM","PTEN","MSH2", "MLH1",
                     "MSH6", "PMS2","APC","TP53")) &
             (.data$CPG_MOI == "AR" |
                is.na(.data$CPG_MOI)))),
          TRUE, FALSE, FALSE
        )
      ) |>
      ## gene-specific BS1 thresholds for core cancer predisposition genes
      dplyr::mutate(
        ACMG_BS1 = dplyr::case_when(
          .data$SYMBOL == "TP53" & .data$gnomADj_FAF_GRPMAX > 0.0003 &
            .data$gnomADj_FAF_GRPMAX < 0.001 &
            .data$ACMG_BA1 == FALSE ~ TRUE,
          .data$SYMBOL == "BRCA1" & .data$gnomAD_NC_FAF_GRPMAX > 0.0001 &
            .data$ACMG_BA1 == FALSE ~ TRUE,
          .data$SYMBOL == "BRCA2" & .data$gnomAD_NC_FAF_GRPMAX > 0.0001 &
            .data$ACMG_BA1 == FALSE ~ TRUE,
          .data$SYMBOL == "PALB2" & .data$gnomAD_NC_FAF_GRPMAX > 0.0001 &
            .data$ACMG_BA1 == FALSE ~ TRUE,
          .data$SYMBOL == "ATM" & .data$gnomADj_FAF_GRPMAX > 0.0005 &
            .data$ACMG_BA1 == FALSE ~ TRUE,
          .data$SYMBOL == "APC" & .data$gnomAD_NC_FAF_GRPMAX >= 0.0001 &
            .data$ACMG_BA1 == FALSE ~ TRUE,
          .data$SYMBOL == "MLH1" & .data$gnomADj_FAF_GRPMAX >= 0.0001 &
            .data$gnomADj_FAF_GRPMAX < 0.001 ~ TRUE,
          .data$SYMBOL == "MSH2" & .data$gnomADj_FAF_GRPMAX >= 0.0001 &
            .data$gnomADj_FAF_GRPMAX < 0.001 ~ TRUE,
          .data$SYMBOL == "PTEN" & .data$gnomADj_FAF_GRPMAX > 0.000043 &
            .data$gnomADj_FAF_GRPMAX <= 0.00056 ~ TRUE,
          .data$SYMBOL == "MSH6" & .data$gnomADj_FAF_GRPMAX >= 0.00022 &
            .data$gnomADj_FAF_GRPMAX < 0.0022 ~ TRUE,
          .data$SYMBOL == "PMS2" & .data$gnomADj_FAF_GRPMAX >= 0.0001 &
            .data$gnomADj_FAF_GRPMAX < 0.001 ~ TRUE,
          TRUE ~ as.logical(.data$ACMG_BS1)
        )) |>
      dplyr::mutate(
        ACMG_BS1_SUPP = dplyr::case_when(
          .data$SYMBOL == "BRCA1" &
            .data$gnomAD_NC_FAF_GRPMAX > 0.00002 &
            .data$gnomAD_NC_FAF_GRPMAX <= 0.0001 &
            .data$ACMG_BS1 == FALSE ~ TRUE,
          .data$SYMBOL == "BRCA2" &
            .data$gnomAD_NC_FAF_GRPMAX > 0.00002 &
            .data$gnomAD_NC_FAF_GRPMAX <= 0.0001 &
            .data$ACMG_BS1 == FALSE ~ TRUE,
          TRUE ~ as.logical(.data$ACMG_BS1_SUPP)
        ))
  }
  return(var_calls)


}

#' Function that assigns ACMG evidence indicators for
#' BP4
#'
#' @param var_calls variants in cancer predisposition genes
#' @return data.frame with ACMG evidence indicators
#' @export
#'
assign_BP4_evidence <- function(var_calls = NULL){

  ## Assign logical ACMG evidence indicators
  #
  # ACMG_BP4 - Multiple lines (>=8) of insilico evidence support a benign effect.
  #
  # Computational evidence for deleterious/benign effect is taken from
  # invidual algorithm predictions in dbNSFP: SIFT,Provean,MutationTaster,
  # MutationAssessor,M_CAP,MutPred,FATHMM_XF, DBNSFP_RNN,dbscSNV_RF,
  # dbscSNV_AdaBoost, AlphaMissense, ClinPred, phactBOOST, PrimateAI, REVEL,
  # Tolerated: Among all possible protein variant effect predictions, at
  #    least six algorithms must have made a call,
  #    with at least 8 predicted as tolerated, and at most 2
  #    predicted as damaging (BP4)
  #    - 0 predictions of splice site affected

  dbnsfp_min_majority <- cpsr::acmg[["insilico_pred_min_majority"]]
  dbnsfp_max_minority <- cpsr::acmg[["insilico_pred_max_minority"]]
  dbnsfp_min_called <- dbnsfp_min_majority

  var_calls$ACMG_BP <- NULL
  var_calls_bp4 <- data.frame()

  if(is.data.frame(var_calls) &
     "N_INSILICO_CALLED" %in% colnames(var_calls) &
     "N_INSILICO_DAMAGING" %in% colnames(var_calls) &
     "N_INSILICO_TOLERATED" %in% colnames(var_calls) &
     "N_INSILICO_SPLICING_AFFECTED" %in% colnames(var_calls) &
     "N_INSILICO_SPLICING_NEUTRAL" %in% colnames(var_calls)){

    var_calls_bp4 <- var_calls |>
      dplyr::mutate(
        ACMG_BP4 = dplyr::if_else(
          .data$N_INSILICO_CALLED >= dbnsfp_min_called &
            .data$N_INSILICO_TOLERATED >= dbnsfp_min_majority &
            .data$N_INSILICO_DAMAGING <= dbnsfp_max_minority &
            .data$N_INSILICO_SPLICING_AFFECTED == 0,
          TRUE,
          FALSE
        )
      ) |>
      dplyr::mutate(
        ACMG_BP4 = dplyr::if_else(
          .data$N_INSILICO_SPLICING_NEUTRAL == 2 &
            .data$N_INSILICO_SPLICING_AFFECTED == 0 &
            .data$N_INSILICO_CALLED == 2,
          TRUE,
          as.logical(.data$ACMG_BP4)
        )
      ) |>

      ## Remove BP4 if MES_PERCENT_CHANGE indicates
      ## splice site effect (lost or gained)
      dplyr::mutate(
        ACMG_BP4 = dplyr::if_else(
          !is.na(.data$MES_PERCENT_CHANGE) &
            (.data$MES_PERCENT_CHANGE <= -20 |
               .data$MES_PERCENT_CHANGE >= 60),
          FALSE,
          as.logical(.data$ACMG_BP4)
        )
      ) |>
      dplyr::select(
        c("VAR_ID","ACMG_BP4")
      )
  }
  if(nrow(var_calls_bp4) > 0){
    var_calls <- dplyr::left_join(
      var_calls,
      var_calls_bp4,
      by = "VAR_ID") |>
      dplyr::mutate(
        ACMG_BP4 = dplyr::if_else(
          is.na(.data$ACMG_BP4),
          FALSE,
          as.logical(.data$ACMG_BP4)
        )
      ) |>
    dplyr::distinct()
  } else {
    var_calls$ACMG_BP4 <- FALSE
  }
  return(var_calls)
}

#' Function that assigns ACMG evidence indicators for
#' BP7
#' @param var_calls variants in cancer predisposition genes
#' @return data.frame with ACMG evidence indicators
#' @export
#'
#'
assign_BP7_evidence <- function(var_calls = NULL){

  var_calls$ACMG_BP7 <- FALSE

  ## Assign logical ACMG evidence indicator
  # ACMG_BP7 - Silent/intronic variant outside of the splice site consensus
  #.          also 3_prime and 5_prime UTR variants included
  if(is.data.frame(var_calls) &
     "INTRON_POSITION" %in% colnames(var_calls) &
     "EXON_POSITION" %in% colnames(var_calls) &
     "MES_PERCENT_CHANGE" %in% colnames(var_calls) &
     "CONSEQUENCE" %in% colnames(var_calls)){

    var_calls <- var_calls |>
      dplyr::mutate(
        ACMG_BP7 =
          dplyr::if_else((
            (as.integer(.data$INTRON_POSITION) < 0 &
               as.integer(.data$INTRON_POSITION) < -3) |
              (as.integer(.data$INTRON_POSITION) > 0 &
                 as.integer(.data$INTRON_POSITION) > 6) |
              (as.integer(.data$EXON_POSITION) < 0 &
                 as.integer(.data$EXON_POSITION) < -2) |
              (as.integer(.data$EXON_POSITION) > 0 &
                 as.integer(.data$EXON_POSITION) > 1)) &
              (is.na(MES_PERCENT_CHANGE) |
                 (!is.na(MES_PERCENT_CHANGE) &
                    MES_PERCENT_CHANGE > -20 &
                    MES_PERCENT_CHANGE < 30)) &
              !stringr::str_detect(
                .data$CONSEQUENCE,
                "splice_(acceptor|donor)|stop_gained|5th"
              ) &
              stringr::str_detect(
                .data$CONSEQUENCE,
                paste0(
                  "^(synonymous_variant|intron_variant|",
                  "splice_region_variant)")
              ),
            TRUE, FALSE, FALSE
          )
      ) |>
      dplyr::mutate(
        ACMG_BP7 = dplyr::if_else(
          !is.na(.data$INTRON_POSITION) &
            .data$INTRON_POSITION == 0 &
            !is.na(.data$EXON_POSITION) &
            .data$EXON_POSITION == 0 &
            stringr::str_detect(
              .data$CONSEQUENCE,
              paste0(
                "^(3_prime_UTR_variant|
                  5_prime_UTR_variant)")
            ),
          TRUE, as.logical(.data$ACMG_BP7), FALSE
        )
      )
  }
  return(var_calls)
}

#' Function that assigns ACMG evidence indicators for
#' PVS1 ('ACMG_PVS1', 'ACMG_PVS1_STR', 'ACMG_PVS1_MOD')
#' @param var_calls variants in cancer predisposition genes
#' @return data.frame with ACMG evidence indicators
#' @export
#'
assign_PVS1_evidence <- function(
    var_calls = NULL){

  var_calls$ACMG_PVS1 <- FALSE
  var_calls$ACMG_PVS1_STR <- FALSE
  var_calls$ACMG_PVS1_MOD <- FALSE

  ## Assign logical ACMG evidence indicator
  # PVS1 - Null variant in a gene where loss-of-function
  # is a known mechanism of disease
  if(is.data.frame(var_calls) &
     "LOSS_OF_FUNCTION" %in% colnames(var_calls) &
     "LOF_FILTER" %in% colnames(var_calls) &
     "REFSEQ_TRANSCRIPT_ID" %in% colnames(var_calls) &
     "PROTEIN_RELATIVE_POSITION" %in% colnames(var_calls) &
     "NMD" %in% colnames(var_calls) &
     "EXONIC_STATUS" %in% colnames(var_calls) &
     "NULL_VARIANT" %in% colnames(var_calls) &
     "INTRON_POSITION" %in% colnames(var_calls) &
     "LAST_INTRON" %in% colnames(var_calls) &
     "CPG_MOD" %in% colnames(var_calls) &
     "MANE_SELECT" %in% colnames(var_calls) &
     "MANE_SELECT2" %in% colnames(var_calls) &
     "MANE_PLUS_CLINICAL" %in% colnames(var_calls) &
     "MANE_PLUS_CLINICAL2" %in% colnames(var_calls) &
     "PFAM_DOMAIN_NAME" %in% colnames(var_calls) &
     "PROTEIN_RELATIVE_POSITION" %in% colnames(var_calls) &
     "MAXENTSCAN" %in% colnames(var_calls) &
     "CONSEQUENCE" %in% colnames(var_calls)){


    var_calls <- var_calls |>
      dplyr::mutate(
        ## Frameshift and stop-gain variants (null variants)
        ## - Classified as loss-of-function
        ## - Predicted NOT to escape NMD
        ## - Loss-of-function known mechanism of disease for gene
        ## - biologically relevant transcript (MANE Select)
        ## -  or in curated set missing in MANE for grch37 (e.g. MEN1)
        ## ---> VeryStrong (PVS1)
        ACMG_PVS1 = dplyr::if_else(
          .data$NULL_VARIANT == TRUE &
            is.na(.data$NMD) &
            (.data$LOSS_OF_FUNCTION == TRUE |
               (.data$LOSS_OF_FUNCTION == FALSE &
                  (!is.na(.data$LOF_FILTER) &
                     .data$LOF_FILTER == "END_TRUNCATION") &
                  (!is.na(.data$PROTEIN_RELATIVE_POSITION) &
                     as.numeric(
                       .data$PROTEIN_RELATIVE_POSITION) < 0.975))) &
            (!is.na(.data$CPG_MOD) &
               .data$CPG_MOD == "LoF") &
            .data$EXONIC_STATUS == "exonic" &
            (!is.na(.data$MANE_SELECT) |
               !is.na(.data$MANE_SELECT2) |
               !is.na(.data$MANE_PLUS_CLINICAL) |
               !is.na(.data$MANE_PLUS_CLINICAL2) |
               .data$REFSEQ_TRANSCRIPT_ID %in%
               cpsr::curated_transcripts$id),
          TRUE,
          as.logical(.data$ACMG_PVS1)
        )) |>
      dplyr::mutate(
        ACMG_PVS1 = dplyr::if_else(

          ## Splice-donor/acceptor variants (+/- 2bp)
          ## - Exonic location (includes splice siste)
          ## - Not last intron
          ## - Not a null variant (stop-gain/frameshift)
          ## - Intron position +/- 1 or 2
          ## - Classified as loss-of-function
          ## - Loss-of-function known mechanism of disease for gene
          ## - Biologically relevant transcript (MANE Select)
          ## ---> VeryStrong (PVS1)
          .data$EXONIC_STATUS == "exonic" &
            .data$LAST_INTRON == FALSE &
            .data$NULL_VARIANT == FALSE &
            (.data$INTRON_POSITION == -2 |
               .data$INTRON_POSITION == -1 |
               .data$INTRON_POSITION == 1 |
               .data$INTRON_POSITION == 2) &
            (.data$CONSEQUENCE == "splice_acceptor_variant" |
               .data$CONSEQUENCE == "splice_donor_variant") &
            .data$LOSS_OF_FUNCTION == TRUE &
            (!is.na(.data$CPG_MOD) &
               .data$CPG_MOD == "LoF") &
            (!is.na(.data$MANE_SELECT) |
               !is.na(.data$MANE_SELECT2) |
               !is.na(.data$MANE_PLUS_CLINICAL) |
               !is.na(.data$MANE_PLUS_CLINICAL2) |
               .data$REFSEQ_TRANSCRIPT_ID %in%
               cpsr::curated_transcripts$id),
          TRUE,
          as.logical(.data$ACMG_PVS1)
        )) |>

      dplyr::mutate(
        ## Frameshift and stop-gain variants (null variants)
        ## - Predicted to escape NMD
        ## - Loss-of-function known mechanism of disease for gene
        ## - Biologically relevant transcript (MANE Select)
        ## - Position overlapping protein domain OR
        ##     removing > 10% of protein
        ## ---> Strong (PVS1_STR)
        ACMG_PVS1_STR = dplyr::case_when(
          .data$NULL_VARIANT == TRUE &
            (!is.na(.data$NMD) &
            .data$NMD == "NMD_escaping_variant") &
            (!is.na(.data$CPG_MOD) &
               .data$CPG_MOD == "LoF") &
            (!is.na(.data$MANE_SELECT) |
               !is.na(.data$MANE_SELECT2) |
               !is.na(.data$MANE_PLUS_CLINICAL) |
               !is.na(.data$MANE_PLUS_CLINICAL2) |
               .data$REFSEQ_TRANSCRIPT_ID %in%
               cpsr::curated_transcripts$id) &
            (!is.na(.data$PFAM_DOMAIN_NAME) |
               (!is.na(.data$PROTEIN_RELATIVE_POSITION) &
                  as.numeric(.data$PROTEIN_RELATIVE_POSITION) < 0.9)) ~ TRUE,

          ## Splice-donor/acceptor variants (+/- 2bp)
          ## - Exonic location (includes splice site)
          ## - Not a null variant (stop-gain/frameshift)
          ## - Classified as loss-of-function
          ## - Intron position +/- 1 or 2
          ## - Last intron (proxy for NMD escaping)
          ## - Loss-of-function known mechanism of disease for gene
          ## - Biologically relevant transcript (MANE Select)
          ## - Position overlapping protein domain OR
          ##     removing > 10% of protein
          ## ---> Strong (PVS1_STR)
          .data$NULL_VARIANT == FALSE &
            .data$LOSS_OF_FUNCTION == TRUE &
            ((.data$INTRON_POSITION == -2 |
               .data$INTRON_POSITION == -1 |
               .data$INTRON_POSITION == 1 |
               .data$INTRON_POSITION == 2) &
               (.data$CONSEQUENCE == "splice_acceptor_variant" |
                  .data$CONSEQUENCE == "splice_donor_variant")) &
            .data$EXONIC_STATUS == "exonic" &
            .data$LAST_INTRON == TRUE &
            (!is.na(.data$CPG_MOD) &
               .data$CPG_MOD == "LoF") &
            (!is.na(.data$MANE_SELECT) |
               !is.na(.data$MANE_SELECT2) |
               !is.na(.data$MANE_PLUS_CLINICAL) |
               !is.na(.data$MANE_PLUS_CLINICAL2) |
               .data$REFSEQ_TRANSCRIPT_ID %in%
               cpsr::curated_transcripts$id) ~ TRUE,
            #(!is.na(.data$PFAM_DOMAIN_NAME) |
               #(!is.na(.data$PROTEIN_RELATIVE_POSITION) &
                  #as.numeric(.data$PROTEIN_RELATIVE_POSITION) < 0.9)) ~ TRUE,
          TRUE ~ as.logical(.data$ACMG_PVS1_STR)
        )
      ) |>

      dplyr::mutate(
        ## Frameshift and stop-gain variants
        ## - Predicted to escape NMD
        ## - Loss-of-function known mechanism of disease for gene
        ## - Biologically relevant transcript (MANE Select)
        ## - Position NOT overlapping protein domain AND
        ##     removing > 10% of protein
        ## ---> Moderate (PVS1_MOD)
        ACMG_PVS1_MOD = dplyr::case_when(
          .data$NULL_VARIANT == TRUE &
            !is.na(.data$NMD) &
            .data$NMD == "NMD_escaping_variant" &
            (!is.na(.data$CPG_MOD) &
               .data$CPG_MOD == "LoF") &
            (!is.na(.data$MANE_SELECT) |
               !is.na(.data$MANE_SELECT2) |
               !is.na(.data$MANE_PLUS_CLINICAL) |
               !is.na(.data$MANE_PLUS_CLINICAL2) |
               .data$REFSEQ_TRANSCRIPT_ID %in%
               cpsr::curated_transcripts$id) &
            is.na(.data$PFAM_DOMAIN_NAME) &
            (!is.na(.data$PROTEIN_RELATIVE_POSITION) &
               as.numeric(.data$PROTEIN_RELATIVE_POSITION) > 0.9) ~ TRUE,

          ## Splice-donor/acceptor variants (+/- 2bp)
          ## - Exonic location (includes splice site)
          ## - Not a null variant (stop-gain/frameshift)
          ## - Intron position +/- 1 or 2
          ## - Last intron (proxy for NMD escaping)
          ## - Loss-of-function known mechanism of disease for gene
          ## - Biologically relevant transcript (MANE Select)
          ## - Position NOT overlapping protein domain AND
          ##     removing < 10% of protein
          ## ---> Moderate (PVS1_MOD)
          # .data$NULL_VARIANT == FALSE &
          #   ((.data$INTRON_POSITION == -2 |
          #       .data$INTRON_POSITION == -1 |
          #       .data$INTRON_POSITION == 1 |
          #       .data$INTRON_POSITION == 2) |
          #      (.data$CONSEQUENCE == "splice_acceptor_variant" |
          #         .data$CONSEQUENCE == "splice_donor_variant")) &
          #   .data$EXONIC_STATUS == "exonic" &
          #   .data$LAST_INTRON == TRUE &
          #   (!is.na(.data$CPG_MOD) &
          #      .data$CPG_MOD == "LoF") &
          #   (!is.na(.data$MANE_SELECT) |
          #      !is.na(.data$MANE_SELECT2) |
          #      !is.na(.data$MANE_PLUS_CLINICAL) |
          #      !is.na(.data$MANE_PLUS_CLINICAL2) |
          #      .data$REFSEQ_TRANSCRIPT_ID %in%
          #      cpsr::curated_transcripts$id) ~ TRUE,
          #   #is.na(.data$PFAM_DOMAIN_NAME) &
          #   #(!is.na(.data$PROTEIN_RELATIVE_POSITION) &
          #      #as.numeric(.data$PROTEIN_RELATIVE_POSITION) > 0.9) ~ TRUE,
          TRUE ~ as.logical(.data$ACMG_PVS1_MOD)
        )
      ) |>
      dplyr::mutate(
        ACMG_PVS1_MOD = dplyr::if_else(
          stringr::str_detect(
            .data$CONSEQUENCE,
            "^(start_lost|intron_variant&splice_donor_5th_base_variant)$"
          ),
          TRUE,
          as.logical(.data$ACMG_PVS1_MOD)
        )
      )

  }
  return(var_calls)

}

#' Function that assigns ACMG evidence indicators for
#' PS1 ('ACMG_PS1')
#' @param var_calls variants in cancer predisposition genes
#' @param pathogenic_codons data.frame with known pathogenic codons
#' pathogenic variants have been observed (from ClinVar)
#' @param pathogenic_nucleotides data.frame with pathogenic
#' nucleotides (non-coding, e.g. splice)
#' @return data.frame with ACMG evidence indicators
#' @export
#'
assign_PS1_evidence <- function(
    var_calls = NULL,
    pathogenic_codons = NULL,
    pathogenic_nucleotides = NULL){


  var_calls$ACMG_PS1 <- NULL

  var_calls_ps1 <- data.frame()
  ## Assign logical ACMG evidence indicator
  # PS1 - Same amino acid change as a previously established
  # pathogenic variant regardless of nucleotide change
  if(is.data.frame(var_calls) &
     NROW(var_calls) > 0 &
     "ENTREZGENE" %in% colnames(var_calls) &
     "CONSEQUENCE" %in% colnames(var_calls) &
     "HGVSP_query" %in% colnames(var_calls) &
     "GRANTHAM_DISTANCE" %in% colnames(var_calls) &
     "VAR_ID" %in% colnames(var_calls) &
     "SYMBOL" %in% colnames(var_calls) &
     "CODON" %in% colnames(var_calls) &

     is.data.frame(pathogenic_codons) &
     "ENTREZGENE" %in% colnames(pathogenic_codons) &
     "CODON" %in% colnames(pathogenic_codons) &
     "NUM_PATHOGENIC_AA" %in% colnames(pathogenic_codons) &
     "MIN_GRANTHAM_DISTANCE" %in% colnames(pathogenic_codons) &
     "PATH_HGVSP" %in% colnames(pathogenic_codons) &
     "PATH_VAR_ID" %in% colnames(pathogenic_codons) &
     NROW(pathogenic_codons) > 0 &

     is.data.frame(pathogenic_nucleotides) &
     "ENTREZGENE" %in% colnames(pathogenic_nucleotides) &
     "PATH_VAR_ID_NUC" %in% colnames(pathogenic_nucleotides) &
     "PATH_NUC_SITE" %in% colnames(pathogenic_nucleotides) &
     NROW(pathogenic_nucleotides) > 0){

    var_calls <-
      var_calls[, !(colnames(var_calls) %in%
                      c("PATH_VAR_ID_NUC",
                        "PATH_NUC_SITE",
                        "NUM_PATHOGENIC_AA",
                        "MIN_GRANTHAM_DISTANCE",
                        "PATH_HGVSP",
                        "PATH_VAR_ID"))]

    var_calls_ps1 <- var_calls |>
      dplyr::select(
        c("CONSEQUENCE","ENTREZGENE","GRANTHAM_DISTANCE",
          "SYMBOL","CODON","VAR_ID","HGVSP_query")
      ) |>
      dplyr::left_join(
        pathogenic_codons,
        by = c("ENTREZGENE","CODON")) |>
      tidyr::separate(
        "VAR_ID", c("CHROM","POS","REF","ALT"),
        sep="_", remove = F) |>
      dplyr::mutate(
        PATH_NUC_SITE = paste(
          .data$CHROM, .data$POS,
          .data$REF, sep="_"
        )
      ) |>
      dplyr::left_join(
        pathogenic_nucleotides,
        by = c("ENTREZGENE",
               "PATH_NUC_SITE")
      ) |>

      ## Check for variants at pathogenic codon/sites
      ## ensuring that
      ## i) amino acid change for query variant is known, AND
      ## iii) not the same genomic variant (VAR_ID)
      dplyr::mutate(ACMG_PS1 = dplyr::if_else(
        !is.na(.data$NUM_PATHOGENIC_AA) &
          .data$NUM_PATHOGENIC_AA >= 1 &
          !is.na(.data$PATH_HGVSP) &
          !is.na(.data$HGVSP_query) &
          ## amino acid change same as known pathogenic
          stringr::str_detect(
            .data$PATH_HGVSP, .data$HGVSP_query) &
          !is.na(.data$PATH_VAR_ID) &
          ## not same variant (VAR_ID)
          !stringr::str_detect(
            .data$PATH_VAR_ID,
            paste0("(",.data$VAR_ID,"(;|$))")),
        TRUE, FALSE, FALSE
      )) |>

      ## Check for variants at essential splice nucleotide sites
      ## - ensure that PS1 has not already been assigned
      ## - ensure that genomic variant is not the same as
      ##   those recorded (VAR_ID)
      dplyr::mutate(ACMG_PS1 = dplyr::if_else(
        .data$ACMG_PS1 == FALSE &
          !is.na(.data$PATH_NUC_SITE) &
          !stringr::str_detect(
            .data$PATH_VAR_ID_NUC,
            paste0("(",.data$VAR_ID,"(;|$))")),
        TRUE, FALSE, FALSE
      )) |>
      dplyr::select(
        c("VAR_ID","ACMG_PS1")
      ) |>
      dplyr::distinct()

  }

  if(NROW(var_calls_ps1) > 0){
    var_calls <- var_calls |>
      dplyr::left_join(
        var_calls_ps1,
        by = c("VAR_ID")
      ) |>
      dplyr::mutate(
        ACMG_PS1 = dplyr::if_else(
          is.na(.data$ACMG_PS1),
          FALSE,
          as.logical(.data$ACMG_PS1)
        )
      )
  }else{
    var_calls$ACMG_PS1 <- FALSE
  }

  return(var_calls)
}


#' Function that assigns ACMG evidence indicators for
#' PM1 ('ACMG_PM1' and 'ACMG_PM1_SUPP')
#'
#' @param var_calls variants in cancer predisposition genes
#' @return data.frame with ACMG evidence indicators
#' @export
#'
assign_PM1_evidence <- function(var_calls = NULL){

  var_calls_hotspot <- data.frame()
  var_calls_motifs <- data.frame()

  if(!is.null(var_calls) & is.data.frame(var_calls)){
    var_calls$ACMG_PM1 <- FALSE
    var_calls$ACMG_PM1_SUPP <- FALSE

    if("MUTATION_HOTSPOT" %in% colnames(var_calls) &
       "MUTATION_HOTSPOT_CANCERTYPE" %in% colnames(var_calls) &
       "AMINO_ACID_START" %in% colnames(var_calls) &
       "CODON" %in% colnames(var_calls) &
       "SYMBOL" %in% colnames(var_calls) &
       "CONSEQUENCE" %in% colnames(var_calls) &
       "HGVSP" %in% colnames(var_calls) &
       "VAR_ID" %in% colnames(var_calls) &
       "AMINO_ACID_END" %in% colnames(var_calls) &
       "ENTREZGENE" %in% colnames(var_calls)){


      ## Make data.frame with sites overlapping with
      ## somatic mutation hotspots (cancerhotspots.org)
      ## foundation for PM1
      var_calls_hotspot <- var_calls |>
        dplyr::filter(
          !is.na(MUTATION_HOTSPOT) &
            stringr::str_detect(
              .data$MUTATION_HOTSPOT,
              "^exonic") &
            stringr::str_detect(
              .data$CONSEQUENCE,
              "^missense_variant"
            ))

      if(NROW(var_calls_hotspot) > 0){
        var_calls_hotspot <- var_calls_hotspot |>
          dplyr::select(
            c("VAR_ID",
              "SYMBOL",
              "AMINO_ACID_END",
              "ENTREZGENE",
              "CONSEQUENCE",
              "MUTATION_HOTSPOT",
              "MUTATION_HOTSPOT_CANCERTYPE",
              "AMINO_ACID_START",
              "HGVSP",
              "CODON")
          ) |>
          tidyr::separate_rows(
            "MUTATION_HOTSPOT_CANCERTYPE",
            sep=","
          ) |>
          tidyr::separate(
            .data$MUTATION_HOTSPOT_CANCERTYPE,
            c("hotspot_site",
              "hotspot_site_num",
              "hotspot_aa_num"),
            sep = "\\|",
            remove = F,
            extra = "drop"
          ) |>
          dplyr::mutate(
            hotspot_aa_num = dplyr::case_when(
              nchar(hotspot_aa_num) == 0 ~ NA_integer_,
              TRUE ~ as.integer(.data$hotspot_aa_num)
            )
          ) |>
          dplyr::group_by(
            .data$ENTREZGENE,
            .data$SYMBOL,
            .data$CONSEQUENCE,
            .data$MUTATION_HOTSPOT,
            .data$VAR_ID,
            .data$HGVSP,
            .data$CODON,
            .data$AMINO_ACID_START,
            .data$AMINO_ACID_END
          ) |>
          dplyr::summarise(
            MUTATION_HOTSPOT_AA_NUM = sum(.data$hotspot_aa_num),
            .groups = "drop"
          ) |>
          dplyr::filter(
            .data$MUTATION_HOTSPOT_AA_NUM >= 2
          )

        if(NROW(var_calls_hotspot) > 0){
          var_calls_hotspot <- var_calls_hotspot |>
            tidyr::separate(
              .data$MUTATION_HOTSPOT,
              c("hotspot_region",
                "hotspot_symbol",
                "hotspot_entrezgene",
                "hotspot_codon",
                "hotspot_aa",
                "hotspot_pvalue"),
              sep = "\\|", remove = F, extra = "drop"
            ) |>
            dplyr::filter(
              !stringr::str_detect(
                .data$hotspot_codon, "-"
              ))

          if(NROW(var_calls_hotspot) > 0){
            var_calls_hotspot <- var_calls_hotspot |>
              dplyr::mutate(
                hotspot_entrezgene = as.character(
                  .data$hotspot_entrezgene)
              ) |>
              dplyr::mutate(
                hotspot_codon =
                  dplyr::if_else(
                    !is.na(.data$hotspot_codon),
                    paste0("p.", .data$hotspot_codon),
                    as.character(NA)
                  )
              ) |>
              dplyr::filter(
                !is.na(.data$hotspot_codon) &
                  !is.na(.data$hotspot_entrezgene) &
                  !is.na(.data$CODON) &
                  !is.na(.data$ENTREZGENE) &
                  .data$hotspot_entrezgene == .data$ENTREZGENE &
                  .data$hotspot_codon == .data$CODON
              )

            if(NROW(var_calls_hotspot) > 0){
              var_calls_hotspot <- var_calls_hotspot |>
                dplyr::select(
                  c("VAR_ID",
                    "SYMBOL",
                    "MUTATION_HOTSPOT_AA_NUM"))
            }
          }
        }
      }
    }


    if("AMINO_ACID_START" %in% colnames(var_calls) &
       "SYMBOL" %in% colnames(var_calls) &
       "CONSEQUENCE" %in% colnames(var_calls) &
       "AMINO_ACID_END" %in% colnames(var_calls) &
       "VAR_ID" %in% colnames(var_calls)){

      ## Make data.frame with sites overlapping with
      ## functional motifs in PTEN and TP53 - PM1
      var_calls_motifs <- var_calls |>
        dplyr::filter(
          !is.na(.data$AMINO_ACID_START) &
            !is.na(.data$AMINO_ACID_END) &
            stringr::str_detect(
              .data$CONSEQUENCE,
              "^missense_variant"
            ))

      if(NROW(var_calls_motifs) > 0){
        var_calls_motifs <- var_calls_motifs |>

          ## PTEN catalytic motifs: 90-94, 123-130, 166-168
          ## TP53 codons: 175, 245, 248, 249, 273, 282
          dplyr::filter(
            (.data$SYMBOL == "PTEN" &
               ((.data$AMINO_ACID_START >= 90 &
                  .data$AMINO_ACID_START <= 94) |
                  (.data$AMINO_ACID_START >= 123 &
               .data$AMINO_ACID_START <= 130) |
                  (.data$AMINO_ACID_START >= 166 &
                     .data$AMINO_ACID_START <= 168))) |
              (.data$SYMBOL == "TP53" &
                 (.data$AMINO_ACID_START %in% c(175, 245, 248, 249, 273, 282)))
          )

        if(NROW(var_calls_motifs) > 0){
          var_calls_motifs <- var_calls_motifs |>
            dplyr::mutate(VAR_FUNCTIONAL_MOTIF = T) |>
            dplyr::select(
              c("VAR_ID",
                "SYMBOL",
                "VAR_FUNCTIONAL_MOTIF",
              ))
        }
      }
    }

    var_calls_pm1 <- var_calls |>
      dplyr::select(
        c("VAR_ID","SYMBOL")
      )

    if(NROW(var_calls_hotspot) > 0){
      var_calls_pm1 <- var_calls |>
        dplyr::select(
          c("VAR_ID","SYMBOL")
        ) |>
        dplyr::left_join(
          var_calls_hotspot,
          by = c("VAR_ID", "SYMBOL")
        )
    }else{
      var_calls_pm1$MUTATION_HOTSPOT_AA_NUM <- NA_integer_
    }

    if(NROW(var_calls_motifs) > 0){
      var_calls_pm1 <- var_calls_pm1 |>
        dplyr::left_join(
          var_calls_motifs,
          by = c("VAR_ID", "SYMBOL")
        ) |>
        dplyr::mutate(
          VAR_FUNCTIONAL_MOTIF = dplyr::if_else(
            is.na(.data$VAR_FUNCTIONAL_MOTIF),
            FALSE,
            as.logical(.data$VAR_FUNCTIONAL_MOTIF)
          )
        )
    }else{
      var_calls_pm1$VAR_FUNCTIONAL_MOTIF <- FALSE
    }

    if(NROW(var_calls_pm1) > 0){
      var_calls <- var_calls |>
        dplyr::left_join(
          var_calls_pm1,
          by = c("VAR_ID", "SYMBOL")
        ) |>
        dplyr::mutate(
          ACMG_PM1 = dplyr::if_else(
            (!(.data$SYMBOL %in%
                 c("BRCA1", "BRCA2",
                   "PALB2","ATM",
                   "MSH2", "MLH1", "MSH6",
                   "PMS2","APC")) &
               ((!is.na(.data$MUTATION_HOTSPOT_AA_NUM) &
                   .data$MUTATION_HOTSPOT_AA_NUM >= 10) |
                  .data$VAR_FUNCTIONAL_MOTIF == TRUE)),
            TRUE,
            FALSE
          )
        ) |>
        dplyr::mutate(
          ACMG_PM1_SUPP = dplyr::if_else(
            (!(.data$SYMBOL %in%
                 c("BRCA1", "BRCA2",
                   "PALB2", "ATM",
                   "MSH2", "MLH1", "MSH6",
                   "PMS2","APC")) &
               .data$ACMG_PM1 == FALSE &
               ((!is.na(.data$MUTATION_HOTSPOT_AA_NUM) &
                   .data$MUTATION_HOTSPOT_AA_NUM >= 2 &
                   .data$MUTATION_HOTSPOT_AA_NUM <= 9))),
            TRUE,
            FALSE
          )
        ) |>
        pcgrr::remove_cols_from_df(
          cnames = c("MUTATION_HOTSPOT_AA_NUM",
                 "VAR_FUNCTIONAL_MOTIF")
        )
    }
  }


  return(var_calls)

}

#' Function that assigns ACMG evidence indicators for
#' PM2_supporting
#'
#' @param var_calls variants in cancer predisposition genes
#' @return data.frame with ACMG evidence indicators
#'
#' @export
#'
assign_PM2_evidence <- function(
    var_calls = NULL){

  var_calls$ACMG_PM2_SUPP <- FALSE

  ## Assign logical ACMG evidence indicator
  # PM2 - Absence of variant in population databases
  if(is.data.frame(var_calls) &
     "gnomAD_NC_AF_POPMAX" %in% colnames(var_calls) &
     "gnomAD_NC_AC_POPMAX" %in% colnames(var_calls) &
     "gnomAD_AF_POPMAX" %in% colnames(var_calls)){

    var_calls <-
      var_calls |>
      dplyr::mutate(
        ACMG_PM2_SUPP =
          dplyr::if_else(
            ## safeguard with gnomAD total population (< 1e-6)
            ## - check that allele frequency is very low
            !(.data$SYMBOL %in%
                c("PALB2","ATM","PTEN",
                  "MSH2", "MLH1", "MSH6",
                  "PMS2","APC","TP53")) &
            ((is.na(.data$gnomAD_NC_AC_POPMAX) |
              (!is.na(.data$gnomAD_NC_AC_POPMAX) &
                 (.data$gnomAD_NC_AC_POPMAX == 0 |
                    .data$gnomAD_NC_AC_POPMAX > 0 &
                    .data$gnomAD_NC_AF_POPMAX < 0.00001))) &
                 (is.na(.data$gnomAD_AF_POPMAX) |
                    (!is.na(.data$gnomAD_AF_POPMAX) &
                       .data$gnomAD_AF_POPMAX < 0.000001))),
            TRUE, FALSE, FALSE
          )
      ) |>
      dplyr::mutate(
        ACMG_PM2_SUPP = dplyr::case_when(
          .data$SYMBOL == "TP53" &
            (.data$gnomAD_AF_POPMAX < 0.00003 |
               is.na(.data$gnomAD_AF_POPMAX)) ~ TRUE,
          .data$SYMBOL == "PALB2" &
            (.data$gnomAD_AF_POPMAX < 0.00000333 |
               is.na(.data$gnomAD_AF_POPMAX)) ~ TRUE,
          .data$SYMBOL == "ATM" &
            (.data$gnomAD_AF_POPMAX <= 0.00001 |
               is.na(.data$gnomAD_AF_POPMAX)) ~ TRUE,
          .data$SYMBOL == "APC" &
            ((.data$gnomAD_NC_AF_POPMAX < 0.000003 &
               .data$gnomAD_NC_AC_POPMAX > 1) |
            (.data$gnomAD_NC_AF_POPMAX < 0.00001 &
               .data$gnomAD_NC_AC_POPMAX <= 1) |
               is.na(.data$gnomAD_NC_AF_POPMAX)) ~ TRUE,
          (.data$SYMBOL == "MLH1" |
             .data$SYMBOL == "MSH2" |
             .data$SYMBOL == "MSH6" |
             .data$SYMBOL == "PMS2") &
            (.data$gnomAD_AF_POPMAX < 0.00002 |
               is.na(.data$gnomAD_AF_POPMAX)) ~ TRUE,
          .data$SYMBOL == "PTEN" &
            (.data$gnomAD_AF_POPMAX < 0.00001 |
               is.na(.data$gnomAD_AF_POPMAX)) ~ TRUE,
          TRUE ~ as.logical(.data$ACMG_PM2_SUPP)
        )
      )
  }
  return(var_calls)
}


#' Function that assigns ACMG evidence indicators for
#' PM4
#'
#' @param var_calls variants in cancer predisposition genes
#' @return data.frame with ACMG evidence indicators
#' #' @export
#'
assign_PM4_evidence <- function(
    var_calls = NULL){

  var_calls$ACMG_PM4 <- FALSE

  ## Assign logical ACMG evidence indicator
  # PM4 - Protein length changes due to in-frame deletions/insertions
  # in a non-repeat region or stop-loss variants
  if(is.data.frame(var_calls) &
     "CONSEQUENCE" %in% colnames(var_calls) &
     "RMSK_HIT" %in% colnames(var_calls) &
     "PFAM_DOMAIN_NAME" %in% colnames(var_calls)){

    var_calls <- var_calls |>
      dplyr::mutate(
        ACMG_PM4 = dplyr::case_when(
          (((.data$CONSEQUENCE == "inframe_deletion" |
             .data$CONSEQUENCE == "inframe_insertion" &
             is.na(.data$RMSKHIT))
          | .data$CONSEQUENCE == "stop_lost") &
            !is.na(.data$PFAM_DOMAIN_NAME)) ~ TRUE,
          TRUE ~ as.logical(.data$ACMG_PM4)
        )
      )
  }
  return(var_calls)

}



#' Function that assigns ACMG evidence indicators for
#' PM5
#' @param var_calls variants in cancer predisposition genes
#' @param pathogenic_codons data.frame with codons where
#' pathogenic variants have been observed (from ClinVar)
#' @return data.frame with ACMG evidence indicators
#' @export
#'
assign_PM5_evidence <- function(
    var_calls = NULL,
    pathogenic_codons = NULL){

  var_calls_pm5 <- data.frame()
  var_calls$ACMG_PM5 <- NULL

  ## Assign logical ACMG evidence indicator
  # PM5 - Novel missense change at an amino acid residue
  # where a different missense change determined to be
  # pathogenic has been seen before
  if(is.data.frame(var_calls) &
     NROW(var_calls) > 0 &
     "ENTREZGENE" %in% colnames(var_calls) &
     "CONSEQUENCE" %in% colnames(var_calls) &
     "HGVSP_query" %in% colnames(var_calls) &
     "GRANTHAM_DISTANCE" %in% colnames(var_calls) &
     "VAR_ID" %in% colnames(var_calls) &
     "SYMBOL" %in% colnames(var_calls) &

     is.data.frame(pathogenic_codons) &
     "ENTREZGENE" %in% colnames(pathogenic_codons) &
     "CODON" %in% colnames(pathogenic_codons) &
     "NUM_PATHOGENIC_AA" %in% colnames(pathogenic_codons) &
     "MIN_GRANTHAM_DISTANCE" %in% colnames(pathogenic_codons) &
     "PATH_HGVSP" %in% colnames(pathogenic_codons) &
     "PATH_VAR_ID" %in% colnames(pathogenic_codons) &
     NROW(pathogenic_codons) > 0){

    var_calls_pm5 <- var_calls |>
      dplyr::select(
        c("CONSEQUENCE","ENTREZGENE","GRANTHAM_DISTANCE",
          "SYMBOL","CODON","VAR_ID","HGVSP_query")
      ) |>
      dplyr::left_join(
        pathogenic_codons,
        by = c("ENTREZGENE","CODON")) |>

      ## Check for variants at pathogenic codon/sites
      ## ensuring that
      ## i) amino acid change for query variant is novel, AND
      ## ii) Grantham distance >= minimum for what is found
      ##     among pathogenic variants at the given codon
      ## iii) not the same variant (VAR_ID)
      dplyr::mutate(ACMG_PM5 = dplyr::if_else(
        !is.na(.data$NUM_PATHOGENIC_AA) &
          .data$NUM_PATHOGENIC_AA >= 1 &
          !is.na(.data$PATH_HGVSP) &
          !is.na(.data$HGVSP_query) &
          ## not same aa change
          !stringr::str_detect(
            .data$PATH_HGVSP, .data$HGVSP_query) &
          ## Grantham distance check
          .data$GRANTHAM_DISTANCE >= .data$MIN_GRANTHAM_DISTANCE &
          !is.na(.data$PATH_VAR_ID) &
          ## not same variant (VAR_ID)
          !stringr::str_detect(
            .data$PATH_VAR_ID, .data$VAR_ID),
        TRUE, FALSE, FALSE
      )) |>
      dplyr::select(
        c("VAR_ID","ACMG_PM5")
      )

  }

  if(NROW(var_calls_pm5) > 0){
    var_calls <- var_calls |>
      dplyr::left_join(
        var_calls_pm5,
        by = c("VAR_ID")
      ) |>
      dplyr::mutate(
        ACMG_PM5 = dplyr::if_else(
          is.na(.data$ACMG_PM5),
          FALSE,
          as.logical(.data$ACMG_PM5)
        )
      )
  }else{
    var_calls$ACMG_PM5 <- FALSE
  }

  return(var_calls)

}

#' Function that assigns ACMG evidence indicators for
#' PP3
#'
#' @param var_calls variants in cancer predisposition genes
#' @return data.frame with ACMG evidence indicators
#' @export
#'
assign_PP3_evidence <- function(
    var_calls = NULL){

  ## Assign logical ACMG evidence indicators
  #
  # ACMG_PP3 - Multiple lines (>=8) of insilico evidence support a
  #             deleterious effect on the gene or gene product
  #          (conservation, evolutionary, splicing impact, etc.)
  #
  # Computational evidence for deleterious/benign effect is taken from
  # invidual algorithm predictions in dbNSFP: SIFT,Provean,MutationTaster,
  # MutationAssessor,M_CAP,MutPred,FATHMM_XF, DBNSFP_RNN,dbscSNV_RF,
  # dbscSNV_AdaBoost, AlphaMissense, ClinPred, phactBOOST, PrimateAI, REVEL,
  # Default scheme (from default TOML file):
  # 1) Damaging: Among all possible protein variant effect predictions, at
  #              least six algorithms must have made a call,
  #              with at least 8 predicted as damaging/D
  #              (possibly_damaging/PD), and at most two
  #              predicted as tolerated/T (PP3)
  #       - at most 1 prediction for a splicing neutral effect
  #    Exception: if both splice site predictions indicate damaging effects;
  #    ignore other criteria

  dbnsfp_min_majority <- cpsr::acmg[["insilico_pred_min_majority"]]
  dbnsfp_max_minority <- cpsr::acmg[["insilico_pred_max_minority"]]
  dbnsfp_min_called <- dbnsfp_min_majority

  var_calls_pp3 <- data.frame()
  var_calls$ACMG_PP3 <- NULL

  if(is.data.frame(var_calls) &
     "N_INSILICO_CALLED" %in% colnames(var_calls) &
     "N_INSILICO_DAMAGING" %in% colnames(var_calls) &
     "N_INSILICO_TOLERATED" %in% colnames(var_calls) &
     "N_INSILICO_SPLICING_AFFECTED" %in% colnames(var_calls) &
     "N_INSILICO_SPLICING_NEUTRAL" %in% colnames(var_calls)){

    var_calls_pp3 <- var_calls |>
      ## Missense predictions
      dplyr::mutate(
        ACMG_PP3 =
          dplyr::if_else(
            .data$N_INSILICO_CALLED >= dbnsfp_min_called &
              .data$N_INSILICO_DAMAGING >= dbnsfp_min_majority &
              .data$N_INSILICO_TOLERATED <= dbnsfp_max_minority &
              .data$N_INSILICO_SPLICING_NEUTRAL <= 1, TRUE,
            FALSE, FALSE
          )
      ) |>

      ## Splice site effect only
      dplyr::mutate(
        ACMG_PP3 = dplyr::case_when(
        .data$N_INSILICO_SPLICING_AFFECTED == 2 &
          .data$N_INSILICO_CALLED >= 2 &
          .data$N_INSILICO_TOLERATED == 0 ~ TRUE,
        TRUE ~ as.logical(.data$ACMG_PP3)
      )) |>
      dplyr::select(c("VAR_ID","ACMG_PP3")) |>
      dplyr::distinct()
  }

  if(NROW(var_calls_pp3) > 0){
    var_calls <- var_calls |>
      dplyr::left_join(
        var_calls_pp3,
        by = c("VAR_ID")
      ) |>
      dplyr::mutate(
        ACMG_PP3 = dplyr::if_else(
          is.na(.data$ACMG_PP3),
          FALSE,
          as.logical(.data$ACMG_PP3)
        )
      )
  }else{
    var_calls$ACMG_PP3 <- FALSE
  }

  return(var_calls)
}

#' Function that calculates the closest distance
#' to a known pathogenic variant in the same gene
#'
#' @param var_calls variants in cancer predisposition genes
#' @param ref_data reference data object with ClinVar
#' annotations
#' @return data.frame with closest distance to known
#' pathogenic variant in same gene
#'
#' @export
#'
min_distance_to_pathogenic <- function(
    var_calls = NULL,
    ref_data = NULL){

  distance_df <- data.frame()
  var_calls$MIN_DISTANCE_TO_PATHOGENIC <- NULL
  ## Calculate minimum distance to known pathogenic

  if(is.data.frame(var_calls) &
     NROW(var_calls) > 0 &
     "ENTREZGENE" %in% colnames(var_calls) &
     "CONSEQUENCE" %in% colnames(var_calls) &
     "CODON" %in% colnames(var_calls) &
     "VAR_ID" %in% colnames(var_calls) &

     is.list(ref_data) &
     "variant" %in% names(ref_data) &
     is.list(ref_data[['variant']]) &
     "clinvar_aa_sites" %in% names(ref_data[['variant']]) &
     is.data.frame(ref_data[['variant']][['clinvar_aa_sites']]) &
     "ENTREZGENE" %in% colnames(ref_data[['variant']][['clinvar_aa_sites']]) &
     "CODON" %in% colnames(ref_data[['variant']][['clinvar_aa_sites']]) &
     "PATHOGENIC" %in% colnames(ref_data[['variant']][['clinvar_aa_sites']]) &
     NROW(ref_data[['variant']][['clinvar_aa_sites']]) > 0){

    pathogenic_df <- ref_data[['variant']][['clinvar_aa_sites']] |>
      dplyr::filter(.data$PATHOGENIC == 1) |>
      dplyr::select(c("ENTREZGENE", "CODON")) |>
      dplyr::filter(!is.na(.data$ENTREZGENE) & !is.na(.data$CODON)) |>
      dplyr::distinct()

    query_df <- var_calls |>
      dplyr::select(c("ENTREZGENE", "CODON","VAR_ID","CONSEQUENCE")) |>
      dplyr::filter(!is.na(.data$ENTREZGENE) & !is.na(.data$CODON)) |>
      dplyr::distinct()

    if(NROW(query_df) == 0){
      var_calls$MIN_DISTANCE_TO_PATHOGENIC <- NA_integer_
      return(var_calls)
    }

    query_df <- query_df |>
      dplyr::filter(
        .data$CONSEQUENCE == "missense_variant"
      )

    if(NROW(query_df) == 0){
      var_calls$MIN_DISTANCE_TO_PATHOGENIC <- NA_integer_
      return(var_calls)
    }

    query_df <- query_df |>
      dplyr::mutate(
        CODON_POSITION = as.numeric(
          gsub("[^0-9]", "", .data$CODON)
      ))

    pathogenic_df <- pathogenic_df |>
      dplyr::mutate(
        CODON_POSITION_PATH = as.numeric(
          gsub("[^0-9]", "", .data$CODON)
      ))


    # Join on gene
    distance_df <- query_df |>
      dplyr::select(
        c("VAR_ID", "ENTREZGENE", "CODON_POSITION")
      ) |>
      dplyr::inner_join(pathogenic_df, by = "ENTREZGENE",
                        relationship = "many-to-many")

    if(NROW(distance_df) > 0){
      distance_df <- distance_df |>
        dplyr::mutate(distance_to_path = abs(
          .data$CODON_POSITION - .data$CODON_POSITION_PATH
        )) |>
        dplyr::group_by(VAR_ID, ENTREZGENE) |>
        dplyr::summarise(
          MIN_DISTANCE_TO_PATHOGENIC = min(distance_to_path), .groups = "drop") |>
        dplyr::filter(.data$MIN_DISTANCE_TO_PATHOGENIC > 0 &
                        .data$MIN_DISTANCE_TO_PATHOGENIC <= 50) |>
        dplyr::distinct() |>
        dplyr::select(c("VAR_ID", "MIN_DISTANCE_TO_PATHOGENIC"))
    }

  }

  if(NROW(distance_df) > 0){
    var_calls <- var_calls |>
      dplyr::left_join(
        distance_df,
        by = c("VAR_ID")
      )
  }
  return(var_calls)

}

#' Function that extracts known benign peptide changes
#' from ClinVar (>= 2 gold stars), as well as known
#' pathogenic nucleotide and codon sites
#'
#' @param ref_data reference data object with ClinVar
#' annotations
#' @return list with data.frames for benign peptide changes,
#' pathogenic nucleotide sites, and pathogenic codon sites
#' @export
#'
known_path_benign_sites <- function(
    ref_data = NULL){

  known_sites <- list()
  known_sites[['benign_peptide']] <- data.frame()
  known_sites[['pathogenic_nucleotide']] <- data.frame()
  known_sites[['pathogenic_codon']] <- data.frame()

  if(is.null(ref_data) |
     !is.list(ref_data) |
     !("variant" %in% names(ref_data)) |
     !is.list(ref_data[['variant']]) |
     !("clinvar_aa_sites" %in% names(ref_data[['variant']])) |
     !is.data.frame(ref_data[['variant']][['clinvar_aa_sites']]) |
     !("clinvar_nuc_sites" %in% names(ref_data[['variant']])) |
     !is.data.frame(ref_data[['variant']][['clinvar_nuc_sites']])){
    return(known_sites)
  }

  known_sites[['benign_peptide']] <-
    ref_data[['variant']][['clinvar_aa_sites']] |>
    dplyr::filter(.data$GOLD_STARS >= 2 & .data$BENIGN == 1) |>
    dplyr::select(c("ENTREZGENE", "HGVSP", "BENIGN", "VAR_ID")) |>
    dplyr::filter(!is.na(.data$ENTREZGENE) & !is.na(.data$HGVSP)) |>
    dplyr::rename(BENIGN_PEPTIDE_CHANGE = "BENIGN",
                  BENIGN_VAR_ID = "VAR_ID") |>
    dplyr::distinct()

  known_sites[['pathogenic_nucleotide']] <-
    ref_data[['variant']][['clinvar_nuc_sites']] |>
    dplyr::filter(.data$GOLD_STARS >= 2) |>
    dplyr::select(c("ENTREZGENE", "PATH_NUC_SITE",
                    "PATH_VAR_ID_NUC"))

  known_sites[['pathogenic_codon']] <-
    ref_data[['variant']][['clinvar_aa_sites']] |>
    dplyr::filter(.data$SYMBOL == "MSH2") |>
    dplyr::filter(.data$PATHOGENIC == 1) |>
    dplyr::select(
      c("ENTREZGENE", "CODON",
        "PATHOGENIC","GOLD_STARS",
        "VAR_ID","GRANTHAM_DISTANCE","HGVSP")) |>
    tidyr::separate_rows("VAR_ID", sep=";") |>
    dplyr::filter(
      !is.na(.data$ENTREZGENE) &
        !is.na(.data$CODON) &
        !is.na(GRANTHAM_DISTANCE)) |>
    dplyr::distinct() |>
    dplyr::group_by(.data$ENTREZGENE,  .data$CODON) |>
    dplyr::reframe(
      MIN_GRANTHAM_DISTANCE = min(.data$GRANTHAM_DISTANCE, na.rm=TRUE),
      MAX_GOLD_STARS = max(.data$GOLD_STARS, na.rm=TRUE),
      PATH_HGVSP = paste(unique(.data$HGVSP), collapse=";"),
      PATH_VAR_ID = paste(unique(.data$VAR_ID), collapse=";")) |>
    dplyr::mutate(
      NUM_PATHOGENIC_AA = stringr::str_count(.data$PATH_HGVSP,";") + 1
    ) |>
    dplyr::filter(MAX_GOLD_STARS >= 2) |>
    dplyr::distinct()

  return(known_sites)
}


#' Function that filters novel (non-ClinVar) variants based
#' on a MAF threshold in gnomAD
#'
#' @param var_calls variants in cancer predisposition genes
#' @param conf CPSR configuration object with run settings
#'
#' @export
exclude_vars_by_maf <- function(
    var_calls = NULL,
    conf = NULL){

  if(is.data.frame(var_calls) &
     "CPSR_CLASSIFICATION_SOURCE" %in% colnames(var_calls) &
     "POPMAX_AF_GNOMAD" %in% colnames(var_calls) &
     "FINAL_CLASSIFICATION" %in% colnames(var_calls) &
     "CPSR_PATHOGENICITY_SCORE" %in% colnames(var_calls) &
     "CANCER_PHENOTYPE" %in% colnames(var_calls) &
     NROW(var_calls) > 0){

      non_clinvar_vars <-
        var_calls |>
        dplyr::filter(
          .data$CPSR_CLASSIFICATION_SOURCE == "CPSR_ACMG")
      clinvar_vars <-
        var_calls |>
        dplyr::filter(
          .data$CPSR_CLASSIFICATION_SOURCE == "ClinVar")

      if(NROW(non_clinvar_vars) > 0){
        non_clinvar_vars_maf_filtered <-
          non_clinvar_vars |>
          dplyr::filter(
            is.na(.data$POPMAX_AF_GNOMAD) |
              (!is.na(.data$POPMAX_AF_GNOMAD) &
                 .data$POPMAX_AF_GNOMAD <=
              conf[["variant_classification"]][["max_af_gnomad"]]))

        pcgrr::log4r_info(
          paste0(
            "Ignoring n = ",
            NROW(non_clinvar_vars) - NROW(non_clinvar_vars_maf_filtered),
            " unclassified variants with a global MAF frequency above ",
            conf[["variant_classification"]][["max_af_gnomad"]]
          )
        )
      }

      var_calls <-
        dplyr::bind_rows(
          non_clinvar_vars, clinvar_vars) |>
        dplyr::arrange(
          dplyr::desc(
            .data$FINAL_CLASSIFICATION),
          dplyr::desc(.data$CANCER_PHENOTYPE),
          dplyr::desc(.data$CPSR_PATHOGENICITY_SCORE))
  }

  return(var_calls)

}

#' Function that filters novel (non-ClinVar) variants based
#' on a MAF threshold in gnomAD
#'
#' @param var_calls variants in cancer predisposition genes
#' @param conf CPSR configuration object with run settings
#'
#' @export
exclude_vars_noncoding <- function(
    var_calls = NULL,
    conf = NULL){

  if(conf[['other']][['show_noncoding']] == T){
    return(var_calls)
  }
  if("CODING_STATUS" %in% colnames(var_calls) &
     NROW(var_calls) > 0){
    n_noncoding_vars <-
      var_calls |>
      dplyr::filter(
        .data$CODING_STATUS == "noncoding") |>
      NROW()
    pcgrr::log4r_info(
      paste0(
        "Excluding n = ",
        n_noncoding_vars,
        " variants for classification (option --ignore_noncoding)"
      )
    )
    var_calls <-
      var_calls |>
      dplyr::filter(
        .data$CODING_STATUS == "coding")

    # if (NROW(var_calls) == 0) {
    #   pcgrr::log4r_warn(paste0(
    #     "There are zero remaining protein-coding ",
    #     "calls in input file - no report will be produced"
    #   ))
    #   return(NULL)
    # }
  }

  return(var_calls)

}

#' Function that ignores/exclude variants attributed to non-cancer
#' related phenotypes in ClinVar
#'
#' @param var_calls variants in cancer predisposition genes
#' @param conf CPSR configuration object with run settings
#'
#' @export
exclude_vars_non_cancer <- function(
    var_calls = NULL,
    conf = NULL){

  if (conf$variant_classification$clinvar_report_noncancer == 0 &
      "VAR_ID" %in% colnames(
        var_calls) &
      "CLINVAR_MSID" %in% colnames(
        var_calls) &
      "CANCER_PHENOTYPE" %in% colnames(
        var_calls)) {

    non_cancer_calls_to_exclude <-
      var_calls |>
      dplyr::filter(!is.na(.data$CLINVAR_MSID)) |>
      dplyr::filter(
        !is.na(.data$CANCER_PHENOTYPE) &
          .data$CANCER_PHENOTYPE == 0)

    if(NROW(non_cancer_calls_to_exclude) > 0){
      non_cancer_calls_to_exclude <-
        non_cancer_calls_to_exclude |>
        dplyr::select("VAR_ID") |>
        dplyr::distinct()
    }

    pcgrr::log4r_info(
      paste0(
        "ClinVar variants related to non-cancer conditions excluded",
        " from report: N = ", NROW(non_cancer_calls_to_exclude)
      )
    )

    if (NROW(non_cancer_calls_to_exclude) > 0) {
      var_calls <-
        var_calls |>
        dplyr::anti_join(
          non_cancer_calls_to_exclude, by = "VAR_ID")

      pcgrr::log4r_info(
        paste0(
          "Variants remaining after exclusion of non-cancer related",
          " ClinVar variants: N = ",
          NROW(var_calls)
        )
      )
    }
  }

  return(var_calls)

}
