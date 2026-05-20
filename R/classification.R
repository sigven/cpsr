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
    "Assigning variant classification codes according to ACMG/AMP criteria")
  classification_settings <-
    settings$conf$variant_classification

  acmg_ev_codes <-
    c(
      "ACMG_PM1_MOD",
      "ACMG_PM1_SUPP",
      "ACMG_PM2_SUPP",
      "ACMG_BA1",
      "ACMG_BP1",
      "ACMG_BP4",
      "ACMG_BP7",
      "ACMG_BS1",
      "ACMG_BS1_SUPP",
      "ACMG_PVS1",
      "ACMG_PVS1_STR",
      "ACMG_PVS1_MOD",
      "ACMG_PS1",
      "ACMG_PP3",
      "ACMG_PM5",
      "ACMG_PM4",
      "ACMG_PM4_SUPP",
      "ACMG_PP2",
      "HGVSP_query"
    )

  var_calls <-
    var_calls[, !(colnames(var_calls) %in% acmg_ev_codes)]

  ## Load known pathogenic and benign sites data from ClinVar
  known_sites <- cpsr::known_path_benign_sites(
    ref_data = ref_data
  )

  ## Assign ACMG evidence codes
  var_calls <- var_calls |>
    tidyr::separate(
      "HGVSp",c("tmpENSP","HGVSp_long"),
      sep=":", remove = F) |>
    dplyr::rename("HGVSP_query" = "HGVSp_long") |>
    dplyr::select(-c("tmpENSP")) |>
    cpsr::assign_BA1_evidence() |>
    cpsr::assign_BS1_evidence() |>
    cpsr::assign_BP1_evidence() |>
    cpsr::assign_BP4_evidence() |>
    cpsr::assign_BP7_evidence() |>
    cpsr::assign_PM2_evidence() |>
    cpsr::assign_PM1_evidence() |>
    cpsr::assign_PVS1_evidence() |>
    cpsr::assign_PP3_evidence() |>
    cpsr::assign_PM4_evidence() |>
    cpsr::assign_PM5_evidence(
      pathogenic_codons =
        known_sites[['pathogenic_codon']]) |>
    cpsr::assign_PS1_evidence(
      pathogenic_codons =
        known_sites[['pathogenic_codon']],
      pathogenic_nucleotides =
        known_sites[['pathogenic_nucleotide']]
    ) |>
    cpsr::assign_PP2_evidence()

  return(var_calls)
}

#' Function that assigns a summary string of
#' all ACMG evidence codes met for a given variant
#' (e.g. "PVS1|PS1|PM2_supporting|PP3")
#'
#' @param var_calls data frame with variant calls in predisposition genes
#'
#' @return var_calls data frame with ACMG_CODE column appended
#' @export
#'
assign_acmg_consensus <- function(
    var_calls = NULL){

  if (!is.data.frame(var_calls)) {
    return(var_calls)
  }

  p_lower_limit <- cpsr::acmg[['score_thresholds']][['p_lower']]
  lp_upper_limit <- cpsr::acmg[['score_thresholds']][['lp_upper']]
  lp_lower_limit <- cpsr::acmg[['score_thresholds']][['lp_lower']]
  vus_upper_limit <- cpsr::acmg[['score_thresholds']][['vus_upper']]
  vus_lower_limit <- cpsr::acmg[['score_thresholds']][['vus_lower']]
  lb_upper_limit <- cpsr::acmg[['score_thresholds']][['lb_upper']]
  lb_lower_limit <- cpsr::acmg[['score_thresholds']][['lb_lower']]
  b_upper_limit <- cpsr::acmg[['score_thresholds']][['b_upper']]

  ## Restrict to evidence codes that are present as columns in var_calls
  ev_codes <- cpsr::acmg[["evidence_codes"]]
  active_codes <- ev_codes[
    ev_codes$cpsr_evidence_code %in% colnames(var_calls), ]

  if (nrow(active_codes) == 0) {
    return(var_calls)
  }

  ## Derive human-readable display names from internal code names:
  ##   ACMG_PVS1_STR  -> PVS1_strong
  ##   ACMG_PVS1_MOD  -> PVS1_moderate
  ##   ACMG_PM2_SUPP  -> PM2_supporting   (etc.)
  active_codes$display_code <- active_codes$cpsr_evidence_code |>
    stringr::str_remove("^ACMG_") |>
    stringr::str_replace("_STR$", "_strong") |>
    stringr::str_replace("_MOD$", "_moderate") |>
    stringr::str_replace("_SUPP$", "_supporting")

  code_names <- active_codes$cpsr_evidence_code
  disp_names <- active_codes$display_code
  score_map  <- stats::setNames(active_codes$path_score, code_names)

  var_calls <- var_calls |>
    dplyr::mutate(
      ACMG_CODE = "",
      ACMG_DOC = "",
      CPSR_PATHOGENICITY_SCORE = 0
    )

  ## Build ACMG_CODE: pipe-separated labels for all TRUE criteria.
  ## Construct a character matrix (nrow x n_codes): each cell is the display
  ## name when the criterion is TRUE, "" otherwise. Then collapse row-wise.
  code_flags <- do.call(cbind, lapply(seq_along(code_names), function(i) {
    col <- var_calls[[code_names[i]]]
    ifelse(!is.na(col) & col, disp_names[i], "")
  }))
  var_calls$ACMG_CODE <- apply(code_flags, 1, function(row) {
    paste(row[nzchar(row)], collapse = "|")
  })

  ## Accumulate CPSR_PATHOGENICITY_SCORE: vectorized column-by-column
  for (code in code_names) {
    score <- score_map[[code]]
    if (score != 0) {
      var_calls$CPSR_PATHOGENICITY_SCORE <-
        var_calls$CPSR_PATHOGENICITY_SCORE +
        score * as.numeric(!is.na(var_calls[[code]]) & var_calls[[code]])
    }
  }

  ## Adjustment: a score driven solely by the single supporting-benign
  ## criterion BP4 should not tip into Likely Benign
  bp4_score <- ev_codes$path_score[ev_codes$cpsr_evidence_code == "ACMG_BP4"]
  if (length(bp4_score) == 1) {
    var_calls <- var_calls |>
      dplyr::mutate(CPSR_PATHOGENICITY_SCORE = dplyr::if_else(
        .data$CPSR_PATHOGENICITY_SCORE == bp4_score &
          stringr::str_detect(.data$ACMG_CODE, "BP4"),
        0,
        as.numeric(.data$CPSR_PATHOGENICITY_SCORE)
      ))
  }

  ## Safeguard: LP/P classification requires converging evidence from at least
  ## two independent pathogenic criteria, unless a strong/very-strong criterion
  ## (path_score >= 4: PVS1 family, PS1) is present. A single moderate or
  ## supporting criterion reaching the LP threshold is a scoring artefact, not
  ## genuine multi-criterion support (e.g. PM4 alone should remain VUS).
  pathogenic_codes <- active_codes$cpsr_evidence_code[
    active_codes$path_score > 0]
  strong_codes <- active_codes$cpsr_evidence_code[
    active_codes$path_score >= 4]

  if (length(pathogenic_codes) > 0) {
    path_matrix <- do.call(cbind, lapply(pathogenic_codes, function(code) {
      as.numeric(!is.na(var_calls[[code]]) & var_calls[[code]])
    }))
    var_calls$.n_pathogenic_criteria <- rowSums(path_matrix)
  } else {
    var_calls$.n_pathogenic_criteria <- 0L
  }

  if (length(strong_codes) > 0) {
    strong_matrix <- do.call(cbind, lapply(strong_codes, function(code) {
      as.numeric(!is.na(var_calls[[code]]) & var_calls[[code]])
    }))
    var_calls$.has_strong_criterion <- rowSums(strong_matrix) > 0
  } else {
    var_calls$.has_strong_criterion <- FALSE
  }

  var_calls <- var_calls |>
    dplyr::mutate(
      CPSR_CLASSIFICATION =
        dplyr::case_when(
          ## Safeguard: only moderate/supporting evidence, single criterion
          .data$CPSR_PATHOGENICITY_SCORE >= lp_lower_limit &
            !.data$.has_strong_criterion &
            .data$.n_pathogenic_criteria < 2 ~ "VUS",
          .data$CPSR_PATHOGENICITY_SCORE <= lb_upper_limit &
            .data$CPSR_PATHOGENICITY_SCORE >= lb_lower_limit ~ "Likely Benign",
          .data$CPSR_PATHOGENICITY_SCORE <= b_upper_limit ~ "Benign",
          .data$CPSR_PATHOGENICITY_SCORE <= vus_upper_limit &
            .data$CPSR_PATHOGENICITY_SCORE >= vus_lower_limit ~ "VUS",
          .data$CPSR_PATHOGENICITY_SCORE >= p_lower_limit ~ "Pathogenic",
          .data$CPSR_PATHOGENICITY_SCORE >= lp_lower_limit &
            .data$CPSR_PATHOGENICITY_SCORE <= lp_upper_limit ~ "Likely Pathogenic",
          TRUE ~ as.character("VUS")
        )
    ) |>
    dplyr::select(-c(".n_pathogenic_criteria", ".has_strong_criterion"))

  return(var_calls)

}


#' Function that combines classifications of novel and
#' pre-classified variants
#'
#' @param var_calls variants in cancer predisposition genes
#' @param conf CPSR configuration object with run settings
#'
#' @export
assign_classification_authority <-
  function(var_calls,
           conf = NULL){

    assertable::assert_colnames(
      var_calls,
      c("CPSR_CLASSIFICATION",
        "CPSR_PATHOGENICITY_SCORE",
        "ACMG_CODE",
        "ACMG_DOC",
        "gnomADe_AF",
        "CLINVAR_PHENOTYPE_CANCER",
        "CLINVAR_CONFLICTED",
        "CLINVAR_GOLD_STARS",
        "CLINVAR_CLASSIFICATION"),
      only_colnames = F,quiet = T
    )

    clinvar_trust_level <- 0
    if (!is.null(conf) &&
        !is.null(conf$variant_classification) &&
        !is.null(conf$variant_classification$clinvar_trust_level)) {
      clinvar_trust_level <- conf$variant_classification$clinvar_trust_level
    }

    stopifnot(
      is.numeric(clinvar_trust_level),
      clinvar_trust_level %in% 0:4
    )

    trust_level_desc <- c(
      "0" = "ClinVar trusted (override conflicted records only)",
      "1" = "Override zero-star ClinVar records",
      "2" = "Override zero- and single-star ClinVar records",
      "3" = "Override low-star and non-cancer-phenotype records",
      "4" = "CPSR always classifies"
    )

    pcgrr::log4r_info(
      paste0(
        "Classification authority assignment - clinvar_trust_level = ",
        clinvar_trust_level,
        " (", trust_level_desc[as.character(clinvar_trust_level)], ")"
      )
    )

    ## Step 1: Determine whether the ClinVar record meets
    ## the user's trust threshold. Build override reasons
    ## as a list of flags, then collapse.

    var_calls <- var_calls |>
      dplyr::mutate(
        ## -- individual override reason flags --
        .clinvar_novel = is.na(.data$CLINVAR_CLASSIFICATION),
        .clinvar_conflicted = (!is.na(.data$CLINVAR_CONFLICTED) &
                                 .data$CLINVAR_CONFLICTED == TRUE) |
          (!is.na(.data$CLINVAR_CLASSIFICATION) &
             .data$CLINVAR_CLASSIFICATION == "CCP Unknown"),
        .clinvar_zero_stars = (!is.na(.data$CLINVAR_GOLD_STARS) &
                                 .data$CLINVAR_GOLD_STARS == 0),
        .clinvar_single_star = (!is.na(.data$CLINVAR_GOLD_STARS) &
                                  .data$CLINVAR_GOLD_STARS == 1),
        .clinvar_no_cancer_pheno = (!is.na(.data$CLINVAR_PHENOTYPE_CANCER) &
                                      .data$CLINVAR_PHENOTYPE_CANCER == 0)
      ) |>
      dplyr::mutate(
        ## -- which flags are active at the chosen trust level --
        .override_novel = .data$.clinvar_novel,
        .override_conflicted = .data$.clinvar_conflicted,
        .override_stars = if (clinvar_trust_level >= 2) {
          .data$.clinvar_zero_stars | .data$.clinvar_single_star
        } else if (clinvar_trust_level >= 1) {
          .data$.clinvar_zero_stars
        } else {
          FALSE
        },
        .override_phenotype = if (clinvar_trust_level >= 3) {
          .data$.clinvar_no_cancer_pheno
        } else {
          FALSE
        },
        .override_always = clinvar_trust_level == 4,

        ## -- should CPSR take authority for this record? --
        .cpsr_overrides = (
          .data$.override_novel |
            .data$.override_conflicted |
            .data$.override_stars |
            .data$.override_phenotype |
            .data$.override_always
        )
      )

    ## Step 2: Build human-readable rationale strings
    var_calls <- var_calls |>
      dplyr::rowwise() |>
      dplyr::mutate(
        ASSERTION_RATIONALE = {
          reasons <- character(0)
          if (.data$.clinvar_novel)
            reasons <- c(reasons, "Novel variant (absent from ClinVar)")
          if (.data$.clinvar_conflicted)
            reasons <- c(reasons, "Conflicting ClinVar interpretations")
          if (!.data$.clinvar_novel & !.data$.clinvar_conflicted) {
            if (.data$.cpsr_overrides) {
              star_label <- paste0(
                "ClinVar review status: ",
                .data$CLINVAR_GOLD_STARS,
                ifelse(.data$CLINVAR_GOLD_STARS == 1,
                       " gold star", " gold stars"))
              reasons <- c(reasons, star_label)
              if (.data$.override_phenotype)
                reasons <- c(reasons,
                             "No cancer-related ClinVar phenotype(s)")
              if (.data$.override_always & length(reasons) == 1)
                reasons <- c(reasons,
                             "CPSR override (trust level 4)")
            } else {
              reasons <- paste0(
                "ClinVar review status: ",
                .data$CLINVAR_GOLD_STARS, " gold stars")
            }
          }
          paste(reasons, collapse = "; ")
        }
      ) |>
      dplyr::ungroup()

    ## Step 3: Assign authority and final classification
    var_calls <- var_calls |>
      dplyr::mutate(
        ASSERTION_AUTHORITY = dplyr::case_when(
          .data$.cpsr_overrides & !is.na(.data$CPSR_CLASSIFICATION) ~ "CPSR",
          !.data$.clinvar_novel & !.data$.cpsr_overrides ~ "ClinVar",
          TRUE ~ NA_character_
        ),
        CLASSIFICATION = dplyr::case_when(
          .data$ASSERTION_AUTHORITY == "ClinVar" ~
            as.character(.data$CLINVAR_CLASSIFICATION),
          .data$ASSERTION_AUTHORITY == "CPSR" ~
            as.character(.data$CPSR_CLASSIFICATION),
          TRUE ~ NA_character_
        )
      )

    ## Step 4: Warn on unresolved records
    n_unresolved <- sum(is.na(var_calls$ASSERTION_AUTHORITY))
    if (n_unresolved > 0) {
      pcgrr::log4r_warn(
        paste0(
          n_unresolved,
          " variant(s) have no assigned ASSERTION_AUTHORITY. ",
          "These lack both ClinVar and CPSR classifications."
        )
      )
    }

    ## Step 5: Log summary of authority assignments
    n_clinvar <- sum(
      var_calls$ASSERTION_AUTHORITY == "ClinVar", na.rm = TRUE)
    n_cpsr <- sum(
      var_calls$ASSERTION_AUTHORITY == "CPSR", na.rm = TRUE)
    pcgrr::log4r_info(
      paste0(
        "Classification authority summary: ",
        n_clinvar, " ClinVar, ",
        n_cpsr, " CPSR, ",
        n_unresolved, " unresolved"
      )
    )

    ## Clean up internal columns
    var_calls <- var_calls |>
      dplyr::select(-dplyr::starts_with(".clinvar_"),
                    -dplyr::starts_with(".override_"),
                    -dplyr::all_of(".cpsr_overrides")) |>
      dplyr::mutate(CLASSIFICATION = dplyr::if_else(
        .data$CLASSIFICATION != "Benign" &
          .data$CLASSIFICATION != "Likely Benign" &
          .data$CLASSIFICATION != "VUS" &
          .data$CLASSIFICATION != "Likely Pathogenic" &
          .data$CLASSIFICATION != "Pathogenic" &
          .data$CLASSIFICATION != "Risk Factor" &
          .data$CLASSIFICATION != "Drug Response",
        "VUS",
        as.character(.data$CLASSIFICATION)
      )) |>
      dplyr::mutate(CLASSIFICATION = dplyr::if_else(
        is.na(.data$CLASSIFICATION),"VUS",
        as.character(.data$CLASSIFICATION)
      )) |>
      dplyr::mutate(CLASSIFICATION = factor(
        .data$CLASSIFICATION,
        levels = c("Risk Factor",
                   "Drug Response",
                   "Benign",
                   "Likely Benign",
                   "VUS",
                   "Likely Pathogenic",
                   "Pathogenic")
      )) |>
      dplyr::arrange(
        dplyr::desc(
          .data$CLASSIFICATION),
        dplyr::desc(.data$CLINVAR_PHENOTYPE_CANCER),
        dplyr::desc(.data$CPSR_PATHOGENICITY_SCORE))

    var_calls$ACMG_DOC <-
      stringr::str_replace_all(
        var_calls$ACMG_DOC,
        "<br>-", ","
      )
    var_calls$ACMG_DOC <-
      stringr::str_replace_all(
        var_calls$ACMG_DOC,
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
          .data$SYMBOL == "TP53" &
            .data$gnomADj_FAF_GRPMAX > 0.0003 &
            .data$gnomADj_FAF_GRPMAX < 0.001 &
            .data$ACMG_BA1 == FALSE ~ TRUE,
          .data$SYMBOL == "BRCA1" &
            .data$gnomAD_NC_FAF_GRPMAX > 0.0001 &
            .data$ACMG_BA1 == FALSE ~ TRUE,
          .data$SYMBOL == "BRCA2" &
            .data$gnomAD_NC_FAF_GRPMAX > 0.0001 &
            .data$ACMG_BA1 == FALSE ~ TRUE,
          .data$SYMBOL == "PALB2" &
            .data$gnomAD_NC_FAF_GRPMAX > 0.0001 &
            .data$ACMG_BA1 == FALSE ~ TRUE,
          .data$SYMBOL == "ATM" &
            .data$gnomADj_FAF_GRPMAX > 0.0005 &
            .data$ACMG_BA1 == FALSE ~ TRUE,
          .data$SYMBOL == "APC" &
            .data$gnomAD_NC_FAF_GRPMAX >= 0.0001 &
            .data$ACMG_BA1 == FALSE ~ TRUE,
          .data$SYMBOL == "MLH1" &
            .data$gnomADj_FAF_GRPMAX >= 0.0001 &
            .data$gnomADj_FAF_GRPMAX < 0.001 ~ TRUE,
          .data$SYMBOL == "MSH2" &
            .data$gnomADj_FAF_GRPMAX >= 0.0001 &
            .data$gnomADj_FAF_GRPMAX < 0.001 ~ TRUE,
          .data$SYMBOL == "PTEN" &
            .data$gnomADj_FAF_GRPMAX > 0.000043 &
            .data$gnomADj_FAF_GRPMAX <= 0.00056 ~ TRUE,
          .data$SYMBOL == "MSH6" &
            .data$gnomADj_FAF_GRPMAX >= 0.00022 &
            .data$gnomADj_FAF_GRPMAX < 0.0022 ~ TRUE,
          .data$SYMBOL == "PMS2" &
            .data$gnomADj_FAF_GRPMAX >= 0.0001 &
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
#' BP1
#'
#' @param var_calls variants in cancer predisposition genes
#' @return data.frame with ACMG evidence indicators
#' @export
#'
assign_BP1_evidence <- function(var_calls = NULL){

  ## Assign a logical ACMG evidence indicator
  # ACMG_BP1 - Missense variant in a gene for which primarily truncating
  # variants (> 90%, as given in Maxwell et al.) are known to cause disease

  var_calls$ACMG_BP1 <- FALSE

  if(is.data.frame(var_calls) &
     "CONSEQUENCE" %in% colnames(var_calls) &
     "PATH_TRUNC_FRAC" %in% colnames(var_calls)){

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

  var_calls$ACMG_BP4 <- NULL
  var_calls_bp4 <- data.frame()

  if(is.data.frame(var_calls) &
     "N_INSILICO_CALLED" %in% colnames(var_calls) &
     "N_INSILICO_DAMAGING" %in% colnames(var_calls) &
     "N_INSILICO_TOLERATED" %in% colnames(var_calls) &
     "MAXENTSCAN" %in% colnames(var_calls) &
     "TF_BINDING_SITE_VARIANT_INFO" %in% colnames(var_calls) &
     "REGULATORY_ANNOTATION" %in% colnames(var_calls) &
     "EFFECT_PREDICTIONS" %in% colnames(var_calls) &
     "N_INSILICO_SPLICING_AFFECTED" %in% colnames(var_calls) &
     "N_INSILICO_SPLICING_NEUTRAL" %in% colnames(var_calls)){

    var_calls_bp4 <- var_calls |>
      dplyr::mutate(
        ACMG_BP4 = dplyr::if_else(
          (.data$N_INSILICO_CALLED >= dbnsfp_min_called &
            .data$N_INSILICO_TOLERATED >= dbnsfp_min_majority &
            .data$N_INSILICO_DAMAGING <= dbnsfp_max_minority &
            .data$N_INSILICO_SPLICING_AFFECTED == 0) |
            (.data$N_INSILICO_CALLED >= dbnsfp_min_called &
               .data$N_INSILICO_TOLERATED >= dbnsfp_min_majority &
               .data$N_INSILICO_DAMAGING <= dbnsfp_max_minority + 1 &
               .data$N_INSILICO_SPLICING_AFFECTED == 0 &
               stringr::str_detect(
                 .data$EFFECT_PREDICTIONS, "alphamissense:T")),
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

      ## add also upstream/downstream gene variant
      ## with no predicted regulatory effect to BP4
       dplyr::mutate(
         ACMG_BP4 = dplyr::if_else(
           stringr::str_detect(
             .data$CONSEQUENCE,
             "^(upstream_gene_variant|downstream_gene_variant)$") &

             ## ensure no predicted regulatory effect (i.e.
             ## no TF binding site or regulatory annotation)
             is.na(.data$REGULATORY_ANNOTATION) &
             is.na(.data$TF_BINDING_SITE_VARIANT_INFO) &

             ## ensure site is not conserved (GERP < 2) and not predicted to affect splicing
             (is.na(.data$GERP_SCORE) |
                (!is.na(.data$GERP_SCORE) &
                   .data$GERP_SCORE < 2)) &
             .data$N_INSILICO_SPLICING_AFFECTED == 0,
           TRUE,
           as.logical(.data$ACMG_BP4)
         )
       ) |>

      ## Remove BP4 if MAXENTSCAN_PCT_CHANGE indicates
      ## splice site effect (lost or gained)
      dplyr::mutate(
        ACMG_BP4 = dplyr::if_else(
          !is.na(.data$MAXENTSCAN) &
           stringr::str_detect(
             .data$MAXENTSCAN,
             "Strong|Moderate|Supporting"),
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
     "MAXENTSCAN" %in% colnames(var_calls) &
     "CONSEQUENCE" %in% colnames(var_calls)){

    var_calls <- var_calls |>
      tidyr::separate(
        "MAXENTSCAN",
        c("tmp_MES","MES_STRATUM","MES_TIER"),
        sep = "\\|", remove = FALSE
      ) |>
      dplyr::mutate(
        MES_STRATUM = dplyr::if_else(
          .data$MES_STRATUM == ".",
          as.character(NA),
          as.character(.data$MES_STRATUM)
        ),
        MES_TIER = dplyr::if_else(
          .data$MES_TIER == ".",
          as.character(NA),
          as.character(.data$MES_TIER)
        )
      ) |>
      dplyr::mutate(
        ACMG_BP7 = dplyr::case_when(
          ((is.na(.data$MES_STRATUM) &
            is.na(.data$MES_TIER)) |
             (.data$MES_STRATUM == "Acceptor_PPT_Distant" &
                .data$MES_TIER == "No_Call") |
             (.data$MES_STRATUM == "Acceptor_PPT_Close" &
                .data$MES_TIER == "No_Call")) &
            !stringr::str_detect(
              .data$CONSEQUENCE,
              "(splice_(acceptor|donor))|(stop_(gained|lost))|5th"
            ) &
            stringr::str_detect(
              .data$CONSEQUENCE,
              paste0(
                "(synonymous_variant|intron_variant)")
            ) ~ TRUE,
          TRUE ~ as.logical(.data$ACMG_BP7)
        )
      ) |>
      dplyr::mutate(
        ACMG_BP7 = dplyr::if_else(
          !is.na(.data$INTRON_POSITION) &
            .data$INTRON_POSITION == 0 &
            !is.na(.data$EXON_POSITION) &
            .data$EXON_POSITION == 0 &
            !stringr::str_detect(
              .data$CONSEQUENCE,
              "(splice_(acceptor|donor))|(stop_(gained|lost))|5th"
            ) &
            stringr::str_detect(
              .data$CONSEQUENCE,
              paste0(
                "^(3_prime_UTR_variant|",
                  "5_prime_UTR_variant)")
            ),
          TRUE,
          as.logical(.data$ACMG_BP7),
          FALSE
        )
      ) |>
      dplyr::select(
        -dplyr::any_of(
          c("tmp_MES","MES_STRATUM","MES_TIER")
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

  hgvsc_values <- rep(NA_character_, nrow(var_calls))
  if ("HGVSc" %in% colnames(var_calls)) {
    hgvsc_values <- as.character(var_calls$HGVSc)

    var_calls$PVS1_PURELY_INTRONIC <-
      cpsr::is_purely_intronic_deletion(var_calls$HGVSc)

  }

  ## Assign logical ACMG evidence indicator
  # PVS1 - Null variant in a gene where loss-of-function
  # is a known mechanism of disease
  if(is.data.frame(var_calls) &
     "LOSS_OF_FUNCTION" %in% colnames(var_calls) &
     "LOF_FILTER" %in% colnames(var_calls) &
     "REFSEQ_TRANSCRIPT_ID" %in% colnames(var_calls) &
     "EXON_INTRON_JUNCTION_SPAN" %in% colnames(var_calls) &
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
      tidyr::separate(
        "MAXENTSCAN",
        c("tmp_MES","MES_STRATUM","MES_TIER"),
        sep = "\\|", remove = FALSE
      ) |>
      dplyr::mutate(
        HGVSC_LOCAL = hgvsc_values,
        PVS1_RELEVANT_TRANSCRIPT =
          !is.na(.data$MANE_SELECT) |
          !is.na(.data$MANE_SELECT2) |
          !is.na(.data$MANE_PLUS_CLINICAL) |
          !is.na(.data$MANE_PLUS_CLINICAL2) |
          .data$REFSEQ_TRANSCRIPT_ID %in% cpsr::curated_transcripts$id,
        PVS1_LOF_GENE = !is.na(.data$CPG_MOD) &
          .data$CPG_MOD == "LoF",
        PVS1_SPLICE_CONSEQUENCE = stringr::str_detect(
          .data$CONSEQUENCE,
          "(^|&)splice_(acceptor|donor)_variant($|&)"
        ),
        PVS1_INFRAME_CONSEQUENCE = stringr::str_detect(
          .data$CONSEQUENCE,
          "(^|&)(inframe_deletion|inframe_insertion)($|&)"
        ),
        PVS1_HGVSC_CANONICAL_SITE = !is.na(.data$HGVSC_LOCAL) &
          stringr::str_detect(
            .data$HGVSC_LOCAL,
            "(\\+|\\-)(1|2)(?![0-9])"
          ) &
          # 2. Variant class is SNV or MNV — not deletion/insertion/indel
          #.data$VARIANT_CLASS %in% c("SNV", "substitution") &

          # 3. Primary consequence is splice acceptor or donor — not frameshift, UTR, inframe
          stringr::str_detect(
            .data$CONSEQUENCE,
            "(^|&)splice_(acceptor|donor)_variant($|&)"
          ) &

          # 4. Consequence does not include UTR, frameshift, or inframe as primary term
          !stringr::str_detect(
            .data$CONSEQUENCE,
            "(UTR|frameshift|inframe)"
          ),
        PVS1_SPLICE_HIGH_IMPACT =
          .data$PVS1_SPLICE_CONSEQUENCE &
          (.data$INTRON_POSITION %in% c(-2, -1, 1, 2) |
            .data$PVS1_HGVSC_CANONICAL_SITE |
            (.data$EXON_INTRON_JUNCTION_SPAN &
               !.data$PVS1_INFRAME_CONSEQUENCE)
          )
      ) |>
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
                     .data$LOF_FILTER == "END_TRUNCATION"))) &
            (!is.na(.data$PROTEIN_RELATIVE_POSITION) &
               as.numeric(
                 .data$PROTEIN_RELATIVE_POSITION) < 0.975) &
            .data$PVS1_LOF_GENE == TRUE &
            .data$EXONIC_STATUS == "exonic" &
            .data$PVS1_RELEVANT_TRANSCRIPT == TRUE,
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
            .data$PVS1_SPLICE_HIGH_IMPACT == TRUE &
            .data$LOSS_OF_FUNCTION == TRUE &
            .data$PVS1_LOF_GENE == TRUE &
            .data$PVS1_RELEVANT_TRANSCRIPT == TRUE,
          TRUE,
          as.logical(.data$ACMG_PVS1)
        )) |>

      dplyr::mutate(
        ## Frameshift and stop-gain variants (null variants)
        ## - Predicted to escape NMD
        ## - Loss-of-function known mechanism of disease for gene
        ## - Biologically relevant transcript (MANE Select)
        ## - Truncated/altered region in critical region (domain) OR
        ##.  removing greater than 10% of protein
        ## ---> Strong (PVS1_STR)
        ACMG_PVS1_STR = dplyr::case_when(
          .data$NULL_VARIANT == TRUE &
            (!is.na(.data$NMD) &
            .data$NMD == "NMD_escaping_variant") &
            .data$PVS1_LOF_GENE == TRUE &
            .data$PVS1_RELEVANT_TRANSCRIPT == TRUE &
            (!is.na(.data$PFAM_DOMAIN_NAME) |
               (is.na(.data$PFAM_DOMAIN_NAME) &
                  .data$PROTEIN_RELATIVE_POSITION < 0.9)) ~ TRUE,

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
            .data$PVS1_SPLICE_HIGH_IMPACT == TRUE &
            .data$EXONIC_STATUS == "exonic" &
            .data$LAST_INTRON == TRUE &
            .data$PVS1_LOF_GENE == TRUE &
            .data$PVS1_RELEVANT_TRANSCRIPT ~ TRUE,
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
            .data$PVS1_LOF_GENE == TRUE &
            .data$PVS1_RELEVANT_TRANSCRIPT == TRUE &
            is.na(.data$PFAM_DOMAIN_NAME) &
            (!is.na(.data$PROTEIN_RELATIVE_POSITION) &
               as.numeric(.data$PROTEIN_RELATIVE_POSITION) > 0.9) ~ TRUE,
          TRUE ~ as.logical(.data$ACMG_PVS1_MOD)
        )
      ) |>
      dplyr::mutate(
        ACMG_PVS1_MOD = dplyr::case_when(
          stringr::str_detect(
            .data$CONSEQUENCE,
            "^(start_lost)$") ~ TRUE,
          (stringr::str_detect(
            .data$CONSEQUENCE,
            "intron_variant") &
          ((.data$MES_STRATUM == "Donor_+5" |
            .data$MES_STRATUM == "Donor_+3") &
            (.data$MES_TIER == "Strong" |
               .data$MES_TIER == "Moderate"))) ~ TRUE,
          TRUE ~ as.logical(.data$ACMG_PVS1_MOD)
        )
      ) |>
      dplyr::select(
        -dplyr::any_of(c(
          "tmp_MES",
          "MES_STRATUM",
          "MES_TIER",
          "HGVSC_LOCAL"
          #"PVS1_RELEVANT_TRANSCRIPT",
          #"PVS1_LOF_GENE"
          #"PVS1_PURELY_INTRONIC",
          #"PVS1_SPLICE_CONSEQUENCE",
          #"PVS1_INFRAME_CONSEQUENCE",
          #"PVS1_HGVSC_CANONICAL_SITE",
          #"PVS1_SPLICE_HIGH_IMPACT"
        ))
      )

  }

  var_calls <- var_calls |>
    dplyr::mutate(
      ACMG_PVS1 = dplyr::if_else(
        is.na(.data$ACMG_PVS1), FALSE, .data$ACMG_PVS1),
      ACMG_PVS1_STR = dplyr::if_else(
        is.na(.data$ACMG_PVS1_STR), FALSE, .data$ACMG_PVS1_STR),
      ACMG_PVS1_MOD = dplyr::if_else(
        is.na(.data$ACMG_PVS1_MOD), FALSE, .data$ACMG_PVS1_MOD)
    )

  return(var_calls)

}

#' Function that assigns ACMG evidence indicators for
#' PS1 ('ACMG_PS1') - Same amino acid change as a previously
#' established pathogenic variant regardless of nucleotide
#' change
#'
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
      dplyr::filter(
        !is.na(.data$HGVSP_query) &
          !is.na(.data$VAR_ID)
      ) |>

      ## Check for variants at pathogenic codons, ensuring that
      ## i) amino acid change for query variant is known (as pathogenic), AND
      ## iii) is not the same genomic variant (VAR_ID)
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
      #
      # IGNORE FOR NOW - ONLY AMINO ACID CHANGES
      #
      # dplyr::mutate(ACMG_PS1 = dplyr::if_else(
      #   .data$ACMG_PS1 == FALSE &
      #     !is.na(.data$PATH_NUC_SITE) &
      #     !stringr::str_detect(
      #       .data$PATH_VAR_ID_NUC,
      #       paste0("(",.data$VAR_ID,"(;|$))")),
      #   TRUE, FALSE, FALSE
      # )) |>
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
          !is.na(.data$MUTATION_HOTSPOT) &
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
                "VAR_FUNCTIONAL_MOTIF"
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
            ((!is.na(.data$MUTATION_HOTSPOT_AA_NUM) &
                .data$MUTATION_HOTSPOT_AA_NUM >= 10) |
               .data$VAR_FUNCTIONAL_MOTIF == TRUE),
            TRUE,
            FALSE
          )
        ) |>
        dplyr::mutate(
          ACMG_PM1_SUPP = dplyr::if_else(
            .data$ACMG_PM1 == FALSE &
              ((!is.na(.data$MUTATION_HOTSPOT_AA_NUM) &
                  .data$MUTATION_HOTSPOT_AA_NUM >= 2 &
                  .data$MUTATION_HOTSPOT_AA_NUM <= 9)),
            TRUE,
            FALSE
          )
        ) |>

        ## Exclude cancer predisposition genes
        ## that explicitly indicate through their
        ## Expert Panel Specifications that PM1
        ## is not applicable
        dplyr::mutate(
          ACMG_PM1 = dplyr::if_else(
            .data$ACMG_PM1 == TRUE &
              .data$SYMBOL %in%
              c("PALB2","ATM","MSH2","MLH1",
                "MSH6","PMS2","APC","BRCA1","BRCA2"),
            FALSE,
            as.logical(.data$ACMG_PM1)
          ),
          ACMG_PM1_SUPP = dplyr::if_else(
            .data$ACMG_PM1_SUPP == TRUE &
              .data$SYMBOL %in%
              c("PALB2","ATM","MSH2","MLH1",
                "MSH6","PMS2","APC","BRCA1","BRCA2"),
            FALSE,
            as.logical(.data$ACMG_PM1_SUPP)
          ),
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
  # PM2 - Absence of variant in population databases (gnomAD non-cancer)
  #.      (or at very rare frequency: < 1e-5)
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
            ((!is.na(.data$gnomAD_NC_AF_POPMAX) &
                .data$gnomAD_NC_AF_POPMAX < 0.000003 &
               .data$gnomAD_NC_AC_POPMAX > 1) |
            (!is.na(.data$gnomAD_NC_AF_POPMAX) &
               .data$gnomAD_NC_AF_POPMAX < 0.00001 &
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
#' @export
#'
assign_PM4_evidence <- function(
    var_calls = NULL){

  var_calls$ACMG_PM4 <- FALSE
  var_calls$ACMG_PM4_SUPP <- FALSE

  ## Assign logical ACMG evidence indicator
  # PM4 - Protein length changes due to in-frame deletions/insertions
  # in a non-repeat region or stop-loss variants in
  # a functional protein domain
  if(is.data.frame(var_calls) &
     "CONSEQUENCE" %in% colnames(var_calls) &
     "SYMBOL" %in% colnames(var_calls) &
     "ACMG_PVS1" %in% colnames(var_calls) &
     "ACMG_PVS1_STR" %in% colnames(var_calls) &
     "ACMG_PVS1_MOD" %in% colnames(var_calls) &
     "REF" %in% colnames(var_calls) &
     "ALT" %in% colnames(var_calls) &
     "RMSK_HIT" %in% colnames(var_calls)){

    has_pfam <- "PFAM_DOMAIN_NAME" %in% colnames(var_calls)

    var_calls <- var_calls |>
      dplyr::mutate(
        ACMG_PM4 = dplyr::case_when(
          .data$ACMG_PVS1 == TRUE |
            .data$ACMG_PVS1_STR == TRUE |
            .data$ACMG_PVS1_MOD == TRUE ~ FALSE,
          .data$SYMBOL == "TP53" ~ FALSE,
          (((.data$CONSEQUENCE == "inframe_deletion" |
             .data$CONSEQUENCE == "inframe_insertion") &
               .data$SYMBOL != "TP53" &
               .data$SYMBOL != "ATM" &
             is.na(.data$RMSK_HIT)) |
           (.data$CONSEQUENCE == "stop_lost" &
            .data$SYMBOL != "TP53" &
            has_pfam & !is.na(.data$PFAM_DOMAIN_NAME))) &
            abs(nchar(.data$REF) - nchar(.data$ALT)) != 3 ~ TRUE,
          TRUE ~ as.logical(.data$ACMG_PM4)
        ),
        ACMG_PM4_SUPP = dplyr::case_when(
          .data$ACMG_PM4 == TRUE ~ FALSE,
          .data$SYMBOL == "TP53" ~ FALSE,
          (((.data$CONSEQUENCE == "inframe_deletion" |
              .data$CONSEQUENCE == "inframe_insertion") &
              .data$SYMBOL != "TP53" &
              .data$SYMBOL != "ATM" &
              is.na(.data$RMSK_HIT)) |
             (.data$CONSEQUENCE == "stop_lost" &
                .data$SYMBOL != "TP53" &
                has_pfam & !is.na(.data$PFAM_DOMAIN_NAME))) &
            abs(nchar(.data$REF) - nchar(.data$ALT)) == 3 ~ TRUE,
          TRUE ~ as.logical(.data$ACMG_PM4_SUPP)
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
      dplyr::filter(
        !is.na(.data$HGVSP_query) &
          !is.na(.data$VAR_ID)) |>

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
          .data$GRANTHAM_DISTANCE >=
          .data$MIN_GRANTHAM_DISTANCE &
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
#' PP2
#'
#' @param var_calls variants in cancer predisposition genes
#' @param max_benign_missense_frac maximum tolerated
#' fraction of benign missense variants for gene
#' @param max_truncating_frac maximum tolerated
#' fraction of pathogenic truncating variants for gene
#' @return data.frame with ACMG evidence indicators
#' @export
#'
assign_PP2_evidence <- function(
    var_calls = NULL,
    max_benign_missense_frac = 0.1,
    max_truncating_frac = 0.5){

  ## Assign logical ACMG evidence indicator
  # ACMG_PP2 - Missense variant in a gene that has a relatively low rate
  # of benign missense variation and where missense variants are a
  # common mechanism of disease
  var_calls$ACMG_PP2 <- FALSE

  if(is.data.frame(var_calls) &
     "BENIGN_MISSENSE_FRAC" %in% colnames(var_calls) &
     "PATH_TRUNC_FRAC" %in% colnames(var_calls) &
     "CONSEQUENCE" %in% colnames(var_calls)){

    var_calls <- var_calls |>
      dplyr::mutate(
        ACMG_PP2 =
          dplyr::if_else(
            (is.na(.data$BENIGN_MISSENSE_FRAC) |
               .data$BENIGN_MISSENSE_FRAC <= max_benign_missense_frac) &
              (is.na(.data$PATH_TRUNC_FRAC) |
                 .data$PATH_TRUNC_FRAC < max_truncating_frac) &
              stringr::str_detect(
                .data$CONSEQUENCE, "^missense_variant"),
            TRUE, FALSE, FALSE
          )
      )
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
     "EFFECT_PREDICTIONS" %in% colnames(var_calls) &
     "TF_BINDING_SITE_VARIANT_INFO" %in% colnames(var_calls) &
     "TF_BINDING_SITE_VARIANT" %in% colnames(var_calls) &
     "N_INSILICO_TOLERATED" %in% colnames(var_calls) &
     "N_INSILICO_SPLICING_AFFECTED" %in% colnames(var_calls) &
     "N_INSILICO_SPLICING_NEUTRAL" %in% colnames(var_calls)){

    var_calls_pp3 <- var_calls |>
      ## Missense predictions
      dplyr::mutate(
        ACMG_PP3 =
          dplyr::if_else(
            (.data$N_INSILICO_CALLED >= dbnsfp_min_called &
              .data$N_INSILICO_DAMAGING >= dbnsfp_min_majority &
              .data$N_INSILICO_TOLERATED <= dbnsfp_max_minority &
              .data$N_INSILICO_SPLICING_NEUTRAL <= 1) |
              (.data$N_INSILICO_CALLED >= dbnsfp_min_called &
                 .data$N_INSILICO_DAMAGING >= dbnsfp_min_majority &
                 .data$N_INSILICO_TOLERATED <= dbnsfp_max_minority + 1 &
                 stringr::str_detect(
                   .data$EFFECT_PREDICTIONS, "alphamissense:D") &
                 .data$N_INSILICO_SPLICING_NEUTRAL <= 1)
              ,
            TRUE,
            FALSE,
            FALSE
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

      ## TF binding site disruption
      dplyr::mutate(
        ACMG_PP3 = dplyr::case_when(
          !is.na(TF_BINDING_SITE_VARIANT_INFO) &
            !is.na(.data$TF_BINDING_SITE_VARIANT) &
            stringr::str_detect(
              .data$TF_BINDING_SITE_VARIANT,
              "Overlap: critical") ~ TRUE,
          TRUE ~ as.logical(.data$ACMG_PP3)
        )
      ) |>
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

#' Function that checks if a variant is a purely intronic deletion
#' (i.e. deletion of intronic sequence that does not extend into exonic sequence
#' at either end)
#'
#' @param hgvsc HGVSc notation for variant
#' @return logical indicating whether variant is a purely intronic deletion
#'
#' @export
#'
is_purely_intronic_deletion <- function(hgvsc) {

  # Regex pattern breakdown:
  # c\.          - literal "c."
  # \d+          - exon position (required anchor)
  # [+-]\d+      - intronic offset (+ or - followed by digits)
  # _            - range separator
  # \d+          - exon position (required anchor)
  # [+-]\d+      - intronic offset (+ or - followed by digits)
  # del          - deletion type

  pattern <- "^[A-Z0-9_\\.]+:c\\.\\d+[+-]\\d+_\\d+[+-]\\d+del"

  grepl(pattern, hgvsc, perl = TRUE)
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
        dplyr::group_by(.data$VAR_ID, .data$ENTREZGENE) |>
        dplyr::summarise(
          MIN_DISTANCE_TO_PATHOGENIC = min(.data$distance_to_path), .groups = "drop") |>
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
#' @param min_gold_stars minimum number of gold stars for
#' considering a variant as known benign/pathogenic
#' @return list with data.frames for benign peptide changes,
#' pathogenic nucleotide sites, and pathogenic codon sites
#' @export
#'
known_path_benign_sites <- function(
    ref_data = NULL,
    min_gold_stars = 2){

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
    dplyr::filter(.data$GOLD_STARS >= min_gold_stars & .data$BENIGN == 1) |>
    dplyr::select(c("ENTREZGENE", "HGVSP", "BENIGN", "VAR_ID")) |>
    dplyr::filter(!is.na(.data$ENTREZGENE) & !is.na(.data$HGVSP)) |>
    dplyr::rename(BENIGN_PEPTIDE_CHANGE = "BENIGN",
                  BENIGN_VAR_ID = "VAR_ID") |>
    dplyr::distinct()

  known_sites[['pathogenic_nucleotide']] <-
    ref_data[['variant']][['clinvar_nuc_sites']] |>
    dplyr::filter(.data$GOLD_STARS >= min_gold_stars) |>
    dplyr::select(c("ENTREZGENE", "PATH_NUC_SITE",
                    "PATH_VAR_ID_NUC"))

  known_sites[['pathogenic_codon']] <-
    ref_data[['variant']][['clinvar_aa_sites']] |>
    dplyr::filter(.data$PATHOGENIC == 1) |>
    dplyr::select(
      c("ENTREZGENE", "CODON",
        "PATHOGENIC","GOLD_STARS",
        "VAR_ID","GRANTHAM_DISTANCE","HGVSP")) |>
    tidyr::separate_rows("VAR_ID", sep=";") |>
    dplyr::filter(
      !is.na(.data$ENTREZGENE) &
        !is.na(.data$CODON) &
        !is.na(.data$GRANTHAM_DISTANCE)) |>
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
    dplyr::filter(.data$MAX_GOLD_STARS >= min_gold_stars) |>
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

  non_clinvar_vars_maf_filtered <- data.frame()
  clinvar_vars <- data.frame()

  if(is.data.frame(var_calls) &
     "ASSERTION_AUTHORITY" %in% colnames(var_calls) &
     "POPMAX_AF_GNOMAD" %in% colnames(var_calls) &
     "CLASSIFICATION" %in% colnames(var_calls) &
     "CPSR_PATHOGENICITY_SCORE" %in% colnames(var_calls) &
     "CLINVAR_PHENOTYPE_CANCER" %in% colnames(var_calls) &
     NROW(var_calls) > 0){

      non_clinvar_vars <-
        var_calls |>
        dplyr::filter(
          .data$ASSERTION_AUTHORITY == "CPSR")

      non_clinvar_vars_maf_filtered <- non_clinvar_vars
      clinvar_vars <-
        var_calls |>
        dplyr::filter(
          .data$ASSERTION_AUTHORITY == "ClinVar")

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
          non_clinvar_vars_maf_filtered, clinvar_vars) |>
        dplyr::arrange(
          dplyr::desc(
            .data$CLASSIFICATION),
          dplyr::desc(.data$CLINVAR_PHENOTYPE_CANCER),
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

  if(isTRUE(conf[['other']][['show_noncoding']])){
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
      "CLINVAR_PHENOTYPE_CANCER" %in% colnames(
        var_calls)) {

    non_cancer_calls_to_exclude <-
      var_calls |>
      dplyr::filter(!is.na(.data$CLINVAR_MSID)) |>
      dplyr::filter(
        !is.na(.data$CLINVAR_PHENOTYPE_CANCER) &
          .data$CLINVAR_PHENOTYPE_CANCER == 0)

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
