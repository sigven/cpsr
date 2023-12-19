#' Function that gets the maximum number of rows across different
#' tier data frames in CPSR report
#'
#' @param cps_report CPSR report structure with tier data frames
#'
#' @return max_row_nr maximum number of rows
#'
#' @export
get_max_rows_pr_datatable <- function(cps_report) {
  max_row_nr <- 0
  if (!is.null(cps_report[["content"]][["snv_indel"]][["disp"]])) {
    for (c in c(
      "class1", "class2", "class3",
      "class4", "class5", "sf", "gwas"
    )) {
      if (NROW(cps_report[["content"]][["snv_indel"]][["disp"]][[c]]) == 0) {
        next
      }
      t1 <- cps_report[["content"]][["snv_indel"]][["disp"]][[c]]
      if (c == "sf" | c == "gwas") {
        if (nrow(t1) > max_row_nr) {
          max_row_nr <- nrow(t1)
        }
      } else {
        t1 <- cps_report[["content"]][["snv_indel"]][["disp"]][[c]]
        num_rows_clinvar <- t1 |>
          dplyr::filter(.data$CPSR_CLASSIFICATION_SOURCE == "ClinVar") |>
          nrow()
        num_rows_other <- t1 |>
          dplyr::filter(.data$CPSR_CLASSIFICATION_SOURCE == "CPSR_ACMG") |>
          nrow()
        if (num_rows_other > max_row_nr) {
          max_row_nr <- num_rows_other
        }
        if (num_rows_clinvar > max_row_nr) {
          max_row_nr <- num_rows_clinvar
        }
      }
    }
  }
  return(max_row_nr)
}

#' Function that counts insilico predictions of variant effects
#' (i.e. damaging/tolerated) from dbNSFP
#'
#' @param cpg_calls sample calls with dbNSFP annotations
#'
#' @return cpg_calls
#'
#' @export
get_insilico_prediction_statistics <- function(cpg_calls) {

  pcgrr::log4r_info(
    "Summarising insilico variant effect predictions (dbNSFP)")


  insilico_pathogenicity_pred_algos <-
    c(
      "DBNSFP_SIFT",
      "DBNSFP_PROVEAN",
      "DBNSFP_META_RNN",
      "DBNSFP_FATHMM",
      "DBNSFP_MUTATIONTASTER",
      "DBNSFP_DEOGEN2",
      "DBNSFP_PRIMATEAI",
      "DBNSFP_MUTATIONASSESSOR",
      "DBNSFP_FATHMM_MKL",
      "DBNSFP_M_CAP",
      "DBNSFP_LIST_S2",
      "DBNSFP_BAYESDEL_ADDAF",
      "DBNSFP_SPLICE_SITE_ADA",
      "DBNSFP_SPLICE_SITE_RF"
    )
  for (v in c(
    "CALLED",
    "DAMAGING",
    "TOLERATED",
    "SPLICING_NEUTRAL",
    "SPLICING_AFFECTED"
  )) {
    cpg_calls[, paste0("N_INSILICO_", v)] <- 0
  }

  for (algo in insilico_pathogenicity_pred_algos) {
    if (algo %in% colnames(cpg_calls)) {
      cpg_calls[!is.na(cpg_calls[, algo]) &
        cpg_calls[, algo] != ".", "N_INSILICO_CALLED"] <-
        cpg_calls[!is.na(cpg_calls[, algo]) &
          cpg_calls[, algo] != ".", "N_INSILICO_CALLED"] + 1
      cpg_calls[
        !is.na(cpg_calls[, algo]) &
          (cpg_calls[, algo] == "D" |
            cpg_calls[, algo] == "PD"),
        "N_INSILICO_DAMAGING"
      ] <-
        cpg_calls[
          !is.na(cpg_calls[, algo]) &
            (cpg_calls[, algo] == "D" |
              cpg_calls[, algo] == "PD"),
          "N_INSILICO_DAMAGING"
        ] + 1
      cpg_calls[
        !is.na(cpg_calls[, algo]) &
          cpg_calls[, algo] == "T",
        "N_INSILICO_TOLERATED"
      ] <-
        cpg_calls[
          !is.na(cpg_calls[, algo]) &
            cpg_calls[, algo] == "T",
          "N_INSILICO_TOLERATED"
        ] + 1
      cpg_calls[
        !is.na(cpg_calls[, algo]) &
          cpg_calls[, algo] == "AS",
        "N_INSILICO_SPLICING_AFFECTED"
      ] <-
        cpg_calls[
          !is.na(cpg_calls[, algo]) &
            cpg_calls[, algo] == "AS",
          "N_INSILICO_SPLICING_AFFECTED"
        ] + 1
      cpg_calls[
        !is.na(cpg_calls[, algo]) &
          cpg_calls[, algo] == "SN",
        "N_INSILICO_SPLICING_NEUTRAL"
      ] <-
        cpg_calls[
          !is.na(cpg_calls[, algo]) &
            cpg_calls[, algo] == "SN",
          "N_INSILICO_SPLICING_NEUTRAL"
        ] + 1
    }
  }
  return(cpg_calls)
}

#' Function that retrieves variants in genes recommended for secondary
#' findings
#'
#' @param calls data frame with variants in genes recommended for
#' secondary findings reporting
#' @param umls_map data frame with UMLS phenotype terms
#'
#' @export
retrieve_secondary_calls <- function(calls, umls_map) {
  assertable::assert_colnames(
    calls,
    colnames = c(
      "CPG_SOURCE",
      "CLINVAR_CLASSIFICATION",
      "CLINVAR_UMLS_CUI",
      "GENOTYPE",
      "SYMBOL",
      "VAR_ID",
      "LOSS_OF_FUNCTION",
      "PROTEIN_CHANGE"
    ),
    only_colnames = F, quiet = T
  )

  secondary_calls <- calls |>
    ## do not consider secondary variant findings if
    ## genotypes have not been retrieved properly
    dplyr::filter(
      !is.na(.data$GENOTYPE) &
        !is.na(.data$SYMBOL) &
        !is.na(.data$CPG_SOURCE) &
        stringr::str_detect(.data$CPG_SOURCE,"ACMG_SF") &
        !is.na(.data$CLINVAR_CLASSIFICATION) &
        stringr::str_detect(
          .data$CLINVAR_CLASSIFICATION,"Pathogenic")) |>
    dplyr::filter(
      .data$GENOTYPE != "undefined"
    )


  if (NROW(secondary_calls) == 0) {
    return(secondary_calls)
  }

  secondary_calls <- secondary_calls |>
    ## only LOF for TTN
    dplyr::filter(.data$SYMBOL != "TTN" |
      (.data$SYMBOL == "TTN" &
        .data$LOSS_OF_FUNCTION == T)) |>
    ## only homozygotes p.Cys282Tyr HFE carriers
    dplyr::filter(.data$SYMBOL != "HFE" |
      (.data$SYMBOL == "HFE" &
        (!is.na(.data$PROTEIN_CHANGE) &
          stringr::str_detect(.data$PROTEIN_CHANGE, "Cys282Tyr")) &
        .data$GENOTYPE == "homozygous"))

  ## AR genes
  min_two_variants_required <-
    data.frame(SYMBOL = "MUTYH", stringsAsFactors = F) |>
    dplyr::bind_rows(
      data.frame(SYMBOL = "CASQ2", stringsAsFactors = F),
      data.frame(SYMBOL = "TRDN", stringsAsFactors = F),
      data.frame(SYMBOL = "RPE65", stringsAsFactors = F),
      data.frame(SYMBOL = "GAA", stringsAsFactors = F),
      data.frame(SYMBOL = "BTD", stringsAsFactors = F),
      data.frame(SYMBOL = "ATP7B", stringsAsFactors = F)
    )

  if (NROW(secondary_calls) > 0) {
    ## Skip MUTYH/RPE65/GAA/BTD/ATP7B if they only occur with one variant
    genes_lacking_twohit_evidence <- secondary_calls |>
      dplyr::group_by(.data$SYMBOL) |>
      dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
      dplyr::filter(.data$n == 1)

    if (nrow(genes_lacking_twohit_evidence) > 0) {
      genes_lacking_twohit_evidence <- genes_lacking_twohit_evidence |>
        dplyr::inner_join(min_two_variants_required, by = "SYMBOL")

      if (nrow(genes_lacking_twohit_evidence) > 0) {
        secondary_calls <- secondary_calls |>
          dplyr::anti_join(
            genes_lacking_twohit_evidence, by = "SYMBOL")
      }
    }
  }

  return(secondary_calls)
}

#' Function that retrieves variants in cancer predisposition genes linked
#' to cancer-related conditions according to ClinVar
#'
#' @param cpg_calls data frame with variants in predisposition_genes
#' @param ref_data object with PCGR/CPSR reference data
#'
#' @export
check_variant2cancer_phenotype <- function(cpg_calls, ref_data) {

  oncotree <- ref_data[['phenotype']][['oncotree']]
  umls_map <- ref_data[['phenotype']][['umls']] |>
    dplyr::filter(.data$MAIN_TERM == TRUE) |>
    dplyr::select(-c("SOURCE")) |>
    dplyr::distinct()

  if (nrow(cpg_calls) > 0 &
    "CLINVAR_UMLS_CUI" %in% colnames(cpg_calls) &
    "VAR_ID" %in% colnames(cpg_calls)) {
    oncotree <- oncotree |>
      dplyr::select("CUI") |>
      dplyr::mutate(CANCER_PHENOTYPE = 1) |>
      dplyr::distinct()

    n_clinvar <- cpg_calls |>
      dplyr::filter(!is.na(.data$CLINVAR_UMLS_CUI)) |>
      nrow()

    if (n_clinvar > 0) {
      cpg_calls_traits <- as.data.frame(
        tidyr::separate_rows(
          cpg_calls, "CLINVAR_UMLS_CUI", sep = ",") |>
          dplyr::select(c("VAR_ID", "CLINVAR_UMLS_CUI")) |>
          dplyr::left_join(
            umls_map, by = c("CLINVAR_UMLS_CUI" = "CUI"),
            relationship = "many-to-many") |>
          dplyr::distinct() |>
          dplyr::filter(!is.na(.data$CUI_NAME)) |>
          dplyr::left_join(
            oncotree, by = c("CLINVAR_UMLS_CUI" = "CUI"),
            relationship = "many-to-many") |>
          dplyr::mutate(
            CANCER_PHENOTYPE = dplyr::if_else(
              is.na(.data$CANCER_PHENOTYPE),
              as.integer(0),
              as.integer(.data$CANCER_PHENOTYPE)
            )
          ) |>
          dplyr::mutate(
            CANCER_PHENOTYPE =
              dplyr::if_else(
                !is.na(.data$CUI_NAME) &
                  stringr::str_detect(
                tolower(.data$CUI_NAME),
                pcgrr::cancer_phenotypes_regex
              ),
              as.integer(1),
              as.integer(.data$CANCER_PHENOTYPE)
              )
          ) |>
          dplyr::group_by(.data$VAR_ID) |>
          dplyr::summarise(
            CLINVAR_PHENOTYPE = paste(
              unique(.data$CUI_NAME), collapse = "; "),
            CANCER_PHENOTYPE = max(.data$CANCER_PHENOTYPE),
            .groups = "drop"
          ) |>
          dplyr::mutate(
            CANCER_PHENOTYPE =
              dplyr::if_else(
                stringr::str_detect(
                  .data$CLINVAR_PHENOTYPE,
                  "^(not specified; not provided|not specified|not provided)"
                ),
                as.integer(1),
                as.integer(.data$CANCER_PHENOTYPE)
              )
          )
      )

      cpg_calls <- cpg_calls |>
        dplyr::left_join(cpg_calls_traits, by = "VAR_ID")
    } else {
      cpg_calls$CLINVAR_PHENOTYPE <- NA
      cpg_calls$CANCER_PHENOTYPE <- NA
    }
  }
  return(cpg_calls)
}

#' Create directory
#'
#' @param d Directory to create.
#'
#' @export
mkdir <- function(d) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
  }
  TRUE
}

