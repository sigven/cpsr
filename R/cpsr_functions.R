#' Function that generates variant predisposition report - CPSR
#'
#' @param yaml_fname YAML file with paths to pre-processed
#' molecular input files, configuration settings, metadata for
#' reference data bundle etc.
#' @export

generate_cpsr_report <- function(yaml_fname = NULL) {
  invisible(assertthat::assert_that(
    !is.null(yaml_fname),
    msg = "Object 'yaml_fname' cannot be NULL"
  ))
  pcgrr::check_file_exists(yaml_fname)

  cps_report <- pcgrr::init_report(
    yaml_fname = yaml_fname,
    report_mode = "CPSR")

  callset_cpsr <-
    pcgrr::load_germline_snv_indel(
      fname = cps_report$settings$molecular_data$fname_mut_tsv,
      ref_data = cps_report$ref_data,
      settings = cps_report$settings
    )

  col_format_output <- cpsr::col_format_output
  for(f in c('html_tier','tsv')){
    col_format_output[[f]] <- c(
      col_format_output[[f]],
      cps_report[["settings"]][["conf"]][["variant_classification"]][["vcftag_gnomad_AF"]])
  }

  if(callset_cpsr$retained_info_tags != ""){
    tags <- stringr::str_split(
      callset_cpsr$retained_info_tags, pattern=",")[[1]]
    for(t in tags){
      if(t %in% colnames(callset_cpsr$variant)){
        col_format_output[['tsv']] <- c(
          col_format_output[['tsv']],
          t
        )
        col_format_output[['html_tier']] <- c(
          col_format_output[['html_tier']],
          t
        )
      }else{
        #issue warning
      }
    }
  }

  primary_cpg_targets <-
    cps_report$settings$conf$gene_panel$panel_genes |>
    dplyr::filter(.data$PRIMARY_TARGET == T) |>
    dplyr::mutate(
      ENTREZGENE = as.character(.data$ENTREZGENE))

  pcgrr::log4r_info("------")
  pcgrr::log4r_info(
    paste0(
      "Considering variants in the targeted predisposition genes: ",
      paste(unique(sort(primary_cpg_targets$SYMBOL)),
            collapse = ", "
      )
    )
  )

  variant_calls <- callset_cpsr[['variant']]
  if (cps_report$settings$conf$other$show_noncoding == F) {
    n_noncoding_vars <- variant_calls |>
      dplyr::filter(.data$CODING_STATUS == "noncoding") |>
      NROW()
    pcgrr::log4r_info(
      paste0(
        "Excluding n = ",
        n_noncoding_vars,
        " variants for classification (option --ignore_noncoding)"
      )
    )
    variant_calls <- variant_calls |>
      dplyr::filter(.data$CODING_STATUS == "coding")
    if (NROW(variant_calls) == 0) {
      pcgrr::log4r_warn(paste0(
        "There are zero remaining protein-coding ",
        "calls in input file - no report will be produced"
      ))
      return(NULL)
    }
  }
  ## get overall call statistics
  call_stats <-
    pcgrr::variant_stats_report(variant_calls, name = "v_stat")

  cpg_calls <-
    dplyr::semi_join(
      variant_calls,
      primary_cpg_targets,
      by = c("ENTREZGENE" = "ENTREZGENE")
    )
  cpg_call_stats <-
    pcgrr::variant_stats_report(
      cpg_calls,
      name = "v_stat_cpg"
    )

  pcgrr::log4r_info(
    paste0(
      "Total number of variants in target cancer predisposition genes (for TIER output): ",
      cpg_call_stats[["v_stat_cpg"]][["n"]]
    )
  )
  pcgrr::log4r_info(
    paste0(
      "Number of coding variants in target cancer predisposition genes (for TIER output): ",
      cpg_call_stats[["v_stat_cpg"]][["n_coding"]]
    )
  )
  pcgrr::log4r_info(
    paste0(
      "Number of non-coding variants in cancer predisposition genes (for TIER output): ",
      cpg_call_stats[["v_stat_cpg"]][["n_noncoding"]]
    )
  )

  if (NROW(cpg_calls) == 0) {
    cps_report$content$snv_indel$v_stat_cpg <-
      cpg_call_stats$v_stat_cpg
    return(cps_report)
  }

  if (cps_report$settings$conf$variant_classification$clinvar_ignore_noncancer == 1 &
      "VAR_ID" %in% colnames(cpg_calls) &
      "CLINVAR_MSID" %in% colnames(cpg_calls) &
      "CANCER_PHENOTYPE" %in% colnames(cpg_calls)) {
    n_clinvar_noncancer <- cpg_calls |>
      dplyr::filter(!is.na(.data$CLINVAR_MSID)) |>
      dplyr::filter(.data$CANCER_PHENOTYPE == 0) |>
      NROW()

    pcgrr::log4r_info(
      paste0(
        "ClinVar variants related to non-cancer conditions excluded",
        " from report: ", n_clinvar_noncancer
      )
    )

    if (n_clinvar_noncancer > 0) {
      cpg_calls_exclude <- cpg_calls |>
        dplyr::filter(!is.na(.data$CLINVAR_MSID)) |>
        dplyr::filter(.data$CANCER_PHENOTYPE == 0) |>
        dplyr::select("VAR_ID")

      cpg_calls <- cpg_calls |>
        dplyr::anti_join(cpg_calls_exclude, by = "VAR_ID")

      pcgrr::log4r_info(
        paste0(
          "Variants remaining after exclusion of non-cancer related",
          " ClinVar variants: ", NROW(cpg_calls)
        )
      )
    }
  }

  if (NROW(cpg_calls) == 0) {
    cps_report$content$snv_indel$v_stat_cpg <-
      cpg_call_stats$v_stat_cpg
    return(cps_report)
  }

  ## Assign calls to tiers (ClinVar calls + CPSR classification
  ## for novel, non-ClinVar variants)
  snv_indel_report <-
    cpsr::assign_variant_tiers(
      cpg_calls,
      config = cps_report[["settings"]][["conf"]],
      col_format_output = col_format_output
    )

  snv_indel_report$v_stat <-
    call_stats$v_stat
  snv_indel_report$v_stat_cpg <-
    cpg_call_stats$v_stat_cpg

  #snv_indel_report$clin_eitem <-
  #  callset_cpsr$biomarker_evidence
  #  cps_report$content$snv_indel$clin_eitem

  if(!is.null(snv_indel_report[['clin_eitem']])){
    for(type in names(snv_indel_report[['clin_eitem']])){
      if(!(type %in% names(callset_cpsr[['biomarker_evidence']]))){
        next
      }
      for(level in names(snv_indel_report[['clin_eitem']][[type]])){
        if(!(level %in% names(callset_cpsr[['biomarker_evidence']][[type]]))){
          next
        }
        if(NROW(callset_cpsr[['biomarker_evidence']][[type]][[level]]) > 0){
          snv_indel_report[['clin_eitem']][[type]][[level]] <-
            callset_cpsr[['biomarker_evidence']][[type]][[level]] |>
            dplyr::select(-c("PRIMARY_SITE")) |>
            dplyr::select(c(
              "SYMBOL",
              "CONSEQUENCE",
              "PROTEIN_CHANGE",
              "EVIDENCE_LEVEL",
              "CANCER_TYPE",
              "CLINICAL_SIGNIFICANCE",
              "RATING",
              "EVIDENCE_DIRECTION",
              "CITATION",
              "THERAPEUTIC_CONTEXT",
              "EVIDENCE_TYPE",
              "EVIDENCE_DESCRIPTION",
              "BIOMARKER_MAPPING",
              "OFFICIAL_GENENAME",
              "CDS_CHANGE",
              "LOSS_OF_FUNCTION",
              "GENOMIC_CHANGE"),
              dplyr::everything()
            ) |>
            dplyr::distinct()
        }
      }
    }
  }

  cps_report <-
    pcgrr::update_report(cps_report, report_data = snv_indel_report)

  gene_hits <- paste(
    unique(
      cps_report[["content"]][["snv_indel"]][["variant_set"]][["tsv"]]$SYMBOL
    ),
    collapse = ", "
  )
  pcgrr::log4r_info(paste0(
    "Variants were found in the following cancer ",
    "predisposition genes: ", gene_hits
  ))

  ## secondary findings
  if (cps_report$settings$conf$variant_classification$secondary_findings == TRUE) {

    secondary_calls <-
      cpsr::retrieve_secondary_calls(
        variant_calls,
        umls_map = cps_report$ref_data$phenotype$umls)

    if(NROW(secondary_calls) > 0){
      secondary_calls <- secondary_calls |>
        dplyr::anti_join(
          cpg_calls, by = "VAR_ID"
        )
    }
    secondary_call_stats <-
      pcgrr::variant_stats_report(
        secondary_calls,
        name = "v_stat_secondary"
      )
    pcgrr::log4r_info(paste0(
      "Assignment of other variants in genes ",
      "recommended for reporting as secondary ",
      "findings (ACMG SF v3.0)"
    ))
    if (NROW(secondary_calls) > 0) {
      cps_report[["content"]][["snv_indel"]][["disp"]][["secondary"]] <-
        secondary_calls |>
        dplyr::arrange(
          .data$LOSS_OF_FUNCTION, .data$CODING_STATUS) |>
        dplyr::select(
          dplyr::one_of(col_format_output[['html_sf']]))
    }
    pcgrr::log4r_info(paste0(
      "Number of pathogenic variants in the ACMG secondary findings list - other ",
      "genes of clinical significance: ",
      cps_report[["content"]][["snv_indel"]][["v_stat_secondary"]][["n_coding"]]
    ))
  }

  cps_report[["content"]][["snv_indel"]][["eval"]] <- TRUE

  if (cps_report$settings$conf$variant_classification$gwas_findings == 1) {
    pcgrr::log4r_info(paste0(
      "Assignment of other variants to hits ",
      "from genome-wide association studies"
    ))

    ## Assign GWAS hits to cps_report object
    cps_report[["content"]][["snv_indel"]][["disp"]][["gwas"]] <-
      dplyr::filter(
        variant_calls,
        !is.na(.data$GWAS_HIT) &
          !is.na(.data$GWAS_CITATION))
    if (NROW(cps_report[["content"]][["snv_indel"]][["disp"]][["gwas"]]) > 0) {
      if (NROW(cps_report[["content"]][["snv_indel"]][["variant_set"]][["tsv"]]) > 0) {
        ## Omit variants that is already present in TIER 1-5
        cps_report[["content"]][["snv_indel"]][["disp"]][["gwas"]] <-
          cps_report[["content"]][["snv_indel"]][["disp"]][["gwas"]] |>
          dplyr::anti_join(
            cps_report[["content"]][["snv_indel"]][["variant_set"]][["tsv"]],
            by = c("GENOMIC_CHANGE")
          ) |>
          dplyr::select(
            dplyr::one_of(col_format_output[['html_gwas']]))
      }
      ## Select variables to include for GWAS hits and arrange results
      if (NROW(cps_report[["content"]][["snv_indel"]][["disp"]][["gwas"]]) > 0) {
        cps_report[["content"]][["snv_indel"]][["disp"]][["gwas"]] <-
          cps_report[["content"]][["snv_indel"]][["disp"]][["gwas"]] |>
          dplyr::arrange(
            .data$LOSS_OF_FUNCTION,
            .data$CODING_STATUS) |>
          dplyr::select(
            dplyr::one_of(col_format_output[['html_gwas']]))
      }
    }
  }

  cps_report[["content"]][["snv_indel"]][["max_dt_rows"]] <-
    get_max_rows_pr_datatable(cps_report)

  return(cps_report)
}


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
          dplyr::filter(.data$CPSR_CLASSIFICATION_SOURCE == "Other") |>
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

#' Function that assign variants to different tiers for
#' prioritization of germline variants
#'
#' @param cpg_calls data frame with variants in predisposition_genes
#' @param config CPSR configuration object with run settings
#' @param col_format_output object with output column formats (html, tsv)
#' @export
assign_variant_tiers <-
  function(cpg_calls,
           config = NULL,
           col_format_output = NULL) {
    # dot_args <- list(...)

    # evidence_codes <- cpsr::acmg[["evidence_codes"]] |>
    #   dplyr::filter(
    #     .data$cpsr_evidence_code != "ACMG_BS2_1" &
    #     .data$cpsr_evidence_code != "ACMG_BS2_2" &
    #     .data$cpsr_evidence_code != "ACMG_BS2_3")
    col_format_output[['tsv']] <-
      c(
        col_format_output[['tsv']],
        #evidence_codes$cpsr_evidence_code,
        c(
          "FINAL_CLASSIFICATION",
          "CPSR_CLASSIFICATION",
          "CPSR_PATHOGENICITY_SCORE",
          "CPSR_CLASSIFICATION_CODE",
          "CPSR_CLASSIFICATION_DOC",
          "CPSR_CLASSIFICATION_SOURCE"
        )
      )

    pcgrr::log4r_info(paste0(
      "Generating tiered set of result variants for ",
      "output in tab-separated values (TSV) file"
    ))

    snv_indel_report <- pcgrr::init_germline_content()

    snv_indel_report[["variant_set"]][["class5"]] <- cpg_calls |>
      dplyr::filter(!is.na(.data$CLINVAR_CLASSIFICATION) &
        .data$CLINVAR_CLASSIFICATION == "Pathogenic") |>
      dplyr::mutate(CPSR_CLASSIFICATION_SOURCE = "ClinVar")

    snv_indel_report[["variant_set"]][["class4"]] <- cpg_calls |>
      dplyr::filter(!is.na(.data$CLINVAR_CLASSIFICATION) &
        .data$CLINVAR_CLASSIFICATION == "Likely_Pathogenic") |>
      dplyr::mutate(CPSR_CLASSIFICATION_SOURCE = "ClinVar")

    snv_indel_report[["variant_set"]][["class3"]] <- cpg_calls |>
      dplyr::filter(!is.na(.data$CLINVAR_CLASSIFICATION) &
        .data$CLINVAR_CLASSIFICATION == "VUS") |>
      dplyr::mutate(CPSR_CLASSIFICATION_SOURCE = "ClinVar")

    snv_indel_report[["variant_set"]][["class2"]] <- cpg_calls |>
      dplyr::filter(!is.na(.data$CLINVAR_CLASSIFICATION) &
        .data$CLINVAR_CLASSIFICATION == "Likely_Benign") |>
      dplyr::mutate(CPSR_CLASSIFICATION_SOURCE = "ClinVar")

    snv_indel_report[["variant_set"]][["class1"]] <- cpg_calls |>
      dplyr::filter(!is.na(.data$CLINVAR_CLASSIFICATION) &
        .data$CLINVAR_CLASSIFICATION == "Benign") |>
      dplyr::mutate(CPSR_CLASSIFICATION_SOURCE = "ClinVar")

    ## identify remaining calls not registered in ClinVar
    all_clinvar_calls <- data.frame()
    for (c in c("class1", "class2", "class3", "class4", "class5")) {
      all_clinvar_calls <- all_clinvar_calls |>
        dplyr::bind_rows(
          dplyr::select(
            snv_indel_report[["variant_set"]][[c]],
            "VAR_ID"
          )
        )
    }
    cpg_calls_non_clinvar <- cpg_calls |>
      dplyr::anti_join(all_clinvar_calls, by = c("VAR_ID"))

    n_nonclinvar <- NROW(cpg_calls_non_clinvar)

    cpg_calls_non_clinvar <- cpg_calls_non_clinvar |>
      dplyr::filter(is.na(.data$gnomAD_AF) |
        .data$gnomAD_AF <=
          config[["variant_classification"]][["maf_upper_threshold"]])
    n_maf_filtered <- n_nonclinvar - nrow(cpg_calls_non_clinvar)
    pcgrr::log4r_info(
      paste0(
        "Ignoring n = ", n_maf_filtered,
        " unclassified variants with a global MAF frequency above ",
        config[["variant_classification"]][["maf_upper_threshold"]]
      )
    )

    n_after_maf_filtering <- nrow(cpg_calls_non_clinvar)

    non_clinvar_calls <- list()
    non_clinvar_calls[["class5"]] <- cpg_calls_non_clinvar |>
      dplyr::filter(.data$CPSR_CLASSIFICATION == "Pathogenic") |>
      dplyr::mutate(CPSR_CLASSIFICATION_SOURCE = "Other")

    non_clinvar_calls[["class4"]] <- cpg_calls_non_clinvar |>
      dplyr::filter(.data$CPSR_CLASSIFICATION == "Likely_Pathogenic") |>
      dplyr::mutate(CPSR_CLASSIFICATION_SOURCE = "Other")

    non_clinvar_calls[["class3"]] <- cpg_calls_non_clinvar |>
      dplyr::filter(.data$CPSR_CLASSIFICATION == "VUS") |>
      dplyr::mutate(CPSR_CLASSIFICATION_SOURCE = "Other")

    non_clinvar_calls[["class2"]] <- cpg_calls_non_clinvar |>
      dplyr::filter(.data$CPSR_CLASSIFICATION == "Likely_Benign") |>
      dplyr::mutate(CPSR_CLASSIFICATION_SOURCE = "Other")

    non_clinvar_calls[["class1"]] <- cpg_calls_non_clinvar |>
      dplyr::filter(.data$CPSR_CLASSIFICATION == "Benign") |>
      dplyr::mutate(CPSR_CLASSIFICATION_SOURCE = "Other")

    for (c in c("class1", "class2", "class3", "class4", "class5")) {
      pcgrr::log4r_info(paste0("Merging ClinVar-classified variants and CPSR-classified (novel) variants - ", c))
      snv_indel_report[["variant_set"]][[c]] <-
        dplyr::bind_rows(
          non_clinvar_calls[[c]],
          snv_indel_report[["variant_set"]][[c]]
        )


      if (nrow(snv_indel_report[["variant_set"]][[c]]) == 0) {
        pcgrr::log4r_info(paste0("Zero variants found - ", c))
        next
      }

      ## set FINAL_CLASSIFICATION col
      snv_indel_report[["variant_set"]][[c]] <-
        snv_indel_report[["variant_set"]][[c]] |>
        dplyr::mutate(
          FINAL_CLASSIFICATION = dplyr::case_when(
            !is.na(.data$CLINVAR_CLASSIFICATION) ~
              as.character(.data$CLINVAR_CLASSIFICATION),
            is.na(.data$CLINVAR_CLASSIFICATION) ~
              as.character(.data$CPSR_CLASSIFICATION),
            TRUE ~ as.character(NA)
          )
        )


      ## If not 'classify_all' is turned on,
      ## remove CPSR classifications for existing
      ## ClinVar classifications
      if (config[["variant_classification"]][["classify_all"]] == 1) {
        snv_indel_report[["variant_set"]][[c]] <-
          snv_indel_report[["variant_set"]][[c]] |>
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
      snv_indel_report[["variant_set"]][[c]] <-
        snv_indel_report[["variant_set"]][[c]] |>
        dplyr::arrange(
          .data$CPSR_CLASSIFICATION_SOURCE,
          dplyr::desc(.data$CANCER_PHENOTYPE),
          dplyr::desc(.data$CPSR_PATHOGENICITY_SCORE)
        )

      # if (config[["visual_reporting"]][["table_display"]] == "full") {
      #   cols_in_tsv_not_display_cols <-
      #     setdiff(cpsr_tsv_cols, cpsr_display_cols)
      #
      #   cpsr_display_cols <-
      #     c(
      #       cpsr_display_cols,
      #       cols_in_tsv_not_display_cols
      #     )
      #
      #   #cpsr_display_cols <-
      #     #cpsr_display_cols[!cpsr_display_cols %in% c("AMINO_ACID_START", "AMINO_ACID_END", "EXON", "CIVIC_ID", "CIVIC_ID_SEGMENT")]
      # }

      snv_indel_report[["disp"]][[c]] <-
        dplyr::select(
          snv_indel_report[["variant_set"]][[c]],
          dplyr::one_of(col_format_output[['html_tier']])
        )
      snv_indel_report[["variant_set"]][[c]] <-
        dplyr::select(
          snv_indel_report[["variant_set"]][[c]],
          dplyr::one_of(col_format_output[['tsv']])
        )

      snv_indel_report[["variant_set"]][[c]]$DBSNP <-
        unlist(lapply(
          stringr::str_match_all(
            snv_indel_report[["variant_set"]][[c]]$DBSNP,
            ">rs[0-9]{1,}<"
          ), paste,
          collapse = ","
        ))
      snv_indel_report[["variant_set"]][[c]]$DBSNP <-
        stringr::str_replace_all(
          snv_indel_report[["variant_set"]][[c]]$DBSNP, ">|<", ""
        )

      snv_indel_report[["variant_set"]][[c]]$CPSR_CLASSIFICATION_DOC <-
        stringr::str_replace_all(
          snv_indel_report[["variant_set"]][[c]]$CPSR_CLASSIFICATION_DOC,
          "<br>-", ","
        )
      snv_indel_report[["variant_set"]][[c]]$CPSR_CLASSIFICATION_DOC <-
        stringr::str_replace_all(
          snv_indel_report[["variant_set"]][[c]]$CPSR_CLASSIFICATION_DOC,
          "^, ", ""
        )

      snv_indel_report[["variant_set"]][[c]] <-
        snv_indel_report[["variant_set"]][[c]] |>
        dplyr::select(
          c(
            "GENOMIC_CHANGE",
            "VAR_ID",
            "GENOTYPE",
            "CPSR_CLASSIFICATION_SOURCE",
            "GENOME_VERSION",
            "VCF_SAMPLE_ID",
            "VARIANT_CLASS",
            "CODING_STATUS",
            "SYMBOL",
            "GENENAME",
            "CCDS",
            "ENTREZGENE",
            "UNIPROT_ID",
            "ENSEMBL_GENE_ID",
            "ENSEMBL_TRANSCRIPT_ID",
            "REFSEQ_TRANSCRIPT_ID",
            "ONCOGENE",
            "TUMOR_SUPPRESSOR",
            "CONSEQUENCE",
            "VEP_ALL_CSQ",
            "REGULATORY_ANNOTATION",
            "PROTEIN_CHANGE",
            "PFAM_DOMAIN_NAME",
            "DBSNP",
            "HGVSp",
            "HGVSc",
            "LAST_EXON",
            "EXON_POSITION",
            "INTRON_POSITION",
            "CDS_CHANGE",
            "MUTATION_HOTSPOT",
            "RMSK_HIT",
            "EFFECT_PREDICTIONS",
            "LOSS_OF_FUNCTION",
            "DBSNP",
            "CANCER_PHENOTYPE",
            "CLINVAR_CLASSIFICATION",
            "CLINVAR_MSID",
            "CLINVAR_VARIANT_ORIGIN",
            "CLINVAR_CONFLICTED",
            "CLINVAR_PHENOTYPE",
            "CLINVAR_REVIEW_STATUS_STARS"
          ),
          dplyr::everything()
        )

      for (col in colnames(snv_indel_report[["variant_set"]][[c]])) {
        if (nrow(snv_indel_report[["variant_set"]][[c]][!is.na(
          snv_indel_report[["variant_set"]][[c]][, col]
        ) &
          snv_indel_report[["variant_set"]][[c]][, col] == "", ]) > 0) {
          snv_indel_report[["variant_set"]][[c]][!is.na(
            snv_indel_report[["variant_set"]][[c]][, col]
          ) &
            snv_indel_report[["variant_set"]][[c]][, col] == "", col] <- NA
        }
      }

      population_tags <- unique(
        c("gnomAD_AF",
          config[["variant_classification"]][["vcftag_gnomad_AF"]])
      )
      for (tag in population_tags) {
        if (tag %in% colnames(snv_indel_report[["disp"]][[c]])) {
          if (nrow(snv_indel_report[["disp"]][[c]][is.na(
            snv_indel_report[["disp"]][[c]][, tag]
          ), ]) > 0) {
            snv_indel_report[["disp"]][[c]][is.na(
              snv_indel_report[["disp"]][[c]][, tag]
            ), tag] <- 0.00
          }
        }
      }
    }

    snv_indel_report[["variant_set"]][["tsv"]] <-
      dplyr::bind_rows(
        snv_indel_report[["variant_set"]][["class5"]],
        snv_indel_report[["variant_set"]][["class4"]],
        snv_indel_report[["variant_set"]][["class3"]],
        snv_indel_report[["variant_set"]][["class2"]],
        snv_indel_report[["variant_set"]][["class1"]]
      )


    return(snv_indel_report)
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
    dplyr::filter(.data$MAIN_TERM == TRUE)

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
        tidyr::separate_rows(cpg_calls, .data$CLINVAR_UMLS_CUI, sep = ",") |>
          dplyr::select(.data$VAR_ID, .data$CLINVAR_UMLS_CUI) |>
          dplyr::left_join(umls_map, by = c("CLINVAR_UMLS_CUI" = "CUI")) |>
          dplyr::distinct() |>
          dplyr::filter(!is.na(.data$CUI_NAME)) |>
          dplyr::left_join(oncotree, by = c("CLINVAR_UMLS_CUI" = "CUI")) |>
          dplyr::mutate(
            CANCER_PHENOTYPE = dplyr::if_else(
              is.na(.data$CANCER_PHENOTYPE),
              as.integer(0),
              as.integer(.data$CANCER_PHENOTYPE)
            )
          ) |>
          dplyr::mutate(
            CANCER_PHENOTYPE =
              dplyr::if_else(stringr::str_detect(
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

#' Function that makes a piechart showing the number of variants at
#' each significance level
#'
#' @param variants_tsv data frame with variants in predisposition_genes
#' @param plot_type ClinVar or Other
#'
#' @export
summary_donut_chart <- function(variants_tsv, plot_type = "ClinVar") {
  title <- "ClinVar variants"
  p <- NULL

  if (nrow(variants_tsv) > 0) {
    set_clinvar <- variants_tsv |>
      dplyr::filter(!is.na(.data$CLINVAR_CLASSIFICATION) &
        !(.data$CLINVAR_CLASSIFICATION == "NA"))
    set_other <- variants_tsv |>
      dplyr::filter(nchar(.data$CPSR_CLASSIFICATION) > 0 &
        (is.na(.data$CLINVAR_CLASSIFICATION) |
          .data$CLINVAR_CLASSIFICATION == "NA"))

    if ((plot_type == "ClinVar" & nrow(set_clinvar) > 0) |
      (plot_type != "ClinVar" & nrow(set_other) > 0)) {
      m <- data.frame()

      if (plot_type == "ClinVar") {
        if (nrow(set_clinvar) > 0) {
          t <- paste0("n = ", nrow(set_clinvar))
          title <- bquote("ClinVar variants, " ~ bold(.(t)))
          m <- as.data.frame(set_clinvar |>
            dplyr::group_by(.data$CLINVAR_CLASSIFICATION) |>
            dplyr::summarise(n = dplyr::n()) |>
            dplyr::rename(level = .data$CLINVAR_CLASSIFICATION)) |>
            dplyr::mutate(
              level =
                factor(
                  .data$level,
                  levels =
                    pcgrr::color_palette[["pathogenicity"]][["levels"]]
                )
            ) |>
            dplyr::arrange(.data$level) |>
            dplyr::mutate(prop = as.numeric(.data$n / sum(.data$n))) |>
            dplyr::mutate(lab.ypos = cumsum(.data$prop) - 0.5 * .data$prop) |>
            dplyr::mutate(n = as.character(.data$n))
        }
      } else {
        if (nrow(set_other) > 0) {
          t <- paste0("n = ", nrow(set_other))
          title <- bquote("Other variants, CPSR-classified, " ~ bold(.(t)))
          m <- as.data.frame(set_other |>
            dplyr::group_by(.data$CPSR_CLASSIFICATION) |>
            dplyr::summarise(n = dplyr::n()) |>
            dplyr::rename(level = .data$CPSR_CLASSIFICATION)) |>
            dplyr::mutate(
              level = factor(
                .data$level,
                levels = pcgrr::color_palette[["pathogenicity"]][["levels"]]
              )
            ) |>
            dplyr::arrange(.data$level) |>
            dplyr::mutate(prop = as.numeric(.data$n / sum(.data$n))) |>
            dplyr::mutate(lab.ypos = cumsum(.data$prop) - 0.5 * .data$prop) |>
            dplyr::mutate(n = as.character(.data$n))
        }
      }


      p <- ggplot2::ggplot(m, ggplot2::aes(x = 2, y = .data$prop, fill = .data$level)) +
        ggplot2::geom_bar(stat = "identity", color = "white") +
        ggplot2::coord_polar(theta = "y", start = 0) +
        ggplot2::geom_text(ggplot2::aes(y = 1 - .data$lab.ypos, label = .data$n),
          color = "white", family = "Helvetica", size = 6
        ) +
        ggplot2::scale_fill_manual(
          values = pcgrr::color_palette[["pathogenicity"]][["values"]],
          labels = pcgrr::color_palette[["pathogenicity"]][["levels"]],
          drop = F
        ) +
        ggplot2::theme_void() +
        ggplot2::xlim(0.5, 2.5) +
        ggplot2::ggtitle(title) +
        ggplot2::theme(
          plot.title =
            ggplot2::element_text(
              family = "Helvetica",
              size = 16, vjust = -1,
              hjust = 0.5
            ),
          legend.title = ggplot2::element_blank(),
          plot.margin = ggplot2::margin(0, 0, 0, 0, "cm"),
          legend.text = ggplot2::element_text(
            family = "Helvetica", size = 14
          )
        )
    }
  }
  return(p)
}


#' Function that retrieves clinical evidence items (CIVIC, CGI) for
#' somatic cancer variants
#'
#' @param sample_calls data frame with germline variant callset from
#' query sample
#' @param colset vector with column names to display for each report element
#' @param eitems all clinical evidence items linking germline variants with
#' impact on therapeutic response/diagnosis/prognosis etc
#'
#' @return list
#' @export
get_germline_biomarkers <- function(sample_calls,
                                    colset = NULL,
                                    eitems = NULL) {
  pcgrr::log4r_info(paste0(
    "Matching variant set with existing genomic ",
    "biomarkers from CIViC (germline)"
  ))

  clin_eitems_list <- list()
  for (type in c("diagnostic", "prognostic", "predictive", "predisposing")) {
    clin_eitems_list[[type]] <- data.frame()
  }

  all_var_evidence_items <- data.frame()
  var_eitems <- list()
  for (m in c("codon", "exon", "gene", "exact")) {
    var_eitems[[m]] <- data.frame()
  }

  sample_calls_p_lp <- sample_calls |>
    dplyr::filter(!is.na(.data$CPSR_CLASSIFICATION) |
      !is.na(.data$CLINVAR_CLASSIFICATION)) |>
    dplyr::filter(
      stringr::str_detect(
        .data$CPSR_CLASSIFICATION, "Pathogenic") |
      stringr::str_detect(
        .data$CLINVAR_CLASSIFICATION, "Pathogenic"))


  ## match clinical evidence items against
  ## query variants (non-regional - exact), civic + cgi
  var_eitems[["exact"]] <-
    pcgrr::match_eitems_to_var(
      sample_calls,
      db = "civic",
      colset = colset,
      eitems = eitems,
      region_marker = F
    )


  ## For evidence items reported at "regional level" - codon/exon/gene
  ## consider only pathogenic/likely pathogenic variants in the sample
  ## as candidates for match
  var_eitems_regional <- data.frame()
  if (nrow(sample_calls_p_lp) > 0) {
    var_eitems_regional <-
      pcgrr::match_eitems_to_var(
        sample_calls_p_lp,
        db = "civic",
        colset = colset,
        eitems = eitems,
        region_marker = T
      )
  }

  ## for regional biomarkers - perform additional quality checks
  ## (making sure variants are of correct consequence,
  ## at the correct codon/exon etc), and is loss-of-function if
  ## this is specified
  for (m in c("codon", "exon", "gene")) {
    if (nrow(var_eitems_regional) > 0) {
      var_eitems[[m]] <-
        pcgrr::qc_var_eitems(
          var_eitems = var_eitems_regional,
          marker_type = m
        )
    }
  }

  var_eitems <- pcgrr::deduplicate_eitems(
    var_eitems = var_eitems,
    target_type = "exact",
    target_other =
      c("codon", "exon", "gene")
  )

  var_eitems <- pcgrr::deduplicate_eitems(
    var_eitems = var_eitems,
    target_type = "codon",
    target_other =
      c("exon", "gene")
  )

  ## log the types and number of clinical
  ## evidence items found (exact / codon / exon)
  pcgrr::log_var_eitem_stats(
    var_eitems = var_eitems,
    target_type = "exact"
  )
  pcgrr::log_var_eitem_stats(
    var_eitems = var_eitems,
    target_type = "codon"
  )
  pcgrr::log_var_eitem_stats(
    var_eitems = var_eitems,
    target_type = "exon"
  )
  pcgrr::log_var_eitem_stats(
    var_eitems = var_eitems,
    target_type = "gene"
  )

  all_var_evidence_items <- all_var_evidence_items |>
    dplyr::bind_rows(var_eitems[["exact"]]) |>
    dplyr::bind_rows(var_eitems[["codon"]]) |>
    dplyr::bind_rows(var_eitems[["exon"]]) |>
    dplyr::bind_rows(var_eitems[["gene"]])


  ## Organize all variants in a list object 'clin_items', organized through
  ## evidence type (diagnostic|prognostic|predictive|predisposing)

  if (nrow(all_var_evidence_items) > 0) {
    for (type in c("prognostic", "diagnostic", "predictive", "predisposing")) {
      clin_eitems_list[[type]] <- all_var_evidence_items |>
        dplyr::filter(.data$EVIDENCE_TYPE == stringr::str_to_title(type)) |>
        dplyr::arrange(.data$EVIDENCE_LEVEL, .data$RATING) |>
        dplyr::select(
          .data$SYMBOL,
          .data$GENE_NAME,
          .data$CANCER_TYPE,
          .data$CLINICAL_SIGNIFICANCE,
          .data$EVIDENCE_LEVEL,
          .data$RATING,
          .data$EVIDENCE_DIRECTION,
          .data$CITATION,
          .data$THERAPEUTIC_CONTEXT,
          .data$EVIDENCE_TYPE,
          .data$DESCRIPTION,
          .data$BIOMARKER_MAPPING,
          .data$CDS_CHANGE,
          .data$LOSS_OF_FUNCTION,
          .data$GENOMIC_CHANGE
        ) |>
        dplyr::distinct()
    }
  }

  return(clin_eitems_list)
}

#' Function that makes a HTML display of virtual gene panel
#'
#' @param gene_df data frame with genes targeted in virtual panel
#'
#'
#' @export
virtual_panel_display_html <- function(gene_df) {

  i <- 1
  gene_df <- gene_df |>
    dplyr::arrange(
      dplyr::desc(.data$CONFIDENCE_LEVEL), .data$SYMBOL) |>
    dplyr::filter(.data$PRIMARY_TARGET == TRUE)

  if(length(unique(gene_df$CONFIDENCE_LEVEL)) == 1){
    if(unique(gene_df$CONFIDENCE_LEVEL) == 5){
      gene_df <- gene_df |>
        dplyr::select(
          c("SYMBOL", "ENTREZGENE", "CONFIDENCE_LEVEL")) |>
        dplyr::distinct() |>
        dplyr::mutate(PANEL_ID = as.character(NA))
    }
  }

  html_string <- "<div id=\"container\">"
  while(i <= nrow(gene_df)) {
    CONFIDENCE_LEVEL <- gene_df[i,"CONFIDENCE_LEVEL"]
    css_class <- "exploratory"
    if (CONFIDENCE_LEVEL == 3) {
      css_class <- "green"
    }
    if (CONFIDENCE_LEVEL == 2) {
      css_class <- "amber"
    }
    if (CONFIDENCE_LEVEL == 1) {
      css_class <- "red"
    }
    if (CONFIDENCE_LEVEL == -1) {
      css_class <- "custom"
    }
    if (CONFIDENCE_LEVEL == 0) {
      css_class <- "nolist"
    }
    if (CONFIDENCE_LEVEL == 5) {
      css_class <- "app_combo"
    }
    symbol <- gene_df[i, "SYMBOL"]
    entrezgene <- gene_df[i, "ENTREZGENE"]
    panel_id <- gene_df[i, "PANEL_ID"]

    gene_url <- paste0("https://www.ncbi.nlm.nih.gov/gene/", entrezgene)
    if (!is.na(panel_id)) {
      gene_url <- paste0("https://panelapp.genomicsengland.co.uk/panels/",
                         panel_id, "/", symbol)
    }

    entry_string <- paste0("  <div class=\"", css_class, "\"><a href=\"",
                           gene_url, "\" target=\"_blank\" title=\"",
                           symbol, "\">", symbol, "</a></div>")
    html_string <- paste0(html_string, entry_string)
    if (i %% 7 == 0) {
      html_string <- paste0(html_string, "</div>  <div id=\"container\">")
    }
    i <- i + 1
  }
  html_string <- paste0(html_string, "</div>")
  return(html_string)
}

