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
    cpsr::load_germline_snv_indel(
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

  if (cps_report$settings$conf$variant_classification$clinvar_report_noncancer == 0 &
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

  cpg_call_stats <-
    pcgrr::variant_stats_report(
      cpg_calls,
      name = "v_stat_cpg"
    )

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
              "BIOMARKER_RESOLUTION",
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
      sort(cps_report[["content"]][["snv_indel"]][["variant_set"]][["tsv"]]$SYMBOL
    )),
    collapse = ", "
  )
  pcgrr::log4r_info(paste0(
    "Variants were found in the following cancer ",
    "predisposition genes: ", gene_hits
  ))

  ## secondary findings
  if (cps_report$settings$conf$variant_classification$secondary_findings == TRUE) {

    if(as.logical(
      cps_report$settings$conf$sample_properties$genotypes_available) == F){
      pcgrr::log4r_warn(paste0(
        "Assessement of secondary variant findings (ACMG SF v3.2) ",
        "NOT possible - variant genotype information unavailable"
      ))
    }else{
      secondary_calls <-
        cpsr::retrieve_secondary_calls(
          variant_calls,
          umls_map = cps_report$ref_data$phenotype$umls)

      ## do not report a secondary variant finding if this is
      ## already reported among the main hits (e.g. BRCA)
      if(NROW(secondary_calls) > 0){
        secondary_calls <- secondary_calls |>
          dplyr::anti_join(
            cpg_calls, by = "VAR_ID"
          )
        cps_report[['content']][['snv_indel']][['variant_set']][['secondary']] <-
          secondary_calls
      }
      secondary_call_stats <-
        pcgrr::variant_stats_report(
          secondary_calls,
          name = "v_stat_secondary"
        )
      cps_report[['content']][['snv_indel']][['v_stat_secondary']] <-
        secondary_call_stats$v_stat_secondary

      pcgrr::log4r_info(paste0(
        "Assessement of secondary variant findings (ACMG SF v3.2)"
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
        "Number of pathogenic ClinVar variants in the ACMG secondary findings list - other ",
        "genes of clinical significance: ",
        cps_report[["content"]][["snv_indel"]][["v_stat_secondary"]][["n_coding"]]
      ))
    }
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

#' Function that writes contents of CPSR report object to various output formats
#' (quarto HTML reports, TSV, XLSX workbooks etc)
#'
#' @param report List object with all report data (CPSR), settings etc.
#' @param output_format file format of output
#' (html//tsv/xlsx etc)

#' @export
write_cpsr_output <- function(report,
                              output_format = "html") {

  settings <- report[['settings']]
  output_dir <- settings[['output_dir']]
  sample_name <- settings[['sample_id']]
  genome_assembly <- settings[['genome_assembly']]

  sample_fname_pattern <-
    paste(sample_name, 'cpsr', genome_assembly, sep = ".")

  fnames <- list()
  fnames[["tsv"]] <-
    file.path(output_dir,
              paste0(sample_fname_pattern,
                     ".snvs_indels.tiers.tsv.gz"))
  fnames[["xlsx"]] <-
    file.path(output_dir,
              paste0(sample_fname_pattern,
                     ".xlsx"))
  fnames[["html"]] <-
    file.path(output_dir,
              paste0(sample_fname_pattern, ".html"))

  ## Path to CPSR reporting templates
  cpsr_rep_template_path <-
    system.file("templates", package = "cpsr")
  quarto_input <- file.path(
    cpsr_rep_template_path, "cpsr_report.qmd")
  report_theme <-
    settings[["conf"]][["visual_reporting"]][["visual_theme"]]

  if (output_format == "html") {
    if(file.exists(quarto_input)){

      ## make temporary directory for quarto report rendering
      tmp_quarto_dir <- file.path(
        output_dir,
        paste0('quarto_', stringi::stri_rand_strings(1, 15))
      )
      quarto_main_template <-
        glue::glue("{tmp_quarto_dir}{.Platform$file.sep}cpsr_report.qmd")
      quarto_main_template_sample <-
        glue::glue("{tmp_quarto_dir}{.Platform$file.sep}cpsr_report_sample.qmd")
      quarto_html <-
        glue::glue("{tmp_quarto_dir}{.Platform$file.sep}cpsr_report_sample.html")

      ## Copy all CPSR quarto reporting templates, bibliography, css etc to
      ## the temporary directory for quarto report rendering
      invisible(cpsr::mkdir(tmp_quarto_dir))
      system(glue::glue("cp -r {cpsr_rep_template_path}{.Platform$file.sep}* {tmp_quarto_dir}"))

      ## Save sample CPSR report object in temporary quarto rendering directory
      rds_report_path <- file.path(
        tmp_quarto_dir, "cps_report.rds")
      saveRDS(report, file = rds_report_path)

      ## Substitute rds object in main quarto template with path to sample rds
      readLines(quarto_main_template) |>
        stringr::str_replace(
          pattern = "<CPSR_REPORT_OBJECT.rds>",
          replacement = rds_report_path) |>
        writeLines(con = quarto_main_template_sample)

      ## Render report (quietly)
      pcgrr::log4r_info("------")
      pcgrr::log4r_info(
        paste0(
        "Generating quarto-based interactive HTML report (.html) with variant findings",
        "- ('",output_format, "')"))

      quarto::quarto_render(
        input = quarto_main_template_sample,
        execute_dir = tmp_quarto_dir,
        quiet = T)

      ## Copy output HTML report from temporary rendering directory
      ## to designated HTML file in output directory
      if(file.exists(quarto_html)){
        system(
          glue::glue(paste0(
            "cp -f {quarto_html} ",
            "{fnames[['html']]}")))
      }else{
        cat("WARNING\n")
      }

      ## remove temporary quarto directory (if debugging is switched off)
      if(!(settings$conf$debug)){
        system(glue::glue("rm -rf {tmp_quarto_dir}"))
      }
      pcgrr::log4r_info("------")
    }
  }

  if (output_format == "tsv") {
    if (NROW(
      report[["content"]][["snv_indel"]][["variant_set"]][[output_format]]) > 0) {
      pcgrr::log4r_info("------")
      pcgrr::log4r_info(
        paste0("Generating SNV/InDel tab-separated values file (.tsv) ",
               "with variant findings - ('",
               output_format, "')"))

      readr::write_tsv(
        report[["content"]][["snv_indel"]][["variant_set"]][[output_format]],
        file = fnames[[output_format]],
        col_names = T,
        quote = "none",
        na = ".")
    }
  }

  if (output_format == "xlsx") {
    pcgrr::log4r_info("------")
    pcgrr::log4r_info(
      paste0("Generating Excel workbook (.xlsx) with ",
             "variant findings - ('",
             output_format, "')"))
    workbook <- openxlsx2::wb_workbook() |>
      openxlsx2::wb_add_worksheet(sheet = "VIRTUAL_PANEL") |>
      openxlsx2::wb_add_worksheet(sheet = "CLASSIFICATION") |>
      openxlsx2::wb_add_worksheet(sheet = "BIOMARKER_EVIDENCE") |>
      openxlsx2::wb_add_worksheet(sheet = "SECONDARY_FINDINGS") |>
      openxlsx2::wb_add_worksheet(sheet = "GWAS_HITS") |>
      openxlsx2::wb_add_data_table(
        sheet = "CLASSIFICATION",
        x = dplyr::select(
          report[["content"]][["snv_indel"]][["variant_set"]][['tsv']],
          cpsr::col_format_output[['xlsx_classification']]),
        start_row = 1,
        start_col = 1,
        col_names = TRUE,
        na.strings = "",
        table_style = "TableStyleMedium15") |>
      openxlsx2::wb_add_data_table(
        sheet = "VIRTUAL_PANEL",
        x = report$settings$conf$gene_panel$panel_genes,
        start_row = 1,
        start_col = 1,
        col_names = TRUE,
        na.strings = "",
        table_style = "TableStyleMedium16") |>
      openxlsx2::wb_set_col_widths(
        sheet = "CLASSIFICATION",
        cols = 1:length(cpsr::col_format_output[['xlsx_classification']]),
        widths = "auto") |>
      openxlsx2::wb_set_col_widths(
        sheet = "VIRTUAL_PANEL",
        cols = 1:ncol(report$settings$conf$gene_panel$panel_genes),
        widths = "auto")

    if(NROW(report$content$snv_indel$variant_set$secondary) > 0){
      workbook <- workbook |>
        openxlsx2::wb_add_data_table(
        sheet = "SECONDARY_FINDINGS",
        x = dplyr::select(
          report$content$snv_indel$variant_set$secondary,
          cpsr::col_format_output[['xlsx_secondary']]),
        start_row = 1,
        start_col = 1,
        col_names = TRUE,
        na.strings = "",
        table_style = "TableStyleMedium17") |>
      openxlsx2::wb_set_col_widths(
        sheet = "SECONDARY_FINDINGS",
        cols = 1:length(cpsr::col_format_output[['xlsx_secondary']]),
        widths = "auto")
    }

    if(NROW(report$content$snv_indel$clin_eitem$all$any) > 0){
      workbook <- workbook |>
        openxlsx2::wb_add_data_table(
          sheet = "BIOMARKER_EVIDENCE",
          x = dplyr::select(
            report$content$snv_indel$clin_eitem$all$any,
            cpsr::col_format_output[['xlsx_biomarker']]),
          start_row = 1,
          start_col = 1,
          col_names = TRUE,
          na.strings = "",
          table_style = "TableStyleMedium18") |>
        openxlsx2::wb_set_col_widths(
          sheet = "BIOMARKER_EVIDENCE",
          cols = 1:length(cpsr::col_format_output[['xlsx_biomarker']]),
          widths = "auto")
    }

    workbook <- workbook |>
      openxlsx2::wb_save(
        fnames[['xlsx']],
        overwrite = TRUE)
  }


}


