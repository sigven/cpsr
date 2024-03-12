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

  settings <- cps_report[['settings']]
  conf <- cps_report[['settings']][['conf']]
  ref_data <- cps_report[['ref_data']]

  callset_cpsr <-
    cpsr::load_germline_snv_indel(
      fname = cps_report$settings$molecular_data$fname_mut_tsv,
      ref_data = cps_report$ref_data,
      settings = cps_report$settings
    )

  cps_report$content$snv_indel[['callset']]$variant$all <-
    callset_cpsr[['variant']][['all']]
  cps_report$content$snv_indel[['callset']]$variant$sf <-
    callset_cpsr[['variant']][['sf']]
  cps_report$content$snv_indel[['callset']]$variant$cpg_non_sf <-
    callset_cpsr$variant$cpg_non_sf
  cps_report$content$snv_indel[['callset']]$variant$gwas <-
    callset_cpsr$variant$gwas
  cps_report$content$snv_indel[['callset']]$retained_info_tags <-
    callset_cpsr[['retained_info_tags']]
  cps_report$content$snv_indel[['callset']]$biomarker_evidence <-
    callset_cpsr[['biomarker_evidence']]
  cps_report$content$snv_indel[['callset']]$variant$bm <-
    callset_cpsr[['biomarker_evidence']]$items

  col_format_output <- cpsr::col_format_output
  for(f in c('html_tier','tsv')){
    col_format_output[[f]] <- c(
      col_format_output[[f]],
      conf[["variant_classification"]][["vcftag_gnomad_AF"]])
  }

  if(cps_report$content$snv_indel[['callset']]$retained_info_tags != ""){
    tags <- stringr::str_split(
      cps_report$content$snv_indel[['callset']]$retained_info_tags, pattern=",")[[1]]
    for(t in tags){
      if(t %in% colnames(cps_report$content$snv_indel[['callset']]$variant$cpg_non_sf)){
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

  ## Exclude novel variants with MAF above user-defined threshold
  cps_report$content$snv_indel[['callset']]$variant$cpg_non_sf <-
    cps_report$content$snv_indel[['callset']]$variant$cpg_non_sf |>
    cpsr::exclude_vars_by_maf(
      conf = conf
    )

  ## Exclude non-coding variants (if option set by user)
  if (conf$other$show_noncoding == F) {
    cps_report$content$snv_indel[['callset']]$variant$cpg_non_sf <-
      cps_report$content$snv_indel[['callset']]$variant$cpg_non_sf |>
      cpsr::exclude_vars_noncoding(
        conf = conf
      )

    if (NROW(cps_report$content$snv_indel[['callset']]$variant$cpg_non_sf) == 0 &
        NROW(cps_report$content$snv_indel[['callset']]$variant$sf)) {
      pcgrr::log4r_warn(paste0(
        "There are zero remaining variants (target genes/secondary findings)",
        "- no report will be produced"
      ))
      return(NULL)
    }
  }


  ## Exclude variants related to non-cancer conditions (if option set by user)
  if (conf$variant_classification$clinvar_report_noncancer == 0){
    cps_report$content$snv_indel[['callset']]$variant$cpg_non_sf <-
      cps_report$content$snv_indel[['callset']]$variant$cpg_non_sf |>
      cpsr::exclude_vars_non_cancer(
        conf = conf
      )

    if (NROW(cps_report$content$snv_indel[['callset']]$variant$cpg_non_sf) == 0 &
        NROW(cps_report$content$snv_indel[['callset']]$variant$sf)) {
      pcgrr::log4r_warn(paste0(
        "There are zero remaining variants (target genes/secondary findings)",
        "- no report will be produced"
      ))
      return(NULL)
    }
  }


  ## get overall call statistics
  callset_all <- list()
  callset_all[['variant']] <-
    cps_report$content$snv_indel[['callset']]$variant$all
  if(NROW(callset_all[['variant']]) > 0){
    cps_report$content$snv_indel$v_stat <-
      pcgrr::variant_stats_report(
        callset = callset_all,
        name = "v_stat")$v_stat
  }

  ## get overall call statistics (cpg targets only)
  callset_cpg <- list()
  callset_cpg[['variant']] <-
    cps_report$content$snv_indel[['callset']]$variant$cpg_non_sf
  if(NROW(callset_cpg[['variant']]) > 0){
    cps_report$content$snv_indel$v_stat_cpg <-
      pcgrr::variant_stats_report(
        callset = callset_cpg,
        name = "v_stat_cpg")$v_stat_cpg
  }

  ## get overall call statistics (sf targets only)
  callset_sf <- list()
  callset_sf[['variant']] <-
    cps_report$content$snv_indel[['callset']]$variant$sf
  if(NROW(callset_sf[['variant']]) > 0){
    cps_report$content$snv_indel$v_stat_sf <-
      pcgrr::variant_stats_report(
        callset = callset_sf,
        name = "v_stat_sf")$v_stat_sf
  }

  callset_bm <- list()
  callset_bm[['variant']] <-
    cps_report$content$snv_indel[['callset']]$variant$bm
  if(NROW(callset_bm[['variant']]) > 0){
    cps_report$content$snv_indel$v_stat_bm <-
      pcgrr::variant_stats_report(
        callset = callset_bm,
        name = "v_stat_bm")$v_stat_bm
  }

  pcgrr::log4r_info(
    paste0(
      "Total number of variants in target cancer predisposition genes: ",
      "N = ", cps_report$content$snv_indel[["v_stat_cpg"]][["n"]]
    )
  )
  pcgrr::log4r_info(
    paste0(
      "Number of coding variants in target cancer predisposition genes: ",
      "N = ", cps_report$content$snv_indel[["v_stat_cpg"]][["n_coding"]]
    )
  )
  pcgrr::log4r_info(
    paste0(
      "Number of non-coding variants in cancer predisposition genes: ",
      "N = ", cps_report$content$snv_indel[["v_stat_cpg"]][["n_noncoding"]]
    )
  )

  gene_hits <- paste(
    unique(
      sort(cps_report[["content"]][["snv_indel"]][["callset"]]$variant$cpg_non_sf$SYMBOL
      )),
    collapse = ", "
  )
  pcgrr::log4r_info(paste0(
    "Variants were found in the following cancer ",
    "predisposition genes: ", gene_hits
  ))

  if(cps_report$content$snv_indel$v_stat_sf$n > 0){
    sf_hits <- paste(
      unique(
        sort(cps_report[["content"]][["snv_indel"]][["callset"]]$variant$sf$SYMBOL
        )),
      collapse = ", "
    )
    pcgrr::log4r_info(paste0(
      "Secondary variant findings were found in the following genes: ", sf_hits
    ))
  }

  cps_report$content$snv_indel[["callset"]][["tsv"]] <-
    dplyr::select(
      cps_report$content$snv_indel[["callset"]][["variant"]]$cpg_non_sf,
      dplyr::any_of(col_format_output[["tsv"]])
    )

  pcgrr::log4r_info(
    "Generating hyperlinked annotations for output data frames"
  )
  for(c in c('sf','cpg_non_sf','gwas','bm')){
    if(NROW(
      cps_report$content$snv_indel[["callset"]]$variant[[c]]) > 0){
      cps_report$content$snv_indel[["callset"]]$variant_display[[c]] <-
        cps_report$content$snv_indel[["callset"]]$variant[[c]] |>
        pcgrr::append_cancer_gene_evidence(
          ref_data = ref_data) |>
        pcgrr::append_dbnsfp_var_link() |>
        pcgrr::append_annotation_links() |>
        pcgrr::append_dbmts_var_link()

      if(c != "all"){
        col_format <- col_format_output$html_tier
        if(c == "sf"){
          col_format <- col_format_output$html_sf
        }
        if(c == "gwas"){
          col_format <- col_format_output$html_gwas
        }
        if(c == "bm"){
          col_format <- col_format_output$html_bm
        }
        cps_report$content$snv_indel[["callset"]]$variant_display[[c]] <-
          dplyr::select(
            cps_report$content$snv_indel[["callset"]]$variant_display[[c]],
            dplyr::any_of(col_format)
          )
      }

    }
  }

  cps_report[["content"]][["snv_indel"]][["eval"]] <- TRUE

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
                     ".classification.tsv.gz"))
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
    if(report$content$snv_indel$v_stat_cpg$n < 2000){
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
        report$ref_data <- NULL
        saveRDS(report, file = rds_report_path)

        ## Substitute rds object in main quarto template with path to sample rds
        readLines(quarto_main_template) |>
          stringr::str_replace(
            pattern = "<CPSR_REPORT_OBJECT.rds>",
            replacement = rds_report_path) |>
          stringr::str_replace(
            pattern = "<SAMPLE_NAME>",
            replacement = sample_name
          ) |>
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
          quiet = !report$settings$conf$debug)

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
    }else{
      pcgrr::log4r_warn("------")
      pcgrr::log4r_warn(
        paste0("Too large variant set (n = ",
               report$content$snv_indel$v_stat_cpg$n,
               ") for display in HTML report - ",
               "skipping report generation"))
      pcgrr::log4r_warn("------")
    }
  }

  if (output_format == "tsv") {
    if (NROW(
      report[["content"]][["snv_indel"]]$callset$tsv) > 0) {
      pcgrr::log4r_info("------")
      pcgrr::log4r_info(
        paste0("Generating SNV/InDel tab-separated values file (.tsv) ",
               "with variant findings - ('",
               output_format, "')"))

      readr::write_tsv(
        report[["content"]][["snv_indel"]]$callset$tsv,
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
      openxlsx2::wb_add_data_table(
        sheet = "CLASSIFICATION",
        x = dplyr::select(
          report[["content"]]$snv_indel$callset$variant$cpg_non_sf,
          dplyr::any_of(
            cpsr::col_format_output[['xlsx_classification']])),
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

    if(NROW(report[["content"]]$snv_indel$callset$variant$sf) > 0){
      workbook <- workbook |>
        openxlsx2::wb_add_data_table(
        sheet = "SECONDARY_FINDINGS",
        x = dplyr::select(
          report[["content"]]$snv_indel$callset$variant$sf,
          dplyr::any_of(
            cpsr::col_format_output[['xlsx_secondary']])),
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

    if(NROW(report$content$snv_indel$callset$variant$bm) > 0){
      workbook <- workbook |>
        openxlsx2::wb_add_data_table(
          sheet = "BIOMARKER_EVIDENCE",
          x = dplyr::select(
            report$content$snv_indel$callset$variant$bm,
            dplyr::any_of(
              cpsr::col_format_output[['xlsx_biomarker']])),
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


