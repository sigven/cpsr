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
    report_mode = "CPSR"
  )

  settings <- cps_report[["settings"]]
  conf <- cps_report[["settings"]][["conf"]]
  ref_data <- cps_report[["ref_data"]]

  callset_cpsr <-
    cpsr::load_germline_snv_indel(
      fname = cps_report$settings$molecular_data$fname_mut_tsv,
      ref_data = cps_report$ref_data,
      settings = cps_report$settings
    )

  for(cat in c('all','sf','cpg_non_sf','gwas','pgx')){
    cps_report$content$snv_indel[["callset"]]$variant[[cat]] <-
      callset_cpsr[["variant"]][[cat]]
  }

  cps_report$content$snv_indel[["callset"]]$retained_info_tags <-
    callset_cpsr[["retained_info_tags"]]
  cps_report$content$snv_indel[["callset"]]$biomarker_evidence <-
    callset_cpsr[["biomarker_evidence"]]
  cps_report$content$snv_indel[["callset"]]$variant$bm <-
    callset_cpsr[["biomarker_evidence"]]$items

  col_format_output <- cpsr::col_format_output
  for (f in c("html_tier", "tsv")) {
    col_format_output[[f]] <- c(
      col_format_output[[f]],
      conf[["variant_classification"]][["vcftag_gnomad_AF"]]
    )
  }

  if (cps_report$content$snv_indel[["callset"]]$retained_info_tags != "") {
    tags <- stringr::str_split(
      cps_report$content$snv_indel[["callset"]]$retained_info_tags,
      pattern = ",")[[1]]
    for (t in tags) {
      if (t %in% colnames(cps_report$content$snv_indel[["callset"]]$variant$cpg_non_sf)) {
        col_format_output[["tsv"]] <- c(
          col_format_output[["tsv"]],
          t
        )
        col_format_output[["html_tier"]] <- c(
          col_format_output[["html_tier"]],
          t
        )
      } else {
        # issue warning
      }
    }
  }

  ## Exclude novel variants with MAF above user-defined threshold
  cps_report$content$snv_indel[["callset"]]$variant$cpg_non_sf <-
    cps_report$content$snv_indel[["callset"]]$variant$cpg_non_sf |>
    cpsr::exclude_vars_by_maf(
      conf = conf
    )

  ## Exclude non-coding variants (if option set by user)
  if (conf$other$show_noncoding == F) {
    cps_report$content$snv_indel[["callset"]]$variant$cpg_non_sf <-
      cps_report$content$snv_indel[["callset"]]$variant$cpg_non_sf |>
      cpsr::exclude_vars_noncoding(
        conf = conf
      )

    if (NROW(cps_report$content$snv_indel[["callset"]]$variant$cpg_non_sf) == 0 &
        NROW(cps_report$content$snv_indel[["callset"]]$variant$sf) == 0 &
        NROW(cps_report$content$snv_indel[["callset"]]$variant$pgx) == 0) {
      pcgrr::log4r_warn(paste0(
        "There are zero remaining variants (main/secondary/pgx targets)",
        "- no report will be produced"
      ))
      return(NULL)
    }
  }


  ## Exclude variants related to non-cancer conditions (if option set by user)
  if (conf$variant_classification$clinvar_report_noncancer == 0) {
    cps_report$content$snv_indel[["callset"]]$variant$cpg_non_sf <-
      cps_report$content$snv_indel[["callset"]]$variant$cpg_non_sf |>
      cpsr::exclude_vars_non_cancer(
        conf = conf
      )

    if (NROW(cps_report$content$snv_indel[["callset"]]$variant$cpg_non_sf) == 0 &
      NROW(cps_report$content$snv_indel[["callset"]]$variant$sf) == 0 &
      NROW(cps_report$content$snv_indel[["callset"]]$variant$pgx) == 0) {
      pcgrr::log4r_warn(paste0(
        "There are zero remaining variants (main/secondary/pgx targets)",
        "- no report will be produced"
      ))
      return(NULL)
    }
  }

  ## get call statistics for different variant categories
  for(cat in c('all','cpg_non_sf','sf','pgx','bm')){
    stat_name <- paste0('v_stat_',cat)
    if (NROW(cps_report$content$snv_indel[["callset"]]$variant[[cat]]) > 0) {
      target_callset <- list()
      target_callset[['variant']] <-
        cps_report$content$snv_indel[["callset"]]$variant[[cat]]
      cps_report$content$snv_indel[[stat_name]] <-
        pcgrr::variant_stats_report(
          callset = target_callset,
          name = stat_name
        )[[stat_name]]
    }
  }

  pcgrr::log4r_info(
    paste0(
      "Total number of variants in target cancer predisposition genes: ",
      "N = ", cps_report$content$snv_indel[["v_stat_cpg_non_sf"]][["n"]]
    )
  )
  pcgrr::log4r_info(
    paste0(
      "Number of coding variants in target cancer predisposition genes: ",
      "N = ", cps_report$content$snv_indel[["v_stat_cpg_non_sf"]][["n_coding"]]
    )
  )
  pcgrr::log4r_info(
    paste0(
      "Number of non-coding variants in cancer predisposition genes: ",
      "N = ", cps_report$content$snv_indel[["v_stat_cpg_non_sf"]][["n_noncoding"]]
    )
  )

  gene_hits <- paste(
    unique(
      sort(cps_report[["content"]][["snv_indel"]][["callset"]]$variant$cpg_non_sf$SYMBOL)
    ),
    collapse = ", "
  )
  if(nchar(gene_hits) > 0){
    pcgrr::log4r_info(paste0(
      "Variants were found in the following cancer ",
      "predisposition genes: ", gene_hits
    ))
  }

  if (cps_report$content$snv_indel$v_stat_sf$n > 0) {
    sf_hits <- paste(
      unique(
        sort(cps_report[["content"]][["snv_indel"]][["callset"]]$variant$sf$SYMBOL)
      ),
      collapse = ", "
    )
    pcgrr::log4r_info(paste0(
      "Secondary variant findings were found in the following genes: ", sf_hits
    ))
  }

  if (cps_report$content$snv_indel$v_stat_pgx$n > 0) {
    pgx_hits <- paste(
      unique(
        sort(cps_report[["content"]][["snv_indel"]][["callset"]]$variant$pgx$SYMBOL)
      ),
      collapse = ", "
    )
    pcgrr::log4r_info(paste0(
      "Potential pharmagenomics findings were found in the following genes: ", pgx_hits
    ))
  }

  cps_report$content$snv_indel[["callset"]][["tsv"]] <-
    dplyr::select(
      cps_report$content$snv_indel[["callset"]][["variant"]]$cpg_non_sf,
      dplyr::any_of(col_format_output[["tsv"]])
    )

  pcgrr::log4r_info(
    "Generating hyperlinked annotations for output data tables"
  )
  for (c in c("sf", "cpg_non_sf", "gwas", "bm", "pgx")) {
    if (NROW(
      cps_report$content$snv_indel[["callset"]]$variant[[c]]) > 0) {
      cps_report$content$snv_indel[["callset"]]$variant_display[[c]] <-
        cps_report$content$snv_indel[["callset"]]$variant[[c]] |>
        pcgrr::append_cancer_gene_evidence(
          ref_data = ref_data
        ) |>
        pcgrr::append_dbnsfp_var_link() |>
        pcgrr::append_annotation_links() |>
        pcgrr::append_dbmts_var_link()

      if (c != "all") {
        col_format <- col_format_output$html_tier
        if (c == "sf") {
          col_format <- col_format_output$html_sf
        }
        if (c == "gwas") {
          col_format <- col_format_output$html_gwas
        }
        if (c == "bm") {
          col_format <- col_format_output$html_bm
        }
        if (c == "pgx") {
          col_format <- col_format_output$html_pgx
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
  settings <- report[["settings"]]
  output_dir <- settings[["output_dir"]]
  sample_name <- settings[["sample_id"]]
  genome_assembly <- settings[["genome_assembly"]]

  sample_fname_pattern <-
    paste(sample_name, "cpsr", genome_assembly, sep = ".")

  fnames <- list()
  fnames[["tsv"]] <-
    file.path(
      output_dir,
      paste0(
        sample_fname_pattern,
        ".classification.tsv.gz"
      )
    )
  fnames[["xlsx"]] <-
    file.path(
      output_dir,
      paste0(
        sample_fname_pattern,
        ".xlsx"
      )
    )
  fnames[["html"]] <-
    file.path(
      output_dir,
      paste0(sample_fname_pattern, ".html")
    )

  ## Path to CPSR reporting templates
  templates_dir <- "templates"
  cpsr_rep_template_path <-
    system.file(templates_dir, package = "cpsr")
  quarto_input <- file.path(
    cpsr_rep_template_path, "cpsr_report.qmd"
  )

  if (output_format == "html") {
    if (report$content$snv_indel$v_stat_cpg_non_sf$n < 15000) {
      if (file.exists(quarto_input)) {
        ## make temporary directory for quarto report rendering
        tmp_quarto_dir1 <- file.path(
          output_dir,
          paste0("quarto_", stringi::stri_rand_strings(1, 15))
        )
        pcgrr::mkdir(tmp_quarto_dir1)
        # files get copied under tmp/templates/
        # see https://github.com/sigven/cpsr/issues/61
        # file.copy(cpsr_rep_template_path, tmp_quarto_dir1, recursive = TRUE, overwrite = TRUE)
        system2("cp", args = c("-r", shQuote(cpsr_rep_template_path), shQuote(tmp_quarto_dir1)))
        # so now overwrite the variable
        tmp_quarto_dir <- file.path(tmp_quarto_dir1, templates_dir)

        quarto_main_template <-
          file.path(tmp_quarto_dir, "cpsr_report.qmd")
        quarto_main_template_sample <-
          file.path(tmp_quarto_dir, "cpsr_report_sample.qmd")
        quarto_html <-
          file.path(tmp_quarto_dir, "cpsr_report_sample.html")

        ## Copy all CPSR quarto reporting templates, bibliography, css etc to
        ## the temporary directory for quarto report rendering

        ## Save sample CPSR report object in temporary quarto rendering directory
        rds_report_path <- file.path(
          tmp_quarto_dir, "cps_report.rds"
        )
        report$ref_data <- NULL
        saveRDS(report, file = rds_report_path)

        ## Substitute rds object in main quarto template with path to sample rds
        readLines(quarto_main_template) |>
          stringr::str_replace(
            pattern = "<CPSR_REPORT_OBJECT.rds>",
            replacement = rds_report_path
          ) |>
          stringr::str_replace(
            pattern = "<SAMPLE_NAME>",
            replacement = sample_name
          ) |>
          writeLines(con = quarto_main_template_sample)

        ## Render report (quietly)
        pcgrr::log4r_info("------")
        pcgrr::log4r_info(
          paste0(
            "Generating quarto-based interactive HTML report ",
            "(.html) with variant findings"
          )
        )

        quarto::quarto_render(
          input = quarto_main_template_sample,
          execute_dir = tmp_quarto_dir,
          quiet = !report$settings$conf$debug
        )

        ## Copy output HTML report from temporary rendering directory
        ## to designated HTML file in output directory
        if (file.exists(quarto_html)) {
          file.copy(quarto_html, fnames[["html"]], overwrite = TRUE)
        } else {
          cat("WARNING\n")
        }

        ## remove temporary quarto directory (if debugging is switched off)
        if (!(settings$conf$debug)) {
          unlink(c(tmp_quarto_dir, tmp_quarto_dir1), force = TRUE, recursive = TRUE)
        }
        pcgrr::log4r_info("------")
      }
    } else {
      pcgrr::log4r_warn("------")
      pcgrr::log4r_warn(
        paste0(
          "Too large variant set (n = ",
          report$content$snv_indel$v_stat_cpg$n,
          ") for display in HTML report - ",
          "skipping report generation"
        )
      )
      pcgrr::log4r_warn("------")
    }
  }

  if (output_format == "tsv") {
    if (NROW(
      report[["content"]][["snv_indel"]]$callset$tsv
    ) > 0) {
      pcgrr::log4r_info("------")
      pcgrr::log4r_info(
        paste0(
          "Generating SNV/InDel tab-separated values file (.tsv.gz) ",
          "with variant findings"
        )
      )
      readr::write_tsv(
        report[["content"]][["snv_indel"]]$callset$tsv,
        file = fnames[[output_format]],
        col_names = T,
        quote = "none",
        na = "."
      )
    }
  }

  if (output_format == "xlsx") {
    pcgrr::log4r_info("------")
    pcgrr::log4r_info(
      paste0(
        "Generating Excel workbook (.xlsx) with ",
        "variant findings"
      )
    )
    workbook <- openxlsx2::wb_workbook() |>
      openxlsx2::wb_add_worksheet(sheet = "VIRTUAL_PANEL") |>
      openxlsx2::wb_add_data_table(
        sheet = "VIRTUAL_PANEL",
        x = report$settings$conf$gene_panel$panel_genes,
        start_row = 1,
        start_col = 1,
        col_names = TRUE,
        na.strings = "",
        table_style = "TableStyleMedium16"
      ) |>
      openxlsx2::wb_set_col_widths(
        sheet = "VIRTUAL_PANEL",
        cols = 1:ncol(report$settings$conf$gene_panel$panel_genes),
        widths = "auto"
      )

    if (NROW(report[["content"]]$snv_indel$callset$variant$cpg_non_sf) > 0) {
      workbook <- workbook |>
        openxlsx2::wb_add_worksheet(sheet = "CLASSIFICATION") |>
        openxlsx2::wb_add_data_table(
          sheet = "CLASSIFICATION",
          x = dplyr::select(
            report[["content"]]$snv_indel$callset$variant$cpg_non_sf,
            dplyr::any_of(
              cpsr::col_format_output[["xlsx_classification"]]
            )
          ),
          start_row = 1,
          start_col = 1,
          col_names = TRUE,
          na.strings = "",
          table_style = "TableStyleMedium15"
        ) |>
        openxlsx2::wb_set_col_widths(
          sheet = "CLASSIFICATION",
          cols = 1:length(cpsr::col_format_output[["xlsx_classification"]]),
          widths = "auto"
        )
    }


    if (NROW(report[["content"]]$snv_indel$callset$variant$sf) > 0) {
      workbook <- workbook |>
        openxlsx2::wb_add_worksheet(sheet = "SECONDARY_FINDINGS") |>
        openxlsx2::wb_add_data_table(
          sheet = "SECONDARY_FINDINGS",
          x = dplyr::select(
            report[["content"]]$snv_indel$callset$variant$sf,
            dplyr::any_of(
              cpsr::col_format_output[["xlsx_secondary"]]
            )
          ),
          start_row = 1,
          start_col = 1,
          col_names = TRUE,
          na.strings = "",
          table_style = "TableStyleMedium17"
        ) |>
        openxlsx2::wb_set_col_widths(
          sheet = "SECONDARY_FINDINGS",
          cols = 1:length(cpsr::col_format_output[["xlsx_secondary"]]),
          widths = "auto"
        )
    }

    if (NROW(report$content$snv_indel$callset$variant$bm) > 0) {
      bm_excel <- report$content$snv_indel$callset$variant$bm |>
        dplyr::mutate(
          BM_MOLECULAR_PROFILE = pcgrr::strip_html(
            .data$BM_MOLECULAR_PROFILE
          )
        )

      workbook <- workbook |>
        openxlsx2::wb_add_worksheet(sheet = "BIOMARKER_EVIDENCE") |>
        openxlsx2::wb_add_data_table(
          sheet = "BIOMARKER_EVIDENCE",
          x = dplyr::select(
            bm_excel,
            dplyr::any_of(
              cpsr::col_format_output[["xlsx_biomarker"]]
            )
          ),
          start_row = 1,
          start_col = 1,
          col_names = TRUE,
          na.strings = "",
          table_style = "TableStyleMedium18"
        ) |>
        openxlsx2::wb_set_col_widths(
          sheet = "BIOMARKER_EVIDENCE",
          cols = 1:length(cpsr::col_format_output[["xlsx_biomarker"]]),
          widths = "auto"
        )
    }

    if (NROW(report[["content"]]$snv_indel$callset$variant$pgx) > 0) {
      workbook <- workbook |>
        openxlsx2::wb_add_worksheet(sheet = "PHARMACOGENETIC_FINDINGS") |>
        openxlsx2::wb_add_data_table(
          sheet = "PHARMACOGENETIC_FINDINGS",
          x = dplyr::select(
            report[["content"]]$snv_indel$callset$variant$pgx,
            dplyr::any_of(
              cpsr::col_format_output[["xlsx_pgx"]]
            )
          ),
          start_row = 1,
          start_col = 1,
          col_names = TRUE,
          na.strings = "",
          table_style = "TableStyleMedium19"
        ) |>
        openxlsx2::wb_set_col_widths(
          sheet = "PHARMACOGENETIC_FINDINGS",
          cols = 1:length(cpsr::col_format_output[["xlsx_pgx"]]),
          widths = "auto"
        )
    }


    openxlsx2::wb_save(
      wb = workbook,
      fnames[["xlsx"]],
      overwrite = TRUE
    )
  }
}
