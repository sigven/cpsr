#' Function that reads and validates an annotated germline SNV/InDel
#' file from CPSR pre-reporting pipeline
#'
#' @param fname Path to file name
#' @param ref_data Object with PCGR/CPSR reference data
#' @param settings Object with CPSR report configuration
#'
#' @export
load_germline_snv_indel <- function(
    fname = NA,
    ref_data = NULL,
    settings = NULL){

  invisible(
    assertthat::assert_that(!is.na(fname)))
  invisible(
    assertthat::assert_that(!is.null(ref_data)))

  conf <- settings[['conf']]
  pcgrr::log4r_info(paste0(
    "Reading annotated molecular dataset (DNA) - germline SNV/InDels"))

  callset <- pcgrr::load_dna_variants(
    fname = fname,
    cols = pcgrr::data_coltype_defs$snv_indel_germline_raw,
    ref_data = ref_data,
    vartype = "snv_indel",
    primary_site = "Any",
    retained_info_tags =
      conf[['other']]$retained_vcf_info_tags,
    variant_origin = "Germline")

  primary_targets <-
    conf[['gene_panel']][['panel_genes']] |>
    dplyr::filter(.data$PRIMARY_TARGET == T) |>
    dplyr::mutate(
      ENTREZGENE = as.character(.data$ENTREZGENE)) |>
    dplyr::select(
      c("ENTREZGENE", "PRIMARY_TARGET")
    )

  if(NROW(callset[['variant']]) > 0){
    callset[['variant']] <- callset[['variant']] |>
      cpsr::append_cpg_properties(ref_data = ref_data) |>
      cpsr::get_insilico_prediction_statistics() |>
      pcgrr::append_tfbs_annotation() |>
      pcgrr::append_gwas_citation_phenotype(
       ref_data = ref_data) |>
      cpsr::assign_pathogenicity_evidence(
        ref_data = ref_data,
        settings = settings) |>
      cpsr::assign_classification() |>
      cpsr::check_variant2cancer_phenotype(
        ref_data = ref_data
      ) |>
      cpsr::combine_novel_and_preclassified(
        conf = conf
      )


    population_tags <- unique(
      c("gnomADe_AF",
        "gnomADg_AF",
        conf$variant_classification$vcftag_gnomad_AF)
    )

    ## set missing population frequences to zero
    for (tag in population_tags) {
      if (tag %in% colnames(callset[['variant']])) {
        if (nrow(callset[['variant']][is.na(
          callset[['variant']][, tag]
        ), ]) > 0) {
          callset[['variant']][is.na(
            callset[['variant']][, tag]
          ), tag] <- 0.00
        }
      }
    }
  }

  ## Only consider biomarker evidence for variants found
  ## in primary targets (virtual gene panel)
  if(!is.null(callset$biomarker_evidence)){
    if(NROW(callset$biomarker_evidence$items) > 0){
      callset$biomarker_evidence$items <-
        callset$biomarker_evidence$items |>
        dplyr::semi_join(
          primary_targets, by = "ENTREZGENE")

      if(NROW(callset$biomarker_evidence$items) > 0){
        callset$biomarker_evidence$items <-
          callset$biomarker_evidence$items |>
          dplyr::group_by(
            dplyr::across(-c("BM_PRIMARY_SITE"))
          ) |>
          dplyr::summarise(
            BM_PRIMARY_SITE = paste(
              sort(.data$BM_PRIMARY_SITE), collapse=", "),
            .groups = "drop"
          )

        key_cols <- c("VARIANT_CLASS","ENTREZGENE","VAR_ID")
        if(all(key_cols %in% colnames(callset$variant)) == TRUE &
           all(key_cols %in%
               colnames(callset$biomarker_evidence$items) == TRUE)){

          callset$biomarker_evidence$items <-
            callset$biomarker_evidence$items |>
            dplyr::left_join(
              callset$variant,
              by = key_cols
            )

          if("FINAL_CLASSIFICATION" %in%
             colnames(callset$biomarker_evidence$items) &
             "CPSR_CLASSIFICATION_SOURCE" %in%
             colnames(callset$variant)){

            callset$biomarker_evidence$items <-
              callset$biomarker_evidence$items |>
              dplyr::filter(
                .data$CPSR_CLASSIFICATION_SOURCE == "ClinVar" &
                (.data$FINAL_CLASSIFICATION == "Pathogenic" |
                  .data$FINAL_CLASSIFICATION == "Likely_Pathogenic"))
          }

        }

      }
    }
  }

  cpsr_callset <- list()
  cpsr_callset[['variant']] <- list()
  cpsr_callset[['variant']][['all']] <- callset$variant |>
    dplyr::left_join(primary_targets, by = "ENTREZGENE") |>
    dplyr::mutate(PRIMARY_TARGET = dplyr::if_else(
      is.na(.data$PRIMARY_TARGET),
      FALSE,
      as.logical(.data$PRIMARY_TARGET)
    ))
  ## Make variant set for tier reporting (virtual panel genes only)
  cpsr_callset[['variant']][['cpg_non_sf']] <- callset$variant |>
    dplyr::semi_join(primary_targets, by = "ENTREZGENE")
  cpsr_callset[['variant']][['sf']] <- data.frame()
  cpsr_callset[['variant']][['gwas']] <- data.frame()
  cpsr_callset[['variant']][['pgx']] <- data.frame()
  cpsr_callset[['retained_info_tags']] <-
    callset$retained_info_tags
  cpsr_callset[['biomarker_evidence']] <-
    callset$biomarker_evidence

  ## Fetch secondary findings (ACMG recommendations)
  if (as.logical(
    conf$variant_classification$secondary_findings) == TRUE) {
    if(as.logical(
      conf$sample_properties$gt_detected) == F){
      pcgrr::log4r_warn(paste0(
        "Assessment of secondary variant findings (ACMG SF v3.2) ",
        "NOT possible - variant genotype information unavailable"
      ))
    }else{
      if(NROW(cpsr_callset$variant$all) > 0){
        cpsr_callset[['variant']][['sf']] <-
          cpsr::retrieve_secondary_calls(
            cpsr_callset[['variant']][['all']]
          )
      }
    }
  }

  ## Fetch chemotherapeutic toxicity variants (DPYD - CPIC)
  if (as.logical(
    conf$variant_classification$pgx_findings) == TRUE) {
    if(as.logical(
      conf$sample_properties$gt_detected) == F){
      pcgrr::log4r_warn(paste0(
        "Assessment of pharmacogenetic variants (Chemotherapy toxicity) ",
        "NOT possible - variant genotype information unavailable"
      ))
    }else{
      if(NROW(cpsr_callset$variant$all) > 0){
        cpsr_callset[['variant']][['pgx']] <-
          cpsr::retrieve_pgx_calls(
            cpsr_callset[['variant']][['all']]
          )
      }
    }
  }

  ## Fetch variants that overlap with GWAS tag SNPs
  if (as.logical(
    conf$variant_classification$gwas_findings) == TRUE &
    NROW(cpsr_callset$variant$all) > 0) {
    cpsr_callset[['variant']][['gwas']] <-
      cpsr_callset[['variant']][['all']] |>
      dplyr::filter(
        .data$PRIMARY_TARGET == FALSE &
          !is.na(.data$GWAS_CITATION)
      )
    if(NROW(cpsr_callset$variant$gwas) > 0){
      cpsr_callset[['variant']][['gwas']] <-
        cpsr_callset[['variant']][['gwas']] |>
        dplyr::arrange(
          .data$LOSS_OF_FUNCTION,
          .data$CODING_STATUS)
    }
  }

  return(cpsr_callset)

}
