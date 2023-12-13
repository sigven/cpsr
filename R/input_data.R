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

  pcgrr::log4r_info(paste0(
    "Reading annotated molecular dataset (DNA) - germline SNV/InDels"))

  callset <- pcgrr::load_dna_variants(
    fname = fname,
    cols = pcgrr::data_coltype_defs$snv_indel_germline_raw,
    ref_data = ref_data,
    retained_info_tags =
      settings[['conf']][['other']]$retained_vcf_info_tags,
    variant_origin = "Germline")

  primary_targets <-
    settings[['conf']][['gene_panel']][['panel_genes']] |>
    dplyr::filter(.data$PRIMARY_TARGET == T) |>
    dplyr::mutate(
      ENTREZGENE = as.character(.data$ENTREZGENE))

  callset[['variant']] <- callset[['variant']] |>
    cpsr::append_cpg_properties(ref_data = ref_data) |>
    pcgrr::append_cancer_gene_evidence(ref_data = ref_data) |>
    cpsr::get_insilico_prediction_statistics() |>
    pcgrr::append_dbnsfp_var_link() |>
    pcgrr::append_annotation_links() |>
    pcgrr::append_tfbs_annotation() |>
    pcgrr::append_gwas_citation_phenotype(
      ref_data = ref_data) |>
    pcgrr::append_dbmts_var_link() |>
    cpsr::assign_pathogenicity_evidence(
      ref_data = ref_data,
      settings = settings) |>
    cpsr::assign_classification() |>
    cpsr::check_variant2cancer_phenotype(
      ref_data = ref_data
    )

  for (type in pcgrr::evidence_types) {
    for (elevel in pcgrr::evidence_levels) {
      if(NROW(callset[['biomarker_evidence']][[type]][[elevel]]) > 0){
        callset[['biomarker_evidence']][[type]][[elevel]] <-
          callset[['biomarker_evidence']][[type]][[elevel]] |>
          dplyr::left_join(
            dplyr::select(
              callset[['variant']],
              c("VAR_ID",
                "GENOTYPE",
                "GENOMIC_CHANGE",
                "GENOME_VERSION",
                "SYMBOL",
                "ENTREZGENE",
                "CONSEQUENCE",
                "PROTEIN_CHANGE",
                "MUTATION_HOTSPOT",
                "CDS_CHANGE",
                "LOSS_OF_FUNCTION",
                "HGVSc",
                "HGVSp",
                "REFSEQ",
                "OFFICIAL_GENENAME",
                "PREDICTED_EFFECT",
                "PROTEIN_DOMAIN",
                "DBSNP",
                "CLINVAR",
                "COSMIC",
                "CLINVAR_CLASSIFICATION",
                "CPSR_CLASSIFICATION",
                "VEP_ALL_CSQ")
            ),
            by = c("VAR_ID")
          ) |>
          dplyr::arrange(
            .data$EVIDENCE_LEVEL,
            .data$PROTEIN_CHANGE,
            dplyr::desc(
              .data$RATING)) |>

          ## only report biomarkers with a pathogenic impact in CPSR
          dplyr::filter(
            (!is.na(CLINVAR_CLASSIFICATION) &
               stringr::str_detect(tolower(CLINVAR_CLASSIFICATION), "pathogenic")) |
              (is.na(CLINVAR_CLASSIFICATION) &
                 !is.na(CPSR_CLASSIFICATION) &
                 stringr::str_detect(tolower(CPSR_CLASSIFICATION), "pathogenic"))
          )

        if(NROW(callset[['biomarker_evidence']][[type]][[elevel]]) > 0){
          callset[['biomarker_evidence']][[type]][[elevel]] <-
            callset[['biomarker_evidence']][[type]][[elevel]] |>
            dplyr::semi_join(primary_targets, by = "ENTREZGENE")
        }
      }
    }
  }

  return(callset)

}
