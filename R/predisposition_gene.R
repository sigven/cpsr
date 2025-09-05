#' Function that appends mechanism-of-inheritance (MOI) and mechanism
#' of disease (MOD) annotations wrt cancer predisposition genes (cpg's),
#' as well as estimated fractions of truncation vs. non-truncating variants
#' etc per predisposition gene (from ClinVar), and pathogenic AF ranges
#' based on gnomAD data, to the cpg_calls data frame.
#'
#' @param cpg_calls data frame with variant calls in predisposition genes
#' @param ref_data PCGR/CPSR reference data object
#'
#' @return cpg_calls data frame with cancer predisposition gene properties
#' appended (mechanism of disease, inheritance patterns etc)
#'
#' @export
append_cpg_properties <- function(cpg_calls, ref_data = NULL){

  pcgrr::log4r_info(
    "Appending cancer predisposition gene properties")


  ## if the cpg_calls data frame contains the column
  ## 'MAX_AF_PATH_<pop>' for any of the populations in gnomAD,
  ## remove it, since it will be annotated here
  for(pop in cpsr::acmg[['gnomAD_pops']]){
    col <- paste0("MAX_AF_PATH_", pop)
    if(col %in% colnames(cpg_calls)){
      cpg_calls[, col] <- NULL
    }
  }

  cpg_calls[,"POPMAX_AF_GNOMAD_PATH"] <- NULL

  gene_pathogenic_af_range_properties <-
    ref_data[['gene']][['cpg_pathogenic_range']] |>
    dplyr::select(
      c("ENTREZGENE",
        dplyr::starts_with("POPMAX"))
        #dplyr::starts_with("MAX_AF_PATH_"))
    ) |>
    dplyr::distinct()

  gene_truncation_properties <-
    ref_data[['variant']][['clinvar_gene_stats']] |>
    dplyr::filter(.data$CONFIDENCE == "min2goldstars") |>
    dplyr::select(
      c("ENTREZGENE",
      "BENIGN_MISSENSE_FRAC",
      "PATH_TRUNC_FRAC",
      "ACMG_PP2_CANDIDATE")
    ) |>
    dplyr::distinct()

  gene_moi_mod_properties <-
    ref_data[['gene']][['cpg']] |>
    dplyr::select(
      c("ENTREZGENE",
      "CPG_MOD",
      "CPG_MOI")
    ) |>
    dplyr::mutate(
      CPG_MOI = dplyr::if_else(
        stringr::str_detect(
      .data$CPG_MOI,
      "AD|AD/AR"),
    "AD",
    as.character(.data$CPG_MOI))) |>
    dplyr::mutate(
      CPG_MOI =
        dplyr::if_else(
          !is.na(.data$CPG_MOI) &
          !stringr::str_detect(
          .data$CPG_MOI,
          "AD|Mosaic"),
        "AR",
        .data$CPG_MOI,
        ))

  cpg_calls <- cpg_calls |>
    dplyr::left_join(
      gene_truncation_properties, by = "ENTREZGENE"
    ) |>
    dplyr::left_join(
      gene_moi_mod_properties, by = "ENTREZGENE"
    ) |>
    dplyr::left_join(
      gene_pathogenic_af_range_properties, by = "ENTREZGENE"
    ) |>
    dplyr::distinct()

  return(cpg_calls)

}
