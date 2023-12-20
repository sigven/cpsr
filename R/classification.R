#' Function that assigns final pathogenicity classification (B, LB, VUS, P, LP)
#' based on accumulated scores from different ACMG criteria and pre-defined
#' cutoffs (calibrated against ClinVar)
#'
#' @param cpg_calls data frame with variant calls in predisposition genes
#'
#' @return cpg_calls data frame with pathogenicity classification appended
#'
#' @export
assign_classification <- function(cpg_calls) {

  evidence_codes <- cpsr::acmg[["evidence_codes"]]

  pcgrr::log4r_info(paste0(
    "Assigning five-tier classifications (P, LP, VUS, LB, B) based on ",
    "aggregated scores from ACMG criteria"))

  path_cols <- c(
    "CPSR_CLASSIFICATION",
    "CPSR_CLASSIFICATION_DOC",
    "CPSR_CLASSIFICATION_CODE",
    "cpsr_score_pathogenic",
    "cpsr_score_benign"
  )
  cpg_calls <- cpg_calls[, !(colnames(cpg_calls) %in% path_cols)]

  cpg_calls$CPSR_CLASSIFICATION <- "VUS"
  cpg_calls$CPSR_CLASSIFICATION_DOC <- ""
  cpg_calls$CPSR_CLASSIFICATION_CODE <- ""
  cpg_calls$cpsr_score_pathogenic <- 0
  cpg_calls$cpsr_score_benign <- 0

  i <- 1
  while (i <= nrow(evidence_codes)) {
    category <- evidence_codes[i, ]$category
    pole <- evidence_codes[i, ]$pathogenicity_pole
    description <- evidence_codes[i, ]$description
    cpsr_evidence_code <- evidence_codes[i, ]$cpsr_evidence_code
    score <- evidence_codes[i, ]$path_score
    if (cpsr_evidence_code %in% colnames(cpg_calls)) {
      cpg_calls <- cpg_calls |>
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
                paste0("- ", description), ""
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

  lb_upper_limit <- -1.5
  lb_lower_limit <- -4.5
  b_upper_limit <- -5.0
  vus_lower_limit <- -1.0
  vus_upper_limit <- 2.0
  lp_lower_limit <- 2.5
  lp_upper_limit <- 4.5
  p_lower_limit <- 5.0


  cpg_calls <- cpg_calls |>
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

  return(cpg_calls)
}


#' Function that assigns variant pathogenicity evidence based on ACMG guidelines
#'
#' @param calls sample calls with dbnsfp annotations
#' @param settings cpsr settings object
#' @param ref_data pcgr data object
#'
#' @return calls
#'
#' @export
assign_pathogenicity_evidence <- function(calls, settings, ref_data) {

  invisible(assertthat::assert_that(!is.null(calls)))
  invisible(assertthat::assert_that(!is.null(settings)))
  invisible(assertthat::assert_that(!is.null(settings$conf)))
  invisible(assertthat::assert_that(!is.null(settings$conf$variant_classification)))
  invisible(assertthat::assert_that(!is.null(ref_data)))
  invisible(assertthat::assert_that(is.data.frame(calls)))
  #pcgrr::validate_settings(settings)

  pcgrr::log4r_info(
    "Assigning pathogenicity codes according to enhanced ACMG criteria/guidelines")
  classification_settings <-
    settings$conf$variant_classification

  gad_population <- toupper(classification_settings[["pop_gnomad"]])
  gad_AN_tag <- classification_settings[['vcftag_gnomad_AN']]
  gad_AF_tag <- classification_settings[['vcftag_gnomad_AF']]
  gad_NHOMALT_tag <- classification_settings[['vcftag_gnomad_NHOMALT']]
  gad_AC_tag <- classification_settings[['vcftag_gnomad_AC']]

  # pathogenic_range_ac <- 20
  pathogenic_range_af <- cpsr::acmg[["pathogenic_range_gnomad"]][["af"]]
  min_an <- cpsr::acmg[["pathogenic_range_gnomad"]][["min_an"]]

  acmg_ev_codes <-
    c(
      "ACMG_BA1_AD",
      ## Very high MAF (> 0.5% in gnomAD non-cancer pop subset) -
      ## min AN = 12,000, - Dominant mechanism of disease
      "ACMG_BS1_1_AD",
      ## High MAF (> 0.1% in gnomAD non-cancer pop subset) -
      ## min AN = 12,000 - Dominant mechanism of disease
      "ACMG_BS1_2_AD",
      ## Somewhat high MAF (> 0.005% in gnomAD non-cancer pop subset) -
      ## min AN = 12,000 - Dominant mechanism of disease
      "ACMG_BA1_AR",
      ## Very high MAF (> 1% in gnomAD non-cancer pop subset) -
      ## min AN = 12,000 - Recessive mechanism of disease
      "ACMG_BS1_1_AR",
      ## High MAF (> 0.3% in gnomAD non-cancer pop subset) -
      ## min AN = 12,000 - Recessive mechanism of disease
      "ACMG_BS1_2_AR",
      ## Somewhat high MAF (> 0.005% in gnomAD non-cancer pop subset) -
      ## min AN = 12,000 - Recessive mechanism of disease
      # "ACMG_BS2_1",
      ## 1 homozygote in gnomAD non-cancer pop subset -
      ## severe, early onset, highly penetrant
      # "ACMG_BS2_2",
      ## 2 homozygotes in gnomAD non-cancer pop subset -
      ## severe, early onset, highly penetrant
      # "ACMG_BS2_3",
      ## 2 homozygotes in gnomAD non-cancer pop subset -
      ## moderate, early onset, variably penetrant
      "ACMG_PM2_1",
      ## Allele count within pathogenic range (MAF < 0.005% in the
      ## population-specific non-cancer gnomAD subset, min AN = 12,000)
      "ACMG_PM2_2",
      ## Alternate allele absent in the population-specific
      ## non-cancer gnomAD subset
      "ACMG_PVS1_1",
      ## Null variant - predicted as LoF by LOFTEE - within pathogenic range
      ## - LoF established for gene
      "ACMG_PVS1_2",
      ## Null variant - not predicted as LoF by LOFTEE -
      ## within pathogenic range - LoF established for gene
      "ACMG_PVS1_3",
      ## Null variant - predicted as LoF by LOFTEE - within pathogenic range -
      ## LoF not established for gene
      "ACMG_PVS1_4",
      ## Null variant - not predicted as LoF by LOFTEE --
      ## within pathogenic range - LoF not established for gene
      "ACMG_PVS1_5",
      ## start lost - within pathogenic range - Lof established for gene
      "ACMG_PVS1_6",
      ## start lost - within pathogenic range - LoF not established for gene
      "ACMG_PVS1_7",
      ## donor/acceptor variant - predicted as LoF by LOFTEE -
      ## within pathogenic range
      ## - not last intron - LoF established for gene
      "ACMG_PVS1_8",
      ## donor/acceptor variant - last intron - within pathogenic range -
      ## LoF established for gene
      "ACMG_PVS1_9",
      ## donor/acceptor variant - not last intron - within pathogenic range
      ## - LoF not established for gene
      "ACMG_PVS1_10",
      ## donor variant at located at the +3, +4 or +5 position of the intron -
      ## within the pathogenic range (i.e. MAF < 0.005% in gnOMAD))
      "ACMG_PS1",
      ## Same amino acid change as a previously established pathogenic
      ## variant (ClinVar) regardless of nucleotide change
      "ACMG_PP2",
      ## Missense variant in a gene that has a relatively low rate of
      ## benign missense variation (<20%) and
      ## where missense variants are a common mechanism of disease
      ## (>50% of high-confidence pathogenic variants (ClinVar))
      "ACMG_PM4",
      ## Protein length changes due to inframe indels or nonstop variant
      ## in non-repetitive regions of genes
      ## that harbor variants with a dominant mode of inheritance.
      "ACMG_PPC1",
      ## Protein length changes due to inframe indels or nonstop variant
      ## in non-repetitive regions of genes
      ## that harbor variants with a recessive mode of inheritance.
      "ACMG_PM5",
      ## Novel missense change at an amino acid residue where a different
      ## missense change determined to be pathogenic
      ## has been seen before (ClinVar)
      "ACMG_PP3",
      ## Multiple lines of computational evidence support a
      ## deleterious effect on the gene or gene product
      ## (conservation, evolutionary, splicing impact, etc. - from dbNSFP
      "ACMG_BP4",
      ## Multiple lines of computational evidence support a benign
      ## effect on the gene or gene product
      ## (conservation, evolutionary, splicing impact, etc. - from dbNSFP
      "ACMG_BMC1",
      ## Peptide change is at the same location of a
      ## known benign change (ClinVar)
      "ACMG_BSC1",
      ## Peptide change is reported as benign (ClinVar),
      "ACMG_BP3",
      ## Variants in promoter or untranslated regions
      "ACMG_BP7",
      ## Silent/intronic variant outside of the splice site consensus
      "ACMG_BP1"
    )
  ## Missense variant in a gene for which primarily truncating
  ## variants are known to cause disease (ClinVar)


  path_columns <-
    c(
      acmg_ev_codes,
      "CODON",
      "PATHOGENIC_CODON",
      "PATHOGENIC_PEPTIDE_CHANGE",
      "BENIGN_CODON",
      "BENIGN_PEPTIDE_CHANGE",
      "hotspot_region",
      "hotspot_symbol",
      "hotspot_entrezgene",
      "hotspot_codon",
      "hotspot_aa",
      "hotspot_pvalue"
    )
  calls <- calls[, !(colnames(calls) %in% path_columns)]

  benign_peptide_changes <-
    ref_data[['variant']][['clinvar_sites']] |>
    dplyr::filter(.data$GOLD_STARS >= 2 & .data$BENIGN == 1) |>
    dplyr::select(c("ENTREZGENE", "HGVSP", "BENIGN", "VAR_ID")) |>
    dplyr::filter(!is.na(.data$ENTREZGENE) & !is.na(.data$HGVSP)) |>
    dplyr::rename(BENIGN_PEPTIDE_CHANGE = .data$BENIGN,
                  VAR_ID_BENIGN_CHANGE = .data$VAR_ID) |>
    dplyr::distinct()

  pathogenic_peptide_changes <-
    ref_data[['variant']][['clinvar_sites']] |>
    dplyr::filter(.data$GOLD_STARS >= 2 & .data$PATHOGENIC == 1) |>
    dplyr::select(c("ENTREZGENE", "HGVSP", "PATHOGENIC", "VAR_ID")) |>
    dplyr::filter(!is.na(.data$ENTREZGENE) & !is.na(.data$HGVSP)) |>
    dplyr::rename(PATHOGENIC_PEPTIDE_CHANGE = .data$PATHOGENIC,
                  VAR_ID_PATH_CHANGE = .data$VAR_ID) |>
    dplyr::distinct()

  benign_codons <- as.data.frame(
    ref_data[['variant']][['clinvar_sites']] |>
    dplyr::filter(.data$GOLD_STARS >= 2 & .data$BENIGN == 1) |>
    dplyr::select(c("ENTREZGENE", "CODON", "BENIGN", "VAR_ID")) |>
    dplyr::filter(!is.na(.data$ENTREZGENE) & !is.na(.data$CODON)) |>
    dplyr::group_by(.data$ENTREZGENE, .data$CODON) |>
    dplyr::reframe(VAR_ID = paste(.data$VAR_ID, collapse=";"),
                   BENIGN = paste(unique(.data$BENIGN), collapse=";")) |>
    dplyr::rename(BENIGN_CODON = .data$BENIGN,
                  VAR_ID_BENIGN_CODON = .data$VAR_ID) |>
    dplyr::distinct()
  )

  pathogenic_codons <-
    ref_data[['variant']][['clinvar_sites']] |>
    dplyr::filter(.data$GOLD_STARS >= 2 & .data$PATHOGENIC == 1) |>
    dplyr::select(c("ENTREZGENE", "CODON", "PATHOGENIC", "VAR_ID")) |>
    dplyr::filter(!is.na(.data$ENTREZGENE) & !is.na(.data$CODON)) |>
    dplyr::group_by(.data$ENTREZGENE, .data$CODON) |>
    dplyr::reframe(VAR_ID = paste(.data$VAR_ID, collapse=";"),
                   PATHOGENIC = paste(unique(.data$PATHOGENIC), collapse=";")) |>
    dplyr::rename(PATHOGENIC_CODON = .data$PATHOGENIC,
                  VAR_ID_PATH_CODON = .data$VAR_ID) |>
    dplyr::distinct()

  ## Assign logical ACMG evidence indicators
  #
  #
  # ACMG_PP3 - Multiple lines (>=5) of insilico evidence support a
  #             deleterious effect on the gene or gene product
  ##           (conservation, evolutionary, splicing impact, etc.)
  # ACMG_BP4 - Multiple lines (>=5) of insilico evidence support a benign effect.
  #
  # Computational evidence for deleterious/benign effect is taken from
  # invidual algorithm predictions in dbNSFP: SIFT,Provean,MutationTaster,
  # MutationAssessor,M_CAP,MutPred,FATHMM,FATHMM-mkl,DBNSFP_RNN,dbscSNV_RF,
  # dbscSNV_AdaBoost
  # Default scheme (from default TOML file):
  # 1) Damaging: Among all possible protein variant effect predictions, at
  #              least six algorithms must have made a call,
  #              with at least 8 predicted as damaging/D
  #              (possibly_damaging/PD), and at most two
  #              predicted as tolerated/T (PP3)
  #       - at most 1 prediction for a splicing neutral effect
  #    Exception: if both splice site predictions indicate damaging effects;
  #    ignore other criteria
  # 2) Tolerated: Among all possible protein variant effect predictions, at
  #    least six algorithms must have made a call,
  #    with at least 8 predicted as tolerated, and at most 2
  #    predicted as damaging (BP4)
  #    - 0 predictions of splice site affected

  dbnsfp_min_majority <- cpsr::acmg[["insilico_pred_min_majority"]]
  dbnsfp_max_minority <- cpsr::acmg[["insilico_pred_max_minority"]]
  dbnsfp_min_called <- dbnsfp_min_majority

  calls <- calls |>
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
    dplyr::mutate(
      ACMG_BP4 = dplyr::if_else(
        .data$N_INSILICO_CALLED >= dbnsfp_min_called &
          .data$N_INSILICO_TOLERATED >= dbnsfp_min_majority &
          .data$N_INSILICO_DAMAGING <= dbnsfp_max_minority &
          .data$N_INSILICO_SPLICING_AFFECTED == 0, TRUE,
        FALSE, FALSE
      )
    ) |>
    dplyr::mutate(ACMG_PP3 = dplyr::case_when(
      .data$N_INSILICO_SPLICING_AFFECTED == 2 ~ TRUE,
      TRUE ~ as.logical(.data$ACMG_PP3)
    ))

  ## Assign logical ACMG evidence indicators based on population frequency
  ## data in non-cancer samples from gnomAD (Dominant vs. recessive
  ## modes of inheritance)
  # 'ACMG_BA1_AD'   - Very high MAF (> 0.5% in gnomAD non-cancer pop subset) -
  #                   min AN = 12,000 - Dominant mechanism of disease
  # 'ACMG_BS1_1_AD' - High MAF (> 0.1% in gnomAD non-cancer pop subset) -
  #                   min AN = 12,000 - Dominant mechanism of disease
  # 'ACMG_BS1_2_AD' - Somewhat high MAF (> 0.005% in gnomAD non-cancer pop
  #                   subset) - Dominant mechanism of disease
  # 'ACMG_BA1_AR'   - Very high MAF (> 1% in gnomAD non-cancer pop subset) -
  #                   min AN = 12,000 - Recessive mechanism of disease
  # 'ACMG_BS1_1_AR' - High MAF (> 0.3% in gnomAD non-cancer pop subset) -
  #                   min AN = 12,000 - Recessive mechanism of disease
  # 'ACMG_BS1_2_AR' - Somewhat high MAF (> 0.005% in gnomAD non-cancer pop
  #                   subset) - Recessive mechanism of disease
  # 'ACMG_PM2_1'    - Allele count within pathogenic range (MAF <= 0.005%
  #                   in the population-specific non-cancer gnomAD subset,
  #                   min AN = 12,000)
  # 'ACMG_PM2_2'    - Alternate allele absent in the population-specific
  #                   non-cancer gnomAD subset
  if (gad_AN_tag %in% colnames(calls) &
      gad_AC_tag %in% colnames(calls) &
      gad_NHOMALT_tag %in% colnames(calls)) {

    calls <- calls |>
      dplyr::mutate(
        gad_af =
          dplyr::if_else(
            !!rlang::sym(gad_AN_tag) >= min_an,
            as.double(!!rlang::sym(gad_AC_tag) /
                        !!rlang::sym(gad_AN_tag)),
            as.double(NA), as.double(NA)
          )
      ) |>
      dplyr::mutate(
        ACMG_PM2_1 =
          dplyr::if_else(
            !!rlang::sym(gad_AN_tag) >= min_an &
              !is.na(!!rlang::sym(gad_AC_tag)) &
              .data$gad_af <= pathogenic_range_af,
            TRUE, FALSE, FALSE
          )
      ) |>
      dplyr::mutate(
        ACMG_PM2_2 = dplyr::if_else(
          is.na(!!rlang::sym(gad_AC_tag)),
                                    TRUE, FALSE, FALSE
        )
      ) |>
      dplyr::mutate(
        ACMG_BA1_AD = dplyr::if_else(
          .data$ACMG_PM2_2 == FALSE &
            .data$gad_af >= 0.005 &
            .data$CPG_MOI == "AD",
          TRUE, FALSE, FALSE
        )
      ) |>
      dplyr::mutate(
        ACMG_BS1_1_AD = dplyr::if_else(
          .data$ACMG_BA1_AD == FALSE &
            .data$ACMG_PM2_2 == FALSE &
            .data$gad_af >= 0.001 &
            .data$CPG_MOI == "AD",
          TRUE, FALSE, FALSE
        )
      ) |>
      dplyr::mutate(
        ACMG_BS1_2_AD = dplyr::if_else(
          .data$ACMG_BS1_1_AD == FALSE &
            .data$ACMG_BA1_AD == FALSE &
            .data$ACMG_PM2_2 == FALSE &
            .data$gad_af > pathogenic_range_af &
            .data$CPG_MOI == "AD",
          TRUE, FALSE, FALSE
        )
      ) |>
      dplyr::mutate(
        ACMG_BA1_AR = dplyr::if_else(
          .data$ACMG_PM2_2 == FALSE &
            .data$gad_af >= 0.01 &
            (.data$CPG_MOI == "AR" |
               is.na(.data$CPG_MOI)),
          TRUE, FALSE, FALSE
        )
      ) |>
      dplyr::mutate(
        ACMG_BS1_1_AR = dplyr::if_else(
          .data$ACMG_BA1_AR == FALSE &
            .data$ACMG_PM2_2 == FALSE &
            .data$gad_af >= 0.003 &
            (.data$CPG_MOI == "AR" |
               is.na(.data$CPG_MOI)),
          TRUE, FALSE, FALSE
        )
      ) |>
      dplyr::mutate(
        ACMG_BS1_2_AR = dplyr::if_else(
          .data$ACMG_BA1_AR == FALSE &
            .data$ACMG_BS1_1_AR == FALSE &
            .data$ACMG_PM2_2 == FALSE &
            .data$gad_af > pathogenic_range_af &
            (.data$CPG_MOI == "AR" |
               is.na(.data$CPG_MOI)),
          TRUE, FALSE, FALSE
        )
      )
  }

  ## Assign logical ACMG evidence indicators on NULL variants in known
  # predisposition genes (LoF established as mechanism of disease or not,
  # presumed loss of mRNA/protein (LOFTEE) or not)
  #
  # 'ACMG_PVS1_1' - Null variant (frameshift, nonsense) -
  #  predicted as LoF by LOFTEE - within pathogenic range - LoF established
  # 'ACMG_PVS1_2' - Null variant (frameshift, nonsense) -
  # not predicted as LoF by LOFTEE - within pathogenic range - LoF established
  # 'ACMG_PVS1_3' - Null variant (frameshift, nonsense) -
  # predicted as LoF by LOFTEE - within pathogenic range - LoF not established
  # 'ACMG_PVS1_4' - Null variant (frameshift, nonsense) -
  # not predicted as LoF by LOFTEE -- within pathogenic range - LoF not
  # established for gene
  # 'ACMG_PVS1_5' - start lost - within pathogenic range - Lof established
  # 'ACMG_PVS1_6' - start lost - within pathogenic range - LoF not established
  # 'ACMG_PVS1_7' - splice acceptor/donor variant - predicted as LoF by LOFTEE
  # - not last intron - within pathogenic range - Lof established
  # 'ACMG_PVS1_8' - splice acceptor/donor variant - predicted as LoF by LOFTEE
  # - last intron - within pathogenic range - Lof established
  # 'ACMG_PVS1_9' - splice acceptor/donor variant - predicted as LoF by LOFTEE
  # - not last intron - within pathogenic range - Lof established
  # 'ACMG_PVS1_10' - splice variant involving a donor at +3A/G, +4A or +5G -
  # predicted as damaging by insilico predictions - within pathogenic range

  calls <- calls |>
    dplyr::mutate(
      ACMG_PVS1_1 =
        dplyr::if_else(
          .data$NULL_VARIANT == T &
            .data$LOSS_OF_FUNCTION == T &
            .data$CPG_MOD == "LoF" &
            (.data$ACMG_PM2_1 == TRUE |
               .data$ACMG_PM2_2 == TRUE),
          TRUE, FALSE, FALSE
        )
    ) |>
    dplyr::mutate(
      ACMG_PVS1_3 =
        dplyr::if_else(
          .data$NULL_VARIANT == T &
            .data$LOSS_OF_FUNCTION == T &
            (is.na(.data$CPG_MOD) | .data$CPG_MOD != "LoF") &
            (.data$ACMG_PM2_1 == TRUE | .data$ACMG_PM2_2 == TRUE),
          TRUE, FALSE, FALSE
        )
    ) |>
    dplyr::mutate(
      ACMG_PVS1_2 =
        dplyr::if_else(
          .data$NULL_VARIANT == T &
            .data$LOSS_OF_FUNCTION == F &
            .data$CPG_MOD == "LoF" &
            (.data$ACMG_PM2_1 == TRUE | .data$ACMG_PM2_2 == TRUE),
          TRUE, FALSE, FALSE
        )
    ) |>
    dplyr::mutate(
      ACMG_PVS1_4 =
        dplyr::if_else(
          .data$NULL_VARIANT == T &
            .data$LOSS_OF_FUNCTION == F &
            (is.na(.data$CPG_MOD) | .data$CPG_MOD != "LoF") &
            (.data$ACMG_PM2_1 == TRUE | .data$ACMG_PM2_2 == TRUE),
          TRUE, FALSE, FALSE
        )
    ) |>
    dplyr::mutate(
      ACMG_PVS1_5 =
        dplyr::if_else(
          .data$CONSEQUENCE == "start_lost" &
            .data$CPG_MOD == "LoF" &
            (.data$ACMG_PM2_1 == TRUE | .data$ACMG_PM2_2 == TRUE),
          TRUE, FALSE, FALSE
        )
    ) |>
    dplyr::mutate(
      ACMG_PVS1_6 =
        dplyr::if_else(
          .data$CONSEQUENCE == "start_lost" &
            (is.na(.data$CPG_MOD) |.data$CPG_MOD != "LoF") &
            (.data$ACMG_PM2_1 == TRUE | .data$ACMG_PM2_2 == TRUE),
          TRUE, FALSE, FALSE
        )
    ) |>
    dplyr::mutate(
      ACMG_PVS1_7 =
        dplyr::if_else(
          .data$LOSS_OF_FUNCTION == T &
            stringr::str_detect(.data$CONSEQUENCE, "_donor|_acceptor") &
            .data$LAST_INTRON == F & .data$CPG_MOD == "LoF" &
            (.data$ACMG_PM2_1 == TRUE | .data$ACMG_PM2_2 == TRUE),
          TRUE, FALSE, FALSE
        )
    ) |>
    dplyr::mutate(
      ACMG_PVS1_8 =
        dplyr::if_else(
          .data$LOSS_OF_FUNCTION == T &
            stringr::str_detect(.data$CONSEQUENCE, "_donor|_acceptor") &
            .data$LAST_INTRON == T & .data$CPG_MOD == "LoF" &
            (.data$ACMG_PM2_1 == TRUE | .data$ACMG_PM2_2 == TRUE),
          TRUE, FALSE, FALSE
        )
    ) |>
    dplyr::mutate(
      ACMG_PVS1_9 =
        dplyr::if_else(
          .data$LOSS_OF_FUNCTION == T &
            stringr::str_detect(.data$CONSEQUENCE, "_donor|_acceptor") &
            .data$LAST_INTRON == F & (is.na(.data$CPG_MOD) | .data$CPG_MOD != "LoF") &
            (.data$ACMG_PM2_1 == TRUE | .data$ACMG_PM2_2 == TRUE),
          TRUE, FALSE, FALSE
        )
    ) |>
    dplyr::mutate(
      ACMG_PVS1_10 =
        dplyr::if_else(
          .data$SPLICE_DONOR_RELEVANT == T & .data$ACMG_PP3 == TRUE &
            (.data$ACMG_PM2_1 == TRUE | .data$ACMG_PM2_2 == TRUE),
          TRUE, FALSE, FALSE
        )
    )



  # Assign logical ACMG evidence indicators
  # # TODO - BA1 -  exceptions for high population germline frequency
  #  (gnomAD) - HFE/SERPINA1

  ## Assign logical ACMG evidence indicator
  # PM4 - Protein length changes (in non-repetitive regions) due to
  # inframe indels or nonstop variant of genes that harbor variants with
  # a dominant mode of inheritance
  #
  # PPC1 - Protein length changes (in non-repetitive regions) due to
  # inframe indels or nonstop variant of genes that harbor variants with a
  # recessive mode of inheritance (and unknown CPG_MOI) - PPC1
  if ("RMSK_HIT" %in% colnames(calls)) {
    calls <- calls |>
      dplyr::mutate(
        ACMG_PM4 =
          dplyr::if_else(
            stringr::str_detect(
              .data$CONSEQUENCE, "stop_lost|inframe_deletion|inframe_insertion"
            ) &
              is.na(.data$RMSK_HIT) & .data$CPG_MOI == "AD",
            TRUE, FALSE, FALSE
          )
      ) |>
      dplyr::mutate(
        ACMG_PPC1 =
          dplyr::if_else(
            stringr::str_detect(
              .data$CONSEQUENCE, "stop_lost|inframe_deletion|inframe_insertion"
            ) &
              is.na(.data$RMSK_HIT) & (.data$CPG_MOI == "AR" | is.na(.data$CPG_MOI)),
            TRUE, FALSE, FALSE
          )
      )
  }

  ## Assign logical ACMG evidence indicator
  # ACMG_PP2 - Missense variant in a gene that has a relatively low rate
  # of benign missense variation and where missense variants are a
  # common mechanism of disease
  calls <- calls |>
    dplyr::mutate(
      ACMG_PP2 =
        dplyr::if_else(
          (is.na(.data$BENIGN_MISSENSE_FRAC) | .data$BENIGN_MISSENSE_FRAC <= 0.1) &
            (is.na(.data$PATH_TRUNC_FRAC) | .data$PATH_TRUNC_FRAC < 0.5) &
            stringr::str_detect(.data$CONSEQUENCE, "^missense_variant"),
          TRUE, FALSE, FALSE
        )
    )

  ## Assign a logical ACMG evidence indicator
  # ACMG_BP1 - Missense variant in a gene for which primarily truncating
  # variants (> 90%, as given in Maxwell et al.) are known to cause disease
  calls <- calls |>
    dplyr::mutate(
      ACMG_BP1 =
        dplyr::if_else(.data$PATH_TRUNC_FRAC > 0.90 &
                         stringr::str_detect(.data$CONSEQUENCE, "^missense_variant"),
                       TRUE, FALSE, FALSE
        )
    )

  ## Assign a logical ACMG evidence indicator
  # ACMG_BP7 - Silent/intronic variant outside of the splice site consensus
  calls <- calls |>
    dplyr::mutate(
      ACMG_BP7 =
        dplyr::if_else((
          (as.integer(.data$INTRON_POSITION) < 0 & as.integer(.data$INTRON_POSITION) < -3) |
            (as.integer(.data$INTRON_POSITION) > 0 & as.integer(.data$INTRON_POSITION) > 6) |
            (as.integer(.data$EXON_POSITION) < 0 & as.integer(.data$EXON_POSITION) < -2) |
            (as.integer(.data$EXON_POSITION) > 0 & as.integer(.data$EXON_POSITION) > 1)) &
            stringr::str_detect(
              .data$CONSEQUENCE,
              "^(synonymous_variant|intron_variant|splice_region_variant)"
            ),
          TRUE, FALSE, FALSE
        )
    )

  ## Assign a logical ACMG evidence indicator
  # ACMG_BP3 - Variants in promoter or untranslated regions
  calls <- calls |>
    dplyr::mutate(
      ACMG_BP3 =
        dplyr::if_else(
          stringr::str_detect(
            .data$CONSEQUENCE,
            "^(downstream|upstream|5_prime_UTR_variant|3_prime_UTR_variant)"
          ),
          TRUE, FALSE, FALSE
        )
    )


  ## Assign logical ACMG evidence indicators
  # ACMG_PS1 - coinciding with known pathogenic missense variants
  # (yet with different nucleotide change)
  # ACMG_PM5 - occurs at the same codon as a known pathogenic missense variant
  # ACMG_BSC1 - coinciding with known benign missense variants
  # ACMG_BMC1 - occurs at the same codon as a known benign missense variant

  calls$ACMG_PM5 <- FALSE
  calls$ACMG_BMC1 <- FALSE
  calls$ACMG_PS1 <- FALSE
  calls$ACMG_BSC1 <- FALSE

  calls <- calls |>
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
    dplyr::left_join(
      pathogenic_codons, by = c("ENTREZGENE","CODON")) |>
    dplyr::left_join(
      pathogenic_peptide_changes, by = c("ENTREZGENE","HGVSP")) |>
    dplyr::left_join(
      benign_peptide_changes, by = c("ENTREZGENE","HGVSP")) |>
    dplyr::left_join(
      benign_codons, by = c("ENTREZGENE","CODON")) |>
    dplyr::mutate(ACMG_PM5 = dplyr::if_else(
      !is.na(.data$PATHOGENIC_CODON) &
      !is.na(.data$VAR_ID_PATH_CODON) &
        !stringr::str_detect(
          .data$VAR_ID, .data$VAR_ID_PATH_CODON),
      TRUE, FALSE, FALSE
    )) |>
    dplyr::mutate(ACMG_BMC1 = dplyr::if_else(
      !is.na(.data$BENIGN_CODON) &
        !is.na(.data$VAR_ID_BENIGN_CODON) &
        !stringr::str_detect(
          .data$VAR_ID, .data$VAR_ID_BENIGN_CODON),
      TRUE, FALSE, FALSE
    )) |>
    dplyr::mutate(ACMG_PS1 = dplyr::if_else(
      !is.na(.data$PATHOGENIC_PEPTIDE_CHANGE) &
        !is.na(.data$VAR_ID_PATH_CHANGE) &
        !stringr::str_detect(
          .data$VAR_ID, .data$VAR_ID_PATH_CHANGE),
      TRUE, FALSE, FALSE)) |>
    dplyr::mutate(ACMG_BSC1 = dplyr::if_else(
      !is.na(.data$BENIGN_PEPTIDE_CHANGE) &
        !is.na(.data$VAR_ID_BENIGN_CHANGE) &
        !stringr::str_detect(
          .data$VAR_ID, .data$VAR_ID_BENIGN_CHANGE),
      TRUE, FALSE, FALSE))

  ## if previously found coinciding with pathogenic variant (ACMG_PS1),
  # set ACMG_PM5 to false
  calls <- calls |>
    dplyr::mutate(
      ACMG_PM5 =
        dplyr::case_when(
          .data$ACMG_PM5 == T &
            .data$ACMG_PS1 == T ~ FALSE,
          TRUE ~ as.logical(.data$ACMG_PM5)
        )
    ) |>
    ## if previously found coinciding with benign variant (ACMG_BSC1),
    ##  set ACMG_BMC1 to false
    dplyr::mutate(
      ACMG_BMC1 =
        dplyr::case_when(
          .data$ACMG_BMC1 == T &
            .data$ACMG_BSC1 == T ~ FALSE,
          TRUE ~ as.logical(.data$ACMG_BMC1)
        )
    )

  ## Assign logical ACMG level
  # PM1 - missense variant in a somatic mutation hotspot as
  # determined by cancerhotspots.org (v2)
  calls$ACMG_PM1 <- FALSE
  if(NROW(calls[!is.na(calls$MUTATION_HOTSPOT),]) > 0){
    calls <- calls |>
      tidyr::separate(
        .data$MUTATION_HOTSPOT,
        c("hotspot_region", "hotspot_entrezgene",
          "hotspot_symbol", "hotspot_codon",
          "hotspot_aa", "hotspot_pvalue"),
        sep = "\\|", remove = F, extra = "drop"
      ) |>
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
      dplyr::mutate(
        ACMG_PM1 =
          dplyr::if_else(
            !is.na(.data$hotspot_codon) &
              !is.na(.data$hotspot_entrezgene) &
              !is.na(.data$CODON) &
              !is.na(.data$ENTREZGENE) &
              .data$hotspot_entrezgene == .data$ENTREZGENE &
              .data$hotspot_codon == .data$CODON,
            TRUE, FALSE
          )
      )
  }

  calls <- calls |>
    pcgrr::remove_cols_from_df(
      cnames = c(
        "PATHOGENIC_CODON",
        "BENIGN_CODON",
        "PATHOGENIC_PEPTIDE_CHANGE",
        "BENIGN_PEPTIDE_CHANGE",
        "VAR_ID_PATH_CHANGE",
        "VAR_ID_BENIGN_CHANGE",
        "VAR_ID_PATH_CODON",
        "VAR_ID_BENIGN_CODON",
        "CODON",
        "gad_af",
        "hotspot_region",
        "hotspot_entrezgene",
        "hotspot_symbol",
        "hotspot_codon",
        "hotspot_aa",
        "hotspot_pvalue"
      )
    ) |>
    dplyr::distinct()

  return(calls)
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
      dplyr::filter(is.na(.data$gnomADe_AF) |
                      .data$gnomADe_AF <=
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
      dplyr::mutate(CPSR_CLASSIFICATION_SOURCE = "CPSR_ACMG")

    non_clinvar_calls[["class4"]] <- cpg_calls_non_clinvar |>
      dplyr::filter(.data$CPSR_CLASSIFICATION == "Likely_Pathogenic") |>
      dplyr::mutate(CPSR_CLASSIFICATION_SOURCE = "CPSR_ACMG")

    non_clinvar_calls[["class3"]] <- cpg_calls_non_clinvar |>
      dplyr::filter(.data$CPSR_CLASSIFICATION == "VUS") |>
      dplyr::mutate(CPSR_CLASSIFICATION_SOURCE = "CPSR_ACMG")

    non_clinvar_calls[["class2"]] <- cpg_calls_non_clinvar |>
      dplyr::filter(.data$CPSR_CLASSIFICATION == "Likely_Benign") |>
      dplyr::mutate(CPSR_CLASSIFICATION_SOURCE = "CPSR_ACMG")

    non_clinvar_calls[["class1"]] <- cpg_calls_non_clinvar |>
      dplyr::filter(.data$CPSR_CLASSIFICATION == "Benign") |>
      dplyr::mutate(CPSR_CLASSIFICATION_SOURCE = "CPSR_ACMG")

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
      if (config[["variant_classification"]][["classify_all"]] == 0) {
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
          c("SAMPLE_ID",
            "GENOMIC_CHANGE",
            "VAR_ID",
            "GENOTYPE",
            "CPSR_CLASSIFICATION_SOURCE",
            "GENOME_VERSION",
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
            "NULL_VARIANT",
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
        c("gnomADe_AF",
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

    for(c in c("class5","class4","class3","class2","class1")){
      if(NROW(snv_indel_report[["variant_set"]][[c]]) > 0){
        snv_indel_report[["variant_set"]][["tsv"]] <-
          snv_indel_report[["variant_set"]][["tsv"]] |>
          dplyr::bind_rows(
            snv_indel_report[["variant_set"]][[c]]
          )
      }
    }

    return(snv_indel_report)
  }
