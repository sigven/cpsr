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
          dplyr::filter(
            .data$ASSERTION_AUTHORITY == "ClinVar") |>
          nrow()
        num_rows_other <- t1 |>
          dplyr::filter(
            .data$ASSERTION_AUTHORITY == "CPSR") |>
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
    "Summarising insilico variant effect predictions (dbNSFP)"
  )

  for (v in c(
    "CALLED",
    "DAMAGING",
    "TOLERATED",
    "SPLICING_NEUTRAL",
    "SPLICING_AFFECTED")) {
    cpg_calls[, paste0("N_INSILICO_", v)] <- 0
  }

  for (algo in cpsr::insilico_path_predictors) {
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
#' @param calls data frame with all calls found
#'
#' @export
retrieve_secondary_calls <- function(calls) {
  assertable::assert_colnames(
    calls,
    colnames = c(
      "CPG_SOURCE",
      "CLASSIFICATION",
      "ASSERTION_AUTHORITY",
      "CLINVAR_GOLD_STARS",
      "PRIMARY_TARGET",
      "GENOTYPE",
      "SYMBOL",
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
        stringr::str_detect(.data$CPG_SOURCE, "ACMG_SF") &
        .data$PRIMARY_TARGET == FALSE &
        !is.na(.data$ASSERTION_AUTHORITY) &
        .data$ASSERTION_AUTHORITY == "ClinVar" &
        !is.na(.data$CLASSIFICATION) &
        stringr::str_detect(
          .data$CLASSIFICATION, "Pathogenic"
        )
    ) |>
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
        .data$LOSS_OF_FUNCTION == TRUE)) |>
    ## only homozygotes p.Cys282Tyr HFE carriers
    dplyr::filter(.data$SYMBOL != "HFE" |
      (.data$SYMBOL == "HFE" &
        (!is.na(.data$PROTEIN_CHANGE) &
          stringr::str_detect(.data$PROTEIN_CHANGE, "Cys282Tyr")) &
        .data$GENOTYPE == "hom_alt"))

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
            genes_lacking_twohit_evidence,
            by = "SYMBOL"
          )
      }
    }
  }

  return(secondary_calls)
}

#' Regular expression matching cancer-related ClinVar trait names
#'
#' Applied (case-insensitive, to lower-cased trait names) when determining
#' whether a ClinVar record is attributed to a cancer-related phenotype.
#' Moved from pcgrr (\code{pcgrr::cancer_phenotypes_regex}), where it was
#' used by CPSR only.
#'
#' @keywords internal
cancer_phenotypes_regex <-
  paste0("(cancer|carcinoma|tumor|neoplasm|gangliom|lymphom|leukem|meningiom",
         "|blastoma|melanom|chordom|adenoma|sarcom|mesotheli|ependymom|glioma",
         "|neurofibro|keratoacan|nevus|brca|polyposis|myelodysplastic|cowden",
         "|gardner|noonan|fanconi|carney|bullosa|schwanno|li-fraumeni|xeroderma",
         "|leiomyom|muir-|nijmegen|neoplasia|trichoepithelioma|brooke|turcot",
         "|exostos|lynch|drash|wilm|perlman|fibrofolliculomas|hippel|hamartom",
         "|bloom|werner|peutz|tuberous|angiomyolipoma",
         "|thymoma|dihydropyrimidine dehydrogenase deficiency",
         "|lymphoproliferative|stat3|teratoma|thrombocytop|tp63|weaver",
         "|pheochromo|gorlin|telangiectasia|hemangiom|osteochondro|response",
         "|polg-related|ras-associated|dyskeratosis",
         ## bone marrow failure / myeloid predisposition
         "|bone marrow failure|marrow hypocellularity|cytopenia|aplastic anemia",
         "|neutropenia|telomere biology|myeloid",
         ## Nijmegen breakage syndrome (OMIM name)
         "|microcephaly, normal intelligence and immunodeficiency",
         "|tyrosine kinase inhibitor|fluorouracil|methotrexate|capecitabine",
         "|imatinib|erlotinib|gefitinib|sunitinib|vemurafenib|crizotinib|cetuximab",
         "|dabrafenib|osimertinib|afatinib|neratinib|lapatinib|ponatinib|regorafenib",
         "|vandetanib|bosutinib|axitinib|cabozantinib|brigatinib|ripretinib|lenvatinib",
         "|midostaurin|duvelisib|copanlisib|alpelisib|larotrectinib|entrectinib",
         "|programmed death ligand|adp-ribose|doxorubicin|cisplatin|carboplatin",
         "|familial adenomatous polyposis|hereditary breast and ovarian cancer",
         "|hereditary nonpolyposis colorectal cancer|multiple endocrine neoplasia",
         "|waardenburg|beckwith|birt-hogg|diamond|frasier",
         "|infantile myofibromatosis|proteus|rothmund|russel)")

#' Regular expression matching uninformative ClinVar trait names
#'
#' Trait names that carry no information about the phenotype (e.g.
#' 'not provided'). Records with only such traits are given an unknown
#' (NA) cancer phenotype status, rather than a non-cancer status.
#'
#' @keywords internal
clinvar_uninformative_traits_regex <-
  paste0("^(not provided|not specified|none provided|",
         "inborn genetic diseases|see cases)$")

#' Function that determines whether variants in cancer predisposition genes
#' are linked to cancer-related conditions according to ClinVar
#'
#' Each ClinVar trait of a variant (trait names in \emph{CLINVAR_TRAITS},
#' and UMLS concept names of \emph{CLINVAR_UMLS_CUI}) is given one of three
#' states: cancer-related (1), not cancer-related (0) or unknown (NA).
#' A trait is cancer-related if it matches the cancer phenotype regex, maps
#' to an OncoTree concept, or is a '<GENE>-related disorder' for a primary
#' target gene of the virtual panel. A trait is unknown if it is
#' uninformative (e.g. 'not provided') or if a UMLS concept cannot be mapped
#' to a name. The variant-level state in \emph{CLINVAR_PHENOTYPE_CANCER} is
#' 1 if any trait is cancer-related, otherwise 0 if any trait is not
#' cancer-related, otherwise NA (also NA for variants not found in ClinVar).
#'
#' @param cpg_calls data frame with variants in predisposition_genes
#' @param ref_data object with PCGR/CPSR reference data
#' @param primary_target_entrez character vector with Entrez gene
#' identifiers of the primary targets (virtual panel genes); used to
#' recognize gene-named conditions ('<GENE>-related disorder') as
#' cancer-related
#'
#' @export
is_clinvar_cancer_phenotype <- function(
    cpg_calls,
    ref_data,
    primary_target_entrez = NULL) {

  cpg_calls$CLINVAR_PHENOTYPE <- NULL
  cpg_calls$CLINVAR_PHENOTYPE_CANCER <- NULL

  if (NROW(cpg_calls) == 0 ||
      !all(c("VAR_ID", "CLINVAR_MSID", "CLINVAR_UMLS_CUI",
             "CLINVAR_TRAITS") %in% colnames(cpg_calls))) {
    cpg_calls$CLINVAR_PHENOTYPE <- rep(NA_character_, NROW(cpg_calls))
    cpg_calls$CLINVAR_PHENOTYPE_CANCER <- rep(NA_integer_, NROW(cpg_calls))
    return(cpg_calls)
  }

  oncotree_cuis <- unique(ref_data[["phenotype"]][["oncotree"]]$CUI)
  umls_map <- ref_data[["phenotype"]][["umls"]] |>
    dplyr::filter(.data$MAIN_TERM == TRUE) |>
    dplyr::select(c("CUI", "CUI_NAME")) |>
    dplyr::distinct(.data$CUI, .keep_all = TRUE)

  na_if_empty <- function(x) {
    x <- as.character(x)
    dplyr::if_else(is.na(x) | x %in% c("", "."), NA_character_, x)
  }

  clinvar_calls <- cpg_calls |>
    dplyr::mutate(
      .msid = na_if_empty(.data$CLINVAR_MSID),
      .traits = na_if_empty(.data$CLINVAR_TRAITS),
      .cuis = na_if_empty(.data$CLINVAR_UMLS_CUI)) |>
    dplyr::filter(!is.na(.data$.msid) &
                    (!is.na(.data$.traits) | !is.na(.data$.cuis)))

  if (NROW(clinvar_calls) == 0) {
    cpg_calls$CLINVAR_PHENOTYPE <- rep(NA_character_, NROW(cpg_calls))
    cpg_calls$CLINVAR_PHENOTYPE_CANCER <- rep(NA_integer_, NROW(cpg_calls))
    return(cpg_calls)
  }

  gene_cols <- intersect(c("SYMBOL", "ENTREZGENE"), colnames(clinvar_calls))
  var_genes <- clinvar_calls |>
    dplyr::select(dplyr::all_of(c("VAR_ID", gene_cols))) |>
    dplyr::distinct(.data$VAR_ID, .keep_all = TRUE)

  ## Trait names as provided by ClinVar (semicolon-separated)
  traits_by_name <- clinvar_calls |>
    dplyr::filter(!is.na(.data$.traits)) |>
    dplyr::select(c("VAR_ID", "TRAIT" = ".traits")) |>
    tidyr::separate_rows("TRAIT", sep = ";") |>
    dplyr::mutate(TRAIT = stringr::str_trim(.data$TRAIT),
                  CUI = NA_character_,
                  SOURCE = "trait")

  ## UMLS concepts (comma-separated), mapped to concept names
  traits_by_cui <- clinvar_calls |>
    dplyr::filter(!is.na(.data$.cuis)) |>
    dplyr::select(c("VAR_ID", "CUI" = ".cuis")) |>
    tidyr::separate_rows("CUI", sep = ",") |>
    dplyr::mutate(CUI = stringr::str_trim(.data$CUI)) |>
    dplyr::left_join(umls_map, by = "CUI") |>
    dplyr::rename(TRAIT = "CUI_NAME") |>
    dplyr::mutate(SOURCE = "cui")

  traits <- dplyr::bind_rows(traits_by_name, traits_by_cui) |>
    dplyr::filter(!(is.na(.data$TRAIT) & is.na(.data$CUI)) &
                    (is.na(.data$TRAIT) | .data$TRAIT != "")) |>
    dplyr::distinct() |>
    dplyr::left_join(var_genes, by = "VAR_ID")

  trait_lc <- tolower(traits$TRAIT)

  ## '<GENE>-related disorder(s)' for a primary target gene
  gene_named <- rep(FALSE, NROW(traits))
  if ("SYMBOL" %in% colnames(traits) && !is.null(primary_target_entrez) &&
      "ENTREZGENE" %in% colnames(traits)) {
    is_target <- as.character(traits$ENTREZGENE) %in%
      as.character(primary_target_entrez)
    symbol_lc <- tolower(dplyr::coalesce(as.character(traits$SYMBOL), ""))
    trait_or_empty <- dplyr::coalesce(trait_lc, "")
    gene_named <- is_target & nchar(symbol_lc) > 0 &
      (startsWith(trait_or_empty, paste0(symbol_lc, "-related")) |
         startsWith(trait_or_empty, paste0(symbol_lc, " related")))
  }

  traits$TRAIT_CANCER <- dplyr::case_when(
    !is.na(traits$CUI) & traits$CUI %in% oncotree_cuis ~ 1L,
    is.na(trait_lc) ~ NA_integer_,
    stringr::str_detect(trait_lc, clinvar_uninformative_traits_regex) ~
      NA_integer_,
    stringr::str_detect(trait_lc, cancer_phenotypes_regex) ~ 1L,
    gene_named ~ 1L,
    TRUE ~ 0L
  )

  ## Variant-level state: 1 if any trait is cancer-related, 0 if any
  ## (informative) trait is not cancer-related, NA otherwise
  phenotype_state <- traits |>
    dplyr::group_by(.data$VAR_ID) |>
    dplyr::summarise(
      CLINVAR_PHENOTYPE_CANCER = dplyr::case_when(
        any(.data$TRAIT_CANCER == 1L, na.rm = TRUE) ~ 1L,
        any(.data$TRAIT_CANCER == 0L, na.rm = TRUE) ~ 0L,
        TRUE ~ NA_integer_),
      ## ClinVar trait names, falling back to UMLS concept names
      CLINVAR_PHENOTYPE = {
        nm <- unique(.data$TRAIT[.data$SOURCE == "trait" &
                                   !is.na(.data$TRAIT)])
        if (length(nm) == 0) {
          nm <- unique(.data$TRAIT[!is.na(.data$TRAIT)])
        }
        if (length(nm) == 0) NA_character_ else paste(nm, collapse = "; ")
      },
      .groups = "drop"
    )

  cpg_calls <- cpg_calls |>
    dplyr::left_join(phenotype_state, by = "VAR_ID")

  pcgrr::log4r_info(
    paste0(
      "ClinVar phenotype status of variants - cancer-related: N = ",
      sum(phenotype_state$CLINVAR_PHENOTYPE_CANCER %in% 1L),
      ", not cancer-related: N = ",
      sum(phenotype_state$CLINVAR_PHENOTYPE_CANCER %in% 0L),
      ", unspecified: N = ",
      sum(is.na(phenotype_state$CLINVAR_PHENOTYPE_CANCER))
    )
  )

  return(cpg_calls)
}
