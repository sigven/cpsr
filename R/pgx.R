## Genotype values (GENOTYPE) of carriers of a variant allele, as
## assigned in CPSR pre-processing (cyvcf2, gts012). Haploid calls
## (e.g. GT = 1 on chrX in males) are labelled 'hom_alt'
pgx_carrier_genotypes <- c("het", "hom_alt")

## CPIC allele function labels of normal function alleles
## (G6PD uses WHO classes, i.e. 'IV/Normal')
pgx_normal_function <- c("Normal function", "IV/Normal")

#' Match sample variants against CPIC allele definitions
#'
#' Sample variants (CHROM/POS/REF/ALT) are matched against CPIC allele
#' definitions (pgx_allele, reference data bundle). For each matched
#' variant, all allele definitions that include the variant are returned,
#' together with the number of allele-defining variants that are
#' assessable from VCF (N_ASSESSABLE) and present in the sample
#' (N_PRESENT).
#'
#' @param calls data frame with sample variants
#' @param pgx_allele data frame with CPIC allele definitions
#'
#' @return data frame with one row per sample variant and allele definition
#' @keywords internal
match_pgx_alleles <- function(calls, pgx_allele) {

  pgx_allele <- pgx_allele |>
    dplyr::mutate(
      N_ALLELE_VARIANTS = as.integer(.data$N_ALLELE_VARIANTS)) |>
    dplyr::select(dplyr::any_of(c(
      "SYMBOL", "CHROM", "ALLELE_DEFINITION_ID", "ALLELE_NAME",
      "N_ALLELE_VARIANTS", "CLINICAL_FUNCTION_STATUS",
      "ACTIVITY_VALUE", "VAR_ID", "VARIANT_NAME", "REFERENCE_ALLELE")))

  sample_vars <- calls |>
    ## do not consider pharmacogenomics-related variants if
    ## genotypes have not been retrieved properly
    dplyr::filter(
      !is.na(.data$GENOTYPE) &
        .data$GENOTYPE %in% pgx_carrier_genotypes
    ) |>
    dplyr::mutate(
      PGX_VAR_ID = paste(
        stringr::str_replace(.data$CHROM, "^chr", ""),
        .data$POS, .data$REF, .data$ALT, sep = "_")) |>
    dplyr::filter(.data$PGX_VAR_ID %in% pgx_allele$VAR_ID) |>
    dplyr::select(c("PGX_VAR_ID", "GENOTYPE")) |>
    dplyr::distinct()

  hit_definitions <- unique(pgx_allele$ALLELE_DEFINITION_ID[
    !is.na(pgx_allele$VAR_ID) &
      pgx_allele$VAR_ID %in% sample_vars$PGX_VAR_ID])

  pgx_allele |>
    dplyr::filter(.data$ALLELE_DEFINITION_ID %in% hit_definitions) |>
    dplyr::group_by(.data$ALLELE_DEFINITION_ID) |>
    dplyr::mutate(
      ## defining variants with a VCF representation (VAR_ID), e.g.
      ## excluding repeat alleles such as NUDT15 GAGTCG[2]
      N_ASSESSABLE = sum(!is.na(.data$VAR_ID)),
      N_PRESENT = sum(!is.na(.data$VAR_ID) &
                        .data$VAR_ID %in% sample_vars$PGX_VAR_ID),
      ALLELE_VAR_IDS = paste(sort(unique(stats::na.omit(.data$VAR_ID))),
                             collapse = ",")) |>
    dplyr::ungroup() |>
    dplyr::inner_join(sample_vars, by = c("VAR_ID" = "PGX_VAR_ID")) |>
    dplyr::rename(PGX_VAR_ID = "VAR_ID") |>
    dplyr::mutate(
      ## all assessable defining variants found in sample
      ALLELE_COMPLETE = .data$N_PRESENT == .data$N_ASSESSABLE,
      NORMAL_FUNCTION = !is.na(.data$CLINICAL_FUNCTION_STATUS) &
        .data$CLINICAL_FUNCTION_STATUS %in% pgx_normal_function,
      ALLELE_LABEL = dplyr::case_when(
        !.data$ALLELE_COMPLETE ~ paste0(
          .data$ALLELE_NAME, " (partial - ", .data$N_PRESENT, "/",
          .data$N_ALLELE_VARIANTS, " defining variants)"),
        .data$N_ASSESSABLE < .data$N_ALLELE_VARIANTS ~ paste0(
          .data$ALLELE_NAME, " (", .data$N_ASSESSABLE, "/",
          .data$N_ALLELE_VARIANTS, " defining variants assessable)"),
        TRUE ~ .data$ALLELE_NAME)) |>
    ## partial matches are ignored for variants that completely define a
    ## normal function allele (e.g. G6PD c.376A>G - 'A' (IV/Normal), also
    ## part of A- (deficient) alleles together with c.202G>A)
    dplyr::group_by(.data$PGX_VAR_ID) |>
    dplyr::filter(
      .data$ALLELE_COMPLETE |
        !any(.data$ALLELE_COMPLETE & .data$NORMAL_FUNCTION)) |>
    dplyr::ungroup()
}

#' Function that retrieves pharmacogenomic variants, i.e. variants
#' that define CPIC alleles (DPYD, TPMT, NUDT15, G6PD) with altered
#' function
#'
#' Sample variants are matched (CHROM/POS/REF/ALT) against CPIC allele
#' definitions in the reference data bundle (pgx_allele). Each matched
#' variant is annotated with the CPIC allele(s) it defines, together with
#' allele function and activity value. Alleles defined by multiple variants
#' (e.g. TPMT *3A, DPYD HapB3) are flagged as partial when not all
#' defining variants are found in the sample (phase is not assessed).
#'
#' @param calls data frame with all calls found
#' @param ref_data object with PCGR/CPSR reference data
#' @param include_normal_function logical, report also variants that
#' only define normal function alleles
#'
#' @export
retrieve_pgx_calls <- function(calls,
                               ref_data = NULL,
                               include_normal_function = FALSE) {
  assertable::assert_colnames(
    calls,
    colnames = c(
      "CHROM",
      "POS",
      "REF",
      "ALT",
      "GENOTYPE",
      "SYMBOL"
    ),
    only_colnames = F, quiet = T
  )

  pgx_allele <- ref_data[['variant']][['pgx_allele']]
  if (is.null(pgx_allele) || NROW(pgx_allele) == 0) {
    pcgrr::log4r_warn(paste0(
      "CPIC allele definitions (pgx_allele) not found in reference data ",
      "bundle - pharmacogenomic findings based on ClinVar classifications"
    ))
    return(retrieve_pgx_calls_clinvar(calls))
  }

  allele_hits <- match_pgx_alleles(calls, pgx_allele)
  if (!include_normal_function) {
    allele_hits <- allele_hits |>
      dplyr::filter(.data$NORMAL_FUNCTION == FALSE)
  }

  if (NROW(allele_hits) == 0) {
    return(calls[0, ])
  }

  allele_hits <- allele_hits |>
    dplyr::arrange(.data$PGX_VAR_ID, dplyr::desc(.data$ALLELE_COMPLETE),
                   .data$ALLELE_NAME) |>
    dplyr::group_by(.data$PGX_VAR_ID) |>
    dplyr::summarise(
      PGX_ALLELE = paste(unique(.data$ALLELE_LABEL), collapse = "; "),
      PGX_ALLELE_COMPLETE = any(.data$ALLELE_COMPLETE),
      PGX_ALLELE_FUNCTION = paste(
        unique(stats::na.omit(.data$CLINICAL_FUNCTION_STATUS)),
        collapse = "; "),
      PGX_ACTIVITY_VALUE = paste(
        unique(stats::na.omit(.data$ACTIVITY_VALUE)), collapse = "; "),
      .groups = "drop") |>
    dplyr::mutate(dplyr::across(
      c("PGX_ALLELE_FUNCTION", "PGX_ACTIVITY_VALUE"),
      ~ dplyr::na_if(.x, "")))

  pgx_calls <- calls |>
    dplyr::mutate(
      PGX_VAR_ID = paste(
        stringr::str_replace(.data$CHROM, "^chr", ""),
        .data$POS, .data$REF, .data$ALT, sep = "_")) |>
    dplyr::filter(
      !is.na(.data$GENOTYPE) &
        .data$GENOTYPE %in% pgx_carrier_genotypes) |>
    dplyr::inner_join(allele_hits, by = "PGX_VAR_ID") |>
    dplyr::select(-c("PGX_VAR_ID")) |>
    dplyr::arrange(.data$SYMBOL, .data$POS) |>
    dplyr::distinct()

  return(pgx_calls)
}

#' Infer CPIC phenotypes for pharmacogenes with detected variant alleles
#'
#' For each gene with at least one detected CPIC allele of altered
#' function, a diplotype (the pair of CPIC alleles on the two gene copies,
#' e.g. TPMT *1/*3C) is assembled from the detected alleles, and mapped to
#' a CPIC phenotype (pgx_phenotype, reference data bundle). The following
#' conventions apply:
#'
#' * Genes without detected altered-function alleles are not reported
#' * CPIC allele positions without a called variant are assumed to be
#'   homozygous reference (not confirmed from the input VCF). A gene copy
#'   without a detected allele is thus the CPIC reference allele (e.g.
#'   TPMT *1, DPYD 'Reference')
#' * Two different heterozygous alleles are assumed to be in trans (phase
#'   is not assessed)
#' * Alleles whose defining variants are a subset of another detected
#'   allele (e.g. TPMT *3B/*3C vs. *3A) are collapsed into the larger allele
#' * X-linked genes (G6PD) are hemizygous in males (one allele, any
#'   variant call counted as a single copy), and diploid in females
#' * The phenotype is not determined ('Not determined') when detected
#'   variants only partially define an allele (e.g. DPYD c.1236G>A without
#'   c.1129-5923C>G), when more altered-function alleles are found than
#'   expected from ploidy, or for X-linked genes when sample sex is unknown
#'
#' @param calls data frame with all calls found
#' @param ref_data object with PCGR/CPSR reference data
#' @param sex sample sex ('MALE', 'FEMALE' or 'UNKNOWN')
#'
#' @return data frame with one row per gene
#' @export
assign_pgx_phenotypes <- function(calls, ref_data = NULL,
                                  sex = "UNKNOWN") {

  pgx_allele <- ref_data[['variant']][['pgx_allele']]
  pgx_phenotype <- ref_data[['variant']][['pgx_phenotype']]
  if (NROW(pgx_allele) == 0 || NROW(pgx_phenotype) == 0) {
    return(data.frame())
  }

  allele_hits <- match_pgx_alleles(calls, pgx_allele) |>
    dplyr::filter(.data$NORMAL_FUNCTION == FALSE)
  if (NROW(allele_hits) == 0) {
    return(data.frame())
  }

  ## labels of detected variants (HGVSc from sample annotation, or the
  ## CPIC variant name), e.g. 'c.1905+1G>A'
  var_labels <- calls |>
    dplyr::mutate(
      PGX_VAR_ID = paste(
        stringr::str_replace(.data$CHROM, "^chr", ""),
        .data$POS, .data$REF, .data$ALT, sep = "_"),
      HGVSC_SHORT = if ("HGVSc" %in% colnames(calls))
        stringr::str_replace(.data$HGVSc, "^[^:]*:", "") else NA_character_) |>
    dplyr::filter(.data$PGX_VAR_ID %in% allele_hits$PGX_VAR_ID) |>
    dplyr::select(c("PGX_VAR_ID", "HGVSC_SHORT")) |>
    dplyr::distinct(.data$PGX_VAR_ID, .keep_all = TRUE)
  allele_hits <- allele_hits |>
    dplyr::left_join(var_labels, by = "PGX_VAR_ID") |>
    dplyr::mutate(
      VARIANT_LABEL = dplyr::coalesce(
        dplyr::na_if(.data$HGVSC_SHORT, ""),
        if ("VARIANT_NAME" %in% colnames(allele_hits))
          .data$VARIANT_NAME else NA_character_,
        .data$PGX_VAR_ID))

  gene_phenotypes <- data.frame()
  for (gene in sort(unique(allele_hits$SYMBOL))) {
    gene_phenotypes <- gene_phenotypes |>
      dplyr::bind_rows(
        infer_pgx_gene_phenotype(
          gene_hits = dplyr::filter(allele_hits, .data$SYMBOL == gene),
          gene_phenotypes = dplyr::filter(pgx_phenotype,
                                          .data$SYMBOL == gene),
          sex = sex))
  }

  return(gene_phenotypes)
}

#' Infer CPIC phenotype for a single gene
#'
#' @param gene_hits allele matches (match_pgx_alleles) for the gene,
#' altered function alleles only
#' @param gene_phenotypes CPIC function/activity value to phenotype
#' mappings for the gene
#' @param sex sample sex ('MALE', 'FEMALE' or 'UNKNOWN')
#'
#' @keywords internal
infer_pgx_gene_phenotype <- function(gene_hits, gene_phenotypes,
                                     sex = "UNKNOWN") {

  gene <- unique(gene_hits$SYMBOL)
  notes <- c()
  ## X-linked genes: hemizygous in males (one allele), diploid in females
  x_linked <- all(stringr::str_replace(gene_hits$CHROM, "^chr", "") == "X")
  hemizygous <- x_linked && identical(sex, "MALE")
  ploidy <- if (hemizygous) 1 else 2
  if (!"VARIANT_LABEL" %in% colnames(gene_hits)) {
    gene_hits$VARIANT_LABEL <- gene_hits$PGX_VAR_ID
  }
  ## CPIC name of the reference allele (e.g. TPMT '*1')
  reference_allele <- if ("REFERENCE_ALLELE" %in% colnames(gene_hits))
    stats::na.omit(unique(gene_hits$REFERENCE_ALLELE))[1] else NA
  if (is.na(reference_allele)) reference_allele <- "Reference"

  ## detected variants (e.g. 'c.719A>G') and their genotypes (same order,
  ## e.g. 'heterozygous'), separated by '; '
  variant_genotypes <- gene_hits |>
    dplyr::distinct(.data$PGX_VAR_ID, .data$VARIANT_LABEL, .data$GENOTYPE) |>
    dplyr::mutate(
      GT = dplyr::case_when(
        hemizygous ~ "hemizygous",
        .data$GENOTYPE == "hom_alt" ~ "homozygous",
        .data$GENOTYPE == "het" ~ "heterozygous",
        TRUE ~ .data$GENOTYPE)) |>
    dplyr::distinct(.data$VARIANT_LABEL, .data$GT) |>
    dplyr::arrange(.data$VARIANT_LABEL)

  result <- data.frame(
    SYMBOL = gene,
    PGX_VARIANTS = paste(variant_genotypes$VARIANT_LABEL, collapse = "; "),
    PGX_GENOTYPES = paste(variant_genotypes$GT, collapse = "; "),
    PGX_DIPLOTYPE = NA_character_,
    PGX_PHENOTYPE = "Not determined",
    PGX_ACTIVITY_SCORE = NA_character_,
    EHR_PRIORITY = NA_character_,
    CONSULTATION_TEXT = NA_character_,
    PGX_PHENOTYPE_NOTE = NA_character_)

  ## alleles where all assessable defining variants are present
  alleles <- gene_hits |>
    dplyr::filter(.data$ALLELE_COMPLETE == TRUE) |>
    dplyr::group_by(dplyr::across(c(
      "ALLELE_DEFINITION_ID", "ALLELE_NAME", "ALLELE_LABEL",
      "CLINICAL_FUNCTION_STATUS", "ACTIVITY_VALUE", "ALLELE_VAR_IDS"))) |>
    dplyr::summarise(
      COPIES = dplyr::if_else(
        !hemizygous & all(.data$GENOTYPE == "hom_alt"), 2, 1),
      HET_CALL = any(.data$GENOTYPE == "het"),
      .groups = "drop")

  ## collapse alleles whose defining variants are contained in another
  ## detected allele (identical variant sets: keep first)
  var_sets <- strsplit(alleles$ALLELE_VAR_IDS, ",")
  keep <- vapply(seq_along(var_sets), function(i) {
    !any(vapply(seq_along(var_sets), function(j) {
      j != i && all(var_sets[[i]] %in% var_sets[[j]]) &&
        (length(var_sets[[j]]) > length(var_sets[[i]]) || j < i)
    }, logical(1)))
  }, logical(1))
  ## phase matters when the variants of a collapsed allele could also
  ## represent two separate alleles (e.g. TPMT *3B + *3C vs. *3A)
  if (sum(!keep) >= 2) {
    notes <- c(notes, paste0(
      "Phase unknown - the detected variants are assumed to be on the ",
      "same gene copy; they may alternatively represent ",
      paste(alleles$ALLELE_NAME[!keep], collapse = " and "),
      " on different gene copies"))
  }
  alleles <- alleles[keep, ]

  ## variants that only partially define an allele, and are not
  ## explained by a complete allele
  explained_vars <- unlist(strsplit(alleles$ALLELE_VAR_IDS, ","))
  partial <- gene_hits |>
    dplyr::filter(.data$ALLELE_COMPLETE == FALSE &
                    !(.data$PGX_VAR_ID %in% explained_vars))

  if (NROW(partial) > 0) {
    result$PGX_PHENOTYPE_NOTE <- paste0(
      "Phenotype not determined - variant(s) only partially define ",
      "allele(s): ", paste(unique(partial$ALLELE_LABEL), collapse = "; "))
    return(result)
  }
  if (x_linked && !(sex %in% c("MALE", "FEMALE"))) {
    result$PGX_DIPLOTYPE <- paste(alleles$ALLELE_LABEL, collapse = ", ")
    result$PGX_PHENOTYPE_NOTE <- paste0(
      "Phenotype not determined - ", gene, " is X-linked, and sample ",
      "sex is unknown (set with 'cpsr --sex')")
    return(result)
  }
  if (hemizygous && any(alleles$HET_CALL)) {
    notes <- c(notes, paste0(
      "Heterozygous variant call(s) on chrX in a male sample - ",
      "interpreted as hemizygous"))
  }

  n_copies <- sum(alleles$COPIES)
  if (n_copies > ploidy) {
    result$PGX_DIPLOTYPE <- paste(alleles$ALLELE_LABEL, collapse = ", ")
    result$PGX_PHENOTYPE_NOTE <- paste0(
      "Phenotype not determined - more altered function alleles detected ",
      "than expected (", ploidy, ") - phase unknown")
    return(result)
  }

  ## allele values used in CPIC phenotype lookup keys: activity values
  ## (e.g. DPYD) or allele function (e.g. TPMT, NUDT15)
  activity_based <- any(!is.na(gene_phenotypes$ACTIVITY_VALUE_1) &
                          gene_phenotypes$ACTIVITY_VALUE_1 != "n/a")
  normal_row <- gene_phenotypes |>
    dplyr::filter(.data$FUNCTION_1 %in% pgx_normal_function &
                    .data$FUNCTION_2 %in% pgx_normal_function)
  normal_value <- if (activity_based) normal_row$ACTIVITY_VALUE_1[1] else
    normal_row$FUNCTION_1[1]
  allele_values <- if (activity_based) alleles$ACTIVITY_VALUE else
    alleles$CLINICAL_FUNCTION_STATUS

  diplotype_values <- c(rep(allele_values, alleles$COPIES),
                        rep(normal_value, ploidy - n_copies))
  ## CPIC-style diplotype - reference allele first, e.g. '*1/*3C'
  diplotype_labels <- c(rep(reference_allele, ploidy - n_copies),
                        sort(rep(alleles$ALLELE_LABEL, alleles$COPIES)))
  result$PGX_DIPLOTYPE <- paste(diplotype_labels, collapse = "/")
  if (hemizygous) {
    result$PGX_DIPLOTYPE <- paste0(result$PGX_DIPLOTYPE, " (hemizygous)")
  }

  if (NROW(alleles) == 2 && ploidy == 2) {
    notes <- c(notes, paste0(
      "Phase unknown - the two alleles are assumed to be on different ",
      "gene copies (in trans)"))
  }

  if (any(is.na(diplotype_values))) {
    result$PGX_PHENOTYPE_NOTE <- paste(c(
      "Phenotype not determined - allele activity value/function missing",
      notes), collapse = ". ")
    return(result)
  }

  lookup_key <- canonical_pgx_key(
    paste0(names(table(diplotype_values)), ":", table(diplotype_values)))
  phenotype_match <- gene_phenotypes |>
    dplyr::filter(vapply(
      strsplit(.data$LOOKUP_KEY, ";"), canonical_pgx_key,
      character(1)) == lookup_key)

  if (NROW(phenotype_match) == 0) {
    result$PGX_PHENOTYPE_NOTE <- paste(c(
      paste0("Phenotype not determined - no CPIC phenotype for ",
             "allele combination ", lookup_key),
      notes), collapse = ". ")
    return(result)
  }

  result$PGX_PHENOTYPE <- phenotype_match$PHENOTYPE[1]
  result$PGX_ACTIVITY_SCORE <- phenotype_match$ACTIVITY_SCORE[1]
  result$EHR_PRIORITY <- phenotype_match$EHR_PRIORITY[1]
  result$CONSULTATION_TEXT <- phenotype_match$CONSULTATION_TEXT[1]
  if (length(notes) > 0) {
    result$PGX_PHENOTYPE_NOTE <- paste(notes, collapse = ". ")
  }

  return(result)
}

#' Map inferred pharmacogene phenotypes to CPIC prescribing
#' recommendations
#'
#' Recommendations are retrieved (pgx_guideline, reference data bundle)
#' for all drugs whose CPIC lookup key involves at least one gene with an
#' inferred phenotype (assign_pgx_phenotypes). For multi-gene
#' recommendations (e.g. thiopurines - TPMT and NUDT15), genes without
#' detected altered-function alleles are assumed to be normal
#' metabolizers, which is stated in PGX_NOTE.
#'
#' @param gene_phenotypes data frame with inferred gene phenotypes
#' @param ref_data object with PCGR/CPSR reference data
#'
#' @return data frame with one row per drug (and CPIC population)
#' @export
assign_pgx_recommendations <- function(gene_phenotypes, ref_data = NULL) {

  pgx_guideline <- ref_data[['variant']][['pgx_guideline']]
  pgx_phenotype <- ref_data[['variant']][['pgx_phenotype']]
  if (NROW(gene_phenotypes) == 0 || NROW(pgx_guideline) == 0) {
    return(data.frame())
  }

  pgx_guideline <- pgx_guideline |>
    dplyr::filter(as.logical(.data$CPSR_CALLABLE) == TRUE &
                    !is.na(.data$LOOKUP_KEY)) |>
    dplyr::mutate(
      KEY_GENES = stringr::str_replace_all(.data$LOOKUP_KEY, ":[^;]*", ""))

  recommendations <- data.frame()
  for (drug in sort(unique(pgx_guideline$DRUG_NAME))) {
    drug_recs <- dplyr::filter(pgx_guideline, .data$DRUG_NAME == drug)
    key_genes <- unique(unlist(strsplit(drug_recs$KEY_GENES, ";")))
    if (!any(key_genes %in% gene_phenotypes$SYMBOL)) {
      next
    }

    notes <- c()
    key_values <- c()
    undetermined <- c()
    for (gene in key_genes) {
      ## lookup values in CPIC recommendations: activity score
      ## (e.g. DPYD) or phenotype (e.g. TPMT, NUDT15)
      gene_key_values <- stringr::str_match(
        drug_recs$LOOKUP_KEY, paste0(gene, ":([^;]*)"))[, 2]
      activity_based <- all(stringr::str_detect(
        stats::na.omit(gene_key_values), "^[0-9.]+$"))

      if (gene %in% gene_phenotypes$SYMBOL) {
        gp <- dplyr::filter(gene_phenotypes, .data$SYMBOL == gene)
        if (gp$PGX_PHENOTYPE[1] == "Not determined") {
          undetermined <- c(undetermined, gene)
          next
        }
        value <- if (activity_based) gp$PGX_ACTIVITY_SCORE[1] else
          gp$PGX_PHENOTYPE[1]
      } else {
        normal_row <- pgx_phenotype |>
          dplyr::filter(.data$SYMBOL == gene &
                          .data$FUNCTION_1 %in% pgx_normal_function &
                          .data$FUNCTION_2 %in% pgx_normal_function)
        value <- if (activity_based) normal_row$ACTIVITY_SCORE[1] else
          normal_row$PHENOTYPE[1]
        notes <- c(notes, paste0(
          gene, " assumed ", value, " (no altered function alleles ",
          "detected)"))
      }
      key_values <- c(key_values, paste0(gene, ":", value))
    }

    rec <- drug_recs[1, ]
    if (length(undetermined) > 0) {
      rec[, c("LOOKUP_KEY", "PHENOTYPE", "ACTIVITY_SCORE", "IMPLICATIONS",
              "RECOMMENDATION", "CLASSIFICATION", "COMMENTS",
              "POPULATION")] <- NA
      notes <- c(notes, paste0(
        "No recommendation - phenotype not determined for ",
        paste(undetermined, collapse = ", ")))
    } else {
      lookup_key <- canonical_pgx_key(key_values)
      rec <- drug_recs |>
        dplyr::filter(vapply(
          strsplit(.data$LOOKUP_KEY, ";"), canonical_pgx_key,
          character(1)) == lookup_key)
      if (NROW(rec) == 0) {
        rec <- drug_recs[1, ]
        rec[, c("PHENOTYPE", "ACTIVITY_SCORE", "IMPLICATIONS",
                "RECOMMENDATION", "CLASSIFICATION", "COMMENTS",
                "POPULATION")] <- NA
        rec$LOOKUP_KEY <- lookup_key
        notes <- c(notes, "No matching CPIC recommendation")
      }
    }
    rec$PGX_NOTE <- if (length(notes) > 0) paste(notes, collapse = ". ") else
      NA_character_
    recommendations <- dplyr::bind_rows(recommendations, rec)
  }

  if (NROW(recommendations) == 0) {
    return(data.frame())
  }

  recommendations |>
    dplyr::select(dplyr::any_of(c(
      "DRUG_NAME", "MOLECULE_CHEMBL_ID", "GENES", "LOOKUP_KEY",
      "PHENOTYPE", "ACTIVITY_SCORE", "IMPLICATIONS", "RECOMMENDATION",
      "CLASSIFICATION", "COMMENTS", "POPULATION", "CPIC_LEVEL",
      "GUIDELINE_NAME", "GUIDELINE_URL", "PGX_NOTE", "CPIC_RELEASE"))) |>
    dplyr::arrange(.data$GENES, .data$DRUG_NAME)
}

#' Readable phenotype of a CPIC recommendation
#'
#' Combines gene-keyed CPIC phenotype and activity score strings, e.g.
#' 'DPYD:Intermediate Metabolizer' and 'DPYD:1.0' into
#' 'DPYD: Intermediate Metabolizer (activity score 1.0)'; multi-gene
#' phenotypes are separated by '; '.
#'
#' @param phenotype gene-keyed CPIC phenotype(s)
#' @param activity_score gene-keyed CPIC activity score(s), or NA
#'
#' @keywords internal
format_pgx_phenotype <- function(phenotype, activity_score = NA) {
  mapply(function(ph, as) {
    if (is.na(ph)) return(NA_character_)
    ph <- stringr::str_split(ph, ";")[[1]]
    as <- if (is.na(as)) character(0) else stringr::str_split(as, ";")[[1]]
    paste(vapply(ph, function(x) {
      gene <- sub(":.*", "", x)
      score <- sub("^[^:]*:", "", as[startsWith(as, paste0(gene, ":"))])
      paste0(gene, ": ", sub("^[^:]*:", "", x),
             if (length(score) == 1) paste0(" (activity score ", score, ")") else "")
    }, character(1)), collapse = "; ")
  }, phenotype, activity_score, USE.NAMES = FALSE)
}

#' Readable gene-keyed CPIC text, e.g. implications ('GENE:text;GENE:text')
#'
#' @keywords internal
format_pgx_gene_text <- function(x) {
  stringr::str_replace_all(
    stringr::str_replace_all(x, ";(?=[A-Z0-9]+:)", "; "),
    "(^|; )([A-Z0-9]+):", "\\1\\2: ")
}

#' Canonical (sorted) representation of a CPIC lookup key
#'
#' @param key_elements character vector with lookup key elements, e.g.
#' c('1.0:1', '0.5:1') or c('TPMT:Normal Metabolizer', 'NUDT15:...')
#'
#' @keywords internal
canonical_pgx_key <- function(key_elements) {
  paste(sort(trimws(key_elements)), collapse = ";")
}

#' Function that retrieves pharmacogenomic variants based on ClinVar
#' classifications (legacy approach, used for data bundles without CPIC
#' allele definitions)
#'
#' @param calls data frame with all calls found
#'
#' @keywords internal
retrieve_pgx_calls_clinvar <- function(calls) {
  assertable::assert_colnames(
    calls,
    colnames = c(
      "CPG_SOURCE",
      "CLASSIFICATION",
      "ASSERTION_AUTHORITY",
      "PRIMARY_TARGET",
      "GENOTYPE",
      "SYMBOL",
      "GENOMIC_CHANGE",
      "LOSS_OF_FUNCTION",
      "PROTEIN_CHANGE",
      "CLINVAR_PHENOTYPE",
      "CLINVAR_GOLD_STARS",
      "CLINVAR_NUM_SUBMITTERS"
    ),
    only_colnames = F, quiet = T
  )

  pgx_calls <- calls |>
    ## do not consider pharmacogenomics-related variants if
    ## genotypes have not been retrieved properly
    dplyr::filter(
      !is.na(.data$GENOTYPE) &
        !is.na(.data$SYMBOL) &
        !is.na(.data$CPG_SOURCE) &
        stringr::str_detect(.data$CPG_SOURCE, "CPIC_PGX_ONCOLOGY") &
        .data$PRIMARY_TARGET == FALSE &
        !is.na(.data$ASSERTION_AUTHORITY) &
        .data$ASSERTION_AUTHORITY == "ClinVar" &
        !is.na(.data$CLASSIFICATION) &
        stringr::str_detect(
          tolower(
            .data$CLASSIFICATION), "drug|pathogenic"
        )
    ) |>
    dplyr::filter(
      .data$GENOTYPE != "undefined"
    )


  if (NROW(pgx_calls) == 0) {
    return(pgx_calls)
  }else{
    pgx_calls <- pgx_calls |>
      dplyr::arrange(
        dplyr::desc(.data$CLINVAR_GOLD_STARS),
        dplyr::desc(.data$CLINVAR_NUM_SUBMITTERS))

    clinvar_phenotypes <- pgx_calls |>
      dplyr::select(c("GENOMIC_CHANGE",
                      "CLINVAR_PHENOTYPE")) |>
      tidyr::separate_rows(
        .data$CLINVAR_PHENOTYPE, sep = "; ") |>
      dplyr::filter(
        stringr::str_detect(
          .data$CLINVAR_PHENOTYPE,
          "^(Fluorouracil|Dihydropyrimidine|Thiopurine)") |
          stringr::str_detect(
            .data$CLINVAR_PHENOTYPE,
            "^$"
          )
      ) |>
      dplyr::distinct() |>
      dplyr::group_by(
        .data$GENOMIC_CHANGE
      ) |>
      dplyr::summarise(
        CLINVAR_PHENOTYPE = paste(
          sort(unique(.data$CLINVAR_PHENOTYPE)),
          collapse = "; "),
        .groups = "drop"
      )

    pgx_calls <- pgx_calls |>
      dplyr::select(-c("CLINVAR_PHENOTYPE")) |>
      dplyr::left_join(
        clinvar_phenotypes,
        by = "GENOMIC_CHANGE"
      ) |>
      dplyr::filter(!is.na(.data$CLINVAR_PHENOTYPE)) |>
      dplyr::distinct()


  }

  return(pgx_calls)
}
