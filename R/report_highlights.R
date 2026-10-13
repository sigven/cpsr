## Key findings ("highlights") of a CPSR report, as a short templated text -
## the most critical findings (pathogenic/likely pathogenic variants in
## cancer predisposition genes) first, supported by other findings (variant
## classification overview, secondary findings, genomic biomarkers,
## pharmacogenomics).
## Markdown output, for use in both the PDF and HTML reports (same framework
## as pcgrr::report_highlights()).

#' Escape markdown emphasis characters in plain text (e.g. 'DPYD *2A')
#'
#' @param x character vector
#' @keywords internal
md_escape_hl <- function(x) {
  stringr::str_replace_all(x, "([\\\\*_`\\[\\]<>])", "\\\\\\1")
}

#' Remove HTML tags from (display) text
#'
#' @param x character vector
#' @keywords internal
strip_html_hl <- function(x) {
  stringr::str_trim(stringr::str_replace_all(x, "<[^>]+>", ""))
}

#' Compact list of items ('A, B, C (+k more)')
#'
#' @param x items (character)
#' @param max_show maximum number of items shown
#' @keywords internal
fmt_list_hl <- function(x, max_show = 5) {
  x <- unique(x[!is.na(x) & x != ""])
  n <- length(x)
  if (n == 0) return("")
  s <- paste(utils::head(x, max_show), collapse = ", ")
  if (n > max_show) s <- paste0(s, " (+", n - max_show, " more)")
  s
}

#' Variant labels for the key findings: gene symbol (bold) and protein
#' change (or coding DNA change for non-coding variants), with zygosity for
#' homozygous variants - ordered by gene symbol, then alteration
#'
#' @param d variant data frame (display columns SYMBOL, ALTERATION, and
#' optionally GENOTYPE)
#' @keywords internal
variant_labels_hl <- function(d) {
  if (NROW(d) == 0) return(character())
  alt <- strip_html_hl(dplyr::coalesce(as.character(d$ALTERATION), ""))
  d <- d[order(as.character(d$SYMBOL), alt), , drop = FALSE]
  alt <- strip_html_hl(dplyr::coalesce(as.character(d$ALTERATION), ""))
  hom <- if ("GENOTYPE" %in% colnames(d))
    d$GENOTYPE %in% c("hom_alt", "homozygous") else rep(FALSE, NROW(d))
  paste0("**", md_escape_hl(d$SYMBOL), "**",
         ifelse(nzchar(alt), paste0(" ", md_escape_hl(alt)), ""),
         ifelse(hom, " (homozygous)", ""))
}

#' Key findings of a CPSR report (highlights) as markdown
#'
#' The most critical findings come first: pathogenic and likely pathogenic
#' variants in the cancer predisposition genes of the virtual panel. These are
#' followed by supporting findings, each included only when assessed: an
#' overview of variant classifications (including VUS), secondary findings
#' (ACMG SF), genomic biomarkers (evidence levels A/B), and pharmacogenomic
#' phenotypes (CPIC).
#'
#' @param report CPSR report object
#' @param max_show maximum number of items listed per finding
#' @param pill label 'Pathogenic' / 'Likely Pathogenic' as colored pills
#' (white bold text on pathogenicity color, as in the Clinical significance
#' column): TRUE gives raw Typst (PDF report), "html" gives an HTML span
#' (HTML report), FALSE (default) gives bold text
#' @param vus_coding_only count protein-coding VUS only (PDF report - more
#' condensed; non-coding VUS are noted separately)
#'
#' @return list with 'lead' (markdown on pathogenic/likely pathogenic
#' variants) and 'supporting' (markdown bullet list, "" if none)
#' @export
report_highlights <- function(report, max_show = 5, pill = FALSE,
                              vus_coding_only = FALSE) {
  snv <- report$content$snv_indel
  conf <- report$settings$conf
  vc <- conf$variant_classification
  has_content <- isTRUE(as.logical(snv$eval))
  cpg <- snv$callset$variant_display$cpg_non_sf
  cls <- if (NROW(cpg) > 0 && "CLASSIFICATION" %in% colnames(cpg))
    as.character(cpg$CLASSIFICATION) else character()

  ## --- Lead: pathogenic/likely pathogenic variants ---
  lead <- c()
  plp <- cpg[cls %in% c("Pathogenic", "Likely Pathogenic"), , drop = FALSE]
  if (NROW(plp) > 0) {
    ## one list item per classification (pathogenic first)
    class_label <- function(cl) {
      if (isFALSE(pill)) return(paste0("**", cl, "**:"))
      pal <- cpsr::color_palette$pathogenicity
      idx <- match(tolower(cl), tolower(pal$levels))
      fill <- if (is.na(idx)) "#999999" else pal$values[idx]
      if (identical(pill, "html")) {
        return(paste0("<span class=\"kf-pill\" style=\"background:", fill,
                      "\">", cl, "</span>"))
      }
      paste0("`#cpsr-hl-pill(\"", cl, "\", \"", fill, "\")`{=typst}")
    }
    class_item <- function(cl) {
      lbl <- variant_labels_hl(plp[plp$CLASSIFICATION == cl, , drop = FALSE])
      if (length(lbl) == 0) return(NULL)
      paste0("- ", class_label(cl), " ", fmt_list_hl(lbl, max_show))
    }
    ## classification source (ClinVar or CPSR's ACMG/AMP implementation)
    auth <- table(factor(plp$ASSERTION_AUTHORITY, levels = c("ClinVar", "CPSR")))
    auth_lbl <- c(ClinVar = "ClinVar-asserted", CPSR = "CPSR-asserted (ACMG/AMP)")
    auth_txt <- paste0(auth_lbl[names(auth)[auth > 0]], ": **n = ",
                       auth[auth > 0], "**", collapse = ", ")
    lead <- c(lead, paste0(
      "**Pathogenic / likely pathogenic variants** in cancer predisposition ",
      "genes: **n = ", NROW(plp), "**",
      if (nzchar(auth_txt)) paste0(" - ", auth_txt), "\n\n",
      paste(c(class_item("Pathogenic"), class_item("Likely Pathogenic")),
            collapse = "\n")))
  } else if (has_content) {
    lead <- c(lead, paste0(
      "**No pathogenic or likely pathogenic variants** were identified in ",
      "the cancer predisposition genes of the virtual panel."))
  } else {
    lead <- c(lead, paste0(
      "**No variants** were found in the cancer predisposition genes of ",
      "the virtual panel."))
  }

  ## --- Supporting findings ---
  supporting <- c()

  ## overview of classifications in the virtual panel
  if (has_content && length(cls) > 0) {
    n_panel <- NROW(dplyr::distinct(dplyr::filter(
      conf$gene_panel$panel_genes, .data$PRIMARY_TARGET == TRUE),
      .data$ENTREZGENE))
    ## VUS - optionally protein-coding only (non-coding VUS noted separately)
    is_vus <- cls == "VUS"
    vus_coding <- if ("CODING_STATUS" %in% colnames(cpg))
      cpg$CODING_STATUS %in% "coding" else rep(TRUE, length(cls))
    n_vus_noncoding <- if (vus_coding_only) sum(is_vus & !vus_coding) else 0
    vus_label <- paste0(
      if (vus_coding_only) "protein-coding " else "",
      "variants of uncertain significance (VUS",
      if (n_vus_noncoding > 0)
        paste0("; +", n_vus_noncoding, " non-coding, not listed"), ")")
    n_cls <- c(
      "pathogenic/likely pathogenic" =
        sum(cls %in% c("Pathogenic", "Likely Pathogenic")),
      stats::setNames(sum(is_vus) - n_vus_noncoding, vus_label),
      "benign/likely benign" = sum(cls %in% c("Benign", "Likely Benign")))
    supporting <- c(supporting, paste0(
      "**Variant classification:** ", length(cls), " variant(s) in ",
      length(unique(cpg$SYMBOL)), " of **n = ", n_panel,
      "** evaluated cancer predisposition genes, ",
      paste0("**n = ", n_cls, "** ", names(n_cls), collapse = ", ")))
  }

  ## secondary findings (ACMG SF genes)
  if (isTRUE(as.logical(vc$secondary_findings))) {
    sf <- snv$callset$variant_display$sf
    supporting <- c(supporting, paste0(
      "**Secondary findings (ACMG SF):** ",
      if (NROW(sf) > 0 && "SYMBOL" %in% colnames(sf))
        fmt_list_hl(variant_labels_hl(sf), max_show) else
          "none identified"))
  }

  ## genomic biomarkers - variants with evidence of level A/B (validated or
  ## clinical evidence), with the types of evidence
  bm <- snv$callset$variant_display$biomarker
  if (has_content && NROW(bm) > 0 &&
      all(c("SYMBOL", "BM_EVIDENCE_LEVEL", "BM_EVIDENCE_TYPE") %in%
          colnames(bm))) {
    bm_ab <- bm[substr(bm$BM_EVIDENCE_LEVEL, 1, 1) %in% c("A", "B"), ,
                drop = FALSE]
    if (NROW(bm_ab) > 0) {
      ev_order <- c("Predictive", "Prognostic", "Diagnostic", "Predisposing")
      bm_ab$VARIANT <- strip_html_hl(stringr::str_replace(
        bm_ab$SAMPLE_ALTERATION, paste0("^", bm_ab$SYMBOL, " - "), ""))
      bm_var <- bm_ab |>
        dplyr::group_by(.data$SYMBOL, .data$VARIANT) |>
        dplyr::summarise(
          TYPES = paste(tolower(intersect(ev_order, .data$BM_EVIDENCE_TYPE)),
                        collapse = ", "),
          .groups = "drop") |>
        dplyr::arrange(.data$SYMBOL, .data$VARIANT)
      supporting <- c(supporting, paste0(
        "**Genomic biomarkers** (evidence levels A/B): ",
        fmt_list_hl(paste0(
          "**", md_escape_hl(bm_var$SYMBOL), "** ",
          md_escape_hl(bm_var$VARIANT),
          ifelse(nzchar(bm_var$TYPES), paste0(" (", bm_var$TYPES, ")"), "")),
          max_show)))
    } else {
      supporting <- c(supporting, paste0(
        "**Genomic biomarkers:** ", length(unique(bm$VAR_ID)),
        " variant(s) with lower-level evidence only (levels C-E)"))
    }
  }

  ## pharmacogenomics - CPIC phenotypes (or ClinVar drug response variants)
  if (isTRUE(as.logical(vc$pgx_findings))) {
    cpic <- identical(snv$pgx_method, "cpic")
    if (isFALSE(as.logical(conf$sample_properties$gt_detected))) {
      pgx_txt <- "not assessed (genotype information missing in input VCF)"
    } else if (cpic && NROW(snv$pgx_phenotype) > 0) {
      ph <- snv$pgx_phenotype
      ph <- ph[order(as.character(ph$SYMBOL)), , drop = FALSE]
      rec <- snv$pgx_recommendation
      ## drugs with CPIC recommendations for the gene
      drugs <- vapply(as.character(ph$SYMBOL), function(g) {
        if (NROW(rec) == 0 || !all(c("GENES", "DRUG_NAME") %in% colnames(rec)))
          return("")
        d <- rec$DRUG_NAME[!is.na(rec$GENES) & stringr::str_detect(
          rec$GENES, paste0("\\b", g, "\\b"))]
        fmt_list_hl(sort(unique(d)), max_show)
      }, character(1))
      pgx_txt <- paste0(
        "**", md_escape_hl(ph$SYMBOL), "** ", md_escape_hl(ph$PGX_PHENOTYPE),
        ifelse(is.na(ph$PGX_DIPLOTYPE), "",
               paste0(" (", md_escape_hl(ph$PGX_DIPLOTYPE), ")")),
        ifelse(nzchar(drugs), paste0(" - relevant for ", drugs), ""),
        collapse = "; ")
    } else if (cpic) {
      pgx_txt <- paste0(
        "no CPIC alleles of altered function detected (DPYD, TPMT, NUDT15, ",
        "G6PD)")
    } else {
      pv <- snv$callset$variant_display$pgx
      pgx_txt <- if (NROW(pv) > 0 && "SYMBOL" %in% colnames(pv))
        paste0("variants related to drug toxicity in ", fmt_list_hl(
          paste0("**", md_escape_hl(sort(unique(pv$SYMBOL))), "**"),
          max_show)) else
            "no variants related to toxicity of chemotherapy"
    }
    supporting <- c(supporting, paste0(
      "**Pharmacogenomics", if (cpic) " (CPIC)", ":** ", pgx_txt))
  }

  list(
    lead = paste(lead, collapse = "\n\n"),
    supporting = if (length(supporting) > 0)
      paste0("- ", supporting, collapse = "\n") else "")
}
