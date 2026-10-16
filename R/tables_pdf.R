## Helpers for tables in the PDF report (tinytable, Typst output)

#' Style the header row of a PDF report table
#'
#' pcgrr::pdf_style_table() with the CPSR colors: dark teal header with
#' white, bold text, and light zebra striping of the body rows
#'
#' @param tt tinytable object
#'
#' @return styled tinytable object
#' @keywords internal
style_tt_header_pdf <- function(tt) {
  pcgrr::pdf_style_table(
    tt, header_color = "#00504c", stripe_color = "#f3f7f7")
}

#' Variant classification as colored pills in a PDF report table -
#' pathogenicity colors as in the HTML report
#'
#' @param classifications classification per table row (e.g. 'Pathogenic',
#' 'Likely Pathogenic', 'VUS')
#'
#' @return character vector with Typst markup (see pcgrr::pdf_pill());
#' unknown categories in grey
#' @keywords internal
classification_pill_pdf <- function(classifications) {
  pal <- cpsr::color_palette$pathogenicity
  idx <- match(tolower(classifications), tolower(pal$levels))
  pcgrr::pdf_pill(
    classifications, ifelse(is.na(idx), "#999999", pal$values[idx]))
}

#' Variant table (classified variants or secondary findings) for the PDF
#' report, with column names as in the HTML report
#'
#' @param variants data frame with variants (variant_display)
#' @param classification_col column with the variant classification
#' @param source_col column with the assertion authority (ClinVar/CPSR),
#' or NULL
#'
#' @return tinytable object
#' @keywords internal
variant_table_pdf <- function(variants,
                              classification_col = "CLASSIFICATION",
                              source_col = "ASSERTION_AUTHORITY") {
  cols <- c(
    "Gene" = "SYMBOL",
    "Consequence" = "CONSEQUENCE",
    "Alteration" = "ALTERATION",
    "Genotype" = "GENOTYPE",
    "Clinical significance" = classification_col)
  if (!is.null(source_col)) {
    cols <- c(cols, "Assertion authority" = source_col)
  }
  ## ClinVar variation ID (from the ClinVar link of the display table) -
  ## alteration linked to ClinVar when ClinVar is the assertion authority
  clinvar_id <- if ("CLINVAR" %in% colnames(variants)) {
    stringr::str_match(variants$CLINVAR, "clinvar/variation/([0-9]+)")[, 2]
  } else rep(NA_character_, NROW(variants))
  clinvar_authority <- if (!is.null(source_col)) {
    variants[[source_col]] %in% "ClinVar"
  } else rep(TRUE, NROW(variants))
  tbl <- variants |>
    dplyr::select(dplyr::all_of(cols)) |>
    dplyr::mutate(
      Genotype = genotype_label(.data$Genotype),
      dplyr::across(dplyr::everything(),
                    ~ dplyr::coalesce(as.character(.x), "")))
  ## raw Typst markup (not escaped): clinical significance and genotype as
  ## colored pills (as in the HTML report), gene as Open Targets link,
  ## alteration as ClinVar link
  tbl[["Clinical significance"]] <-
    classification_pill_pdf(tbl[["Clinical significance"]])
  tbl[["Genotype"]] <- pgx_genotype_pill_pdf(tbl[["Genotype"]])
  tbl[["Gene"]] <- pcgrr::pdf_gene_link(tbl[["Gene"]], variants[["ENSEMBL_GENE_ID"]])
  tbl[["Alteration"]] <- pcgrr::pdf_link(
    tbl[["Alteration"]],
    ifelse(clinvar_authority & !is.na(clinvar_id),
           paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/",
                  clinvar_id, "/"), NA_character_))
  j_raw <- which(colnames(tbl) %in%
                   c("Gene", "Clinical significance", "Genotype", "Alteration"))
  ## full page width, as the other tables of the PDF report
  widths <- c(0.8, 1.6, 1.3, 1.1, 1.4, 1.2)[seq_len(ncol(tbl))]
  tinytable::tt(tbl, width = widths) |>
    tinytable::format_tt(j = setdiff(seq_len(ncol(tbl)), j_raw),
                         escape = TRUE) |>
    style_tt_header_pdf() |>
    tinytable::style_tt(j = 1, bold = TRUE)
}

#' Ensembl gene ids of gene symbols, looked up among the variants of the
#' report (all variant sets) and the genes of the virtual panel
#'
#' @param report CPSR report object
#' @param symbol gene symbol(s)
#'
#' @return Ensembl gene id per symbol (NA if not found)
#' @keywords internal
gene_ensembl_id_pdf <- function(report, symbol) {
  sets <- c(report$content$snv_indel$callset$variant,
            list(report$settings$conf$gene_panel$panel_genes))
  ids <- dplyr::bind_rows(lapply(sets, function(d) {
    if (NROW(d) == 0 || !all(c("SYMBOL", "ENSEMBL_GENE_ID") %in% colnames(d)))
      return(NULL)
    data.frame(SYMBOL = as.character(d$SYMBOL),
               ENSEMBL_GENE_ID = as.character(d$ENSEMBL_GENE_ID))
  }))
  ids <- ids[grepl("^ENSG[0-9]+", ids$ENSEMBL_GENE_ID), , drop = FALSE]
  ids$ENSEMBL_GENE_ID[match(as.character(symbol), ids$SYMBOL)]
}

#' URL of the biomarker entry in the evidence source: the CIViC variant page
#' (by CIViC variant id), or the CGI biomarker catalogue (no pages for
#' individual CGI biomarkers)
#'
#' @param source_db evidence source ('civic' or 'cgi')
#' @param variant_id variant identifier in the source (CIViC variant id)
#'
#' @return URL per entry (NA if not available)
#' @keywords internal
bm_entry_url_pdf <- function(source_db, variant_id) {
  civic_id <- stringr::str_extract(as.character(variant_id), "^[0-9]+$")
  dplyr::case_when(
    source_db %in% "civic" & !is.na(civic_id) ~
      paste0("https://civicdb.org/variants/", civic_id, "/summary"),
    source_db %in% "cgi" ~
      "https://www.cancergenomeinterpreter.org/2021/biomarkers",
    .default = NA_character_)
}

#' Convert the HTML elements of a documentation note (markdown with inline
#' HTML, as used in the HTML report) to markdown, for the PDF report
#'
#' @param x documentation note (character)
#'
#' @return documentation note with <b>/<i> as markdown emphasis, line
#' breaks (<br>) removed, and literature links/mentions as citations of the
#' report bibliography (see doc_citations_pdf)
#' @keywords internal
doc_note_pdf <- function(x) {
  x <- x |>
    stringr::str_replace_all("(?s)<b[^>]*>(.*?)</b>", "**\\1**") |>
    stringr::str_replace_all("(?s)<i>(.*?)</i>", "*\\1*") |>
    stringr::str_replace_all("<br\\s*/?>", "")
  for (pattern in names(doc_citations_pdf)) {
    x <- stringr::str_replace_all(x, pattern, doc_citations_pdf[[pattern]])
  }
  x
}

#' Literature references in the documentation notes (as hyperlinks or plain
#' mentions in the HTML report) and their citations (Pandoc, keys of
#' cpsr_bibliography.bib) in the PDF report - regular expression (names)
#' and replacement (values)
#'
#' @keywords internal
doc_citations_pdf <- c(
  "(\\[Genomics England PanelApp\\]\\([^)]+\\))" = "\\1 [@Martin2019-nq]",
  "\\[Maxwell et al\\.[^\\]]*\\]\\([^)]+\\)" = "@Maxwell2016-mr",
  "\\[Huang et al\\.[^\\]]*\\]\\([^)]+\\)" = "@Huang2018-ah",
  "\\*\\*ACMG/AMP criteria\\*\\*" =
    "**ACMG/AMP criteria** [@Richards2015-kj; @Amendola2016-pu]")

#' Biomarker clinical significance as colored pills (Typst markup) for the
#' PDF report - colors as in the HTML report (rt_cell_bm_significance)
#'
#' @param significance clinical significance per table row (e.g.
#' 'Sensitivity/Response', 'Poor Outcome')
#'
#' @return character vector with raw Typst markup (rounded box with white,
#' bold text); unknown categories in grey
#' @keywords internal
bm_significance_pill_pdf <- function(significance) {
  pal <- cpsr::color_palette$biomarker_types
  idx <- match(significance, pal$levels)
  pcgrr::pdf_pill(
    significance, ifelse(is.na(idx), "#999999", pal$values[idx]))
}

#' CPIC phenotypes as colored pills in a PDF report table - blues, darker
#' for more impact (pgx_phenotype_color, as in the HTML report)
#'
#' @param phenotypes CPIC phenotype per table row
#'
#' @return character vector with Typst markup (see pcgrr::pdf_pill())
#' @keywords internal
pgx_phenotype_pill_pdf <- function(phenotypes) {
  out <- rep("", length(phenotypes))
  ok <- !is.na(phenotypes) & phenotypes != ""
  out[ok] <- pcgrr::pdf_pill(
    phenotypes[ok], pgx_phenotype_color(phenotypes[ok]))
  out
}

#' Genotypes of the detected variants of a gene ('; '-separated:
#' Heterozygous, Homozygous, Hemizygous) as colored pills in a PDF report
#' table, one per line - colors as in the HTML report
#'
#' @param genotypes genotypes per table row
#'
#' @return character vector with Typst markup (see pcgrr::pdf_pill())
#' @keywords internal
pgx_genotype_pill_pdf <- function(genotypes) {
  pal <- cpsr::color_palette$genotypes
  vapply(genotypes, function(x) {
    if (is.na(x) || x == "") return("")
    gts <- strsplit(x, "; ")[[1]]
    level <- ifelse(grepl("homozygous|hemizygous", gts, ignore.case = TRUE), "hom_alt", "het")
    paste(pcgrr::pdf_pill(gts, pal$bgcolor_values[match(level, pal$levels)]),
          collapse = " #linebreak() ")
  }, character(1), USE.NAMES = FALSE)
}

#' Cover page of the PDF report: sample metadata (cpsr-cover(), defined in
#' typst/_preamble.typ) and key findings (report_highlights()), followed by
#' the table of contents
#'
#' @param report CPSR report object
#'
#' @return raw Typst blocks around markdown (character), to be printed with
#' cat() in a chunk with output 'asis'
#' @keywords internal
pdf_cover_typst <- function(report) {
  conf <- report$settings$conf

  ## virtual panel description
  gene_panel_desc <- conf$gene_panel$description
  if (!is.null(conf$gene_panel$description_trait) &&
      !conf$gene_panel$panel_id %in% c("-1", "0")) {
    gene_panel_desc <- paste0(
      gene_panel_desc, ": ", conf$gene_panel$description_trait)
  }
  ## sex and age of the subject (age at testing: see Documentation)
  sp <- conf$subject_properties
  subject_sex <- if (is.null(sp$sex)) "UNKNOWN" else as.character(sp$sex)
  subject_age <- if (is.null(sp$age) || identical(sp$age, "NA")) "N/A" else
    paste(sp$age, "years")

  v <- c(
    sample_id = report$settings$sample_id,
    genome_assembly = report$settings$genome_assembly,
    gene_panel = gene_panel_desc,
    sex_age = paste0(subject_sex, " / ", subject_age),
    cpsr_version = report$settings$software$cpsr_version)
  args <- paste0("  ", names(v), ': "', pcgrr::typst_str(v), '",', collapse = "\n")

  ## key findings (markdown, see report_highlights()) in a framed box
  hl <- report_highlights(report, pill = TRUE, vus_coding_only = TRUE,
                          links = TRUE)
  paste0(
    "```{=typst}\n#cpsr-cover(\n", args, "\n)\n",
    "#cpsr-highlights[\n#cpsr-hl-lead[\n```\n\n",
    hl$lead, "\n\n",
    "```{=typst}\n]\n",
    if (nzchar(hl$supporting)) paste0(
      "#cpsr-hl-support[\n```\n\n", hl$supporting, "\n\n",
      "```{=typst}\n]\n"),
    "]\n#cpsr-cover-end()\n```\n")
}
