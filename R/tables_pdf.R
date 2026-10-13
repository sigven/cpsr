## Helpers for tables in the PDF report (tinytable, Typst output)

#' Style the header row of a PDF report table
#'
#' Dark teal header with white, bold text, white separator lines between
#' the header columns, and light zebra striping of the body rows
#'
#' @param tt tinytable object
#'
#' @return styled tinytable object
#' @keywords internal
style_tt_header_pdf <- function(tt) {
  tt <- tinytable::style_tt(
    tt, i = 0, background = "#00504c", color = "white", bold = TRUE)
  n_cols <- ncol(tt@data)
  n_rows <- nrow(tt@data)
  if (n_rows >= 2) {
    tt <- tinytable::style_tt(
      tt, i = seq(2, n_rows, by = 2), background = "#f3f7f7")
  }
  if (n_cols > 1) {
    tt <- tinytable::style_tt(
      tt, i = 0, j = seq_len(n_cols - 1), line = "r",
      line_color = "white", line_width = 0.15)
  }
  tt
}

#' Colored pill (rounded box with white, bold text) in a PDF report table
#' cell - raw Typst markup; the column must not be escaped by format_tt()
#'
#' @param label pill text(s)
#' @param fill background color(s) (hex)
#'
#' @return character vector with Typst markup ("" for empty labels)
#' @keywords internal
pill_pdf <- function(label, fill) {
  esc <- stringr::str_replace_all(label, "([\\\\#*_\\[\\]<>@$`])", "\\\\\\1")
  ifelse(
    is.na(label) | label == "", "",
    paste0('#box(fill: rgb("', fill, '"), radius: 3pt, ',
           'inset: (x: 4pt, y: 2.5pt))[#text(fill: white, weight: "bold", ',
           'size: 0.92em)[', esc, ']]'))
}

#' Variant classification as colored pills in a PDF report table -
#' pathogenicity colors as in the HTML report
#'
#' @param classifications classification per table row (e.g. 'Pathogenic',
#' 'Likely Pathogenic', 'VUS')
#'
#' @return character vector with Typst markup (see pill_pdf); unknown
#' categories in grey
#' @keywords internal
classification_pill_pdf <- function(classifications) {
  pal <- cpsr::color_palette$pathogenicity
  idx <- match(tolower(classifications), tolower(pal$levels))
  pill_pdf(classifications, ifelse(is.na(idx), "#999999", pal$values[idx]))
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
  tbl[["Gene"]] <- gene_link_pdf(tbl[["Gene"]], variants[["ENSEMBL_GENE_ID"]])
  tbl[["Alteration"]] <- link_pdf(
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

#' Escape Typst markup characters in plain text
#'
#' @param x character vector
#' @keywords internal
typst_escape <- function(x) {
  stringr::str_replace_all(x, "([\\\\#*_\\[\\]<>@$`])", "\\\\\\1")
}

#' Text as a link in a PDF report table cell (raw Typst markup; the column
#' must not be escaped by format_tt())
#'
#' @param label link text(s)
#' @param url URL(s); NA for plain (escaped) text
#'
#' @return character vector with Typst markup
#' @keywords internal
link_pdf <- function(label, url) {
  esc <- typst_escape(label)
  ifelse(is.na(url) | is.na(label) | label == "", esc,
         paste0('#link("', url, '")[', esc, ']'))
}

#' Gene symbols as links to the Open Targets Platform target page (by
#' Ensembl gene id) in a PDF report table cell (raw Typst markup; the column
#' must not be escaped by format_tt()). Symbols without a valid Ensembl gene
#' id are shown as plain (escaped) text.
#'
#' @param symbol gene symbol(s)
#' @param ensembl_gene_id Ensembl gene id(s) (e.g. ENSG00000139618)
#'
#' @return character vector with Typst markup
#' @keywords internal
gene_link_pdf <- function(symbol, ensembl_gene_id = NULL) {
  symbol <- dplyr::coalesce(as.character(symbol), "")
  if (is.null(ensembl_gene_id)) {
    ensembl_gene_id <- rep(NA_character_, length(symbol))
  }
  id <- stringr::str_extract(as.character(ensembl_gene_id), "ENSG[0-9]+")
  link_pdf(symbol, ifelse(
    is.na(id), NA_character_,
    paste0("https://platform.opentargets.org/target/", id)))
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

#' Biomarker evidence sources (CIViC, CGI) as logos (Typst markup) - images
#' copied with the templates to the rendering directory (typst/logo_*.png);
#' unknown sources as plain text
#'
#' @param sources source identifiers (e.g. 'civic', 'cgi')
#' @param height logo height (Typst length)
#' @keywords internal
bm_source_logo_pdf <- function(sources, height = "1.4em") {
  ids <- tolower(dplyr::coalesce(as.character(sources), ""))
  ifelse(
    ids %in% c("civic", "cgi"),
    paste0('#box(baseline: 20%, image("typst/logo_', ids, '.png", height: ',
           height, '))'),
    typst_escape(dplyr::coalesce(as.character(sources), "")))
}

#' Legend of the biomarker evidence source logos (raw Typst block)
#'
#' @param sources source identifiers shown in the legend
#' @keywords internal
bm_source_legend_pdf <- function(sources = c("civic", "cgi")) {
  names <- c(civic = "CIViC", cgi = "CGI")
  items <- paste0(bm_source_logo_pdf(sources, height = "1.2em"), " ",
                  names[tolower(sources)])
  paste0("\n```{=typst}\n#block(above: 0.4em, below: 0.8em)[",
         "#set text(size: 8.5pt, fill: cpsr-muted)\n",
         "Evidence sources: ", paste(items, collapse = " #h(0.8em) "),
         "]\n```\n\n")
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
#' @return documentation note with <b>/<i> as markdown emphasis, and line
#' breaks (<br>) removed
#' @keywords internal
doc_note_pdf <- function(x) {
  x |>
    stringr::str_replace_all("(?s)<b[^>]*>(.*?)</b>", "**\\1**") |>
    stringr::str_replace_all("(?s)<i>(.*?)</i>", "*\\1*") |>
    stringr::str_replace_all("<br\\s*/?>", "")
}

#' Message box for a PDF report section without findings
#'
#' @param text message (plain text; Typst markup characters are escaped)
#' @param title short title shown in bold before the message
#'
#' @return raw Typst block (character), to be printed with cat() in a chunk
#' with output 'asis'
#' @keywords internal
no_findings_pdf <- function(text, title = "No findings") {
  esc <- function(x) stringr::str_replace_all(x, "([\\\\#*_\\[\\]<>@$`])", "\\\\\\1")
  paste0(
    "\n```{=typst}\n",
    "#block(width: 100%, inset: (x: 10pt, y: 8pt), radius: 3pt, ",
    "fill: luma(245), stroke: 0.6pt + luma(200), above: 0.8em, below: 1em)[",
    "#text(weight: \"bold\")[", esc(title), "] #h(0.4em) ",
    "#text(fill: luma(70))[", esc(text), "]]\n",
    "```\n\n")
}

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
  pill_pdf(significance, ifelse(is.na(idx), "#999999", pal$values[idx]))
}

#' CPIC phenotypes as colored pills in a PDF report table - blues, darker
#' for more impact (pgx_phenotype_color, as in the HTML report)
#'
#' @param phenotypes CPIC phenotype per table row
#'
#' @return character vector with Typst markup (see pill_pdf)
#' @keywords internal
pgx_phenotype_pill_pdf <- function(phenotypes) {
  out <- rep("", length(phenotypes))
  ok <- !is.na(phenotypes) & phenotypes != ""
  out[ok] <- pill_pdf(phenotypes[ok], pgx_phenotype_color(phenotypes[ok]))
  out
}

#' Genotypes of the detected variants of a gene ('; '-separated:
#' heterozygous, homozygous, hemizygous) as colored pills in a PDF report
#' table, one per line - colors as in the HTML report
#'
#' @param genotypes genotypes per table row
#'
#' @return character vector with Typst markup (see pill_pdf)
#' @keywords internal
pgx_genotype_pill_pdf <- function(genotypes) {
  pal <- cpsr::color_palette$genotypes
  vapply(genotypes, function(x) {
    if (is.na(x) || x == "") return("")
    gts <- strsplit(x, "; ")[[1]]
    level <- ifelse(grepl("homozygous|hemizygous", gts), "hom_alt", "het")
    paste(pill_pdf(gts, pal$bgcolor_values[match(level, pal$levels)]),
          collapse = " #linebreak() ")
  }, character(1), USE.NAMES = FALSE)
}

#' Escape a value for a Typst string literal ("...")
#'
#' @param x character
#' @keywords internal
typst_str <- function(x) {
  x <- dplyr::coalesce(as.character(x), "")
  stringr::str_replace_all(x, '(["\\\\])', "\\\\\\1")
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
  sp <- conf$sample_properties
  subject_sex <- if (is.null(sp$sex)) "UNKNOWN" else as.character(sp$sex)
  subject_age <- if (is.null(sp$age) || identical(sp$age, "NA")) "N/A" else
    paste(sp$age, "years")

  v <- c(
    sample_id = report$settings$sample_id,
    genome_assembly = report$settings$genome_assembly,
    gene_panel = gene_panel_desc,
    sex_age = paste0(subject_sex, " / ", subject_age),
    cpsr_version = report$settings$software$cpsr_version)
  args <- paste0("  ", names(v), ': "', typst_str(v), '",', collapse = "\n")

  ## key findings (markdown, see report_highlights()) in a framed box
  hl <- report_highlights(report, pill = TRUE, vus_coding_only = TRUE)
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
