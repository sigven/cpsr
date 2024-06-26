#' Function that makes a piechart showing the number of variants at
#' each significance level
#'
#' @param variants_tsv data frame with variants in predisposition_genes
#' @param plot_type ClinVar or Other
#'
#' @export
plot_summary_statistics <- function(variants_tsv, plot_type = "ClinVar") {
  title <- "ClinVar variants"
  p <- NULL

  if (nrow(variants_tsv) > 0) {
    set_clinvar <- variants_tsv |>
      dplyr::filter(!is.na(.data$CLINVAR_CLASSIFICATION) &
        !(.data$CLINVAR_CLASSIFICATION == "NA"))
    set_other <- variants_tsv |>
      dplyr::filter(nchar(.data$CPSR_CLASSIFICATION) > 0 &
        (is.na(.data$CLINVAR_CLASSIFICATION) |
          .data$CLINVAR_CLASSIFICATION == "NA"))

    if ((plot_type == "ClinVar" & nrow(set_clinvar) > 0) |
      (plot_type != "ClinVar" & nrow(set_other) > 0)) {
      m <- data.frame()

      if (plot_type == "ClinVar") {
        if (nrow(set_clinvar) > 0) {
          t <- paste0("n = ", nrow(set_clinvar))
          title <- bquote("ClinVar variants, " ~ bold(.(t)))
          m <- as.data.frame(set_clinvar |>
            dplyr::group_by(.data$CLINVAR_CLASSIFICATION) |>
            dplyr::summarise(n = dplyr::n()) |>
            dplyr::rename(level = .data$CLINVAR_CLASSIFICATION)) |>
            dplyr::mutate(
              level =
                factor(
                  .data$level,
                  levels =
                    pcgrr::color_palette[["pathogenicity"]][["levels"]]
                )
            ) |>
            dplyr::arrange(.data$level) |>
            dplyr::mutate(prop = as.numeric(.data$n / sum(.data$n))) |>
            dplyr::mutate(lab.ypos = cumsum(.data$prop) - 0.5 * .data$prop) |>
            dplyr::mutate(n = as.character(.data$n))
        }
      } else {
        if (nrow(set_other) > 0) {
          t <- paste0("n = ", nrow(set_other))
          title <- bquote("Other variants, CPSR-classified, " ~ bold(.(t)))
          m <- as.data.frame(set_other |>
            dplyr::group_by(.data$CPSR_CLASSIFICATION) |>
            dplyr::summarise(n = dplyr::n()) |>
            dplyr::rename(level = .data$CPSR_CLASSIFICATION)) |>
            dplyr::mutate(
              level = factor(
                .data$level,
                levels = pcgrr::color_palette[["pathogenicity"]][["levels"]]
              )
            ) |>
            dplyr::arrange(.data$level) |>
            dplyr::mutate(prop = as.numeric(.data$n / sum(.data$n))) |>
            dplyr::mutate(lab.ypos = cumsum(.data$prop) - 0.5 * .data$prop) |>
            dplyr::mutate(n = as.character(.data$n))
        }
      }


      p <- ggplot2::ggplot(m, ggplot2::aes(x = 2, y = .data$prop, fill = .data$level)) +
        ggplot2::geom_bar(stat = "identity", color = "white") +
        ggplot2::coord_polar(theta = "y", start = 0) +
        ggplot2::geom_text(ggplot2::aes(y = 1 - .data$lab.ypos, label = .data$n),
          color = "white", family = "Helvetica", size = 6
        ) +
        ggplot2::scale_fill_manual(
          values = pcgrr::color_palette[["pathogenicity"]][["values"]],
          labels = pcgrr::color_palette[["pathogenicity"]][["levels"]],
          drop = F
        ) +
        ggplot2::theme_void() +
        ggplot2::xlim(0.5, 2.5) +
        ggplot2::ggtitle(title) +
        ggplot2::theme(
          plot.title =
            ggplot2::element_text(
              family = "Helvetica",
              size = 16, vjust = -1,
              hjust = 0.5
            ),
          legend.title = ggplot2::element_blank(),
          plot.margin = ggplot2::margin(0, 0, 0, 0, "cm"),
          legend.text = ggplot2::element_text(
            family = "Helvetica", size = 14
          )
        )
    }
  }
  return(p)
}

#' Function that makes a HTML display of virtual gene panel
#'
#' @param gene_df data frame with genes targeted in virtual panel
#'
#'
#' @export
plot_virtual_panels <- function(gene_df) {

  i <- 1
  gene_df <- gene_df |>
    dplyr::arrange(
      dplyr::desc(.data$CONFIDENCE_LEVEL), .data$SYMBOL) |>
    dplyr::filter(.data$PRIMARY_TARGET == TRUE)

  if(length(unique(gene_df$CONFIDENCE_LEVEL)) == 1){
    if(unique(gene_df$CONFIDENCE_LEVEL) == 5){
      gene_df <- gene_df |>
        dplyr::select(
          c("SYMBOL", "ENTREZGENE", "CONFIDENCE_LEVEL")) |>
        dplyr::distinct() |>
        dplyr::mutate(PANEL_ID = as.character(NA))
    }
  }

  html_string <- "<div id=\"container\">"
  while(i <= nrow(gene_df)) {
    CONFIDENCE_LEVEL <- gene_df[i,"CONFIDENCE_LEVEL"]
    css_class <- "exploratory"
    #if (CONFIDENCE_LEVEL == 3) {
    #  css_class <- "green"
    #}
    #if (CONFIDENCE_LEVEL == 2) {
    #  css_class <- "amber"
    #}
    #if (CONFIDENCE_LEVEL == 1) {
    #  css_class <- "red"
    #}
    if (CONFIDENCE_LEVEL == -1) {
      css_class <- "custom"
    }
    #if (CONFIDENCE_LEVEL == 0) {
    #  css_class <- "nolist"
    #}
    if (CONFIDENCE_LEVEL == 5) {
      css_class <- "app_combo"
    }
    symbol <- gene_df[i, "SYMBOL"]
    entrezgene <- gene_df[i, "ENTREZGENE"]
    panel_id <- gene_df[i, "PANEL_ID"]

    gene_url <- paste0("https://www.ncbi.nlm.nih.gov/gene/", entrezgene)
    if (!is.na(panel_id)) {
      gene_url <- paste0("https://panelapp.genomicsengland.co.uk/panels/",
                         panel_id, "/", symbol)
    }

    entry_string <- paste0("  <div class=\"", css_class, "\"><a href=\"",
                           gene_url, "\" target=\"_blank\" title=\"",
                           symbol, "\">", symbol, "</a></div>")
    html_string <- paste0(html_string, entry_string)
    if (i %% 9 == 0) {
      html_string <- paste0(html_string, "</div>  <div id=\"container\">")
    }
    i <- i + 1
  }
  html_string <- paste0(html_string, "</div>")
  return(html_string)
}

