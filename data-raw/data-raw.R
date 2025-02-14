
col_format_output <- list()

col_format_output[['html_tier']] <-
  c(
    "SYMBOL",
    "CLINVAR_PHENOTYPE",
    "CONSEQUENCE",
    "ALTERATION",
    "GENOTYPE",
    "GENENAME",
    "PROTEIN_DOMAIN",
    "DP_CONTROL",
    "PROTEIN_CHANGE",
    "HGVSp",
    "HGVSc",
    "HGVSc_RefSeq",
    "ENSEMBL_GENE_ID",
    "REFSEQ_TRANSCRIPT_ID",
    "ENSEMBL_TRANSCRIPT_ID",
    "CDS_CHANGE",
    "MUTATION_HOTSPOT",
    "RMSK_HIT",
    "PREDICTED_EFFECT",
    "miRNA_TARGET_HIT",
    "miRNA_TARGET_HIT_PREDICTION",
    "TF_BINDING_SITE_VARIANT",
    "TF_BINDING_SITE_VARIANT_INFO",
    "GERP_SCORE",
    "LOSS_OF_FUNCTION",
    "LOF_FILTER",
    "DBSNP_RSID",
    "CLINVAR",
    "CLINVAR_CLASSIFICATION",
    "CLINVAR_REVIEW_STATUS_STARS",
    "CLINVAR_CONFLICTED",
    "CLINVAR_VARIANT_ORIGIN",
    "CPSR_CLASSIFICATION_SOURCE",
    "CPSR_CLASSIFICATION",
    "CPSR_PATHOGENICITY_SCORE",
    "CPSR_CLASSIFICATION_DOC",
    "CPSR_CLASSIFICATION_CODE",
    "FINAL_CLASSIFICATION",
    "ONCOGENE",
    "TUMOR_SUPPRESSOR",
    "gnomADe_AF",
    "GENOMIC_CHANGE",
    "GENOME_VERSION"
  )

## define tags/variables to display in data tables (secondary findings)
col_format_output[['html_sf']] <-
  c(
    "SYMBOL",
    "CONSEQUENCE",
    "FINAL_CLASSIFICATION",
    "CLINVAR_PHENOTYPE",
    "ALTERATION",
    "GENOTYPE",
    "GENENAME",
    "PROTEIN_DOMAIN",
    "DP_CONTROL",
    "PROTEIN_CHANGE",
    "HGVSp",
    "HGVSc",
    "HGVSc_RefSeq",
    "ENSEMBL_GENE_ID",
    "REFSEQ_TRANSCRIPT_ID",
    "ENSEMBL_TRANSCRIPT_ID",
    "CDS_CHANGE",
    "PREDICTED_EFFECT",
    "LOSS_OF_FUNCTION",
    "LOF_FILTER",
    "DBSNP_RSID",
    "CLINVAR_CLASSIFICATION",
    "CLINVAR",
    "CLINVAR_REVIEW_STATUS_STARS",
    "CLINVAR_CONFLICTED",
    "CPSR_CLASSIFICATION_SOURCE",
    "CPSR_PATHOGENICITY_CODE",
    "CPSR_PATHOGENICITY_SCORE",
    "gnomADe_AF",
    "GENOMIC_CHANGE",
    "GENOME_VERSION"
  )

## define tags/variables to display in data tables (PGx findings)
col_format_output[['html_pgx']] <-
  c(
    "SYMBOL",
    "ALTERATION",
    "CLINVAR_CLASSIFICATION",
    "CLINVAR_PHENOTYPE",
    "GENOTYPE",
    "GENENAME",
    "CONSEQUENCE",
    "PROTEIN_DOMAIN",
    "DP_CONTROL",
    "PROTEIN_CHANGE",
    "HGVSp",
    "HGVSc",
    "HGVSc_RefSeq",
    "ENSEMBL_GENE_ID",
    "REFSEQ_TRANSCRIPT_ID",
    "ENSEMBL_TRANSCRIPT_ID",
    "CDS_CHANGE",
    "PREDICTED_EFFECT",
    "LOSS_OF_FUNCTION",
    "LOF_FILTER",
    "DBSNP_RSID",
    "CLINVAR",
    "CLINVAR_REVIEW_STATUS_STARS",
    "CLINVAR_CONFLICTED",
    "gnomADe_AF",
    "GENOMIC_CHANGE",
    "GENOME_VERSION"
  )


## define tags/variables to display in data tables (GWAS findings)
col_format_output[['html_gwas']] <-
  c(
    "SYMBOL",
    "CONSEQUENCE",
    "GWAS_CITATION",
    "ALTERATION",
    "GENOTYPE",
    "GENENAME",
    "LOSS_OF_FUNCTION",
    "PROTEIN_DOMAIN",
    "DP_CONTROL",
    "PROTEIN_CHANGE",
    "GWAS_PHENOTYPE",
    "HGVSp",
    "HGVSc",
    "HGVSc_RefSeq",
    "ENSEMBL_GENE_ID",
    "REFSEQ_TRANSCRIPT_ID",
    "ENSEMBL_TRANSCRIPT_ID",
    "CDS_CHANGE",
    "CODING_STATUS",
    "miRNA_TARGET_HIT",
    "miRNA_TARGET_HIT_PREDICTION",
    "TF_BINDING_SITE_VARIANT",
    "TF_BINDING_SITE_VARIANT_INFO",
    "GERP_SCORE",
    "PREDICTED_EFFECT",
    "DBSNP_RSID",
    "gnomADe_AF",
    "GENOMIC_CHANGE",
    "GENOME_VERSION"
  )

## define tags/variables to display in output TSV
col_format_output[['tsv']] <-
  c("SAMPLE_ID",
    "GENOMIC_CHANGE",
    "VAR_ID",
    "GENOME_VERSION",
    "GENOTYPE",
    "DP_CONTROL",
    "CPSR_CLASSIFICATION_SOURCE",
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
    "ALTERATION",
    "PROTEIN_CHANGE",
    "PFAM_DOMAIN",
    "PFAM_DOMAIN_NAME",
    "HGVSp",
    "HGVSc",
    "HGVSc_RefSeq",
    "CDS_CHANGE",
    "LAST_EXON",
    "EXON",
    "EXON_AFFECTED",
    "EXON_POSITION",
    "INTRON_POSITION",
    "VEP_ALL_CSQ",
    "CANCER_PHENOTYPE",
    "MUTATION_HOTSPOT",
    "RMSK_HIT",
    "EFFECT_PREDICTIONS",
    "LOSS_OF_FUNCTION",
    "LOF_FILTER",
    "NULL_VARIANT",
    "DBMTS",
    "REGULATORY_ANNOTATION",
    "TF_BINDING_SITE_VARIANT",
    "TF_BINDING_SITE_VARIANT_INFO",
    "GERP_SCORE",
    "DBSNP_RSID",
    "CLINVAR_CLASSIFICATION",
    "CLINVAR_MSID",
    "CLINVAR_VARIANT_ORIGIN",
    "CLINVAR_CONFLICTED",
    "CLINVAR_PHENOTYPE",
    "CLINVAR_REVIEW_STATUS_STARS",
    "N_INSILICO_CALLED",
    "N_INSILICO_DAMAGING",
    "N_INSILICO_TOLERATED",
    "N_INSILICO_SPLICING_NEUTRAL",
    "N_INSILICO_SPLICING_AFFECTED",
    "gnomADe_AF",
    "FINAL_CLASSIFICATION",
    "CPSR_CLASSIFICATION",
    "CPSR_PATHOGENICITY_SCORE",
    "CPSR_CLASSIFICATION_CODE",
    "CPSR_CLASSIFICATION_DOC",
    "CPSR_CLASSIFICATION_SOURCE"
  )

col_format_output[['html_bm']] <-
  c('SYMBOL',
    'GENENAME',
    'ALTERATION',
    'CONSEQUENCE',
    'BM_EVIDENCE_LEVEL',
    'BM_MOLECULAR_PROFILE',
    'BM_REFERENCE',
    'GENOTYPE',
    'PROTEIN_CHANGE',
    'DP_CONTROL',
    'BM_CANCER_TYPE',
    'BM_DISEASE_ONTOLOGY_ID',
    'BM_PRIMARY_SITE',
    'BM_CLINICAL_SIGNIFICANCE',
    'BM_THERAPEUTIC_CONTEXT',
    'BM_RATING',
    'BM_EVIDENCE_TYPE',
    'BM_EVIDENCE_DIRECTION',
    'BM_EVIDENCE_DESCRIPTION',
    'BM_SOURCE_DB',
    'BM_EVIDENCE_ID',
    'BM_VARIANT_ORIGIN',
    'BM_MATCH',
    'BM_RESOLUTION',
    'PROTEIN_DOMAIN',
    'CDS_CHANGE',
    'HGVSc',
    "HGVSc_RefSeq",
    'HGVSp',
    'FINAL_CLASSIFICATION',
    'PREDICTED_EFFECT',
    'DBSNP_RSID',
    'ENSEMBL_GENE_ID',
    'ENSEMBL_TRANSCRIPT_ID',
    'REFSEQ_TRANSCRIPT_ID',
    'GENOMIC_CHANGE',
    'GENOME_VERSION')

## define tags/variables to display in output Excel
col_format_output[['xlsx_classification']] <-
  c("SAMPLE_ID",
    "GENOMIC_CHANGE",
    "GENOTYPE",
    "DP_CONTROL",
    "GENOME_VERSION",
    "VARIANT_CLASS",
    "SYMBOL",
    "GENENAME",
    "CONSEQUENCE",
    "PROTEIN_CHANGE",
    "FINAL_CLASSIFICATION",
    "CPSR_CLASSIFICATION_SOURCE",
    "CLINVAR_CLASSIFICATION",
    "CPSR_CLASSIFICATION",
    "CPSR_CLASSIFICATION_CODE",
    "CPSR_PATHOGENICITY_SCORE",
    "ENTREZGENE",
    "ENSEMBL_GENE_ID",
    "ENSEMBL_TRANSCRIPT_ID",
    "REFSEQ_TRANSCRIPT_ID",
    "ONCOGENE",
    "TUMOR_SUPPRESSOR",
    "PFAM_DOMAIN_NAME",
    "HGVSp",
    "HGVSc",
    "HGVSc_RefSeq",
    "CDS_CHANGE",
    "CODING_STATUS",
    "MUTATION_HOTSPOT",
    "EFFECT_PREDICTIONS",
    "LOSS_OF_FUNCTION",
    "LOF_FILTER",
    "NULL_VARIANT",
    "DBSNP_RSID",
    "CLINVAR_MSID",
    "CLINVAR_VARIANT_ORIGIN",
    "CLINVAR_CONFLICTED",
    "CLINVAR_PHENOTYPE",
    "CLINVAR_REVIEW_STATUS_STARS",
    "gnomADe_AF"
  )

## define tags/variables to display in output Excel
col_format_output[['xlsx_secondary']] <-
  c("SAMPLE_ID",
    "GENOMIC_CHANGE",
    "GENOTYPE",
    "DP_CONTROL",
    "GENOME_VERSION",
    "VARIANT_CLASS",
    "SYMBOL",
    "GENENAME",
    "CONSEQUENCE",
    "PROTEIN_CHANGE",
    "CPSR_CLASSIFICATION_SOURCE",
    "FINAL_CLASSIFICATION",
    "CPSR_CLASSIFICATION",
    "CPSR_PATHOGENICITY_SCORE",
    "CPSR_CLASSIFICATION_CODE",
    "CLINVAR_CLASSIFICATION",
    "CLINVAR_MSID",
    "CLINVAR_VARIANT_ORIGIN",
    "CLINVAR_CONFLICTED",
    "CLINVAR_PHENOTYPE",
    "CLINVAR_REVIEW_STATUS_STARS",
    "ENSEMBL_GENE_ID",
    "ENSEMBL_TRANSCRIPT_ID",
    "REFSEQ_TRANSCRIPT_ID",
    "PFAM_DOMAIN_NAME",
    "HGVSp",
    "HGVSc",
    "HGVSc_RefSeq",
    "CDS_CHANGE",
    "CODING_STATUS",
    "MUTATION_HOTSPOT",
    "EFFECT_PREDICTIONS",
    "LOSS_OF_FUNCTION",
    "LOF_FILTER",
    "NULL_VARIANT",
    "DBSNP_RSID",
    "gnomADe_AF"
  )

## define tags/variables to display in output Excel
col_format_output[['xlsx_pgx']] <-
  c("SAMPLE_ID",
    "GENOMIC_CHANGE",
    "GENOTYPE",
    "DP_CONTROL",
    "GENOME_VERSION",
    "VARIANT_CLASS",
    "SYMBOL",
    "GENENAME",
    "CONSEQUENCE",
    "PROTEIN_CHANGE",
    "CLINVAR_CLASSIFICATION",
    "CLINVAR_MSID",
    "CLINVAR_VARIANT_ORIGIN",
    "CLINVAR_CONFLICTED",
    "CLINVAR_PHENOTYPE",
    "CLINVAR_REVIEW_STATUS_STARS",
    "ENSEMBL_GENE_ID",
    "ENSEMBL_TRANSCRIPT_ID",
    "REFSEQ_TRANSCRIPT_ID",
    "PFAM_DOMAIN_NAME",
    "HGVSp",
    "HGVSc",
    "HGVSc_RefSeq",
    "CDS_CHANGE",
    "CODING_STATUS",
    "MUTATION_HOTSPOT",
    "EFFECT_PREDICTIONS",
    "LOSS_OF_FUNCTION",
    "LOF_FILTER",
    "NULL_VARIANT",
    "DBSNP_RSID",
    "gnomADe_AF"
  )


## define tags/variables to display in output Excel
col_format_output[['xlsx_biomarker']] <-
  c("SAMPLE_ID",
    "GENOMIC_CHANGE",
    "GENOTYPE",
    "DP_CONTROL",
    "GENOME_VERSION",
    "VARIANT_CLASS",
    "SYMBOL",
    "GENENAME",
    "CONSEQUENCE",
    "PROTEIN_CHANGE",
    "CPSR_CLASSIFICATION_SOURCE",
    "FINAL_CLASSIFICATION",
    "CPSR_CLASSIFICATION",
    "CPSR_PATHOGENICITY_SCORE",
    "CPSR_CLASSIFICATION_CODE",
    "CLINVAR_CLASSIFICATION",
    "BM_CANCER_TYPE",
    "BM_DISEASE_ONTOLOGY_ID",
    "BM_PRIMARY_SITE",
    "BM_CLINICAL_SIGNIFICANCE",
    "BM_THERAPEUTIC_CONTEXT",
    "BM_CITATION",
    "BM_RATING",
    "BM_MOLECULAR_PROFILE",
    "BM_EVIDENCE_TYPE",
    "BM_EVIDENCE_LEVEL",
    "BM_EVIDENCE_DIRECTION",
    "BM_EVIDENCE_DESCRIPTION",
    "BM_SOURCE_DB",
    "BM_EVIDENCE_ID",
    "BM_VARIANT_ORIGIN",
    "BM_MATCH",
    "BM_RESOLUTION"
  )


#---- acmg ----#
acmg <- list()
acmg[["score2tier"]] <- data.frame()
acmg[["evidence_codes"]] <-
  utils::read.table(file = "data-raw/acmg_evidence.tsv",
                    header = T, stringsAsFactors = F,
                    comment.char = "", na.strings = c("NA"),
                    sep = "\t")
acmg[["pathogenic_range_gnomad"]] <- list()
acmg[["pathogenic_range_gnomad"]][["af"]] <- 0.0005
acmg[["pathogenic_range_gnomad"]][["min_an"]] <- 12000
acmg[["insilico_pred_min_majority"]] <- 8
acmg[["insilico_pred_max_minority"]] <- 2

acmg[['score_thresholds']] <- list()
acmg[['score_thresholds']][['p_lower']] <- 4.5
acmg[['score_thresholds']][['lp_upper']] <- 4.0
acmg[['score_thresholds']][['lp_lower']] <- 2.0
acmg[['score_thresholds']][['vus_upper']] <- 1.5
acmg[['score_thresholds']][['vus_lower']] <- -1.0
acmg[['score_thresholds']][['lb_upper']] <- -1.5
acmg[['score_thresholds']][['lb_lower']] <- -2.5
acmg[['score_thresholds']][['b_upper']] <- -3.0

acmg[["score2tier"]] <-
  data.frame("CPSR_CLASSIFICATION" = "Pathogenic",
             "CPSR_PATHOGENICITY_SCORE" = "<b>[4.5, ]</b>")
acmg[["score2tier"]] <-
  dplyr::bind_rows(
    acmg[["score2tier"]],
    data.frame("CPSR_CLASSIFICATION" = "Likely Pathogenic",
               "CPSR_PATHOGENICITY_SCORE" = "<b>[2.0, 4.0]</b>"))
acmg[["score2tier"]] <-
  dplyr::bind_rows(
    acmg[["score2tier"]],
    data.frame("CPSR_CLASSIFICATION" = "VUS",
               "CPSR_PATHOGENICITY_SCORE" = "<b>[-1.0, 1.5]</b>"))
acmg[["score2tier"]] <-
  dplyr::bind_rows(
    acmg[["score2tier"]],
    data.frame("CPSR_CLASSIFICATION" = "Likely Benign",
               "CPSR_PATHOGENICITY_SCORE" = "<b>[-2.5, -1.5]</b>"))
acmg[["score2tier"]] <-
  dplyr::bind_rows(
    acmg[["score2tier"]],
    data.frame("CPSR_CLASSIFICATION" = "Benign",
               "CPSR_PATHOGENICITY_SCORE" = "<b>[, -3.0]</b>"))


color_palette <- list()
color_palette[['report']] <- "#007a74"
color_palette[['none']] <- "#8B8989"
color_palette[['genotypes']] <- list()
color_palette[['genotypes']][['values']] <- c('#bdbdbd','#bdbdbd','#007a74','#252525')
color_palette[['genotypes']][['levels']] <- c('undefined','hom_ref','het','hom_alt')


usethis::use_data(color_palette, overwrite = T)
usethis::use_data(acmg, overwrite = T)
usethis::use_data(col_format_output, overwrite = T)




# my_log4r_layout <- function(level, ...) {
#   paste0(format(Sys.time()), " - cpsr-report-generation - ",
#          level, " - ", ..., "\n", collapse = "")
# }
#
# log4r_logger <-
#   log4r::logger(
#     threshold = "INFO", appenders = log4r::console_appender(my_log4r_layout))
#
# # this gets passed on to all the log4r_* functions inside the pkg
# options("PCGRR_LOG4R_LOGGER" = log4r_logger)
#
# panel_zero <- list()
# for(build in c('grch37','grch38')){
#   ref_data <- pcgrr::load_reference_data(
#     pcgr_db_assembly_dir =
#       file.path(
#         "/Users/sigven/project_data/data/data__pcgrdb/dev/pcgrdb",
#         "20250217/data",
#         build),
#     genome_assembly = build
#   )
#
#   set1 <- ref_data$gene$cpg |>
#     dplyr::filter(CPG_SOURCE != "ACMG_SF" &
#                     CPG_SOURCE != "CPIC_PGX_ONCOLOGY") |>
#     dplyr::filter(!is.na(ENSEMBL_GENE_ID)) |>
#     dplyr::inner_join(
#       dplyr::select(ref_data$gene$gene_xref,
#                     ENTREZGENE,
#                     ENSEMBL_GENE_ID,
#                     GENE_BIOTYPE,
#                     GENENAME,
#                     TSG,
#                     TSG_SUPPORT,
#                     ONCOGENE,
#                     ONCOGENE_SUPPORT),
#       by = c("ENTREZGENE","ENSEMBL_GENE_ID")
#     )
#
#   set2 <- ref_data$gene$cpg |>
#     dplyr::filter(CPG_SOURCE != "ACMG_SF" &
#                     CPG_SOURCE != "CPIC_PGX_ONCOLOGY") |>
#     dplyr::filter(is.na(ENSEMBL_GENE_ID)) |>
#     dplyr::select(-c("ENSEMBL_GENE_ID")) |>
#     dplyr::inner_join(
#       dplyr::select(ref_data$gene$gene_xref,
#                     ENTREZGENE,
#                     ENSEMBL_GENE_ID,
#                     GENE_BIOTYPE,
#                     GENENAME,
#                     TSG,
#                     TSG_SUPPORT,
#                     ONCOGENE,
#                     ONCOGENE_SUPPORT),
#       by = c("ENTREZGENE")
#     )
#
#   panel_zero[[build]] <- dplyr::bind_rows(set1, set2) |>
#     dplyr::mutate(
#       PANEL_NAME = "CPSR superpanel of cancer predisposition genes",
#       PANEL_VERSION = "v2025_02") |>
#     dplyr::rename(
#       TUMOR_SUPPRESSOR = TSG,
#       TUMOR_SUPPRESSOR_SUPPORT = TSG_SUPPORT
#     ) |>
#     dplyr::left_join(
#       dplyr::select(
#         dplyr::filter(
#           ref_data$variant$clinvar_gene_stats,
#           .data$CONFIDENCE == "min2goldstars"),
#         c("ENTREZGENE",
#           "N_TRUNC_PATH",
#           "N_NONTRUNC_PATH",
#           "N_MISSENSE_PATH",
#           "N_MISSENSE_BENIGN",
#           "BENIGN_MISSENSE_FRAC",
#           "PATH_TRUNC_FRAC")
#       ), by = "ENTREZGENE"
#     ) |>
#     dplyr::select(
#       dplyr::any_of(
#         c("ENTREZGENE",
#           "SYMBOL",
#           "GENENAME",
#           "GENE_BIOTYPE",
#           "ENSEMBL_GENE_ID",
#           "TUMOR_SUPPRESSOR",
#           "TUMOR_SUPPRESSOR_SUPPORT",
#           "ONCOGENE",
#           "ONCOGENE_SUPPORT",
#           "CPG_SOURCE",
#           "CPG_MOD",
#           "CPG_MOI",
#           "CPG_PHENOTYPES",
#           "CPG_CANCER_CUI",
#           "CPG_SYNDROME_CUI")
#       ),
#       dplyr::everything()
#     ) |>
#     dplyr::distinct()
# }
# #
# workbook <- openxlsx2::wb_workbook() |>
#   openxlsx2::wb_add_worksheet(sheet = "CPSR_SUPERPANEL.GRCH37") |>
#   openxlsx2::wb_add_worksheet(sheet = "CPSR_SUPERPANEL.GRCH38") |>
#   openxlsx2::wb_add_data_table(
#     sheet = "CPSR_SUPERPANEL.GRCH37",
#     x = panel_zero[['grch37']],
#     start_row = 1,
#     start_col = 1,
#     col_names = TRUE,
#     na.strings = "NA",
#     table_style = "TableStyleMedium15") |>
#   openxlsx2::wb_add_data_table(
#     sheet = "CPSR_SUPERPANEL.GRCH38",
#     x = panel_zero[['grch38']],
#     start_row = 1,
#     start_col = 1,
#     col_names = TRUE,
#     na.strings = "NA",
#     table_style = "TableStyleMedium16")
#
# openxlsx2::wb_save(
#   wb = workbook,
#   "pkgdown/assets/cpsr_superpanel_2025_02.xlsx",
#   overwrite = TRUE)
#
# #
# panel_zero_display <- panel_zero$grch38 |>
#   dplyr::select(
#     c("ENTREZGENE",
#       "SYMBOL",
#       "ENTREZGENE",
#       "ENSEMBL_GENE_ID",
#       "GENENAME",
#       "CPG_PHENOTYPES",
#       "CPG_MOI",
#       "CPG_MOD",
#       "CPG_SOURCE",
#     )
#   ) |>
#   dplyr::mutate(
#     CPG_SOURCE = stringr::str_replace_all(
#       CPG_SOURCE, "&", ", "
#     )) |>
#   dplyr::mutate(
#     CPG_SOURCE = stringr::str_replace_all(
#       CPG_SOURCE, "ACMG_SF", ""
#     )
#   ) |>
#   dplyr::mutate(
#     GENE = paste0(
#       "<a href='https://www.ncbi.nlm.nih.gov/gene/",
#       .data$ENTREZGENE,
#       "' target='_blank'>",
#       .data$SYMBOL, "</a>"
#     )
#   ) |>
#   dplyr::select(
#     c("GENE","ENTREZGENE","ENSEMBL_GENE_ID",
#       "CPG_MOD", "CPG_MOI", "GENENAME",
#       "CPG_SOURCE", "CPG_PHENOTYPES")
#   )
# #
# readr::write_tsv(
#   panel_zero_display, file = "inst/extdata/panel_zero.tsv.gz",
#   na = "NA", col_names = T,quote = "none"
# )

