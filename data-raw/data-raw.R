
col_format_output <- list()

col_format_output[['html_tier']] <-
  c(
    "SYMBOL",
    "CLINVAR_PHENOTYPE",
    "CONSEQUENCE",
    "PROTEIN_CHANGE",
    "GENOTYPE",
    "GENENAME",
    "PROTEIN_DOMAIN",
    "HGVSp",
    "HGVSc",
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
    "PROTEIN_CHANGE",
    "GENOTYPE",
    "GENENAME",
    "PROTEIN_DOMAIN",
    "HGVSp",
    "HGVSc",
    "ENSEMBL_GENE_ID",
    "REFSEQ_TRANSCRIPT_ID",
    "ENSEMBL_TRANSCRIPT_ID",
    "CDS_CHANGE",
    "PREDICTED_EFFECT",
    "LOSS_OF_FUNCTION",
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

## define tags/variables to display in data tables (GWAS findings)
col_format_output[['html_gwas']] <-
  c(
    "SYMBOL",
    "CONSEQUENCE",
    "GWAS_PHENOTYPE",
    "PROTEIN_CHANGE",
    "GENOTYPE",
    "LOSS_OF_FUNCTION",
    "GENENAME",
    "PROTEIN_DOMAIN",
    "GWAS_CITATION",
    "HGVSp",
    "HGVSc",
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
    "PROTEIN_CHANGE",
    "PFAM_DOMAIN_NAME",
    "HGVSp",
    "HGVSc",
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
    'PROTEIN_CHANGE',
    'CONSEQUENCE',
    'BM_EVIDENCE_LEVEL',
    'GENOTYPE',
    'BM_CANCER_TYPE',
    'BM_DISEASE_ONTOLOGY_ID',
    'BM_PRIMARY_SITE',
    'BM_CLINICAL_SIGNIFICANCE',
    'BM_THERAPEUTIC_CONTEXT',
    'BM_REFERENCE',
    'BM_MOLECULAR_PROFILE_NAME',
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
    "CDS_CHANGE",
    "CODING_STATUS",
    "MUTATION_HOTSPOT",
    "EFFECT_PREDICTIONS",
    "LOSS_OF_FUNCTION",
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
    "CDS_CHANGE",
    "CODING_STATUS",
    "MUTATION_HOTSPOT",
    "EFFECT_PREDICTIONS",
    "LOSS_OF_FUNCTION",
    "NULL_VARIANT",
    "DBSNP_RSID",
    "gnomADe_AF"
  )


## define tags/variables to display in output Excel
col_format_output[['xlsx_biomarker']] <-
  c("SAMPLE_ID",
    "GENOMIC_CHANGE",
    "GENOTYPE",
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
    "BM_MOLECULAR_PROFILE_NAME",
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

acmg[["score2tier"]] <-
  data.frame("CPSR_CLASSIFICATION" = "Pathogenic",
             "CPSR_PATHOGENICITY_SCORE" = "<b>[5, ]</b>")
acmg[["score2tier"]] <-
  dplyr::bind_rows(
    acmg[["score2tier"]],
    data.frame("CPSR_CLASSIFICATION" = "Likely Pathogenic",
               "CPSR_PATHOGENICITY_SCORE" = "<b>[2.5, 4.5]</b>"))
acmg[["score2tier"]] <-
  dplyr::bind_rows(
    acmg[["score2tier"]],
    data.frame("CPSR_CLASSIFICATION" = "VUS",
               "CPSR_PATHOGENICITY_SCORE" = "<b>[-1.0, 2.0]</b>"))
acmg[["score2tier"]] <-
  dplyr::bind_rows(
    acmg[["score2tier"]],
    data.frame("CPSR_CLASSIFICATION" = "Likely Benign",
               "CPSR_PATHOGENICITY_SCORE" = "<b>[-4.5, -1.5]</b>"))
acmg[["score2tier"]] <-
  dplyr::bind_rows(
    acmg[["score2tier"]],
    data.frame("CPSR_CLASSIFICATION" = "Benign",
               "CPSR_PATHOGENICITY_SCORE" = "<b>[, -5]</b>"))


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
#         "20240309/data",
#         build)
#   )
#   panel_zero[[build]] <- ref_data$gene$cpg |>
#     dplyr::filter(CPG_SOURCE != "ACMG_SF") |>
#     dplyr::mutate(
#       PANEL_NAME = "CPSR superpanel of cancer predisposition genes",
#       PANEL_VERSION = "v2024_03") |>
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
#     ) |>
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
#
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
#   "pkgdown/assets/cpsr_superpanel_2024_03.xlsx",
#   overwrite = TRUE)
